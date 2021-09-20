/**
 * @file ProteinBindStatus.hpp
 * @author wenyan4work (wenyan4work@gmail.com)
 * @brief
 * @version 0.1
 * @date 2019-02-07
 *
 * @copyright Copyright (c) 2019
 *
 */
#ifndef PROTEINBINDSTATUS_HPP_
#define PROTEINBINDSTATUS_HPP_

// SimToolbox
#include "FDPS/particle_simulator.hpp"
#include "Util/EigenDef.hpp"
#include "Util/GeoUtil.hpp"
#include "Util/IOHelper.hpp"

#include <cassert>
#include <limits>
#include <type_traits>

constexpr int ID_UB = -1; ///< the default Unbind state
constexpr double NAND = std::numeric_limits<double>::quiet_NaN(); ///< nan
constexpr double EPSD = std::numeric_limits<double>::epsilon();   ///< epsilon

/**
 * @brief time-varying ProteinBindStatus
 *
 * This class offers only basic setBind() and setUnbind() functions.
 * Other protein behaviors are defined in ProteinData
 *
 */
struct ProteinBindStatus {
  public:
    struct Link {
        int prev = -1;
        int next = -1;
    };

    // TODO: Test this
    double pos[3] = {NAND, NAND, NAND};
    ///< Case1 both unbind -> protein is a point only
    ///< Case2 one bind -> protein is a point only, pos = posEndBind
    ///< Case3 both bind -> pos = mid-point of posEndBind

    // time varying properties
    bool changeBind[2] = {true, true}; ///< flag for if binding status changes
    int idBind[2] = {ID_UB, ID_UB};    ///< global ID of bind MT
    int indexBind[2] = {ID_UB, ID_UB}; ///< globalIndex of bind MT
    int rankBind[2] = {ID_UB, ID_UB};  ///< mpi rank of bind MT
    double lenBind[2] = {NAND, NAND};  ///< length of bind MT
    Link links[2];                     ///< Linked list to flexible filament

    double distBind[2] = {NAND, NAND};
    ///< the distance to bind MT center,
    ///< [-lenBind/2,lenBind/2], +towards plus end

    double posEndBind[2][3] = {{NAND, NAND, NAND}, {NAND, NAND, NAND}};
    ///< position (in lab frame) of two ends when bind
    ///< when unbind, the position is NAND

    double centerBind[2][3] = {{NAND, NAND, NAND}, {NAND, NAND, NAND}};
    ///< center position of two MTs

    double directionBind[2][3] = {{NAND, NAND, NAND}, {NAND, NAND, NAND}};
    ///< direction of two MTs

    /**
     * @brief clear bind status for both ends.
     *  FDPS interface
     */
    void clear() {
        changeBind[0] = true;
        changeBind[1] = true;
        setUnBind(0);
        setUnBind(1);
        std::fill(pos, pos + 3, NAND);
    }

    /**
     * @brief Set the end to unbind status
     *
     * @param end 0 or 1
     */
    void setUnBind(int end) {
        idBind[end] = ID_UB;
        indexBind[end] = ID_UB;
        rankBind[end] = ID_UB;
        lenBind[end] = NAND;
        distBind[end] = NAND;
        for (int i = 0; i < 3; i++) {
            directionBind[end][i] = NAND;
            centerBind[end][i] = NAND;
            posEndBind[end][i] = NAND;
        }
        // TODO unset Link
    }

    /**
     * @brief Set position for completely unbound protein
     * 
     * @param newPos 
     */
    void setUnBindPos(double newPos[3]) {
        pos[0] = newPos[0];
        pos[1] = newPos[1];
        pos[2] = newPos[2];
    }

    /**
     * @brief Set the bind status
     * 
     * @param end 
     * @param gid 
     * @param index 
     * @param directionLine 
     * @param centerLine 
     * @param centerDist 
     * @param length 
     * @param rank 
     */
    void setBind(int end, const int gid, const int index,
                 const double directionLine[3], const double centerLine[3],
                 const double centerDist, const double length, const int rank) {
        assert(end == 0 || end == 1); // Make sure you choose a viable head
        assert(idBind[end] == ID_UB); // Make sure end is originally unbound
        assert(gid != ID_UB);
        /* TODO: Set links <05-03-21, ARL> */

        idBind[end] = gid;
        indexBind[end] = index;
        distBind[end] = centerDist;
        lenBind[end] = length;
        rankBind[end] = rank;
        std::copy(directionLine, directionLine + 3, directionBind[end]);
        std::copy(centerLine, centerLine + 3, centerBind[end]);
        updatePosEndBind(end);
    }

    /**
     * @brief calculate the position of end with bind information
     *
     * @param end
     */
    void updatePosEndBind(const int end) {
        if (idBind[end] == ID_UB) {
            setUnBind(end);
        } else if (lenBind[end] == 0) {
            std::copy(posEndBind[end], posEndBind[end] + 3, centerBind[end]);
        } else {
            for (int i = 0; i < 3; i++) { // posEnd = direction * dist + center
                posEndBind[end][i] =
                    directionBind[end][i] * distBind[end] + centerBind[end][i];
            }
        }
    }

    /**
     * @brief
     * Case 1: if doubly bound, set pos to center of two ends.
     * Case 2: if singly bound, set pos to that end.
     * Case 3: if unbound, do nothing
     * This must be called when posEndBind is valid
     */
    void updatePosWithEndBind() {
        if (idBind[0] != ID_UB && idBind[1] == ID_UB) { // Case 2
            std::copy(posEndBind[0], posEndBind[0] + 3, pos);
        } else if (idBind[0] == ID_UB && idBind[1] != ID_UB) { // Case 2
            std::copy(posEndBind[1], posEndBind[1] + 3, pos);
        } else if (idBind[0] != ID_UB && idBind[1] != ID_UB) { // Case 1
            for (int i = 0; i < 3; i++) {
                pos[i] = 0.5 * (posEndBind[0][i] + posEndBind[1][i]);
            }
        }
    }

    /**
     * @brief clamp posEnd to MT range
     *
     * @param end
     */
    void updatePosEndClamp(int end) {
        if (idBind[end] != ID_UB) {
            const double lenHalf = lenBind[end] * 0.5;
            if (distBind[end] > lenHalf)
                distBind[end] = lenHalf;
            if (distBind[end] < -lenHalf)
                distBind[end] = -lenHalf;
        }
    }
};

static_assert(std::is_trivially_copyable<ProteinBindStatus>::value, "");
static_assert(std::is_default_constructible<ProteinBindStatus>::value, "");

#endif
