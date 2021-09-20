/**
 * @file TubuleBind.hpp
 * @author Wen Yan (wenyan4work@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2021-03-17
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef TUBULEBIND_HPP_
#define TUBULEBIND_HPP_

#include "SimToolbox/FDPS/particle_simulator.hpp"
#include "SimToolbox/Sylinder/Sylinder.hpp"
#include "SimToolbox/Util/EigenDef.hpp"

/**
 * @brief Essential type for Tubule binding calculations 
 *
 * Essential Particle Class for FDPS
 */
struct TubuleBindEP {
  public:
    int gid;         ///< global unique id
    int globalIndex; ///< sequentially ordered unique index in sylinder map
    int rank;        ///< mpi rank of owning rank
    double radius;   ///< radius
    double length;   ///< length

    double pos[3];       ///< position
    double direction[3]; ///< direction (unit norm vector)

    /**
     * @brief Get gid
     *
     * @return int
     */
    int getGid() const { return gid; }

    /**
     * @brief Get global index (sequentially ordered in sylinder map)
     *
     * @return int
     */
    int getGlobalIndex() const { return globalIndex; }

    /**
     * @brief copy data fields from full type Sylinder
     *
     * interface for FDPS
     * @param fp
     */
    void copyFromFP(const Sylinder &fp) {
        gid = fp.gid;
        globalIndex = fp.globalIndex;
        rank = fp.rank;

        radius = fp.radius;
        length = fp.length;

        std::copy(fp.pos, fp.pos + 3, pos);
        Evec3 q = ECmapq(fp.orientation) * Evec3(0, 0, 1);
        direction[0] = q[0];
        direction[1] = q[1];
        direction[2] = q[2];
    }

    /**
     * @brief Get pos as a PS::F64vec3 object
     *
     * interface for FDPS
     * @return PS::F64vec
     */
    PS::F64vec getPos() const { return PS::F64vec3(pos[0], pos[1], pos[2]); }

    /**
     * @brief get search radius
     *
     * interface for FDPS
     * FDPS does not support search with rI+rJ.
     * Here length*2 ensures contact is detected with Symmetry search mode
     * @return PS::F64
     */
    PS::F64 getRSearch() const { return (length * 0.5 + radius) * 1.2; }

    /**
     * @brief Set pos with a PS::F64vec3 object
     *
     * interface for FDPS
     * @param newPos
     */
    void setPos(const PS::F64vec3 &newPos) {
        pos[0] = newPos.x;
        pos[1] = newPos.y;
        pos[2] = newPos.z;
    }

    bool isSphere() const { return length < radius * 2; }
};

static_assert(std::is_trivially_copyable<TubuleBindEP>::value, "");
static_assert(std::is_default_constructible<TubuleBindEP>::value, "");

#endif