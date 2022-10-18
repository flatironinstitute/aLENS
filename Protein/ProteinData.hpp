/**
 * @file ProteinData.hpp
 * @author wenyan4work (wenyan4work@gmail.com)
 * @brief
 * @version 0.1
 * @date 2019-01-04
 *
 * @copyright Copyright (c) 2019
 *
 */

#ifndef PROTEINDATA_HPP_
#define PROTEINDATA_HPP_

#include "ProteinBindStatus.hpp"
#include "ProteinType.hpp"

// KMC submodule
#include "KMC/helpers.hpp"
#include "KMC/lookup_table.hpp"
#include "KMC/lut_filler.hpp"
#include "KMC/macros.hpp"

// SimToolbox submodule
#include "Collision/DCPQuery.hpp"
#include "FDPS/particle_simulator.hpp"
#include "Util/EigenDef.hpp"
#include "Util/GeoUtil.hpp"
#include "Util/IOHelper.hpp"

#include <cassert>
#include <cmath>
#include <limits>
#include <type_traits>

/**
 * @brief FullParticle type for Protein managed by PS::ParticleSystem
 *
 *  TODO: clear up this class
 * Naming convention:
 * setXXX(), getXXX(): usual set&get
 * updateXXX(): update member state, return void
 * calcXXX() const: calc with some input, return value
 */
struct ProteinData {
  public:
    int gid;                 ///< unique global ID
    ProteinType property;    ///< intrinsic properties
    ProteinBindStatus bind;  ///< time-varying bind status
    double forceBind[2][3];  ///< spring force on bind MTs
    double torqueBind[2][3]; ///< spring torque on bind MTs

    /**
     * @brief Set default values for protein
     *
     * @return double
     */
    void setMockProtein() {
        bind.clear();
        property.setMockProteinType();
        std::fill(bind.pos, bind.pos + 3, 0);
    }

    /**********************************
     *
     * Get() const Functions
     *
     ***********************************/

    /**
     * @brief Get the pointer to doubel pos[3]
     *
     * @return const double*
     */
    const double *getPosPtr() const { return bind.pos; }

    /**
     * @brief Get gid
     *
     * @return int
     */
    int getGid() const { return gid; }

    /**
     * @brief Get BindingFactor from Unbound to Singly bound
     *
     * @param e
     * @return double
     */
    double getBindingFactorUS(int e) const {
        assert(e == 0 || e == 1);
        return property.ko_s[e] * property.Ka[e] * property.eps;
    }

    /**
     * @brief Get UnbindingFactor from Singly bound to Unbound
     *
     * @param e
     * @return double
     */
    double getUnbindingFactorSU(int e) const { return property.ko_s[e]; }

    /**
     * @brief Get BindingFactor from Singly to Doubly bound
     *
     * @param e Unbound end index
     * @param direction Direction of unbound tubule
     * @return double
     */
    double getBindingFactorSD(int e, const double direction[3]) const {
        // If angle between bound and unbound tubule is less than 90 degrees and
        // greater than -90, then use parallel binding affinity
        // Prefact should never exceed 1
        bool parallel = sgn(dot3(direction, bind.directionBind[1 - e])) > 0;
        double prefact = 1.;
        if (parallel && property.PtoAPratio < 1.) {
            prefact = property.PtoAPratio;
        } else if (!parallel && property.PtoAPratio > 1.) {
            prefact /= property.PtoAPratio;
        }
        assert(prefact >= 0 && prefact <= 1.);
        return prefact * property.ko_d[e] * property.Ke[e] * property.eps;
    }

    /**
     * @brief Get UnbindingFactor from Doubly to Singly boun
     * 
     * @param e 
     * @param KBT 
     * @return double 
     */
    double getUnbindingFactorDS(int e, double KBT) const {
        if (property.lambda == 0)
            return property.ko_d[e];
        else {
            double M = property.lambda * .5 * property.kappa *
                       SQR(getProteinForceLength() - property.freeLength) / KBT;
            return property.ko_d[e] * exp(M);
        }
    }

    /**
     * @brief Get cut off calculation range when calculating 0 to 1 head bound
     *
     * @return double
     */
    double getRcutUS() const { return property.rc; }

    /**
     * @brief Get cut off calculation range when calculating 1 to 2 head bound
     * @return double
     */
    double getRcutSD() const {
        // use the cutoff length as LUT construction
        return property.LUTablePtr->getDsbound();
    }

    /*! \brief Return the ID of the rod protein end is bound to
     *
     *
     * \param end Which protein end
     * \return ID of rod protein is bound to
     */
    int getBindID(int end) const { return bind.idBind[end]; }

    /**
     * @brief Get the ProteinForceLength, subtracting tubule diameter
     *
     * @return double the calculated length
     */
    double getProteinForceLength() const {
        double length = 0;
        if (bind.idBind[0] == ID_UB || bind.idBind[1] == ID_UB) {
            length = 0;
        } else {
            // consistent tubuleDiameter as the original LUT construction
            double tubuleD = property.LUTablePtr->getRodDiameter();
            double ell = getProteinEndEndLength();
            if (property.lookupType == 0) {
                length = ell - tubuleD;
            } else if (property.lookupType == 1) {
                Evec3 point = ECmap3(bind.posEndBind[0]);
                Evec3 center = ECmap3(bind.centerBind[1]);
                Evec3 direct = ECmap3(bind.directionBind[1]);
                direct.normalize();
                Evec3 minus = center - (0.5 * bind.lenBind[1]) * direct;
                Evec3 plus = center + (0.5 * bind.lenBind[1]) * direct;
                Evec3 foot;
                double d = DistPointSeg(point, minus, plus, foot);
                // corrections to tubule radius
                length = (1 - tubuleD / d) * ell;
            }
        }
        return length;
    }

    /**
     * @brief Get the ProteinEndEndLength
     *
     * @return double the calculated length
     */
    double getProteinEndEndLength() const {
        if (bind.idBind[0] == ID_UB || bind.idBind[1] == ID_UB) {
            return 0;
        } else {
            Evec3 r = ECmap3(bind.posEndBind[0]) - ECmap3(bind.posEndBind[1]);
            return r.norm();
        }
    }

    /**
     * @brief Get the reference to Protein Property
     *
     * @return const double&
     */
    double getDiffU() const { return property.diffUnbound; }

    /**
     * @brief Get the reference to Protein Property
     *
     * @return const ProteinType&
     */
    const ProteinType &getProperty() const { return property; }

    /**
     * @brief get the ptr to LookupTable object
     *
     * @return LookupTable*
     */
    const LookupTable *getLUTablePtr() const { return property.LUTablePtr; }

    /**
     * @brief Get if the protein is walking along MT or not
     *
     * @return true
     * @return false
     */
    bool getWalkOrNot() const {
        return (bind.idBind[0] != ID_UB || bind.idBind[1] != ID_UB);
    }

    /**********************************
     *
     * Set() Functions
     *
     ***********************************/

    /**
     * @brief Set new position
     * Binding location is also jumped
     * @param newPos
     */
    void setPos(const double newPos[3]) {
        // current rvec
        Evec3 jump = ECmap3(newPos) - Emap3(bind.pos);

        // set new pos
        std::copy(newPos, newPos + 3, bind.pos);

        // move bind center if bind
        for (int e = 0; e < 2; e++) {
            if (bind.idBind[e] != ID_UB) {
                Emap3(bind.centerBind[e]) += jump;
            }
        }
        // update posEndBind and posProtein
        updateGeometryWithBind();
    }

    /**
     * @brief Set protein status from file information
     * 
     * WARNING: after this function, the status of this protein is incomplete. 
     * the bind info must be updated with idBind before simulation steps
     * 
     * @param gid_ 
     * @param tag_ 
     * @param end0 
     * @param end1 
     * @param idBind
     * @param property_ 
     */
    void setFromFileInput(const int gid_, const int tag_, const double end0[3],
                          const double end1[3], const int idBind[2],
                          const ProteinType &property_) {
        bind.clear();
        gid = gid_;
        for (int k = 0; k < 3; k++) {
            // set pos to prepare for bind reconstruction
            bind.posEndBind[0][k] = end0[k];
            bind.posEndBind[1][k] = end1[k];
            bind.pos[k] = 0.5 * (end0[k] + end1[k]);
        }
        bind.idBind[0] = idBind[0];
        bind.idBind[1] = idBind[1];
        property = property_;
        if (property.tag != tag_) {
            printf("protein tag error. dat/yaml file mismatch\n");
            std::exit(1);
        }
    }

    /**********************************
     *
     * Update() Functions
     *
     ***********************************/

    /**
     * @brief update position by walking (S or D)
     * 
     * @param KBT 
     * @param dt 
     * @param U01 
     * @param N01a 
     * @param N01b 
     */
    void updatePosWalk(double KBT, double dt, double U01, double N01a,
                       double N01b) {
        assert(getWalkOrNot());
        if (bind.idBind[0] != ID_UB && bind.idBind[1] != ID_UB) {
            double v0 = calcEndWalkVelocity(0) + calcEndDragVelocity(0, KBT);
            double v1 = calcEndWalkVelocity(1) + calcEndDragVelocity(1, KBT);
            const double diff0 = sqrt(dt * property.diffBoundD[0] * 2);
            const double diff1 = sqrt(dt * property.diffBoundD[1] * 2);
            bind.distBind[0] += v0 * dt + diff0 * N01a;
            bind.distBind[1] += v1 * dt + diff1 * N01b;
        } else if (bind.idBind[0] != ID_UB && bind.idBind[1] == ID_UB) {
            double v0 = property.vmax[0];
            const double diff0 = sqrt(dt * property.diffBoundS[0] * 2);
            bind.distBind[0] += v0 * dt + diff0 * N01a;
        } else if (bind.idBind[0] == ID_UB && bind.idBind[1] != ID_UB) {
            double v1 = property.vmax[1];
            const double diff1 = sqrt(dt * property.diffBoundS[1] * 2);
            bind.distBind[1] += v1 * dt + diff1 * N01b;
        } else {
            printf("walk while unbound is error\n");
            exit(1);
        }
        if (!property.walkOff) {
            updateWalkClampEnd();
        } else {
            updateWalkOffEnd(U01);
        }
        updateGeometryWithBind();
    }

    /**
     * @brief update bind when walk off
     *
     */
    void updateWalkOffEnd(double U01) {

        bool endOff[2] = {false, false};
        for (int end = 0; end < 2; end++) {
            if (bind.distBind[end] > bind.lenBind[end] * 0.5 ||
                bind.distBind[end] < bind.lenBind[end] * (-0.5)) {
                endOff[end] = true;
            }
        }

        if (endOff[0] && endOff[1]) {
            // both off, choose to unbind 1
            if (U01 < 0.5) {
                bind.setUnBind(0);
                bind.updatePosEndClamp(1);
            } else {
                bind.setUnBind(1);
                bind.updatePosEndClamp(0);
            }
        } else if (endOff[0]) {
            bind.setUnBind(0);
        } else if (endOff[1]) {
            bind.setUnBind(1);
        }
    }

    /**
     * @brief clamp distBind to end
     *
     */
    void updateWalkClampEnd() {
        bind.updatePosEndClamp(0);
        bind.updatePosEndClamp(1);
    }

    /**
     * @brief update position by diffusion (U)
     *
     * @param dt timestep size
     * @param N01x rngN01
     * @param N01y rngN01
     * @param N01z rngN01
     */
    void updatePosDiffuse(double dt, double N01x, double N01y, double N01z) {
        assert(!getWalkOrNot());
        const double diff = sqrt(dt * property.diffUnbound * 2);
        bind.pos[0] += diff * N01x;
        bind.pos[1] += diff * N01y;
        bind.pos[2] += diff * N01z;
    }

    /**
     * @brief update bind.posEndBind and bind.pos
     *
     */
    void updateGeometryWithBind() {
        bind.updatePosEndBind(0);
        bind.updatePosEndBind(1);
        bind.updatePosWithEndBind();
    }

    /**
     * @brief compute the spring binding force
     *
     * @return double
     */
    void updateForceTorqueBind() {
        if (bind.idBind[0] != ID_UB && bind.idBind[1] != ID_UB) {
            Evec3 r = Emap3(bind.posEndBind[0]) - Emap3(bind.posEndBind[1]);
            double force = (getProteinForceLength() - property.freeLength) *
                           property.kappa;
            Evec3 f0 = -force * r.normalized();
            Evec3 f1 = -f0;
            // torque = r x f
            Evec3 torque0 =
                bind.distBind[0] * (Emap3(bind.directionBind[0]).cross(f0));
            Evec3 torque1 =
                bind.distBind[1] * (Emap3(bind.directionBind[1]).cross(f1));
            for (int k = 0; k < 3; k++) {
                forceBind[0][k] = f0[k];
                forceBind[1][k] = f1[k];
                torqueBind[0][k] = torque0[k];
                torqueBind[1][k] = torque1[k];
            }
        } else {
            // spring bind force is zero
            forceBind[0][0] = forceBind[0][1] = forceBind[0][2] = 0;
            forceBind[1][0] = forceBind[1][1] = forceBind[1][2] = 0;
            torqueBind[0][0] = torqueBind[0][1] = torqueBind[0][2] = 0;
            torqueBind[1][0] = torqueBind[1][1] = torqueBind[1][2] = 0;
        }
    }

    /**********************************
     *
     * Calc() const {} functions
     *
     ***********************************/

    /**
     * @brief compute the walking velocity of one end
     *
     * positive vmax towards plus end
     *
     * @param end 0 or 1
     * @return double computed walking velocity
     */
    double calcEndWalkVelocity(int end) const {
        assert(bind.idBind[end] != ID_UB);
        // sgn_fac = +1 for moving towards plus end, -1 for moving towards minus end
        double sgn_fac = sgn(property.vmax[end]);
        double fproj =
            sgn_fac *
            ECmap3(forceBind[end]).dot(ECmap3(bind.directionBind[end])) /
            std::abs(property.fstall);
        double vfrac = std::max(0.0, std::min(1.0, 1.0 + fproj));
        if (bind.idBind[1 - end] == ID_UB) { // only this end is bound
            return property.vmax[end] * vfrac;
        } else {
            double vel =
                ECmap3(bind.directionBind[end])
                            .dot(ECmap3(bind.directionBind[1 - end])) > 0
                    ? property.vmax[end]
                    : property.vmaxAP[end];
            return vel * vfrac;
        }
    }

    /**
     * @brief compute the velocity of a passive diffusing head being dragged
     *
     * positive v towards plus end
     *
     * @param end 0 or 1
     * @param KBT Boltzmann constant * temperature
     * @return double computed drag velocity from force
     */
    double calcEndDragVelocity(int end, double KBT) const {
        assert(bind.idBind[end] != ID_UB);
        if (property.vdrag[end]) {
            double fproj =
                ECmap3(forceBind[end]).dot(ECmap3(bind.directionBind[end]));
            return fproj * property.diffBoundD[end] / KBT;
        } else {
            return 0;
        }
    }

    /**********************************
     *
     * FDPS interface Functions
     *
     ***********************************/

    /**
     * @brief Get the Pos data field
     *  interface requried for FDPS
     * @return PS::F64vec3
     */
    PS::F64vec3 getPos() const {
        return PS::F64vec(bind.pos[0], bind.pos[1], bind.pos[2]);
    }

    /**
     * @brief get neighbor search radius interface requried for FDPS
     *
     * @return double
     */
    double getRSearch() const { return property.freeLength * 10; }

    /**
     * @brief copyFromFP()
     *  interface requried for FDPS
     * @param fp
     */
    void copyFromFP(const ProteinData &fp) { *this = fp; }

    /**
     * @brief Set the Pos data field
     *  interface requried for FDPS
     * @param newPos
     */
    void setPos(const PS::F64vec3 &newPos) {
        double newPos_[3] = {newPos.x, newPos.y, newPos.z};
        setPos(newPos_);
    }

    /**
     * @brief writeAscii file
     *  interface requried for FDPS IO
     * both cases can be read with the same routine
     * @param fptr
     */
    void writeAscii(FILE *fptr) const {
        if (bind.idBind[0] != ID_UB && bind.idBind[1] != ID_UB) {
            // protein has finite length
            // this should NOT out put nan
            fprintf(fptr, "P %d %d %.6g %.6g %.6g %.6g %.6g %.6g %d %d\n", //
                    gid, property.tag,                                     //
                    bind.posEndBind[0][0], bind.posEndBind[0][1],
                    bind.posEndBind[0][2], //
                    bind.posEndBind[1][0], bind.posEndBind[1][1],
                    bind.posEndBind[1][2], //
                    bind.idBind[0], bind.idBind[1]);
        } else {
            // protein has zero length
            fprintf(fptr, "P %d %d %.6g %.6g %.6g %.6g %.6g %.6g %d %d\n", //
                    gid, property.tag,                                     //
                    bind.pos[0], bind.pos[1], bind.pos[2],                 //
                    bind.pos[0], bind.pos[1], bind.pos[2],                 //
                    bind.idBind[0], bind.idBind[1]);
        }
    }

    /**********************************
     *
     * VTK format IO Functions
     *
     ***********************************/

    /**
     * @brief write VTK XML Parallel VTP file
     * 
     * @tparam Container 
     * @param protein 
     * @param proteinNumber 
     * @param prefix 
     * @param postfix 
     * @param rank 
     */
    template <class Container>
    static void writeVTP(const Container &protein, const int proteinNumber,
                         const std::string &prefix, const std::string &postfix,
                         int rank) {
        // each protein is a vtk polyline with two vertices and one cell

        // write VTP for basic data
        //  use float to save some space
        // point and point data
        // position always in Float64
        std::vector<double> pos(6 * proteinNumber);

        // point connectivity of line
        std::vector<int32_t> connectivity(2 * proteinNumber);
        std::vector<int32_t> offset(proteinNumber);

        // protein data
        std::vector<int32_t> gid(proteinNumber);
        std::vector<int32_t> tag(proteinNumber);
        std::vector<int32_t> idBind(2 * proteinNumber);

#pragma omp parallel for
        for (int i = 0; i < proteinNumber; i++) {
            const auto &p = protein[i];
            // point and point data
            if (p.bind.idBind[0] != ID_UB && p.bind.idBind[1] != ID_UB) {
                // if both bind, posEndBind is valid
                pos[6 * i + 0] = p.bind.posEndBind[0][0];
                pos[6 * i + 1] = p.bind.posEndBind[0][1];
                pos[6 * i + 2] = p.bind.posEndBind[0][2];
                pos[6 * i + 3] = p.bind.posEndBind[1][0];
                pos[6 * i + 4] = p.bind.posEndBind[1][1];
                pos[6 * i + 5] = p.bind.posEndBind[1][2];
            } else {
                // else, shrink a point
                pos[6 * i + 0] = p.bind.pos[0];
                pos[6 * i + 1] = p.bind.pos[1];
                pos[6 * i + 2] = p.bind.pos[2];
                pos[6 * i + 3] = p.bind.pos[0];
                pos[6 * i + 4] = p.bind.pos[1];
                pos[6 * i + 5] = p.bind.pos[2];
            }

            // connectivity
            connectivity[2 * i] = 2 * i;         // index of point 0 in line
            connectivity[2 * i + 1] = 2 * i + 1; // index of point 1 in line
            // offset is the end of each line. in fortran indexing
            offset[i] = 2 * i + 2;

            // protein data
            // point data
            idBind[2 * i] = p.bind.idBind[0];
            idBind[2 * i + 1] = p.bind.idBind[1];
            // cell data
            gid[i] = p.gid;
            tag[i] = p.property.tag;
        }

        std::ofstream file(prefix + std::string("Protein_") + "r" +
                               std::to_string(rank) + std::string("_") +
                               postfix + std::string(".vtp"),
                           std::ios::out);

        IOHelper::writeHeadVTP(file);

        file << "<Piece NumberOfPoints=\"" << proteinNumber * 2
             << "\" NumberOfLines=\"" << proteinNumber << "\">\n";
        // Points
        file << "<Points>\n";
        IOHelper::writeDataArrayBase64(pos, "position", 3, file);
        file << "</Points>\n";
        // cell definition
        file << "<Lines>\n";
        IOHelper::writeDataArrayBase64(connectivity, "connectivity", 1, file);
        IOHelper::writeDataArrayBase64(offset, "offsets", 1, file);
        file << "</Lines>\n";
        // point data
        file << "<PointData Scalars=\"scalars\">\n";
        IOHelper::writeDataArrayBase64(idBind, "idBind", 1, file);
        file << "</PointData>\n";
        // cell data
        file << "<CellData Scalars=\"scalars\">\n";
        IOHelper::writeDataArrayBase64(gid, "gid", 1, file);
        IOHelper::writeDataArrayBase64(tag, "tag", 1, file);
        file << "</CellData>\n";
        file << "</Piece>\n";

        IOHelper::writeTailVTP(file);
        file.close();
    }

    /**
     * @brief write VTK XML Parallel Header
     *
     * @param prefix
     * @param postfix
     * @param nProcs
     */
    static void writePVTP(const std::string &prefix, const std::string &postfix,
                          const int nProcs) {
        std::vector<std::string> pieceNames;

        std::vector<IOHelper::FieldVTU> pointDataFields;
        pointDataFields.emplace_back(1, IOHelper::IOTYPE::Int32, "idBind");

        std::vector<IOHelper::FieldVTU> cellDataFields;
        cellDataFields.emplace_back(1, IOHelper::IOTYPE::Int32, "gid");
        cellDataFields.emplace_back(1, IOHelper::IOTYPE::Int32, "tag");

        for (int i = 0; i < nProcs; i++) {
            pieceNames.emplace_back(std::string("Protein_") + std::string("r") +
                                    std::to_string(i) + "_" + postfix + ".vtp");
        }

        IOHelper::writePVTPFile(prefix + "Protein_" + postfix + ".pvtp",
                                pointDataFields, cellDataFields, pieceNames);
    }
};

static_assert(std::is_trivially_copyable<ProteinData>::value, "");
static_assert(std::is_default_constructible<ProteinData>::value, "");

/**
 * @brief FDPS writeAscii file header
 */
class ProteinAsciiHeader {
  public:
    int nparticle; ///< global number of particles to write
    double time;   ///< current time
    /**
     * @brief write the dat file header
     * 
     * @param fp 
     */
    void writeAscii(FILE *fp) const {
        fprintf(fp, "%d \n %lf\n", nparticle, time);
    }
};

#endif
