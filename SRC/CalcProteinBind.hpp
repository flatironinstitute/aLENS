/**
 * @file CalcProteinBind.hpp
 * @author Wen Yan (wenyan4work@gmail.com) and Adam Lamson
 * @brief
 * @version 0.1
 * @date 2019-01-04
 *
 * @copyright Copyright (c) 2019
 *
 */

#ifndef CALCPROTEINBIND_HPP_
#define CALCPROTEINBIND_HPP_

#include "Protein/ProteinData.hpp"
#include "TubuleBind.hpp"
#include "kmc_step_funcs.hpp"

// KMC module
#include "KMC/kmc.hpp"
// SimToolbox module
#include "MPI/MixPairInteraction.hpp"
#include "Util/Logger.hpp"

/**
 * @brief The functor class for KMC protein-tubule binding
 *
 */
class CalcProteinBind {
    double dt;                            ///< timestep size
    double KBT;                           ///< KBT for protein KMC calculation
    std::shared_ptr<TRngPool> rngPoolPtr; ///< rng generator

  public:
    /**
     * @brief Construct a new CalcProteinBind object
     *
     * @param dt_ Time step
     * @param KBT_ Temperature used in binding kinetics
     * @param rngPoolPtr_ Random number generator
     */
    CalcProteinBind(double dt_, double KBT_,
                    std::shared_ptr<TRngPool> &rngPoolPtr_) {
        dt = dt_;
        KBT = KBT_;
        rngPoolPtr = rngPoolPtr_;
    }

    /**
     * @brief Functor interface required by MixPairInteraction
     *
     * @param trgPtr Array of objects that will bind to srcs (i.e. Proteins)
     * @param nTrg Number of binding objects
     * @param srcPtr Array of objects that will be bound by trgs (i.e. Spheres and Sylinders)
     * @param nSrc Number of to-be-bound objects
     * @param mixForcePtr Array of current bind states of srcs
     */
    template <class Tubule>
    void operator()(const MixEPI<ProteinData> *const trgPtr, const PS::S32 nTrg,
                    const MixEPJ<Tubule> *const srcPtr, const PS::S32 nSrc,
                    ProteinBindStatus *const mixForcePtr) {
        const int threadID = omp_get_thread_num();

        // Build Sylinder pointer vector used for all protein bindings
        std::vector<const Tubule *> srcPtrArr;
        for (int s = 0; s < nSrc; s++) {
            if (srcPtr[s].srcFlag)
                srcPtrArr.push_back(&(srcPtr[s].epSrc));
        }

        for (int t = 0; t < nTrg; t++) {
            auto &trg = trgPtr[t];
            if (!trg.trgFlag) {
                continue;
            }
            // make a modifiable copy of target pData
            auto pData = trg.epTrg;
            auto &bindStatus = pData.bind;
            // update information of already bound MT
            bool bindFound[2] = {false, false};
            for (const auto &ptr : srcPtrArr) {
                for (int end = 0; end < 2; end++) {
                    if (bindStatus.idBind[end] == ptr->gid) {
                        std::copy(ptr->pos, ptr->pos + 3,
                                  bindStatus.centerBind[end]);
                        std::copy(ptr->direction, ptr->direction + 3,
                                  bindStatus.directionBind[end]);
                        bindStatus.lenBind[end] = ptr->length;
                        bindStatus.rankBind[end] = ptr->rank;
                        bindStatus.indexBind[end] = ptr->globalIndex;
                        bindFound[end] = true;
                    }
                }
            }
            if (bindFound[0] == false && pData.property.fixedEnd0 == false) {
#ifndef NDEBUG
                if (bindStatus.idBind[0] != ID_UB)
                    spdlog::error(
                        "Unbinding end 0 of protein {} because sylinder {} "
                        "could not be found.",
                        pData.gid, bindStatus.idBind[0]);
#endif
                bindStatus.setUnBind(0);
            }
            if (bindFound[1] == false) {
#ifndef NDEBUG
                if (bindStatus.idBind[1] != ID_UB)
                    spdlog::error(
                        "Unbinding end 1 of protein {} because sylinder {} "
                        "could not be found.",
                        pData.gid, bindStatus.idBind[1]);
#endif
                bindStatus.setUnBind(1);
            }
            bindStatus.updatePosEndBind(0);
            bindStatus.updatePosEndBind(1);
            bindStatus.updatePosWithEndBind();

            auto &bindStatusResult = mixForcePtr[t];
            bindStatusResult = bindStatus;
            // Up to this point pData.bind is consistent.
            // ALL of the following functions should use pData.bind
            // and store results in bindStatusResult
            // NOT trg.epTrg
            int stage = -1;
            // trg.epTrg is a ProteinData object
            bool end1isUB = (bindStatus.idBind[0] == ID_UB);
            bool end2isUB = (bindStatus.idBind[1] == ID_UB);
            if (end1isUB && end2isUB)
                stage = 0;
            else if (end1isUB || end2isUB)
                stage = 1;
            else if (!(end1isUB && end2isUB))
                stage = 2;

            assert(stage != -1); // A stage should always be found
            double roll[3] = {}; // rng to be setup

            switch (stage) {
            case 0:
                // Unbound protein -> 1 head bound protein
                // or
                // Unbound protein -> Unbound protein
                roll[0] = rngPoolPtr->getU01(threadID);
                KMC_U(pData, srcPtrArr, dt, roll[0], bindStatusResult);
                break;
            case 1:
                // 1 head bound protein -> Unbound protein
                // or
                // 1 head bound protein -> 2 head bound protein
                // or
                // 1 head bound protein -> 1 head bound protein
                for (int i = 0; i < 3; ++i) {
                    roll[i] = rngPoolPtr->getU01(threadID);
                }
                KMC_S(pData, srcPtrArr, dt, KBT, roll, bindStatusResult);
                break;
            case 2:
                // 2 head bound protein -> 1 head bound protein
                // or
                // 2 head bound protein -> 2 head bound protein
                roll[0] = rngPoolPtr->getU01(threadID);
                KMC_D(pData, srcPtrArr, dt, KBT, roll[0], bindStatusResult);
                break;
            default:
                spdlog::critical(
                    " Could not execute correct stage in CalcProteinBind.");
                exit(1);
            }
        }
    }
};

static_assert(std::is_trivially_copyable<ProteinData>::value, "");
static_assert(std::is_default_constructible<ProteinData>::value, "");

static_assert(std::is_trivially_copyable<ProteinBindStatus>::value, "");
static_assert(std::is_default_constructible<ProteinBindStatus>::value, "");

#endif
