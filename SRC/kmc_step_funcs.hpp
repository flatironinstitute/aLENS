#ifndef KMC_STEP_FUNCS_HPP
#define KMC_STEP_FUNCS_HPP

#include "kmc_step_funcs.hpp"

#include "Protein/ProteinBindStatus.hpp"
#include "Protein/ProteinData.hpp"
#include "Protein/ProteinType.hpp"

// KMC module
#include "KMC/helpers.hpp"
#include "KMC/kmc.hpp"
#include "KMC/kmc_choose.hpp"
#include "KMC/lookup_table.hpp"
#include "KMC/macros.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <vector>

/**
 * @brief: Perform kinetic monte carlo step of protein object when no heads
 * are attached.
 *
 * @param: &pData Refernce for head object
 * @param: Npj Number of rods near head object
 * @param: ep_j Array of rod object pointers near head
 * @param: &uFilterJ Filter of non-unique rods based on periodic
 *         boundary conditions
 *
 * @return: void, Change protein data if head binds
 */
template <class Tubule>
void KMC_U(const ProteinData &pData, const std::vector<const Tubule *> &ep_j,
           double dt, double roll, ProteinBindStatus &pBind) {
    // Assert the heads are not attached
    assert(pData.getBindID(0) == ID_UB && pData.getBindID(1) == ID_UB);

    // Organize tubule objects into spheres and sylinders
    std::vector<const Tubule *> sphPtrArr; //Spheres
    std::vector<const Tubule *> syPtrArr;  // Sylinders
    sphPtrArr.reserve(ep_j.size());        //Set memory for speed
    syPtrArr.reserve(ep_j.size());

    int group_bind = pData.property.group_bind;
    // Split up for loops for performance even if it is ugly
    if (group_bind < 0) { // Don't check tubule group type
        for (const auto &src : ep_j) {
            src->isSphere() ? sphPtrArr.push_back(src)
                            : syPtrArr.push_back(src);
        }
    } else {
        for (const auto &src : ep_j) {
            if (src->group == group_bind) {
                src->isSphere() ? sphPtrArr.push_back(src)
                                : syPtrArr.push_back(src);
            }
        }
    }
    int Nsy = syPtrArr.size();
    int Nsph = sphPtrArr.size();

    double rcut01 = pData.getRcutUS();

    // Create KMC objects
    KMC<Tubule, Tubule> kmc_end0(pData.bind.pos, Nsy, Nsph, rcut01,
                                 pData.getDiffU(), dt);
    KMC<Tubule, Tubule> kmc_end1(pData.bind.pos, Nsy, Nsph, rcut01,
                                 pData.getDiffU(), dt);

    std::vector<double> bindFactors0(Nsy + Nsph, pData.getBindingFactorUS(0));
    std::vector<double> bindFactors1(Nsy + Nsph, pData.getBindingFactorUS(1));

    // Loop over object to bind to and calculate binding probabilities
    kmc_end0.CalcTotProbsUS(syPtrArr, sphPtrArr, bindFactors0);
    kmc_end1.CalcTotProbsUS(syPtrArr, sphPtrArr, bindFactors1);

    // Which protein data is activated to bind
    int activated_end =
        choose_kmc_double(kmc_end0.getTotProb(), kmc_end1.getTotProb(), roll);
    int rj; // Tubule j;
    double bindPos = 0;

    switch (activated_end) {
    case 0:
        rj = kmc_end0.whichObjBindUS(ep_j, bindPos, roll);
        break;
    case 1:
        rj = kmc_end1.whichObjBindUS(ep_j, bindPos, roll);
        break;
    default:
        // No binding occured
        return;
    }

    // Should not give -1 since roll is within the total binding probability
    // If failure, make sure rolls are shifted properly
    assert(rj != -1);

    // Bind head to object
    auto &obj = rj < Nsy ? *(syPtrArr[rj]) : *(sphPtrArr[rj - Nsy]);
    // const double rPos[3] = rod.getPos();
    PS::F64vec3 rVec = obj.getPos();
    double rPos[3] = {rVec[0], rVec[1], rVec[2]};
#ifndef NDEBUG
    printf("U->S Binding\n");
#endif

    pBind.setBind(activated_end, obj.gid, obj.globalIndex, obj.direction, rPos,
                  bindPos, obj.length, obj.rank);
    return;
}

/**
 * @brief: Perform kinetic monte carlo step of protein with 1 head attached.
 *
 * @param: &pData Refernce for protein object
 * @param: Npj Number of rods near protein object
 * @param: ep_j Array of rod objects near protein object
 * @param: &uFilterJ Filter of non-unique rods based on periodic
 *         boundary conditions
 * @param: rollVec
 *
 * @return: void, Change protein data if it binds or unbinds
 */
template <class Tubule>
void KMC_S(const ProteinData &pData, const std::vector<const Tubule *> &ep_j,
           double dt, double KBT, double rollVec[3], ProteinBindStatus &pBind) {
    int Npj = ep_j.size();
    double roll = rollVec[0];
    // Find out which head is bound
    int bound_end = (pData.getBindID(0) == ID_UB) ? 1 : 0;

    const ProteinType *pType = &pData.property;

    std::vector<double> bindFactors(Npj, 0);

    // Create KMC objects for probability and step calculations
    std::vector<const Tubule *> sphPtrArr; //Spheres
    std::vector<const Tubule *> syPtrArr;  // Sylinders
    sphPtrArr.reserve(ep_j.size());        //Set memory for speed
    syPtrArr.reserve(ep_j.size());

    const int group_bind = pData.property.group_bind;
    // Loop over spheres and sylinders, adding to particle arrays
    for (int i = 0; i < Npj; ++i) {
        ep_j[i]->isSphere() ? sphPtrArr.push_back(ep_j[i])
                            : syPtrArr.push_back(ep_j[i]);
        // Specific binding to a group 
        if (ep_j[i]->group == group_bind || group_bind < 0)
            bindFactors[i] =
                pData.getBindingFactorSD(1 - bound_end, ep_j[i]->direction);
    }

    unsigned int Nsy = syPtrArr.size();
    unsigned int Nsph = sphPtrArr.size();

    // Set up KMC objects, then calculate and store probabilities
    KMC<Tubule, Tubule> kmc_unbind(pData.bind.posEndBind[bound_end], dt);
    KMC<Tubule, Tubule> kmc_bind(pData.bind.posEndBind[bound_end], Nsy, Nsph,
                                 dt, pType->LUTablePtr);
    kmc_unbind.CalcProbSU(pData.getUnbindingFactorSU(bound_end));
    kmc_bind.LUCalcTotProbsSD(syPtrArr, sphPtrArr, pData.getBindID(bound_end),
                              bindFactors);
    int activated_end = -1;
    if (bound_end == 1) { // The bound end index is 1
        activated_end = choose_kmc_double(kmc_bind.getTotProb(),
                                          kmc_unbind.getTotProb(), roll);
    } else {
        activated_end = choose_kmc_double(kmc_unbind.getTotProb(),
                                          kmc_bind.getTotProb(), roll);
    }
    // assert(head_activate != -1);
    // Change status of activated head
    if (activated_end == -1) {
        return;
    } else if (bound_end == activated_end) { // Unbind bound head
        rollVec[0] = roll;
        double pos[3] = {};
        for (int i = 0; i < 3; ++i) {
            pos[i] = pData.bind.pos[i];
        }
        kmc_unbind.whereUnbindSU(pData.getRcutUS(), pData.getDiffU(), rollVec,
                                 pos);
        pBind.setUnBind(bound_end);
        pBind.setUnBindPos(pos);
#ifndef NDEBUG
        printf("S->U Unbinding\n");
#endif
        assert(pBind.idBind[0] == ID_UB && pBind.idBind[1] == ID_UB);
    } else if (activated_end == (1 - bound_end)) { // Bind unbound head
        double bindPos; // Position on rod where protein will bind,
                        // passed by reference.
        // Pick rod to bind to
        // Should not give -1 since roll is within the total binding probability
        int rj = kmc_bind.whichObjBindSD(bindPos, roll);
        assert(rj != -1);
        // Bind head to object
        auto &obj = rj < Nsy ? *(syPtrArr[rj]) : *(sphPtrArr[rj - Nsy]);
        // Get position of rod
        PS::F64vec3 rVec = obj.getPos();
        double rPos[3] = {rVec[0], rVec[1], rVec[2]};
#ifndef NDEBUG
        printf("S->D Binding\n");
#endif
        // Bind protein to object
        pBind.setBind(activated_end, obj.gid, obj.globalIndex, obj.direction,
                      rPos, bindPos, obj.length, obj.rank);
        assert(pBind.idBind[0] != ID_UB && pBind.idBind[1] != ID_UB);
    }
    return;
}

/**
 * @brief: Perform kinetic monte carlo step of protein with 2 heads of
 * protein object attached.
 *
 * @param: &pData Refernce for protein object
 * @param: Npj Number of rods near protein object
 * @param: ep_j Array of rod objects near protein object
 * @param: &uFilterJ Filter of non-unique rods based on periodic
 *         boundary conditions
 *
 * @return: void, Change protein data if protein unbinds
 */
template <class Tubule>
void KMC_D(const ProteinData &pData, const std::vector<const Tubule *> &ep_j,
           double dt, double KBT, double roll, ProteinBindStatus &pBind) {
    KMC<Tubule> kmc_UB0(pData.bind.posEndBind[0], dt);
    KMC<Tubule> kmc_UB1(pData.bind.posEndBind[1], dt);
    kmc_UB0.CalcProbDS(pData.getUnbindingFactorDS(0, KBT));
    kmc_UB1.CalcProbDS(pData.getUnbindingFactorDS(1, KBT));
    // double totProb = kmc_UB0.getTotProb() + kmc_UB1.getTotProb();
    int activated_end =
        choose_kmc_double(kmc_UB0.getTotProb(), kmc_UB1.getTotProb(), roll);
    if (activated_end != -1) {
        pBind.setUnBind(activated_end);
#ifndef NDEBUG
        printf("D->S Unbinding\n");
#endif
    }

    return;
}

#endif /* KMC_STEP_FUNCS_HPP */
