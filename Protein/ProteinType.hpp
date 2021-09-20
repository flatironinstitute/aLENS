/**
 * @file ProteinType.hpp
 * @author wenyan4work (wenyan4work@gmail.com)
 * @brief
 * @version 0.1
 * @date 2019-01-04
 *
 * @copyright Copyright (c) 2019
 *
 */

#ifndef PROTEINTYPE_HPP_
#define PROTEINTYPE_HPP_

#include "Util/YamlHelper.hpp"

#include "KMC/kmc.hpp"
#include "KMC/lookup_table.hpp"
#include <cstdio>

/**
 * @brief Specify protein behavior
 *
 */
struct ProteinType {
  public:
    // per protein properties
    int tag = 0;              ///< user-assigned integer tag for different types
    bool walkOff = true;      ///< walf off the end, i.e., no 'end-dewelling'
    double fixedLocation = 0; ///< in [-1,1]
    double freeLength = 0.05; ///< the free length. unit um
    double kappa = 204.7; ///< Spring constant when attached to MT. unit pN/um
    double fstall = 1.0;  ///< stall force. unit pN
    double diffUnbound = 4.5; ///< Unbounded diffusivity, unit um^2/s
    double lambda = 0;        ///< dimensionless unbinding load sensitivity

    double PtoAPratio = 1.; ///< Parallel to anti-parallel binding ratio
    bool fixedEnd0 = false; ///< end0 is fixed to a MT at the given location
    bool useBindVol = false;
    ///< Calculate KMC_s_d factor using effective volume unbound exists in

    // per head properties
    bool vdrag[2] = {false, false}; ///< if including dragged motion
    double vmax[2] = {0, 0};        ///< max velocity for each end.  um/s
    double vmaxAP[2] = {0, 0};
    ///< max velocity for each end, anti parallel. um/s

    double diffBoundS[2] = {0, 0};
    ///< diffusivity for ends when singly bound. unit um^2/s
    double diffBoundD[2] = {0, 0};
    ///< diffusivity for ends when doubly bound. unit um^2/s

    // kMC data field
    double rc;  ///< the capture radius of protein (user set)
    double eps; ///< Number of crosslinker binding sites per um of MT
    double ko_s[2];
    ///< Turnover rate for each end for KMC_u_s and KMC_s_u steps
    double ko_d[2];
    ///< Turnover rate for each end for KMC_s_d and KMC_d_s steps
    double Ka[2];
    ///< Association constant when neither head is bound
    double Ke[2];
    ///< Equilibrium constant when other end is bound

    ///< The type of potential to use when using the binding lookup table
    ///<    0: U = k/2 * (sqrt(d^2 + s^2) - ell0 - D)^2
    ///<    1: U = k/2 * ((d-D)*sqrt(1 +  (s/d)^2) - ell0)^2
    ///< d: perpendicular distance to potentially bound filament
    ///< s: parallel distance along filament from closest point
    int lookupType = 0;
    int lookupGrid = 256;

    // kMC lookup table
    LookupTable *LUTablePtr = nullptr; ///< Pointer to lookup table for binding

    /**
     * @brief read protein type from yaml file
     * 
     * @param p 
     */
    void readFromYaml(const YAML::Node &p) {
        // int
        readConfig(p, VARNAME(tag), tag, "");
        // bool
        readConfig(p, VARNAME(walkOff), walkOff, "");
        readConfig(p, VARNAME(PtoAPratio), PtoAPratio, "");
        readConfig(p, VARNAME(fixedEnd0), fixedEnd0, "");
        useBindVol = false; // optional, default to false
        readConfig(p, VARNAME(useBindVol), useBindVol, "", true);

        // double
        readConfig(p, VARNAME(freeLength), freeLength, "");
        readConfig(p, VARNAME(kappa), kappa, "");
        readConfig(p, VARNAME(fstall), fstall, "");
        readConfig(p, VARNAME(lambda), lambda, "");
        readConfig(p, VARNAME(diffUnbound), diffUnbound, "");
        readConfig(p, VARNAME(rc), rc, "");
        readConfig(p, VARNAME(eps), eps, "");

        // array[2]
        readConfig(p, VARNAME(vmax), vmax, 2, "");
        // vmaxAP is optional
        // if not exist in yaml file, vmaxAP=vmax
        vmaxAP[0] = vmax[0];
        vmaxAP[1] = vmax[1];
        readConfig(p, VARNAME(vmaxAP), vmaxAP, 2, "", true); // optional
        // vdrag is optional
        // if not exist in yaml file, set to false
        vdrag[0] = false;
        vdrag[1] = false;
        readConfig(p, VARNAME(vdrag), vdrag, 2, "", true); // optional
        readConfig(p, VARNAME(diffBoundS), diffBoundS, 2, "");
        readConfig(p, VARNAME(diffBoundD), diffBoundD, 2, "");
        readConfig(p, VARNAME(Ka), Ka, 2, "");
        readConfig(p, VARNAME(Ke), Ke, 2, "");
        readConfig(p, VARNAME(ko_s), ko_s, 2, "");
        readConfig(p, VARNAME(ko_d), ko_d, 2, "");
        // Convert from inverse molarity to inverse number density (um^3)
        Ka[0] /= 602;
        Ka[1] /= 602;
        if (useBindVol) { // Only convert units if useBindVol is not false
            Ke[0] /= 602;
            Ke[1] /= 602;
        }
        lookupType = 0;
        readConfig(p, VARNAME(lookupType), lookupType, "", true); // optional
        lookupGrid = 256;
        readConfig(p, VARNAME(lookupGrid), lookupGrid, "", true); // optional
    }

    /**
     * @brief show protein type
     * 
     */
    void echo() const {
        printf("------Protein Type Properties-----\n");
        printf("tag: %d\n", tag);
        printf("walk off: %d\n", walkOff);
        // printf("bindAntiParallel: %d\n", bindAntiParallel);
        printf("PtoAPratio: %g\n", PtoAPratio);
        printf("fixedEnd0: %d\n", fixedEnd0);
        printf("freeLength: %g um\n", freeLength);
        printf("spring lappa: %g pN/um\n", kappa);
        printf("fstall: %g pN\n", fstall);
        printf("capture radius: %g um\n", rc);
        printf("Unbinding load sensitivity: %g \n", lambda);
        printf("Binding site density along MT: %g um^{-1} \n", eps);
        printf("vmax: %g um/s, %g um/s\n", vmax[0], vmax[1]);
        printf("vmaxAP: %g um/s, %g um/s\n", vmaxAP[0], vmaxAP[1]);
        printf("singly bound diffusivity: %g um^2/s, %g um^2/s\n",
               diffBoundS[0], diffBoundS[1]);
        printf("doubly bound diffusivity: %g um^2/s, %g um^2/s\n",
               diffBoundD[0], diffBoundD[1]);
        printf("singly bound turnover rate: %g s^{-1},%g s^{-1}\n", ko_s[0],
               ko_s[1]);
        printf("doubly bound turnover rate: %g s^{-1},%g s^{-1}\n", ko_d[0],
               ko_d[1]);
        printf("Ka: %g um^3,%g um^3\n", Ka[0], Ka[1]);
        if (useBindVol)
            printf("Ke: %g um^3,%g um^3\n", Ke[0], Ke[1]);
        else
            printf("Ke: %g ,%g \n", Ke[0], Ke[1]);
        printf("useBindVol: %d\n", useBindVol);
        printf("lookupType: %d\n", lookupType);
        if (lookupType)
            printf("WARNING: Current LUT method for lookupType=1 gives large "
                   "errors due to sbound. Do not use this for now\n");
        printf("----------------------------------\n");
    }

    void setMockProteinType() {
        Ka[0] = Ka[1] = 1.0;
        Ke[0] = Ke[1] = 1.0;
        ko_s[0] = ko_s[1] = 1.0;
        ko_d[0] = ko_d[1] = 1.0;
        vmax[0] = vmax[1] = 1.0; // max velocity. unit um/s
        diffBoundS[0] = diffBoundS[1] = 0.0;
        diffBoundD[0] = diffBoundD[1] = 0.0;
        kappa = 1.0; // Spring constant when attached to MT
        eps = 1.0;
        rc = .5;          // capture radius
        freeLength = 1.0; // the 'bind search' radius
        fstall = 1.0;     ///< stall force. unit pN
        tag = 0;          ///< user-assigned integer tag for different types
        walkOff = true;   ///< walf off the end, i.e., no 'end-dewelling'
        lambda = 0;
        PtoAPratio = 1.;
    }

    /**
     * @brief Use diagnostics in KMC to make sure protein parameters are within
     * adequate range to use kinetic Monte-Carlo algorithm.
     *
     * Warning: If binding paradigm changes, this function must also change
     * in order to give good diagnostic information. Might want to change this
     * in the future for better maintainability
     *
     *
     * TODO Right now we over estimate binding because KMC diagnostics does
     * not have the capability to track different head types through
     * binding processes.
     *
     */
    void testKMCStepSize(const double dt, LookupTable *LUT) const {

        // KMC does not need to be templated for diagnostics.
        // "int" is a place holder.
        KMC<int> kmc_diag(rc, diffUnbound, dt, LUT);

        const double prob_thresh = 1e-4;

        double probs[4];
        double ko_s_max = std::max(ko_s[0], ko_s[1]);
        double ko_d_max = std::max(ko_d[0], ko_s[1]);
        double Ka_max = std::max(Ka[0], Ka[1]);
        double Ke_max = std::max(Ke[0], Ke[1]);
        double PPA_ratio = std::max(PtoAPratio, 1.);

        // Constant rate factors for (un)binding. Change if bind model changes.
        double u_s_fact = ko_s_max * Ka_max * eps;
        double s_u_fact = ko_s_max;
        double s_d_fact = PPA_ratio * ko_d_max * Ke_max * eps;
        double d_s_fact = ko_d_max;

        if (lambda == 0) {
            kmc_diag.Diagnostic(u_s_fact, s_u_fact, s_d_fact, d_s_fact, probs);
        } else {
            kmc_diag.DiagnosticUnBindDep(u_s_fact, s_u_fact, s_d_fact, probs);
        }

        // If any probabilities of a double binding event occuring is too high,
        // throw a warning.
        if (probs[0] > prob_thresh || std::isnan(probs[0])) {
            printf(" !!!WARNING: Probability of double event (U->S->U = %f) is "
                   "too high. Try decreasing dt, diffUnbound, or singly "
                   "(un)binding parameters.\n",
                   probs[0]);
        }
        if (probs[1] > prob_thresh || std::isnan(probs[1])) {
            printf(" !!!WARNING: Probability of double event (U->S->D = %f) is "
                   "too high. Try decreasing dt, diffUnbound, or binding "
                   "parameters.\n",
                   probs[1]);
        }
        if (probs[2] > prob_thresh || std::isnan(probs[2])) {
            printf(" !!!WARNING: Probability of double event (S->D->S = %f) is "
                   "too high. Try decreasing dt or doubly (un)binding "
                   "parameters.\n",
                   probs[2]);
        }
        if (probs[3] > prob_thresh || std::isnan(probs[3])) {
            printf(" !!!WARNING: Probability of double event (D->S->U = %f) is "
                   "too high. Try decreasing dt or unbinding parameters.\n",
                   probs[3]);
        }
    }
};

#endif
