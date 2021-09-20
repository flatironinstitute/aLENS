/**
 * @file TubuleSystem.hpp
 * @author wenyan4work (wenyan4work@gmail.com)
 * @brief Mixed system of microtubule and protein
 * @version 0.1
 * @date 2018-12-14
 *
 * @copyright Copyright (c) 2018
 *
 */
#ifndef TUBULESYSTEM_HPP_
#define TUBULESYSTEM_HPP_

// Protein
#include "Protein/ProteinConfig.hpp"
#include "Protein/ProteinData.hpp"

#include "TubuleBind.hpp"

// Sylinder as Tubule
#include "Constraint/ConstraintSolver.hpp"
#include "FDPS/particle_simulator.hpp"
#include "MPI/MixPairInteraction.hpp"
#include "Sylinder/SylinderSystem.hpp"
#include "Trilinos/TpetraUtil.hpp"
#include "Util/TRngPool.hpp"

#include <memory>

#include <mpi.h>
#include <omp.h>

/**
 * @brief Mixed system of microtubule and protein
 *
 */
class TubuleSystem {

  public:
    const ProteinConfig proteinConfig; ///< the ProteinConfig file

    /**
     * @brief Construct a new TubuleSystem object
     *
     * @param configFileSystem
     * @param posFileTubule
     * @param configFileProtein
     * @param posFileProtein
     * @param argc
     * @param argv
     */
    TubuleSystem(const std::string &configFileSystem,
                 const std::string &posFileTubule,
                 const std::string &configFileProtein,
                 const std::string &posFileProtein, int argc, char **argv);

    /**
     * @brief Construct a new Tubule System:: Tubule System object
     * 
     * @param configFileSystem 
     * @param configFileProtein 
     * @param restartFile saved restarting file
     * @param argc 
     * @param argv 
     */
    TubuleSystem(const std::string &configFileSystem,
                 const std::string &configFileProtein,
                 const std::string &restartFile, //
                 int argc, char **argv);

    /**
     * @brief Destroy the TubuleSystem object
     *
     */
    ~TubuleSystem() = default;

    // forbid copy
    TubuleSystem(const TubuleSystem &) = delete;
    TubuleSystem &operator=(const TubuleSystem &) = delete;

    /**
     * @brief one time step forward
     *
     */
    void step();

    /**
     * @brief mark the end of time-stepping
     * 
     * @return true 
     * @return false 
     */
    bool end();

  private:
    int rank;                             ///< mpi rank
    int nProcs;                           ///< mpi size
    std::shared_ptr<TRngPool> rngPoolPtr; ///< point to rodSystem.rngPoolPtr
    std::vector<LookupTable> LUTArr;      ///< protein LookupTable Holder

    SylinderSystem rodSystem; ///< each tubule modeled as a sylinder

    // same DomainInfo for protein and tubule
    PS::ParticleSystem<ProteinData> proteinContainer; ///< all proteins

    MixPairInteraction<ProteinData, Sylinder,     // FPT, FPS
                       ProteinData, TubuleBindEP, // EPT, EPS
                       ProteinBindStatus>         // 'Force'
        bindInteraction;                          ///< mixed interaction

    std::vector<double> proteinForceTorqueOnTubule; ///< bind force and torque

    /**
     * @brief prepare necessary data structures for a step
     *
     */
    void prepareStep();

    /**
     * @brief Set the configuration for sylinder system
     *
     * @param sylinderConfig
     */
    void setSylinderConfig(SylinderConfig &sylinderConfig);

    /**
     * @brief calculate and update protein.bind with bindInteraction tree
     *
     */
    void calcBindInteraction();

    /**
     * @brief calculate protein diffusion or walking
     *
     */
    void updateProteinMotion();

    /**
     * @brief set protein as bilateral constraints in sylinderSystem
     *
     */
    void setProteinConstraints();

    /**
     * @brief update protein.bind with tubule gid
     * used in two cases
     * (1) tubule moved. update protein.bind.centerBind and directionBind only
     * (2) new data. reconstruct full protein.bind
     *
     * if distBind[e] in [-lenMT/2,lenMT/2], reconstruct
     * otherwise, input file mismatch and recalc with pos
     */
    void updateBindWithGid(bool reconstruct = false);

    /**
     * @brief find home rank of each tubule by gid
     *
     */
    void findTubuleRankWithGid();

    /**
     * @brief write snapshot
     *
     */
    void outputProteinData();

    /**
     * @brief Set initial protein data from config file
     *
     */
    void setInitialProteinFromConfig();

    /**
     * @brief Set initial protein data from config and dat files
     *
     * @param posFilename
     */
    void setInitialProteinFromFile(const std::string &posFilename);

    /**
     * @brief Build all lookup tables in LUTArr
     *
     */
    void buildLookupTable();

    /**
     * @brief Make filler object pointer to initialize lookup tables. 
     * Pointer because there are several different types of fillers available
     * inhereted from the abstract class LUTFiller.
     *
     */
    LUTFiller *makeLUTFiller(const ProteinType &ptype);

    /**
     * @brief set the LUT ptr in every protein
     * 
     */
    void setLookupTablePtr();

    /**
     * @brief write protein dat file
     *
     */
    void writeProteinAscii();

    /**
     * @brief write protein XML VTK file
     *
     */
    void writeProteinVTK();

    /**
     * @brief read protein XML VTK file
     * 
     */
    void readProteinVTK(const std::string &pvtpFileName);
};

#endif
