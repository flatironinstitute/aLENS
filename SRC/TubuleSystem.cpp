#include "TubuleSystem.hpp"
#include "CalcProteinBind.hpp"

#include "KMC/kmc.hpp"
#include "Protein/ProteinConfig.hpp"

// SimToolbox module
#include "MPI/CommMPI.hpp"
#include "Trilinos/ZDD.hpp"
#include "Util/GeoUtil.hpp"

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTypeInt32Array.h>
#include <vtkTypeUInt8Array.h>
#include <vtkXMLPPolyDataReader.h>
#include <vtkXMLPolyDataReader.h>

#include <unordered_map>

TubuleSystem::TubuleSystem(const std::string &configFileSystem,
                           const std::string &configFileProtein,
                           const std::string &restartFile, //
                           int argc, char **argv)
    : proteinConfig(configFileProtein) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    if (rank == 0) {
        proteinConfig.echo();
    }

    // step 1 reinitialize tubule system
    // move rods after protein bindings are reconstructed
    // snapID and stepCount both ++ in this reinitialize() function
    rodSystem.reinitialize(configFileSystem, restartFile, argc, argv, false);
    MPI_Barrier(MPI_COMM_WORLD);

    // step 2 initialize shared resource
    rngPoolPtr = rodSystem.getRngPoolPtr();

    // step 3 reinitialize proteins and distribute
    buildLookupTable();
    proteinContainer.initialize();
    // load protein ascii file:
    // snapID has ++ in rodSystem.reinitialize()
    auto snapID = rodSystem.getSnapID() - 1;
    std::string baseFolder = rodSystem.getResultFolderWithID(snapID);
    std::string proteinVTKFile =
        baseFolder + std::string("Protein_") + std::to_string(snapID) + ".pvtp";
    readProteinVTK(proteinVTKFile);

    rodSystem.stepEuler();

    // step 5 setup MixPairInteraction object
    skipBindKinetics_ = proteinConfig.checkSkipBindKinetics();
    if (!skipBindKinetics_)
        bindInteraction.initialize();
    else
        spdlog::warn("All binding kinetics have been turned off.");

    Teuchos::TimeMonitor::zeroOutTimers();

    return;
}

TubuleSystem::TubuleSystem(const std::string &configFileSystem,
                           const std::string &posFileTubule,
                           const std::string &configFileProtein,
                           const std::string &posFileProtein, //
                           int argc, char **argv)
    : proteinConfig(configFileProtein) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    if (rank == 0) {
        proteinConfig.echo();
    }

    // step 1 initialize tubule system
    rodSystem.initialize(configFileSystem, posFileTubule, argc, argv);
    MPI_Barrier(MPI_COMM_WORLD);

    // step 2 initialize shared resource
    rngPoolPtr = rodSystem.getRngPoolPtr();

    // step 3 initialize proteins and distribute
    buildLookupTable();
    proteinContainer.initialize();
    if (IOHelper::fileExist(posFileProtein)) {
        // if posFileProtein file exists, ignore numbers in configFileProtein
        setInitialProteinFromFile(posFileProtein);
    } else {
        setInitialProteinFromConfig();
    }

    // step 4 output initial configuration
    // simBox and Sylinder has been written inside SylinderSystem
    // Here, write protein Info only
    outputProteinData();

    rodSystem.writeResult();

    // step 5 setup MixPairInteraction object
    skipBindKinetics_ = proteinConfig.checkSkipBindKinetics();
    if (!skipBindKinetics_)
        bindInteraction.initialize();
    else
        spdlog::warn("All binding kinetics have been turned off.");

    rodSystem.setTimer(true);
    Teuchos::TimeMonitor::zeroOutTimers();

    return;
}

bool TubuleSystem::end() {
    const auto dt = rodSystem.runConfig.dt;
    const auto timeTotal = rodSystem.runConfig.timeTotal;
    return dt * rodSystem.getStepCount() > timeTotal;
}

void TubuleSystem::prepareStep() {
    // step 1 prepare rodSystem
    // repartitioned if necessary
    rodSystem.prepareStep();
    auto &dinfo = rodSystem.getDomainInfo();

    // protein partition follow tubule partition
    proteinContainer.adjustPositionIntoRootDomain(rodSystem.getDomainInfo());
    proteinContainer.exchangeParticle(rodSystem.getDomainInfoNonConst());

    const auto &runConfig = rodSystem.runConfig;
    if (runConfig.monolayer) {
        const double monoZ =
            (runConfig.simBoxHigh[2] + runConfig.simBoxLow[2]) / 2;
        const int nLocalProtein = proteinContainer.getNumberOfParticleLocal();
#pragma omp parallel for
        for (int i = 0; i < nLocalProtein; i++) {
            auto &pr = proteinContainer[i];
            auto pos = pr.getPosPtr();
            double posNew[3] = {pos[0], pos[1], monoZ};
            pr.setPos(posNew);
        }
    }

    setLookupTablePtr();
}

void TubuleSystem::thermEquil() {
    const double dt = rodSystem.runConfig.dt;
    const double tot_equil_time = rodSystem.runConfig.thermEquilTime;
    spdlog::warn("Thermally equilibrating system for {:8g}", tot_equil_time);
    //Turn off calcBindInteractions when using this loop
    skipBindKinetics_ = true;
    for (double t = 0; t < tot_equil_time; t += dt) {
        thermEquilStep(t);
    }
    skipBindKinetics_ = proteinConfig.checkSkipBindKinetics();
}

void TubuleSystem::thermEquilStep(double t) {
    const double dt = rodSystem.runConfig.dt;
    spdlog::warn("Thermal equilibration time {:8g}", t);

    using Teuchos::Time;
    using Teuchos::TimeMonitor;
    TimeMonitor::zeroOutTimers();

    // Global Timers across all MPI ranks
    // Timers do not accumulate over timesteps
    Teuchos::RCP<Teuchos::Time> therm_loop_timer =
        Teuchos::TimeMonitor::getNewCounter("aLENS thermal equil loop");
    Teuchos::RCP<Time> prepare_step_timer =
        TimeMonitor::getNewCounter("1 prepareStep Time");
    Teuchos::RCP<Time> updateProteinMotionTimer =
        TimeMonitor::getNewCounter("2 updateProteinMotion Time");
    Teuchos::RCP<Time> setProteinConstraintTimer =
        TimeMonitor::getNewCounter("3 setProteinConstraint Time");
    Teuchos::RCP<Time> rodSystemTimer =
        TimeMonitor::getNewCounter("4 rodSystem Time");

    {
        TimeMonitor mon(*therm_loop_timer);
        // step 1 prepare.
        // nothing moves
        {
            TimeMonitor mon(*prepare_step_timer);
            prepareStep();
            // rodSystem.calcOrderParameter();
            spdlog::debug("prepareStep");
        }

        // step 2
        // MTs have moved at the end of the last timestep
        // MT info should be updated and protein move according to this updated MT
        // configuration
        {
            TimeMonitor mon(*updateProteinMotionTimer);
            updateBindWithGid();
            // updateProteinMotion();
            proteinContainer.adjustPositionIntoRootDomain(
                rodSystem.getDomainInfo());
            spdlog::debug("updateProteins");
        }

        // step 3 calculate bilateral constraints with protein binding information
        {
            TimeMonitor mon(*setProteinConstraintTimer);
            setProteinConstraints();
            spdlog::debug("setProteinConstraints");
        }

        // MAJOR STEP:
        // move tubules with binding force and Brownian & collision & bilateral
        // tubule data and protein data written in this step before moving.
        {
            TimeMonitor mon(*rodSystemTimer);
            rodSystem.runStep(false);
            rodSystem.calcConStress();
            spdlog::debug("rodSystemStep");
        }
    }

    rodSystem.printTimingSummary();
}

void TubuleSystem::step() {
    const double dt = rodSystem.runConfig.dt;
    spdlog::warn("CurrentTime {:8g}", rodSystem.getStepCount() * dt);

    using Teuchos::Time;
    using Teuchos::TimeMonitor;
    TimeMonitor::zeroOutTimers();

    // Global Timers across all MPI ranks
    // Timers do not accumulate over timesteps
    Teuchos::RCP<Teuchos::Time> mainLoopTimer =
        Teuchos::TimeMonitor::getNewCounter("aLENS main loop");
    Teuchos::RCP<Time> prepareStepTimer =
        TimeMonitor::getNewCounter("1 prepareStep Time");
    Teuchos::RCP<Time> updateProteinMotionTimer =
        TimeMonitor::getNewCounter("2 updateProteinMotion Time");
    Teuchos::RCP<Time> calcBindInteractionTimer =
        TimeMonitor::getNewCounter("3 calcBindInteraction Time");
    Teuchos::RCP<Time> outputProteinDataTimer =
        TimeMonitor::getNewCounter("4 outputProteinData Time");
    Teuchos::RCP<Time> setProteinConstraintTimer =
        TimeMonitor::getNewCounter("5 setProteinConstraint Time");
    Teuchos::RCP<Time> rodSystemTimer =
        TimeMonitor::getNewCounter("6 rodSystem Time");

    {
        TimeMonitor mon(*mainLoopTimer);
        // step 1 prepare.
        // nothing moves
        {
            TimeMonitor mon(*prepareStepTimer);
            prepareStep();
            rodSystem.calcOrderParameter();
            spdlog::debug("prepareStep");
        }

        // step 2
        // MTs have moved at the end of the last timestep
        // MT info should be updated and protein move according to this updated MT
        // configuration this move includes diffusion and walking
        {
            TimeMonitor mon(*updateProteinMotionTimer);
            updateBindWithGid();
            updateProteinMotion();
            proteinContainer.adjustPositionIntoRootDomain(
                rodSystem.getDomainInfo());
            spdlog::debug("updateProteinMotion");
        }

        // step 3 compute bind interaction.
        // protein ends have moved inside this function
        // this move includes only KMC binding/unbinding kinetics
        if (!skipBindKinetics_) {
            TimeMonitor mon(*calcBindInteractionTimer);
            calcBindInteraction();
            proteinContainer.adjustPositionIntoRootDomain(
                rodSystem.getDomainInfo());
            spdlog::debug("calcBindInteraction");
        }

        // write protein data
        {
            TimeMonitor mon(*outputProteinDataTimer);
            outputProteinData();
            spdlog::debug("outputProteinData");
        }

        // step 4 calculate bilateral constraints with protein binding information
        {
            TimeMonitor mon(*setProteinConstraintTimer);
            setProteinConstraints();
            spdlog::debug("setProteinConstraints");
        }

        // MAJOR STEP:
        // move tubules with binding force and Brownian & collision & bilateral
        // tubule data and protein data written in this step before moving.
        {
            TimeMonitor mon(*rodSystemTimer);
            rodSystem.runStep();
            rodSystem.calcConStress();
            spdlog::debug("rodSystemStep");
        }
    }

    rodSystem.printTimingSummary();
}

void TubuleSystem::calcBindInteraction() {
    auto &dinfo = rodSystem.getDomainInfoNonConst();
    bindInteraction.updateSystem(proteinContainer, rodSystem.getContainer(),
                                 dinfo);
    spdlog::debug("mixSystemUpdated");
    bindInteraction.updateTree();
    spdlog::debug("mixTreeUpdated");
    // bindInteraction.dumpSystem();

    CalcProteinBind interactionFtr(rodSystem.runConfig.dt, proteinConfig.KBT,
                                   rngPoolPtr);
    bindInteraction.computeForce(interactionFtr, dinfo);
    spdlog::debug("forceComputed");

    auto &result = bindInteraction.getForceResult();
    const int nProteinLocal = proteinContainer.getNumberOfParticleLocal();
    assert(result.size() == nProteinLocal);
#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
#pragma omp for
        for (int i = 0; i < nProteinLocal; i++) {
            auto &p = proteinContainer[i];
            p.bind = result[i];
            p.updateGeometryWithBind();
        }
    }
}

void TubuleSystem::updateProteinMotion() {
    const int nProteinLocal = proteinContainer.getNumberOfParticleLocal();
    const double dt = rodSystem.runConfig.dt;
    const double KBT = proteinConfig.KBT;

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
#pragma omp for
        for (int i = 0; i < nProteinLocal; i++) {
            auto &protein = proteinContainer[i];
            if (protein.getWalkOrNot()) { // At least one head is bound
                protein.updatePosWalk(KBT, dt, rngPoolPtr->getU01(threadId),
                                      rngPoolPtr->getN01(threadId),
                                      rngPoolPtr->getN01(threadId));
            } else { // No heads are bound, diffuse in solution
                protein.updatePosDiffuse(dt, rngPoolPtr->getN01(threadId),
                                         rngPoolPtr->getN01(threadId),
                                         rngPoolPtr->getN01(threadId));
            }
        }
    }

    // process protein-boundary motion
    for (const auto &bPtr : rodSystem.runConfig.boundaryPtr) {
#pragma omp parallel
        {
            const int threadId = omp_get_thread_num();
#pragma omp for
            for (int i = 0; i < nProteinLocal; i++) {
                auto &protein = proteinContainer[i];
                if (!protein.getWalkOrNot()) {
                    // No heads are bound, diffuse in solution
                    const auto Query = protein.getPosPtr();
                    double Proj[3] = {0, 0, 0};
                    double delta[3] = {0, 0, 0};
                    bPtr->project(Query, Proj, delta);
                    if (Emap3(delta).dot(ECmap3(Query) - Emap3(Proj)) < 0) {
                        // protein outside of boundary
                        protein.setPos(Proj);
                    }
                }
            }
        }
    }
}

void TubuleSystem::outputProteinData() {
    if (rodSystem.getIfWriteResultCurrentStep()) {
        std::string baseFolder = rodSystem.getCurrentResultFolder();
        IOHelper::makeSubFolder(baseFolder);
        writeProteinAscii();
        writeProteinVTK();
    }
}

void TubuleSystem::writeProteinAscii() {
    // write a single ascii .dat file
    const int nGlobal = proteinContainer.getNumberOfParticleGlobal();
    auto snapID = rodSystem.getSnapID();
    std::string baseFolder = rodSystem.getCurrentResultFolder();
    std::string name = baseFolder + std::string("ProteinAscii_") +
                       std::to_string(snapID) + ".dat";
    ProteinAsciiHeader header;
    header.nparticle = nGlobal;
    header.time = rodSystem.getStepCount() * rodSystem.runConfig.dt;
    proteinContainer.writeParticleAscii(name.c_str(), header);
    MPI_Barrier(MPI_COMM_WORLD);
}

void TubuleSystem::setInitialProteinFromConfig() {
    double boxEdge[3];
    auto initBoxLow = rodSystem.runConfig.initBoxLow;
    auto initBoxHigh = rodSystem.runConfig.initBoxHigh;
    for (int i = 0; i < 3; i++) {
        boxEdge[i] = initBoxHigh[i] - initBoxLow[i];
    }

    // x axis circular cross section
    Evec3 centerCrossSec =
        Evec3(0, (initBoxHigh[1] - initBoxLow[1]) * 0.5 + initBoxLow[1],
              (initBoxHigh[2] - initBoxLow[2]) * 0.5 + initBoxLow[2]);
    double radiusCrossSec = 0.5 * std::min(initBoxHigh[2] - initBoxLow[2],
                                           initBoxHigh[1] - initBoxLow[1]);

    const auto &tubuleContainer = rodSystem.getContainer();
    const int nTubuleLocal = tubuleContainer.getNumberOfParticleLocal();
    const int nTubuleGlobal = tubuleContainer.getNumberOfParticleGlobal();
    const auto &tubuleMap =
        getTMAPFromLocalSize(nTubuleLocal, rodSystem.getCommRcp());

    /**
     * Gid order:
     * 0 Tubule: [0,nTubuleGlobal - 1]
     * 1 Protein Type 1: [0, freeNumberType1 - 1] [0, nFixedPerMT1 *
     * (nTubuleGlobal-1)]
     * 2 Protein Type 2: [0, freeNumberType2 - 1] [0, nFixedPerMT2 *
     * (nTubuleGlobal-1)]
     * ....
     * All gids start from 0 for Tubules and each protein type
     */

    // free protein initialized on rank 0
    // fixedEnd0 protein initialized on all ranks
    const int nType = proteinConfig.types.size();
    if (rank == 0) {
        for (int iType = 0; iType < nType; iType++) {
            // free proteins initialized only on rank 0
            const int freeNumber = proteinConfig.freeNumber[iType];
            if (freeNumber > 0)
                spdlog::debug("initializing free proteins for tag {}",
                              proteinConfig.types[iType].tag);
            for (int i = 0; i < freeNumber; i++) {
                ProteinData newProtein;
                newProtein.gid = i;
                newProtein.property = proteinConfig.types[iType];
                newProtein.property.LUTablePtr = &(LUTArr[iType]);
                newProtein.bind.clear();
                double pos[3] = {0, 0, 0};
                for (int k = 0; k < 3; k++) {
                    pos[k] = rngPoolPtr->getU01(0) * boxEdge[k] + initBoxLow[k];
                }
                if (rodSystem.runConfig.initCircularX) {
                    double y, z = 0;
                    getRandPointInCircle(radiusCrossSec, rngPoolPtr->getU01(0),
                                         rngPoolPtr->getU01(0), y, z);
                    pos[1] = y + centerCrossSec[1];
                    pos[2] = z + centerCrossSec[2];
                }
                newProtein.setPos(pos);
                newProtein.updateGeometryWithBind();
                proteinContainer.addOneParticle(newProtein);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // fixedEnd0 protein initialized on all ranks
    // this includes permanently fixed proteins and initially fixed
    // proteins
    for (int iType = 0; iType < nType; iType++) {
        const int nFixedPerMT = proteinConfig.fixedLocations[iType].size();
        const int gidFixedBase = tubuleMap->getMinGlobalIndex();
        if (nFixedPerMT > 0)
            for (int t = 0; t < nTubuleLocal; t++) {
                for (int k = 0; k < nFixedPerMT; k++) {
                    ProteinData newProtein;
                    newProtein.property = proteinConfig.types[iType];
                    newProtein.bind.clear();
                    newProtein.gid = proteinConfig.freeNumber[iType] +
                                     nFixedPerMT * (gidFixedBase + t) + k;
                    auto &tubule = tubuleContainer[t];
                    newProtein.bind.idBind[0] = tubule.gid;
                    if (proteinConfig.fixedLocations[iType][k] > 1 ||
                        proteinConfig.fixedLocations[iType][k] < -1) {
                        // random location along MT
                        newProtein.bind.distBind[0] =
                            (rngPoolPtr->getU01(0) - 0.5) * tubule.length;
                    } else {
                        newProtein.bind.distBind[0] =
                            proteinConfig.fixedLocations[iType][k] *
                            tubule.length * 0.5;
                    }
                    newProtein.bind.lenBind[0] = tubule.length;
                    newProtein.bind.centerBind[0][0] = tubule.pos[0];
                    newProtein.bind.centerBind[0][1] = tubule.pos[1];
                    newProtein.bind.centerBind[0][2] = tubule.pos[2];
                    Evec3 direction =
                        ECmapq(tubule.orientation) * Evec3(0, 0, 1);
                    newProtein.bind.directionBind[0][0] = direction[0];
                    newProtein.bind.directionBind[0][1] = direction[1];
                    newProtein.bind.directionBind[0][2] = direction[2];
                    newProtein.updateGeometryWithBind();
                    proteinContainer.addOneParticle(newProtein);
                }
            }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    proteinContainer.adjustPositionIntoRootDomain(rodSystem.getDomainInfo());
    proteinContainer.exchangeParticle(rodSystem.getDomainInfoNonConst());
}

void TubuleSystem::updateBindWithGid(bool reconstruct) {
    // both tubule and protein are all distributed on all ranks
    // protein must have valid pos[], adjusted into the root domain
    // this is for reconstruction of bind information from initial
    // configuration newProtein.bind.posEndBind is valid

    // step 1  put data into tubuleDataDirectory
    const auto &tubuleContainer = rodSystem.getContainer();
    const int nTubuleLocal = tubuleContainer.getNumberOfParticleLocal();
    ZDD<SylinderNearEP> tubuleDataDirectory(nTubuleLocal);
    tubuleDataDirectory.gidOnLocal.resize(nTubuleLocal);
    tubuleDataDirectory.dataOnLocal.resize(nTubuleLocal);
#pragma omp parallel for
    for (int t = 0; t < nTubuleLocal; t++) {
        tubuleDataDirectory.gidOnLocal[t] = tubuleContainer[t].gid;
        tubuleDataDirectory.dataOnLocal[t].copyFromFP(tubuleContainer[t]);
    }
    tubuleDataDirectory.buildIndex();

    // step 2 put id to find. two ids per protein
    const int nProteinLocal = proteinContainer.getNumberOfParticleLocal();
    tubuleDataDirectory.gidToFind.resize(2 * nProteinLocal);
    tubuleDataDirectory.dataToFind.resize(2 * nProteinLocal);
#pragma omp parallel for
    for (int p = 0; p < nProteinLocal; p++) {
        // for idBind = ID_UB, ZDD fills findData with invalid data.
        tubuleDataDirectory.gidToFind[2 * p + 0] =
            proteinContainer[p].bind.idBind[0];
        tubuleDataDirectory.gidToFind[2 * p + 1] =
            proteinContainer[p].bind.idBind[1];
    }
    tubuleDataDirectory.find();

    auto simBoxLow = rodSystem.runConfig.simBoxLow;
    auto simBoxHigh = rodSystem.runConfig.simBoxHigh;
    auto simBoxPBC = rodSystem.runConfig.simBoxPBC;

    // step 3 update
#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
#pragma omp for
        for (int p = 0; p < nProteinLocal; p++) {
            auto &protein = proteinContainer[p];

            // printf("%d,%lf,%lf,%lf\n", protein.gid, protein.bind.pos[0],
            // protein.bind.pos[1], protein.bind.pos[2]);

            for (int e = 0; e < 2; e++) { // check both ends
                if (protein.bind.idBind[e] == ID_UB) {
                    protein.bind.setUnBind(e);
                } else { // get tubule data from ZDD
                    auto &tubuleBind =
                        tubuleDataDirectory.dataToFind[2 * p + e];
                    if (protein.bind.idBind[e] != tubuleBind.gid) {
                        printf("gid does not match\n");
                        std::exit(1);
                    }
                    // step 1 find actual bind PBC Image of tubule
                    for (int dim = 0; dim < 3; dim++) {
                        if (simBoxPBC[dim]) {
                            findPBCImage(simBoxLow[dim], simBoxHigh[dim],
                                         tubuleBind.pos[dim],
                                         protein.bind.pos[dim]);
                        }
                    }

                    // step 2 reconstruct binding information here
                    /***
                     * For each end:
                     * 1. update center/direction
                     * 2. if distBind not valid (setup initial data):
                     *      compute distBind with posEnd.
                     *      if this distBind is invalid (data error),
                     *      then generage random valid distBind
                     *    if distBind is valid (during simulation):
                     *      do nothing
                     *
                     * After both ends updated:
                     *    protein.calcPosEndWithDistBind()
                     *
                     * Result: updated data always consistent with
                     *  valid given posEnd or distBind
                     *
                     */
                    protein.bind.lenBind[e] = tubuleBind.length;
                    for (int dim = 0; dim < 3; dim++) {
                        protein.bind.directionBind[e][dim] =
                            tubuleBind.direction[dim];
                        protein.bind.centerBind[e][dim] = tubuleBind.pos[dim];
                    }
                    if (reconstruct || skipBindKinetics_) {
                        // update distBind rebuild with posEndBind
                        Evec3 bindFoot = Emap3(protein.bind.posEndBind[e]);
                        double distBind = (bindFoot - Emap3(tubuleBind.pos))
                                              .dot(Emap3(tubuleBind.direction));
                        if (distBind > -tubuleBind.length * 0.5 &&
                            distBind < tubuleBind.length * 0.5) {
                            // Case 1, valid data
                            protein.bind.distBind[e] = distBind;
                        } else if ((distBind >
                                        -tubuleBind.length * 0.5 * 1.01 &&
                                    distBind <
                                        tubuleBind.length * 0.5 * 1.01) ||
                                   tubuleBind.isSphere()) {
                            // Case 2, valid data with some floating point error
                            protein.bind.distBind[e] = distBind;
                            protein.bind.updatePosEndClamp(e);
                        } else {
                            // Case 3, invalid data, set unbind
                            spdlog::error(
                                "posEnd {} invalid at distBind {} for "
                                "protein {}",
                                e, distBind, protein.gid);
                            spdlog::error("set end {} to unbind", e);
                            protein.bind.setUnBind(e);
                        }
                    }
                }
            }
            if (reconstruct || skipBindKinetics_)
                protein.updateGeometryWithBind();
        }
    }
    proteinContainer.adjustPositionIntoRootDomain(rodSystem.getDomainInfo());

    if (reconstruct) {
        //buildLookupTable();
        setLookupTablePtr();
    }

#pragma omp parallel for
    for (int i = 0; i < nProteinLocal; i++) {
        // update protein force/torque of binding by actual updated binding position.
        auto &p = proteinContainer[i];
        p.updateForceTorqueBind();
    }
}

void TubuleSystem::buildLookupTable() {
    LUTArr.resize(proteinConfig.types.size());
    double D = rodSystem.runConfig.sylinderDiameter;
    for (int i = 0; i < proteinConfig.types.size(); i++) {
        LUTFiller *lut_filler_ptr = makeLUTFiller(proteinConfig.types[i]);
        LUTArr[i] =
            LookupTable(lut_filler_ptr, proteinConfig.types[i].useBindVol);
        if (proteinConfig.types[i].lookupType == 0) {
            proteinConfig.types[i].testKMCStepSize(rodSystem.runConfig.dt,
                                                   &LUTArr[i]);
        }
        if (proteinConfig.types[i].useBindVol)
            spdlog::critical("bind volume: {}", LUTArr[i].getBindVolume());

        // New pointer was created in makeLUTFiller method. Clean it up.
        delete lut_filler_ptr;
    }
}

LUTFiller *TubuleSystem::makeLUTFiller(const ProteinType &ptype) {
    double D = rodSystem.runConfig.sylinderDiameter;
    const int grid_num = ptype.lookupGrid;
    switch (ptype.lookupType) {
    case 0: {
        // Energy dependent lookup table filler object
        LUTFillerEdep *lut_filler_ptr = new LUTFillerEdep(grid_num, grid_num);
        // Exponent pre-factor in Boltzmann factor of lookup table
        double M = .5 * (1. - ptype.lambda) * ptype.kappa / proteinConfig.KBT;
        // Add tubule diameter to freeLength to approximate binding to the
        // surface of sylinder instead of center.
        double ell0 = ptype.freeLength + D;
        lut_filler_ptr->Init(M, ell0, D);
        return lut_filler_ptr;
    }
    case 1: {
        LUTFiller2ndOrder *lut_filler_ptr =
            new LUTFiller2ndOrder(grid_num, grid_num);
        // Exponent pre-factor in Boltzmann factor of lookup table
        double M = .5 * (1. - ptype.lambda) * ptype.kappa / proteinConfig.KBT;
        lut_filler_ptr->Init(M, ptype.freeLength, D);
        return lut_filler_ptr;
    }
    default:
        return nullptr;
        break;
    }
}

void TubuleSystem::setLookupTablePtr() {
    // set LUT ptrs
    const int numProtein = proteinContainer.getNumberOfParticleLocal();
    const auto &tagLookup = proteinConfig.tagLookUp;
#pragma omp parallel for
    for (int i = 0; i < numProtein; i++) {
        auto &protein = proteinContainer[i];
        const int tag = protein.property.tag;
        const auto &index = tagLookup.find(tag);
        if (index == tagLookup.end()) {
            spdlog::critical("protein tagLookup error");
            std::exit(1);
        }
        protein.property.LUTablePtr = &(LUTArr[index->second]);
    }
}

void TubuleSystem::findTubuleRankWithGid() {
    // using ZDD to find distributed data
    // step 1  put data into tubuleDataDirectory
    const auto &tubuleContainer = rodSystem.getContainer();
    const int nTubuleLocal = tubuleContainer.getNumberOfParticleLocal();
    ZDD<int> tubuleDataDirectory(nTubuleLocal);
    tubuleDataDirectory.gidOnLocal.resize(nTubuleLocal);
    tubuleDataDirectory.dataOnLocal.resize(nTubuleLocal);
#pragma omp parallel for
    for (int t = 0; t < nTubuleLocal; t++) {
        tubuleDataDirectory.gidOnLocal[t] = tubuleContainer[t].gid;
        tubuleDataDirectory.dataOnLocal[t] = rank;
    }
    tubuleDataDirectory.buildIndex();

    // step 2 put id to find. two ids per protein
    const int nProteinLocal = proteinContainer.getNumberOfParticleLocal();
    tubuleDataDirectory.gidToFind.resize(2 * nProteinLocal);
    tubuleDataDirectory.dataToFind.resize(2 * nProteinLocal);
#pragma omp parallel for
    for (int p = 0; p < nProteinLocal; p++) {
        // for idBind = ID_UB, ZDD fills findData with invalid data.
        tubuleDataDirectory.gidToFind[2 * p + 0] =
            proteinContainer[p].bind.idBind[0];
        tubuleDataDirectory.gidToFind[2 * p + 1] =
            proteinContainer[p].bind.idBind[1];
    }
    tubuleDataDirectory.find();

    // step 3 update
#pragma omp parallel for
    for (int p = 0; p < nProteinLocal; p++) {
        auto &protein = proteinContainer[p];
        for (int e = 0; e < 2; e++) {
            if (protein.bind.idBind[e] != ID_UB) {
                protein.bind.rankBind[e] =
                    tubuleDataDirectory.dataToFind[2 * p + e];
                assert(protein.bind.rankBind[e] >= 0 &&
                       protein.bind.rankBind[e] < nProcs);
            } else {
                protein.bind.rankBind[e] = -1;
            }
        }
    }
}

void TubuleSystem::setProteinConstraints() {
    const int nLocal = proteinContainer.getNumberOfParticleLocal();

    auto &conPool = rodSystem.getConstraintPoolNonConst();
    const int nThreads = conPool.size();
    const double tubuleDiameter = rodSystem.runConfig.sylinderDiameter;
#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        auto &biQue = conPool[tid];
#pragma omp for
        for (int i = 0; i < nLocal; i++) {
            const auto &pr = proteinContainer[i];
            if (pr.getBindID(0) == ID_UB || pr.getBindID(1) == ID_UB) {
                // not doubly bound, not a constraint
                continue;
            }
            // geometry of MT I and J
            const Evec3 centerI = ECmap3(pr.bind.centerBind[0]);
            const Evec3 directionI = ECmap3(pr.bind.directionBind[0]);
            const Evec3 Ploc = ECmap3(pr.bind.posEndBind[0]);
            const Evec3 centerJ = ECmap3(pr.bind.centerBind[1]);
            const Evec3 directionJ = ECmap3(pr.bind.directionBind[1]);
            const Evec3 Qloc = ECmap3(pr.bind.posEndBind[1]);
            // information of constraint block
            const Evec3 PQvec = Qloc - Ploc;
            const double delta0 =
                pr.getProteinForceLength() - pr.property.freeLength;
            const double gamma = -delta0 * pr.property.kappa;
            const Evec3 normI = (Ploc - Qloc).normalized();
            const Evec3 normJ = -normI;
            const Evec3 posI = Ploc - centerI;
            const Evec3 posJ = Qloc - centerJ;
            biQue.emplace_back(
                delta0, gamma, // current separation, initial guess of gamma
                pr.bind.idBind[0], pr.bind.idBind[1], //
                pr.bind.indexBind[0],                 //
                pr.bind.indexBind[1],                 //
                normI.data(), normJ.data(), // direction of constraint force
                posI.data(),
                posJ.data(), // location relative to particle center
                Ploc.data(), Qloc.data(), // location in lab frame
                false, true, pr.property.kappa, 0.0);
            Emat3 stressIJ;
            CalcSylinderNearForce::collideStress(
                directionI, directionJ, centerI, centerJ, //
                pr.bind.lenBind[0], pr.bind.lenBind[1], tubuleDiameter / 2,
                tubuleDiameter / 2, 1.0, Ploc, Qloc, stressIJ);
            biQue.back().setStress(stressIJ);
        }
    }
}

void TubuleSystem::writeProteinVTK() {
    // write parallel XML VTK files from all ranks
    std::string baseFolder = rodSystem.getCurrentResultFolder();
    auto snapID = rodSystem.getSnapID();
    ProteinData::writeVTP<PS::ParticleSystem<ProteinData>>(
        proteinContainer, proteinContainer.getNumberOfParticleLocal(),
        baseFolder, std::to_string(snapID), rank);

    if (rank == 0) { // write parallel head
        ProteinData::writePVTP(baseFolder, std::to_string(snapID), nProcs);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void TubuleSystem::setInitialProteinFromFile(const std::string &posFilename) {
    // read file all to rank 0
    std::vector<ProteinData> proteinReadFromFile;
    if (rank == 0) {
        proteinConfig.echo();
        spdlog::warn("Read protein position from data file");
        std::ifstream myfile(posFilename);
        std::string line;
        std::getline(myfile, line); // read two header lines
        std::getline(myfile, line);

        while (std::getline(myfile, line)) {
            char typeChar;
            std::istringstream liness(line);
            liness >> typeChar;
            if (typeChar == 'P') {
                int gid, tag, idBind[2];
                double end0[3];
                double end1[3];
                liness >> gid >> tag >> end0[0] >> end0[1] >> end0[2] >>
                    end1[0] >> end1[1] >> end1[2] >> idBind[0] >> idBind[1];
                int iType = proteinConfig.tagLookUp.find(tag)->second;
                ProteinData newProtein;
                newProtein.setFromFileInput(gid, tag, end0, end1, idBind,
                                            proteinConfig.types[iType]);
                proteinReadFromFile.push_back(newProtein);
                typeChar = 'N';
            }
        }
        myfile.close();
    } else { // other rank no protein in the beginning
    }

    MPI_Barrier(MPI_COMM_WORLD);
    const int nProteinInFile = proteinReadFromFile.size();
    proteinContainer.setNumberOfParticleLocal(nProteinInFile);
#pragma omp parallel for
    for (int i = 0; i < nProteinInFile; i++) {
        proteinContainer[i] = proteinReadFromFile[i];
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // protein.setPos() cannot be called before the bound MT data is updated
    // this must be called before adjustPosition and exchange
    updateBindWithGid(true);

    proteinContainer.adjustPositionIntoRootDomain(rodSystem.getDomainInfo());
    proteinContainer.exchangeParticle(rodSystem.getDomainInfoNonConst());
}

void TubuleSystem::readProteinVTK(const std::string &pvtpFileName) {
    auto &commRcp = rodSystem.getCommRcp();
    if (commRcp->getRank() != 0) {
        proteinContainer.setNumberOfParticleLocal(0);
    } else {
        vtkSmartPointer<vtkXMLPPolyDataReader> reader =
            vtkSmartPointer<vtkXMLPPolyDataReader>::New();
        spdlog::warn("Reading " + pvtpFileName);
        reader->SetFileName(pvtpFileName.c_str());
        reader->Update();

        // Extract the polydata (At this point, the polydata is unsorted)
        vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();

        // Extract the point/cell data
        vtkSmartPointer<vtkPoints> posData = polydata->GetPoints();
        // cell data
        vtkSmartPointer<vtkTypeInt32Array> gidData =
            vtkArrayDownCast<vtkTypeInt32Array>(
                polydata->GetCellData()->GetAbstractArray("gid"));
        vtkSmartPointer<vtkTypeInt32Array> tagData =
            vtkArrayDownCast<vtkTypeInt32Array>(
                polydata->GetCellData()->GetAbstractArray("tag"));
        // point data
        vtkSmartPointer<vtkTypeInt32Array> idBindData =
            vtkArrayDownCast<vtkTypeInt32Array>(
                polydata->GetPointData()->GetArray("idBind"));

        // two points per protein
        const int proteinNumberInFile = posData->GetNumberOfPoints() / 2;
        std::vector<ProteinData> proteinReadFromFile(proteinNumberInFile);
        // set local
        proteinContainer.setNumberOfParticleLocal(proteinNumberInFile);
        spdlog::debug("Protein number in file: {}", proteinNumberInFile);
#pragma omp parallel for
        for (int i = 0; i < proteinNumberInFile; i++) {
            double end0[3] = {0, 0, 0};
            double end1[3] = {0, 0, 0};
            posData->GetPoint(i * 2, end0);
            posData->GetPoint(i * 2 + 1, end1);
            int gid = gidData->GetTypedComponent(i, 0);
            int tag = tagData->GetTypedComponent(i, 0);
            int idBind[2];
            idBind[0] = idBindData->GetTypedComponent(2 * i, 0);
            idBind[1] = idBindData->GetTypedComponent(2 * i + 1, 0);
            auto &newProtein = proteinContainer[i];
            int iType = proteinConfig.tagLookUp.find(tag)->second;
            newProtein.setFromFileInput(gid, tag, end0, end1, idBind,
                                        proteinConfig.types[iType]);
            // std::cout << gid << " " << tag << " " << idBind[0] << " "
            //           << idBind[1] << " " << Emap3(end0).transpose() << " "
            //           << Emap3(end1).transpose() << std::endl;
        }
    }
    commRcp->barrier();
    // protein.setPos() cannot be called before the bound MT data is updated
    // this must be called before adjustPosition and exchange
    updateBindWithGid(true);

    proteinContainer.adjustPositionIntoRootDomain(rodSystem.getDomainInfo());
    proteinContainer.exchangeParticle(rodSystem.getDomainInfoNonConst());
}
