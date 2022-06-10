#include "TubuleSystem.hpp"

#include "SimToolbox/Util/Logger.hpp"

#include <string>
#include <vector>

#include <mpi.h>

/*! \brief Main function
 *
 * \param argc number of input arguments
 * \param **argv Command line argument array
 * \return int 0=success, 1=error
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    Eigen::setNbThreads(1); // disable threading in eigen

    MPI_Init(&argc, &argv);
    PS::Initialize(argc, argv); // init FDPS system

    Logger::setup_mpi_spdlog();

    // create a system and distribute it to all ranks
    std::string configFile = "RunConfig.yaml";
    std::string configFileProtein = "ProteinConfig.yaml";
    std::string restartFile = "TimeStepInfo.txt";

    // check if restarting
    if (IOHelper::fileExist(restartFile)) {
        TubuleSystem mySystem(configFile, configFileProtein, restartFile, argc,
                              argv);
        while (!mySystem.end()) {
            mySystem.step();
        }
        MPI_Barrier(MPI_COMM_WORLD); // barrier before destructing mySystem;
    } else {
        std::string posFileTubule = "TubuleInitial.dat";
        std::string posFileProtein = "ProteinInitial.dat";
        TubuleSystem mySystem(configFile, posFileTubule, configFileProtein,
                              posFileProtein, argc, argv);
        //Run thermal equilibration loop (no crosslinking kinetics)
        //   Allows for systems to find start in different  initial configurations
        //   before full simulation begins. Does not count towards total number of steps
        mySystem.thermEquil();

        // main simulation loop
        while (!mySystem.end()) {
            mySystem.step();
        }
        MPI_Barrier(MPI_COMM_WORLD); // barrier before destructing mySystem;
    }

    PS::Finalize();
    MPI_Finalize();
    return 0;
}
