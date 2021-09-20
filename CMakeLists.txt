cmake_minimum_required(VERSION 3.10) # same requirement as trilinos
project(aLENS)

# global settings
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

# compiler flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qno-offload")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qno-offload")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

# required variables
set(SFTPATH
    $ENV{HOME}/local
    CACHE FILEPATH "Path to external dependency libraries")
message(STATUS "SFTPATH: ${SFTPATH}")

option(ENABLE_TEST "enable tests")

# external libraries, static link
find_package(TRNG REQUIRED)
find_package(yaml-cpp REQUIRED PATHS ${SFTPATH} NO_DEFAULT_PATH)
find_package(VTK REQUIRED PATHS ${SFTPATH} NO_DEFAULT_PATH)
# external libraries, header only
find_package(Eigen3 REQUIRED HINTS ${SFTPATH}/share/eigen3/cmake)
find_package(Boost REQUIRED)

# Get Trilinos as one entity
find_package(Trilinos REQUIRED HINTS ${SFTPATH}/lib/cmake)
# clean up the too-long list
list(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
list(REMOVE_DUPLICATES Trilinos_TPL_INCLUDE_DIRS)
list(REVERSE Trilinos_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_LIBRARIES)
list(REVERSE Trilinos_LIBRARIES)
list(REVERSE Trilinos_TPL_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
list(REVERSE Trilinos_TPL_LIBRARIES)
# Echo trilinos build info just for fun
message("\nFound Trilinos!  Here are the details: ")
message("   Trilinos_DIR = ${Trilinos_DIR}")
message("   Trilinos_VERSION = ${Trilinos_VERSION}")
message("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
message("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
message("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
message("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
message("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
message("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
message("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
message("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
message("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
message("End of Trilinos details\n")

# internal submodules
set(KMC ${CMAKE_SOURCE_DIR}/KMC)
set(PROTEIN ${CMAKE_SOURCE_DIR}/Protein)
set(SIMTOOLBOX ${CMAKE_SOURCE_DIR}/SimToolbox)
set(SRC ${CMAKE_SOURCE_DIR}/SRC)

#[[
#################################################
    main executable aLENS.X
#################################################
]]
add_executable(
  aLENS.X
  ${SRC}/main.cpp
  ${SRC}/TubuleSystem.cpp
  ${PROTEIN}/ProteinConfig.cpp
  ${SIMTOOLBOX}/Boundary/Boundary.cpp
  ${SIMTOOLBOX}/Constraint/BCQPSolver.cpp
  ${SIMTOOLBOX}/Constraint/ConstraintCollector.cpp
  ${SIMTOOLBOX}/Constraint/ConstraintOperator.cpp
  ${SIMTOOLBOX}/Constraint/ConstraintSolver.cpp
  ${SIMTOOLBOX}/Sylinder/Sylinder.cpp
  ${SIMTOOLBOX}/Sylinder/SylinderConfig.cpp
  ${SIMTOOLBOX}/Sylinder/SylinderSystem.cpp
  ${SIMTOOLBOX}/Trilinos/TpetraUtil.cpp
  ${SIMTOOLBOX}/Util/Base64.cpp)

target_compile_definitions(
  aLENS.X
  PRIVATE # PARTICLE_SIMULATOR_DEBUG_PRINT # for debug
          PARTICLE_SIMULATOR_MPI_PARALLEL #
          PARTICLE_SIMULATOR_THREAD_PARALLEL)

target_include_directories(
  aLENS.X PRIVATE ${CMAKE_SOURCE_DIR} ${KMC} ${SIMTOOLBOX} ${TRNG_INCLUDE_DIR}
                  ${EIGEN3_INCLUDE_DIR} ${Trilinos_INCLUDE_DIRS})

target_link_libraries(
  aLENS.X
  PRIVATE ${TRNG_LIBRARY}
          ${Trilinos_LIBRARIES}
          ${Trilinos_TPL_LIBRARIES}
          yaml-cpp
          VTK::IOXML
          Boost::boost
          OpenMP::OpenMP_CXX
          MPI::MPI_CXX)

# add git information to installation directory
add_custom_command(
  TARGET aLENS.X
  PRE_BUILD
  COMMAND
    ${CMAKE_COMMAND} -Dlocal_dir="${CMAKE_SOURCE_DIR}"
    -Doutput_dir="${CMAKE_BINARY_DIR}" -P
    "${CMAKE_SOURCE_DIR}/cmake/gitversion.cmake")

install(DIRECTORY ${CMAKE_SOURCE_DIR}/Run/result
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY ${CMAKE_SOURCE_DIR}/Run/scripts
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ${CMAKE_BINARY_DIR}/gitversion.txt
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(TARGETS aLENS.X DESTINATION ${CMAKE_INSTALL_PREFIX})

if(ENABLE_TEST)
  enable_testing()
  add_subdirectory(Protein)
  add_subdirectory(KMC EXCLUDE_FROM_ALL)
endif()
