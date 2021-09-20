# Force recent cmake and C++
cmake_minimum_required(VERSION 3.4)

add_executable(test_protein test_protein.cpp ProteinConfig.cpp)
# Set the include directory for builds (tests, examples, etc) and for clients
target_include_directories(
  test_protein
  PRIVATE ${CMAKE_SOURCE_DIR}
          ${SIMTOOLBOX}
          ${KMC}
          ${EIGEN3_INCLUDE_DIR}
          ${Trilinos_INCLUDE_DIRS}
          ${TRNG_INCLUDE_DIR}
          ${YAML_CPP_INCLUDE_DIR}
          # $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}> $<INSTALL_INTERFACE:>)
)

target_link_libraries(
  test_protein
  PRIVATE
  PRIVATE ${Trilinos_LIBRARIES}
          ${Trilinos_TPL_LIBRARIES}
          ${TRNG_LIBRARY}
          yaml-cpp
          VTK::IOXML
          Boost::boost
          OpenMP::OpenMP_CXX
          MPI::MPI_CXX)

add_test(NAME protein_test COMMAND test_protein)
