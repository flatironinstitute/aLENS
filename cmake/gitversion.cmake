#https://github.com/tikonen/blog/tree/master/cmake/git_version
cmake_minimum_required(VERSION 3.10.0)

find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY "${local_dir}"
    OUTPUT_VARIABLE _build_version
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "GIT hash: ${_build_version}")
else()
  message(STATUS "GIT not found")
endif()

string(TIMESTAMP _time_stamp)

configure_file(${local_dir}/cmake/gitversion.txt.in
               ${output_dir}/gitversion.txt @ONLY)
