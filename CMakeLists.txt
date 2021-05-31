cmake_minimum_required(VERSION 3.16)

#set the project name
project(SpinParser VERSION 1.0 LANGUAGES C CXX)

#define compile options
option(SPINPARSER_BUILD_TESTS "Build tests" ON)
option(SPINPARSER_BUILD_DOCUMENTATION "Build documentation" ON)
option(SPINPARSER_ENABLE_ASSERTIONS "Additional assertions for consistency checks and memory bounds enabled" OFF)
option(SPINPARSER_DISABLE_OMP "Disable OpenMP support" OFF)
option(SPINPARSER_DISABLE_MPI "Disable MPI support" OFF)

#set default build type to release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#locate libraries
#locate boost library
find_package(Boost 1.71 REQUIRED COMPONENTS regex thread system program_options filesystem timer date_time unit_test_framework)
#locate HDF5 library
find_package(HDF5 REQUIRED COMPONENTS C)
#locate OpenMP library
if(NOT SPINPARSER_DISABLE_OMP)
    find_package(OpenMP REQUIRED)
endif()
#locate MPI library
if(NOT SPINPARSER_DISABLE_MPI)
    find_package(MPI REQUIRED)
endif()

#set global compiler flags
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#add subdirectories
add_subdirectory(src)
if(SPINPARSER_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
if(SPINPARSER_BUILD_DOCUMENTATION)
    add_subdirectory(doc)
endif()

#install supplemental files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/res/ DESTINATION res)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/examples/ DESTINATION examples)