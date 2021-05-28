#####################################
#  add target library               #
#####################################

set(SPINPARSER_SOURCE_FILES
    CommandLineOptions.cpp 
    TaskFileParser.cpp 
    FrgCommon.cpp 
    SpinParser.cpp 
    Measurement.cpp 
    LatticeModelFactory.cpp 
    FrgCoreFactory.cpp 
    lib/Log.cpp 
    SU2/SU2FrgCore.cpp 
    SU2/SU2MeasurementCorrelation.cpp 
    XYZ/XYZFrgCore.cpp 
    XYZ/XYZMeasurementCorrelation.cpp 
    TRI/TRIFrgCore.cpp 
    TRI/TRIMeasurementCorrelation.cpp
)
add_library(${CMAKE_PROJECT_NAME}Lib STATIC ${SPINPARSER_SOURCE_FILES} )
target_include_directories(${CMAKE_PROJECT_NAME}Lib PUBLIC ${PROJECT_SOURCE_DIR}/src)

#set compiler flags
if(SPINPARSER_DISABLE_OMP)
    target_compile_definitions(${CMAKE_PROJECT_NAME}Lib PUBLIC DISABLE_OMP)
endif()

if(SPINPARSER_DISABLE_MPI)
    target_compile_definitions(${CMAKE_PROJECT_NAME}Lib PUBLIC DISABLE_MPI)
endif()

if(SPINPARSER_ENABLE_ASSERTIONS)
    target_compile_definitions(${CMAKE_PROJECT_NAME}Lib PUBLIC ENABLE_ASSERTIONS)
endif()

#link libraries
#link boost library
target_link_libraries(${CMAKE_PROJECT_NAME}Lib PUBLIC Boost::regex)
target_link_libraries(${CMAKE_PROJECT_NAME}Lib PUBLIC Boost::thread)
target_link_libraries(${CMAKE_PROJECT_NAME}Lib PUBLIC Boost::system)
target_link_libraries(${CMAKE_PROJECT_NAME}Lib PUBLIC Boost::program_options)
target_link_libraries(${CMAKE_PROJECT_NAME}Lib PUBLIC Boost::filesystem)
target_link_libraries(${CMAKE_PROJECT_NAME}Lib PUBLIC Boost::timer)
target_link_libraries(${CMAKE_PROJECT_NAME}Lib PUBLIC Boost::date_time)
target_link_libraries(${CMAKE_PROJECT_NAME}Lib PUBLIC Boost::unit_test_framework)
#link HDF5 library
target_include_directories(${CMAKE_PROJECT_NAME}Lib PUBLIC ${HDF5_C_INCLUDE_DIRS})
target_link_libraries(${CMAKE_PROJECT_NAME}Lib PUBLIC ${HDF5_C_LIBRARIES})
#link OpenMP library
if(NOT SPINPARSER_DISABLE_OMP)
    target_link_libraries(${CMAKE_PROJECT_NAME}Lib PUBLIC OpenMP::OpenMP_CXX)
endif()
#link MPI library
if(NOT SPINPARSER_DISABLE_MPI)
    target_link_libraries(${CMAKE_PROJECT_NAME}Lib PUBLIC MPI::MPI_CXX)
endif()


#####################################
#  add executable                   #
#####################################

add_executable(${CMAKE_PROJECT_NAME} main.cpp)
target_link_libraries(${CMAKE_PROJECT_NAME} ${CMAKE_PROJECT_NAME}Lib)
install(TARGETS ${CMAKE_PROJECT_NAME} DESTINATION bin)