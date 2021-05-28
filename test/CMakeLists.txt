############################################
#  add tests                               #
############################################

#add unit tests
set(SPINPARSER_UNIT_TEST_FILES
	test_CutoffDiscretization.cpp
	test_FrequencyDiscretization.cpp
	test_Geometry.cpp
	test_InputParser.cpp
	test_Integrator.cpp
	test_Lattice.cpp
	test_SU2VertexSingleParticle.cpp
	test_SU2VertexTwoParticle.cpp
	test_TRIVertexSingleParticle.cpp
	test_TRIVertexTwoParticle.cpp
	test_XYZVertexSingleParticle.cpp
	test_ValueBundle.cpp
	test_XYZVertexTwoParticle.cpp
)

if(SPINPARSER_DISABLE_MPI)
	list(APPEND SPINPARSER_UNIT_TEST_FILES test_LoadManager.cpp )
else()
	set(SPINPARSER_UNIT_TEST_FILES_MPI test_LoadManager.cpp)
endif()

foreach(TEST_SOURCE IN LISTS SPINPARSER_UNIT_TEST_FILES)
	string(REGEX REPLACE "(^test_)|(\\.[ch]pp)" "" TEST_BASE_NAME ${TEST_SOURCE})
	add_executable(${TEST_BASE_NAME}Test ${TEST_SOURCE})
	target_link_libraries(${TEST_BASE_NAME}Test ${CMAKE_PROJECT_NAME}Lib)
	add_test(NAME ${TEST_BASE_NAME}Test COMMAND ${TEST_BASE_NAME}Test)
endforeach()

foreach(TEST_SOURCE IN LISTS SPINPARSER_UNIT_TEST_FILES_MPI)
	string(REGEX REPLACE "(^test_)|(\\.[ch]pp)" "" TEST_BASE_NAME ${TEST_SOURCE})
	add_executable(${TEST_BASE_NAME}Test ${TEST_SOURCE})
	target_link_libraries(${TEST_BASE_NAME}Test ${CMAKE_PROJECT_NAME}Lib)
	set(TEST_PARAMETERS ${MPIEXEC_NUMPROC_FLAG} 4 ${TEST_BASE_NAME}Test)
	add_test(NAME ${TEST_BASE_NAME}Test COMMAND ${MPIEXEC_EXECUTABLE} ${TEST_PARAMETERS})
endforeach()

#add scripted tests
set(SPINPARSER_SCRIPTED_TEST_FILES
	test_reference1.sh
	test_reference2.sh
	test_reference3.sh
	test_checkpoint.sh
	test_defer.sh
	test_pythonObs.sh
)
if(NOT SPINPARSER_DISABLE_MPI)
	list(APPEND SPINPARSER_SCRIPTED_TEST_FILES test_MPI.sh)
endif()

set(SPINPARSER_SCRIPTED_TEST_FAILURE_FILES
	test_referenceFail.sh
)

set(TEST_ENVIRONMENT
	TEST_ROOT_DIR=${PROJECT_SOURCE_DIR}
	TEST_WORK_DIR=${CMAKE_BINARY_DIR}
	TEST_SCRIPT_DIR=${PROJECT_SOURCE_DIR}/test/scripted
	TEST_EXECUTABLE=${CMAKE_BINARY_DIR}/src/${CMAKE_PROJECT_NAME}\ -r\ ${PROJECT_SOURCE_DIR}/res
	TEST_MPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
	TEST_MPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
)

foreach(TEST_SCRIPT IN LISTS SPINPARSER_SCRIPTED_TEST_FAILURE_FILES)
	string(REGEX REPLACE "(^test_)|(\\.sh)" "" TEST_BASE_NAME ${TEST_SCRIPT})
	add_test(NAME ${TEST_BASE_NAME}ScriptedTest COMMAND ${PROJECT_SOURCE_DIR}/test/scripted/${TEST_SCRIPT})
	set_property(TEST ${TEST_BASE_NAME}ScriptedTest PROPERTY ENVIRONMENT ${TEST_ENVIRONMENT})
	set_property(TEST ${TEST_BASE_NAME}ScriptedTest PROPERTY WILL_FAIL TRUE)
endforeach()

foreach(TEST_SCRIPT IN LISTS SPINPARSER_SCRIPTED_TEST_FILES)
	string(REGEX REPLACE "(^test_)|(\\.sh)" "" TEST_BASE_NAME ${TEST_SCRIPT})
	add_test(NAME ${TEST_BASE_NAME}ScriptedTest COMMAND ${PROJECT_SOURCE_DIR}/test/scripted/${TEST_SCRIPT})
	set_property(TEST ${TEST_BASE_NAME}ScriptedTest PROPERTY ENVIRONMENT ${TEST_ENVIRONMENT})
endforeach()