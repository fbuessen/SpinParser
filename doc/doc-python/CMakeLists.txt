find_package(PythonInterp REQUIRED)

add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/html"
    COMMAND "${PYTHON_EXECUTABLE}" "-m" "pdoc" "--html" "-c" "latex_math=True" "-o" "${CMAKE_CURRENT_BINARY_DIR}/html" "${PROJECT_SOURCE_DIR}/opt/python/spinparser"
)

add_custom_target(SpinParserDocumentationPython ALL
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/html"
)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/spinparser/ DESTINATION doc/doc-python)