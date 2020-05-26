add_custom_target(doxygen_HelmholtzTransmissionProblem
COMMAND ${DOXYGEN_EXECUTABLE} O_DIRECTORY=${CMAKE_BINARY_DIR}/doxygen_HelmholtzTransmissionProblem doxygen Doxyfile
COMMENT "Generating HelmholtzTransmissionProblem Documentation"
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doxygen
)
