add_custom_target(doxygen_parametricbem2d
COMMAND O_DIRECTORY=${CMAKE_BINARY_DIR}/doxygen_parametricbem2d doxygen Doxyfile
COMMENT "Generating 2DParametricBEM Documentation"
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doxygen
)
