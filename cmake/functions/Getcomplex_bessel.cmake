#execute_process(COMMAND git clone https://www.github.com/DiegoRenner/complex_bessel
#				COMMAND cd complex_bessel
#				COMMAND bash build.sh
#				COMMAND cd build
#				COMMAND make install
#				COMMAND echo test)
execute_process(COMMAND echo test)
execute_process(COMMAND git clone git@github.com:DiegoRenner/complex_bessel
				WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
			execute_process(COMMAND bash build.sh "${PROJECT_SOURCE_DIR}"
				WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/complex_bessel")
execute_process(COMMAND make install
				WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/complex_bessel/build")
execute_process(COMMAND ls complex_bessel_lib
				WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
execute_process(COMMAND cp -r complex_bessel_lib/include/. include
				WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
execute_process(COMMAND cp -r complex_bessel_lib/lib .
				WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
link_directories(${PROJECT_SOURCE_DIR}/lib)
execute_process(COMMAND rm -rf complex_bessel
			WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
execute_process(COMMAND rm -rf complex_bessel_lib
				WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
