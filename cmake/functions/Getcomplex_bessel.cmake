if(NOT EXISTS /usr/lib/libcomplex_bessel.so)
	if(NOT EXISTS "${PROJECT_SOURCE_DIR}/complex_bessel_lib")
		message(STATUS "No complex_bessel library found, installing locally.")
		execute_process(COMMAND git clone https://github.com/joeydumont/complex_bessel.git
						WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
		execute_process(COMMAND mkdir complex_bessel_lib
						WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
		execute_process(COMMAND bash build.sh "${PROJECT_SOURCE_DIR}/complex_bessel_lib"
						WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/complex_bessel")
		execute_process(COMMAND make install
						WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/complex_bessel/build")
		include_directories(${PROJECT_SOURCE_DIR}/complex_bessel_lib/include)
		link_directories(${PROJECT_SOURCE_DIR}/complex_bessel_lib/lib)
		execute_process(COMMAND rm -rf complex_bessel
					WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
	else()
		include_directories(${PROJECT_SOURCE_DIR}/complex_bessel_lib/include)
		link_directories(${PROJECT_SOURCE_DIR}/complex_bessel_lib/lib)
		message(STATUS "Found complex_bessel library installed locally.")
	endif()
else()
	message(STATUS "Found complex_bessel library installed in /usr/lib directory.") 
	if(NOT EXISTS "${PROJECT_SOURCE_DIR}/include/complex_bessel.h")
		message(STATUS "Headers are missing, grabbing them.")
		execute_process(COMMAND git clone https://github.com/joeydumont/complex_bessel.git
						WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
		execute_process(COMMAND mkdir complex_bessel_lib
						WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
		execute_process(COMMAND bash build.sh "${PROJECT_SOURCE_DIR}/complex_bessel_lib"
						WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/complex_bessel")
		execute_process(COMMAND make install
						WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/complex_bessel/build")
		execute_process(COMMAND cp -a complex_bessel_lib/include/. include
						WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
		execute_process(COMMAND rm -rf complex_bessel
					WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
		execute_process(COMMAND rm -rf complex_bessel_lib
					WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
	endif()
endif()

