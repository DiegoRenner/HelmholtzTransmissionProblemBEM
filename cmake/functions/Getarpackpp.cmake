if(NOT EXISTS "${PROJECT_SOURCE_DIR}/arpackpp_lib")
		message(STATUS "No arpackpp library found, installing from github.")
		execute_process(COMMAND git clone git@github.com:m-reuter/arpackpp.git
						WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
		execute_process(COMMAND mkdir -p arpackpp_lib
						WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
		execute_process(COMMAND cp -r arpackpp/include arpackpp_lib/include
						WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
		include_directories(${PROJECT_SOURCE_DIR}/arpackpp_lib/include)
		execute_process(COMMAND rm -rf arpackpp
					WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
endif()

