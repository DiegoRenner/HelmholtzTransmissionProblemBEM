macro(subdirlist result curdir)
    file(GLOB children RELATIVE ${curdir} ${curdir}/*)
    set(dirlist "")
    foreach(child ${children})
        if(IS_DIRECTORY ${curdir}/${child})
            list(APPEND dirlist ${child})
        endif()
    endforeach()
    set(${result} ${dirlist})
endmacro()

macro(add_all_subdirectories directory)
    subdirlist(subdirs ${directory})

    foreach(subdir ${subdirs})
        if(EXISTS "${directory}/${subdir}/CMakeLists.txt")
        add_subdirectory(${subdir})
        else()
            # check whether the subdirectory contains cpp files
            execute_process(COMMAND bash -c "find ${directory}/${subdir} | grep .cpp"
                RESULT_VARIABLE contains_cpp_code
                OUTPUT_QUIET)

            # if the directory contains cpp files we print a warning
            if(${contains_cpp_code} EQUAL 0)
                message("Skipping ${directory}/$subdir (containing cpp code) because no CMakeLists.txt file was found")
            else() # otherwise just a debug message
#               message("Skipping ${directory}/${subdir} since it did not contain any cpp code (e.g. not ported yet)")
            endif()
        endif()
    endforeach()
endmacro()

