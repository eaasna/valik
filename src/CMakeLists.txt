cmake_minimum_required (VERSION 3.16)

# A interface to reuse common properties.
# You can add more external include paths of other projects that are needed for your project.

# Shared interface_
add_library ("${PROJECT_NAME}_interface" INTERFACE)
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE seqan3::seqan3)
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE sharg::sharg)
target_include_directories ("${PROJECT_NAME}_interface" INTERFACE ../include)
# SYSTEM keyword turns off compiler warnings for subprojects
target_include_directories ("${PROJECT_NAME}_interface" SYSTEM INTERFACE ../lib/stellar3/include)
target_include_directories ("${PROJECT_NAME}_interface" SYSTEM INTERFACE ../lib/seqan/include)
target_compile_options ("${PROJECT_NAME}_interface" INTERFACE "-pedantic" "-Wall" "-Wextra" "-Wno-interference-size")

# IBF size 
add_library ("${PROJECT_NAME}_index_size_lib" STATIC prepare/compute_bin_size.cpp)
target_link_libraries ("${PROJECT_NAME}_index_size_lib" PUBLIC "${PROJECT_NAME}_interface")

# Parameter tuning
add_library ("${PROJECT_NAME}_threshold_lib" STATIC threshold/find.cpp)
target_link_libraries ("${PROJECT_NAME}_threshold_lib" PUBLIC "${PROJECT_NAME}_interface")

# Valik split
add_library ("${PROJECT_NAME}_split_lib" STATIC valik_split.cpp)
target_link_libraries ("${PROJECT_NAME}_split_lib" PUBLIC "${PROJECT_NAME}_interface")

add_library ("${PROJECT_NAME}_split_write_lib" STATIC split/write_seg_sequences.cpp)
target_link_libraries ("${PROJECT_NAME}_split_write_lib" PUBLIC "${PROJECT_NAME}_interface")

# Valik build
add_library ("${PROJECT_NAME}_build_lib" STATIC valik_build.cpp)
target_link_libraries ("${PROJECT_NAME}_build_lib" PUBLIC "${PROJECT_NAME}_interface")

# Valik search
add_library ("raptor_threshold" STATIC
             raptor/threshold/threshold.cpp
             raptor/threshold/one_indirect_error_model.cpp
             raptor/threshold/multiple_error_model.cpp
             raptor/threshold/pascal_row.cpp
             raptor/threshold/one_error_model.cpp
             raptor/threshold/precompute_correction.cpp
             raptor/threshold/precompute_threshold.cpp
)
target_link_libraries ("raptor_threshold" PUBLIC "${PROJECT_NAME}_interface")

add_library ("${PROJECT_NAME}_search_lib" STATIC valik_search.cpp)
target_link_libraries ("${PROJECT_NAME}_search_lib" PUBLIC "raptor_threshold")

# Consolidating matches from split Stellar runs
add_library ("${PROJECT_NAME}_consolidation_lib" STATIC consolidate/consolidate_matches.cpp)
target_link_libraries ("${PROJECT_NAME}_consolidation_lib" PUBLIC "raptor_threshold")

add_library ("${PROJECT_NAME}_merge_processes_lib" STATIC consolidate/merge_processes.cpp)
target_link_libraries ("${PROJECT_NAME}_merge_processes_lib" PUBLIC "raptor_threshold")

# Argument parsing
add_library ("${PROJECT_NAME}_argument_parsing_shared_lib" STATIC argument_parsing/shared.cpp)
target_link_libraries ("${PROJECT_NAME}_argument_parsing_shared_lib" PUBLIC "${PROJECT_NAME}_interface")

add_library ("${PROJECT_NAME}_argument_parsing_build_lib" STATIC argument_parsing/build.cpp)
target_link_libraries ("${PROJECT_NAME}_argument_parsing_build_lib" PUBLIC "${PROJECT_NAME}_argument_parsing_shared_lib")

add_library ("${PROJECT_NAME}_argument_parsing_search_lib" STATIC argument_parsing/search.cpp)
target_link_libraries ("${PROJECT_NAME}_argument_parsing_search_lib" PUBLIC "${PROJECT_NAME}_argument_parsing_shared_lib")

add_library ("${PROJECT_NAME}_argument_parsing_top_level_lib" STATIC argument_parsing/top_level.cpp)
target_link_libraries ("${PROJECT_NAME}_argument_parsing_top_level_lib" PUBLIC "${PROJECT_NAME}_argument_parsing_shared_lib")

# Valik library
add_library ("${PROJECT_NAME}_lib" INTERFACE)
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_argument_parsing_build_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_argument_parsing_search_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_argument_parsing_shared_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_argument_parsing_top_level_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_split_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_split_write_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_build_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_search_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_consolidation_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_merge_processes_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_threshold_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_index_size_lib")

# Valik executable
add_executable ("${PROJECT_NAME}" valik_main.cpp)
target_link_libraries ("${PROJECT_NAME}" PUBLIC "${PROJECT_NAME}_lib")

# Log git commit
execute_process(COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE COMMIT_ID
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_custom_command(TARGET "${PROJECT_NAME}" POST_BUILD
    WORKING_DIRECTORY
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    DEPENDS
        ${COMMIT_ID}
    COMMAND
        ${CMAKE_COMMAND} -E echo ${COMMIT_ID} > version.md
    COMMENT
        "Generating file version.md"
    VERBATIM)
