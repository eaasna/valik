cmake_minimum_required (VERSION 3.8)

# A interface to reuse common properties.
# You can add more external include paths of other projects that are needed for your project.

# Shared interface_
add_library ("${PROJECT_NAME}_interface" INTERFACE)
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE seqan3::seqan3)
target_include_directories ("${PROJECT_NAME}_interface" INTERFACE ../include)
target_include_directories ("${PROJECT_NAME}_interface" INTERFACE ../lib/robin-hood-hashing/src/include)
#target_include_directories ("${PROJECT_NAME}_interface" INTERFACE ${Boost_INCLUDE_DIRS})
target_compile_options ("${PROJECT_NAME}_interface" INTERFACE "-pedantic" "-Wall" "-Wextra")

# Sliding window split
add_library ("${PROJECT_NAME}_split_lib" STATIC valik_split.cpp)
target_link_libraries ("${PROJECT_NAME}_split_lib" PUBLIC "${PROJECT_NAME}_interface")

# Sliding window build
add_library ("${PROJECT_NAME}_build_lib" STATIC valik_build.cpp)
target_link_libraries ("${PROJECT_NAME}_build_lib" PUBLIC "${PROJECT_NAME}_interface")

# Sliding window search
add_library ("${PROJECT_NAME}_search_helper_lib" STATIC search/detail/helper.cpp)
target_link_libraries ("${PROJECT_NAME}_search_helper_lib" PUBLIC "${PROJECT_NAME}_interface")

add_library ("${PROJECT_NAME}_minimiser_model_lib" STATIC search/minimiser_model.cpp)
target_link_libraries ("${PROJECT_NAME}_minimiser_model_lib" PUBLIC "${PROJECT_NAME}_search_helper_lib")

add_library ("${PROJECT_NAME}_simple_model_lib" STATIC search/compute_simple_model.cpp)
target_link_libraries ("${PROJECT_NAME}_simple_model_lib" PUBLIC "${PROJECT_NAME}_minimiser_model_lib")

add_library ("${PROJECT_NAME}_search_lib" STATIC valik_search.cpp)
target_link_libraries ("${PROJECT_NAME}_search_lib" PUBLIC "${PROJECT_NAME}_simple_model_lib")

# Sliding window argument parsing
add_library ("${PROJECT_NAME}_argument_parsing_shared_lib" STATIC argument_parsing/shared.cpp)
target_link_libraries ("${PROJECT_NAME}_argument_parsing_shared_lib" PUBLIC "${PROJECT_NAME}_interface")

add_library ("${PROJECT_NAME}_argument_parsing_split_lib" STATIC argument_parsing/split.cpp)
target_link_libraries ("${PROJECT_NAME}_argument_parsing_split_lib" PUBLIC "${PROJECT_NAME}_argument_parsing_shared_lib")

add_library ("${PROJECT_NAME}_argument_parsing_build_lib" STATIC argument_parsing/build.cpp)
target_link_libraries ("${PROJECT_NAME}_argument_parsing_build_lib" PUBLIC "${PROJECT_NAME}_argument_parsing_shared_lib")

add_library ("${PROJECT_NAME}_argument_parsing_search_lib" STATIC argument_parsing/search.cpp)
target_link_libraries ("${PROJECT_NAME}_argument_parsing_search_lib" PUBLIC "${PROJECT_NAME}_argument_parsing_shared_lib")

add_library ("${PROJECT_NAME}_argument_parsing_top_level_lib" STATIC argument_parsing/top_level.cpp)
target_link_libraries ("${PROJECT_NAME}_argument_parsing_top_level_lib" PUBLIC "${PROJECT_NAME}_argument_parsing_shared_lib")

# Sliding window library
add_library ("${PROJECT_NAME}_lib" INTERFACE)
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_argument_parsing_split_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_argument_parsing_build_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_argument_parsing_search_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_argument_parsing_shared_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_argument_parsing_top_level_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_split_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_build_lib")
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE "${PROJECT_NAME}_search_lib")

# Sliding window executable
add_executable ("${PROJECT_NAME}" valik_main.cpp)
target_link_libraries ("${PROJECT_NAME}" PUBLIC "${PROJECT_NAME}_lib")
