project(common)

include_directories(.)
add_library(common STATIC
        cl_parser.cpp oneline_utils.hpp coverage_utils.cpp)
target_link_libraries (common m)