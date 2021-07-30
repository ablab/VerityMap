project(common)

include_directories(.)
add_library(common STATIC
        cl_parser.cpp oneline_utils.hpp)
target_link_libraries (common m)