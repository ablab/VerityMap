cmake_minimum_required(VERSION 3.12)
project(veritymap)

set(CMAKE_CXX_STANDARD 20)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-ggdb3 -pthread -lz -lrt -fopenmp" )
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg")
set(CMAKE_SHARED_LINKER_FLAGS "-static-libgcc -static-libstdc++ -Wall -Wc++-compat -O2 -msse4.1 -DHAVE_KALLOC -DKSW_CPU_DISPATCH -D_FILE_OFFSET_BITS=64 -ltbb -fsigned-char -fsanitize=address -pg")
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

# Define various dirs
set(MAIN_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
# set(LIB_SRC_DIR ${MAIN_SRC_DIR}/lib)
set(PROJECTS_SRC_DIR ${MAIN_SRC_DIR}/projects)
# set(TESTS_SRC_DIR ${MAIN_SRC_DIR}/tests)
set(TOOLS_SRC_DIR ${MAIN_SRC_DIR}/tools)

set(VM_SRC_DIR ${PROJECTS_SRC_DIR}/veritymap)



find_package(Git)
find_package (ZLIB)
include_directories(SYSTEM "${ZLIB_INCLUDE_DIRS}")

# the commit's SHA1, and whether the building workspace was dirty or not
execute_process(COMMAND
        "${GIT_EXECUTABLE}" describe --match=NeVeRmAtCh --always --abbrev=40 --dirty
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_SHA1
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

# the date of the commit
execute_process(COMMAND
        "${GIT_EXECUTABLE}" log -1 --format=%ad --date=local
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_DATE
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

# the subject of the commit
execute_process(COMMAND
        "${GIT_EXECUTABLE}" log -1 --format=%s
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_COMMIT_SUBJECT
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

# generate version.cpp
configure_file("${TOOLS_SRC_DIR}/version/version.cpp.in" ${TOOLS_SRC_DIR}/version/version.cpp @ONLY)
set(VERSION ${TOOLS_SRC_DIR}/version/version.cpp)

include_directories(${TOOLS_SRC_DIR})
# include_directories(${TESTS_SRC_DIR})
include_directories(${PROJECT_SRC_DIR})

add_subdirectory(src/projects)
add_subdirectory(src/tools)


