include(FetchContent)

cmake_minimum_required (VERSION 3.15)
project (rlz-repair)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(VERSION_MAJOR "1")
SET(VERSION_MINOR "0")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

###############################################################################
# Set install directory if not specified
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "Installation directory")
################################################################################
# Specify build type
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install directory: ${CMAKE_INSTALL_PREFIX}")
################################################################################
# Locating Git
find_package(Git)
if(GIT_FOUND)
    message(STATUS "git found: ${GIT_EXECUTABLE}")
else()
    message(FATAL_ERROR "git not found. Cloning of submodules will not work.")
endif()
################################################################################
# Version number of project
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "Git current branch of project: ${GIT_BRANCH}")
message(STATUS "Git commit hash of project: ${GIT_COMMIT_HASH}")
################################################################################
# Compiler and linker flags
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(BUILD_TYPE_FLAGS "-O3")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(BUILD_TYPE_FLAGS "-g")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(BUILD_TYPE_FLAGS "-O3 -g -DNDEBUG")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BUILD_TYPE_FLAGS}")
#################################################################################
# Check dependencies
message(STATUS "The source directory of the project is ${CMAKE_SOURCE_DIR}")

# include_directories has been legacy since 2012 
# https://stackoverflow.com/questions/13703647/how-to-properly-add-include-directories-with-cmake
set(PROJECT_SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
set(THIRDPARTY_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/thirdparty)

message(STATUS "The include directories of the project is ${INCLUDE_DIRECTORIES}")

# Log directory for external content builds
make_directory(${CMAKE_CURRENT_BINARY_DIR}/logs)
set(INSTALL_FILE  ${CMAKE_CURRENT_BINARY_DIR}/install.cmake)


find_program(RLZ rlz)
if (NOT RLZ)
    message(STATUS "RLZ not found. Building as external content.")
    FetchContent_Declare(
            rlz
            GIT_REPOSITORY "https://github.com/rvarki/RLZ.git"
            GIT_TAG update_build
    )
    if (NOT rlz_POPULATED)
        FetchContent_MakeAvailable(rlz)
    endif()
endif()
