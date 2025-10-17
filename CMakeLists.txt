cmake_minimum_required (VERSION 3.15)
project (rlz-repair)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS "-O3")
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
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
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
message(STATUS "The source directory of the project is ${CMAKE_CURRENT_LIST_DIR}")

# include_directories has been legacy since 2012 
# https://stackoverflow.com/questions/13703647/how-to-properly-add-include-directories-with-cmake
set(PROJECT_SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/src)
set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/include)
set(THIRDPARTY_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty)

message(STATUS "The include directories of the project is ${INCLUDE_DIRECTORIES}")

# Log directory for external content builds
make_directory(${CMAKE_CURRENT_BINARY_DIR}/logs)
set(INSTALL_FILE  ${CMAKE_CURRENT_BINARY_DIR}/install.cmake)

# Set install directories
make_directory(${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")

# Set search paths for CMake
list(APPEND CMAKE_LIBRARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/lib")
list(APPEND CMAKE_INCLUDE_PATH "${CMAKE_CURRENT_BINARY_DIR}/include")
list(APPEND CMAKE_PREFIX_PATH  "${CMAKE_CURRENT_BINARY_DIR}")

# RLZ
find_program(RLZ rlz)
if (NOT RLZ)
    message(STATUS "rlz not found. Building from submodule...")
    # Ensure the submodule is initialized
    if(NOT EXISTS ${THIRDPARTY_INCLUDE_DIR}/RLZ)
        message(FATAL_ERROR "RLZ submodule is not initialized. Run 'git submodule update --init --recursive'")
    endif()
    add_subdirectory(${THIRDPARTY_INCLUDE_DIR}/RLZ ${CMAKE_CURRENT_BINARY_DIR}/rlz_build)
    set(RLZ ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/rlz CACHE FILEPATH "Path to rlz executable" FORCE)
else()
    message(STATUS "rlz found at ${RLZ}.")
endif()

# RePair
add_executable(repair 
               ${PROJECT_SRC_DIR}/repair.cpp 
               ${PROJECT_SRC_DIR}/repair/heap.c
               ${PROJECT_SRC_DIR}/repair/array.c
               ${PROJECT_SRC_DIR}/repair/basics.c
               ${PROJECT_SRC_DIR}/repair/records.c
               ${PROJECT_SRC_DIR}/repair/hash.c)

add_library(malloc_count STATIC ${THIRDPARTY_INCLUDE_DIR}/malloc_count/malloc_count.c)
target_link_libraries(repair PRIVATE malloc_count dl)
target_include_directories(repair PUBLIC ${PROJECT_INCLUDE_DIR} ${PROJECT_INCLUDE_DIR}/repair ${THIRDPARTY_INCLUDE_DIR} ${THIRDPARTY_INCLUDE_DIR}/CLI ${THIRDPARTY_INCLUDE_DIR}/malloc_count)

# Gonzalo's Integer RePair
add_executable(gonzalo_repair
               ${PROJECT_SRC_DIR}/repair_large/irepair.c
               ${PROJECT_SRC_DIR}/repair_large/array.c
               ${PROJECT_SRC_DIR}/repair_large/hash.c
               ${PROJECT_SRC_DIR}/repair_large/heap.c
               ${PROJECT_SRC_DIR}/repair_large/records.c
               ${PROJECT_SRC_DIR}/repair_large/basics.c)

target_include_directories(gonzalo_repair PUBLIC ${PROJECT_INCLUDE_DIR}/repair_large)

# Despair
add_executable(despair
               ${PROJECT_SRC_DIR}/decompress/despair.c
               ${PROJECT_SRC_DIR}/decompress/basics.c)
target_include_directories(despair PUBLIC ${PROJECT_INCLUDE_DIR}/decompress)

# Setting up python pipeline script
configure_file(${CMAKE_CURRENT_LIST_DIR}/pipeline/rlz-repair.in
                ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/rlz-repair
                @ONLY)

file(CHMOD ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/rlz-repair PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

#################### Packaging ##############################
include(GNUInstallDirs)

install(TARGETS repair gonzalo_repair despair
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(PROGRAMS "${RLZ}" DESTINATION ${CMAKE_INSTALL_BINDIR})

install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/rlz-repair" DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES README.md LICENSE
        DESTINATION ${CMAKE_INSTALL_DOCDIR})

set(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
set(CPACK_PACKAGE_VERSION "${VERSION}")


set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "RLZ-RePair: A tool for Relative Lempel-Ziv compression combined with RePair")
set(CPACK_PACKAGE_VENDOR "University of Florida")
set(CPACK_PACKAGE_CONTACT "rvarki@ufl.edu")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/rvarki/RLZ-RePair")

set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CMAKE_SYSTEM_NAME}")

if(APPLE)
    # Use DragNDrop for .dmg files, plus tarballs
    set(CPACK_GENERATOR "DragNDrop;STGZ;TGZ")
else()
    # Use DEB for .deb files, plus tarballs
    set(CPACK_GENERATOR "STGZ;TGZ;DEB")
endif()

# --- Set Debian Package Architecture ---
# Get CMake's idea of the CPU
string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" deb_arch)

# Translate common CMake names to Debian names
if(deb_arch MATCHES "x86_64")
    set(deb_arch "amd64")
elseif(deb_arch MATCHES "aarch64")
    set(deb_arch "arm64")
elseif(deb_arch MATCHES "i[3-6]86")
    set(deb_arch "i386")
endif()

# Set the CPack variable
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${deb_arch})

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Rahul Varki")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

include(CPack)