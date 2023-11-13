#This cmakelist is based on music cmakelist file

cmake_minimum_required(VERSION 3.4)
project(Preequilibrium VERSION 2.0 LANGUAGES CXX)


set(exename "Preequilibrium")

# require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

set(CMAKE_INSTALL_RPATH ${CMAKE_HOME_DIRECTORY})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# Find GSL and openMP for this project:
find_package(GSL)
if (GSL_FOUND)
    message("Found GSL library ${GSL_INCLUDE_DIR}")
    include_directories(${GSL_INCLUDE_DIR})
elseif()
    message("Did not find GSL library")
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
        message(FATAL_ERROR "Insufficient gcc version")
    endif()
endif()


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    message("Using Clang compiler without OpenMP parallelization... ")
else ()
    find_package(OpenMP)
    if (OpenMP_FOUND)
        set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS})
        set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${OpenMP_EXE_LINKER_FLAGS})
    endif ()
endif()


# default build type: Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_HOME_DIRECTORY}/"
        CACHE PATH "Install path prefix, prepended onto install directories."
        FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "-ipo ${OpenMP_CXX_FLAGS} -std=c++11")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -std=c++11")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS "-std=c++11  -stdlib=libc++")
endif()

string(APPEND CMAKE_CXX_FLAGS " -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-g -O3")

add_subdirectory (src)