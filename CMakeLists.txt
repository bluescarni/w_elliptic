project(w_elliptic)

enable_testing()

cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")

message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")

# Set default build type to "Release".
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()

# Build option: enable test set.
option(BUILD_TESTS "Build test set." OFF)
# Build option: build the Python bindings.
option(BUILD_PYTHON_BINDINGS "Build the Python bindings." OFF)

# Setup compiler.
include(w_ellipticCompilerLinkerSettings)

# Tests setup.
if(BUILD_TESTS)
    add_subdirectory("${CMAKE_SOURCE_DIR}/tests")
endif()

# Python bindings setup.
if(BUILD_PYTHON_BINDINGS)
    add_subdirectory("${CMAKE_SOURCE_DIR}/python")
endif()

# Install the headers.
#install(FILES src/arbpp.hpp DESTINATION include/arbpp)
