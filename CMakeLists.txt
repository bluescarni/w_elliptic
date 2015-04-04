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

# Compiler setup.
if(BUILD_TESTS)
    # Setup compiler.
    include(w_ellipticCompilerLinkerSettings)
    # Boost unit test library.
    find_package(Boost 1.48.0 REQUIRED COMPONENTS "unit_test_framework")
    include_directories(${Boost_INCLUDE_DIRS})
    # Assemble all libraries and add the tests subdirectory.
    set(MANDATORY_LIBRARIES ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    add_subdirectory("${CMAKE_SOURCE_DIR}/tests")
endif()

# Install the headers.
#install(FILES src/arbpp.hpp DESTINATION include/arbpp)
