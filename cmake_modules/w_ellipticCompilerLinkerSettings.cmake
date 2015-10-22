INCLUDE(CheckCXXCompilerFlag)

message(STATUS "The C++ compiler ID is: ${CMAKE_CXX_COMPILER_ID}")

# Clang detection:
# http://stackoverflow.com/questions/10046114/in-cmake-how-can-i-test-if-the-compiler-is-clang
# http://www.cmake.org/cmake/help/v2.8.10/cmake.html#variable:CMAKE_LANG_COMPILER_ID
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_COMPILER_IS_CLANGXX 1)
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_COMPILER_IS_INTELXX 1)
endif()

macro(W_ELLIPTIC_CHECK_ENABLE_CXX_FLAG flag)
    set(W_ELLIPTIC_CHECK_CXX_FLAG)
    check_cxx_compiler_flag("${flag}" W_ELLIPTIC_CHECK_CXX_FLAG)
    if(W_ELLIPTIC_CHECK_CXX_FLAG)
        message(STATUS "Enabling the '${flag}' compiler flag.")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
    else()
        message(STATUS "Disabling the '${flag}' compiler flag.")
    endif()
    unset(W_ELLIPTIC_CHECK_CXX_FLAG CACHE)
endmacro()

macro(W_ELLIPTIC_CHECK_ENABLE_DEBUG_CXX_FLAG flag)
    set(W_ELLIPTIC_CHECK_DEBUG_CXX_FLAG)
    check_cxx_compiler_flag("${flag}" W_ELLIPTIC_CHECK_DEBUG_CXX_FLAG)
    if(W_ELLIPTIC_CHECK_DEBUG_CXX_FLAG)
        message(STATUS "Enabling the '${flag}' debug compiler flag.")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${flag}")
    else()
        message(STATUS "Disabling the '${flag}' debug compiler flag.")
    endif()
    unset(W_ELLIPTIC_CHECK_DEBUG_CXX_FLAG CACHE)
endmacro()

# Configuration for GCC.
if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "GNU compiler detected, checking version.")
    try_compile(GCC_VERSION_CHECK ${CMAKE_BINARY_DIR} "${CMAKE_SOURCE_DIR}/cmake_modules/gcc_check_version.cpp")
    if(NOT GCC_VERSION_CHECK)
        MESSAGE(FATAL_ERROR "Unsupported GCC version, please upgrade your compiler.")
    endif()
    message(STATUS "GCC version is ok.")
    # Color diagnostic available since GCC 4.9.
    W_ELLIPTIC_CHECK_ENABLE_CXX_FLAG(-fdiagnostics-color=auto)
endif()

# Configuration for Clang.
if(CMAKE_COMPILER_IS_CLANGXX)
    message(STATUS "Clang compiler detected, checking version.")
    try_compile(CLANG_VERSION_CHECK ${CMAKE_BINARY_DIR} "${CMAKE_SOURCE_DIR}/cmake_modules/clang_check_version.cpp")
    if(NOT CLANG_VERSION_CHECK)
        MESSAGE(FATAL_ERROR "Unsupported Clang version, please upgrade your compiler.")
    endif()
endif()

# Configuration for the Intel compiler.
if(CMAKE_COMPILER_IS_INTELXX)
    message(STATUS "Intel compiler detected, checking version.")
    try_compile(INTEL_VERSION_CHECK ${CMAKE_BINARY_DIR} "${CMAKE_SOURCE_DIR}/cmake_modules/intel_check_version.cpp")
    if(NOT INTEL_VERSION_CHECK)
        MESSAGE(FATAL_ERROR "Unsupported Intel compiler version, please upgrade your compiler.")
    endif()
    # These diagnostic messages from ICC are often wrong. Use another compiler when debugging.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -diag-disable 2304,2305,1682,2259,3373")
endif()

# Common configuration for GCC, Clang and Intel.
if(CMAKE_COMPILER_IS_CLANGXX OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_INTELXX)
    # Enable the C++11 flag. Need to check for either c++0x or c++11.
    check_cxx_compiler_flag("-std=c++11" W_ELLIPTIC_CHECK_CPP11_FLAG)
    if(W_ELLIPTIC_CHECK_CPP11_FLAG)
        message(STATUS "Enabling the '-std=c++11' flag.")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    else()
        message(STATUS "Enabling the '-std=c++0x' flag.")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    endif()
    W_ELLIPTIC_CHECK_ENABLE_CXX_FLAG(-Wall)
    W_ELLIPTIC_CHECK_ENABLE_CXX_FLAG(-Wextra)
    W_ELLIPTIC_CHECK_ENABLE_CXX_FLAG(-Wnon-virtual-dtor)
    W_ELLIPTIC_CHECK_ENABLE_CXX_FLAG(-Wnoexcept)
    W_ELLIPTIC_CHECK_ENABLE_CXX_FLAG(-Wlogical-op)
    W_ELLIPTIC_CHECK_ENABLE_CXX_FLAG(-Wconversion)
    # NOTE: this can be useful, but at the moment it triggers lots of warnings.
    # W_ELLIPTIC_CHECK_ENABLE_CXX_FLAG(-Wold-style-cast)
    # NOTE: disable this for now, as it results in a lot of clutter from Boost.
    # W_ELLIPTIC_CHECK_ENABLE_CXX_FLAG(-Wzero-as-null-pointer-constant)
    W_ELLIPTIC_CHECK_ENABLE_CXX_FLAG(-pedantic-errors)
    W_ELLIPTIC_CHECK_ENABLE_CXX_FLAG(-Wdisabled-optimization)
    W_ELLIPTIC_CHECK_ENABLE_CXX_FLAG(-fvisibility-inlines-hidden)
    W_ELLIPTIC_CHECK_ENABLE_CXX_FLAG(-fvisibility=hidden)
    W_ELLIPTIC_CHECK_ENABLE_DEBUG_CXX_FLAG(-fstack-protector-all)
endif()
