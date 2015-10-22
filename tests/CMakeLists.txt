# Boost unit test library setup.
find_package(Boost 1.48.0 REQUIRED COMPONENTS "unit_test_framework")
include_directories(${Boost_INCLUDE_DIRS})
set(MANDATORY_LIBRARIES ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# NOTE: a bit hackish, as in principle Boost could
# be compiled statically. Proper way would be to check
# this with a try_compile or similar.
add_definitions(-DBOOST_TEST_DYN_LINK)

# http://stackoverflow.com/questions/1346731/optional-argument-in-cmake-macro
macro(ADD_W_ELLIPTIC_TESTCASE arg1)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_executable(${arg1} ${arg1}.cpp ${ARGN})
        target_link_libraries(${arg1} ${MANDATORY_LIBRARIES})
        add_test(${arg1} ${arg1})
    endif()
endmacro()

macro(ADD_W_ELLIPTIC_PERFORMANCE_TESTCASE arg1)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        add_executable("${arg1}_perf" ${arg1}.cpp)
        target_link_libraries("${arg1}_perf" ${MANDATORY_LIBRARIES})
        add_test("${arg1}_perf" "${arg1}_perf")
    endif()
endmacro()

ADD_W_ELLIPTIC_TESTCASE(w_elliptic)

# Add the example separately, as it needs a special -D.
add_executable(example example.cpp)
target_include_directories(example PUBLIC "${CMAKE_SOURCE_DIR}/src")
target_link_libraries(example ${MANDATORY_LIBRARIES})
add_test(example example)
