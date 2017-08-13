# Initial Python setup.
include(w_ellipticPythonSetup)

# Python headers.
include_directories(${PYTHON_INCLUDE_PATH})

# Boost Python library setup.
if(${PYTHON_VERSION_MAJOR} EQUAL 2)
    find_package(Boost 1.48.0 REQUIRED COMPONENTS "python")
    set(MANDATORY_LIBRARIES ${Boost_PYTHON_LIBRARY})
else()
    find_package(Boost 1.48.0 REQUIRED COMPONENTS "python3")
    set(MANDATORY_LIBRARIES ${Boost_PYTHON3_LIBRARY})
endif()
include_directories(${Boost_INCLUDE_DIRS})

add_definitions(-DBOOST_ALL_NO_LIB)

# Install path.
set(W_ELLIPTIC_INSTALL_PATH "${PYTHON_MODULES_PATH}/w_elliptic")

add_library(_core SHARED
    # General.
    core.cpp
    # Python files.
    #__init__.py
)

# Prefix and suffix setup.
set_target_properties(_core PROPERTIES PREFIX "")
if(PYDEXTENSION)
    set_target_properties(_core PROPERTIES SUFFIX ".pyd")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # Apparently on OS X Python expects the .so extension for compiled modules.
    message(STATUS "OS X detected, the Python bindings' module suffix will be '.so'.")
    set_target_properties(_core PROPERTIES SUFFIX ".so")
endif()

# Strict aliasing setup.
# NOTE: Python < 3 is not strict C, as per:
# http://www.python.org/dev/peps/pep-3123/
# When using GCC or Clang, we add -fno-strict-aliasing here and when compiling core.cpp to prevent possible problems.
# Other compilers should probably use similar flags.
if(${PYTHON_LIBRARY_VERSION} LESS 30 AND (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANGXX))
    message(STATUS "Python < 3 version detected: the Python bindings will be compiled with the '-fno-strict-aliasing' flag.")
    set_target_properties(_core PROPERTIES COMPILE_FLAGS "-fno-strict-aliasing" LINK_FLAGS "-fno-strict-aliasing")
endif()

# Link & install.
target_link_libraries(_core ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARIES})
install(TARGETS _core
    RUNTIME DESTINATION ${W_ELLIPTIC_INSTALL_PATH}
    LIBRARY DESTINATION ${W_ELLIPTIC_INSTALL_PATH}
)

# Add the global init file.
install(FILES __init__.py DESTINATION ${W_ELLIPTIC_INSTALL_PATH})
