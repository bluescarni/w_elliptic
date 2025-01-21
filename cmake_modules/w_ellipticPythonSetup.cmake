# Copyright (C) 2009-2011 by Francesco Biscani
# bluescarni@gmail.com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# NOTE: if and when we decide to support OSX, recover some of the quirks in the PaGMO
# version of this file. Let's keep it basic for the moment.

INCLUDE(FindPythonLibs)
# We need the Python interpreter to figure out Python's version in certain cases.
INCLUDE(FindPythonInterp)

# Find Python libraries
FIND_PACKAGE(PythonLibs REQUIRED)
MESSAGE(STATUS "Python libraries: " "${PYTHON_LIBRARIES}")
MESSAGE(STATUS "Python library: " "${PYTHON_LIBRARY}")

# These flags are used to signal the need to override the default extension of the Python modules
# depending on the architecture. Under Windows, for instance, CMake produces shared objects as
# .dll files, but Python from 2.5 onwards requires .pyd files (hence the need to override).
SET(PYDEXTENSION FALSE)

IF(UNIX)
	# We need the Python interpreter in order to detect the appropriate directory of modules installation.
	IF(NOT PYTHONINTERP_FOUND)
		MESSAGE(FATAL_ERROR "Unable to locate the Python interpreter.")
	ENDIF(NOT PYTHONINTERP_FOUND)
	MESSAGE(STATUS "Python interpreter is: ${PYTHON_EXECUTABLE}")
	# Now we must establish if the installation dir for Python modules is named 'site-packages' (as usual)
	# or 'dist-packages' (apparently Ubuntu 9.04 or maybe Python 2.6, it's not clear).
	EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cmake_modules/python_packages_dir.py
		OUTPUT_VARIABLE PY_PACKAGES_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
	MESSAGE(STATUS "Python packages dir is: ${PY_PACKAGES_DIR}")
	# In Unix systems we can fetch the Python version number directly from the library.
	STRING(REGEX MATCH libpython[0-9]*\\.[0-9]* PYTHON_LIBRARY_VERSION_DOT ${PYTHON_LIBRARY})
	# Remove the dot from the Python version.
	STRING(REGEX REPLACE libpython "" PYTHON_LIBRARY_VERSION_DOT ${PYTHON_LIBRARY_VERSION_DOT})
	STRING(REGEX REPLACE \\. "" PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY_VERSION_DOT})
	# Let's use CMAKE_INSTALL_PREFIX, so that if we specify a different install path it will be respected.
	SET(PYTHON_MODULES_PATH lib/python${PYTHON_LIBRARY_VERSION_DOT}/${PY_PACKAGES_DIR})
ELSE(UNIX)
	# On Windows, we will install directly into the install path of the Python interpreter.
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
		OUTPUT_VARIABLE _PYTHON_MODULES_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
	set(PYTHON_MODULES_PATH "${_PYTHON_MODULES_PATH}" CACHE PATH "Install path for Python modules.")
		mark_as_advanced(PYTHON_MODULES_PATH)
	IF(WIN32)
		message(STATUS "Windows platform detected.")
		message(STATUS "Output extension for compiled modules will be '.pyd'.")
		# .pyd extension is default on Windows with supported Python versions.
		SET(PYDEXTENSION TRUE)
		MESSAGE(STATUS "Windows platform detected. Output extension for compiled modules will be '.pyd'.")
		STRING(REGEX MATCH python[0-9]* PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY})
		STRING(REGEX REPLACE python "" PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY_VERSION})
	ENDIF(WIN32)
ENDIF(UNIX)

MESSAGE(STATUS "Python library version: " ${PYTHON_LIBRARY_VERSION})
MESSAGE(STATUS "Python modules install path: " "${PYTHON_MODULES_PATH}")

#IF(${PYTHON_LIBRARY_VERSION} LESS 26)
#	MESSAGE(FATAL_ERROR "Minimum supported Python version is 2.6.")
#ENDIF(${PYTHON_LIBRARY_VERSION} LESS 26)
