# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.18)

include(GNUInstallDirs)

set(_dftd4_python_standalone FALSE)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(_dftd4_python_standalone TRUE)
  project(
    dftd4-python
    LANGUAGES C
    VERSION 4.2.0
    DESCRIPTION "Python API of the DFT-D4 project"
  )
endif()

find_package(Python3 3.7 REQUIRED COMPONENTS "Interpreter" "Development.Module")
execute_process(
  COMMAND "${Python3_EXECUTABLE}" -c "import cffi, setuptools"
  RESULT_VARIABLE _dftd4_cffi_status
  ERROR_QUIET
)
if(NOT _dftd4_cffi_status EQUAL 0)
  message(FATAL_ERROR
    "The Python interpreter does not provide the cffi and setuptools modules"
  )
endif()

set(_dftd4_target "")
set(_dftd4_include_dir "")

# An in-tree build already provides the native target. A standalone Python
# build first tries an installed CMake package and then the source bundled in
# the Python sdist.
if(TARGET dftd4::dftd4)
  set(_dftd4_target "dftd4::dftd4")
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../include/dftd4.h")
    set(_dftd4_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/../include")
  endif()
elseif(TARGET dftd4::dftd4-lib)
  set(_dftd4_target "dftd4::dftd4-lib")
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../include/dftd4.h")
    set(_dftd4_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/../include")
  endif()
elseif(_dftd4_python_standalone)
  find_package(dftd4 "${PROJECT_VERSION}" CONFIG QUIET)
  if(TARGET dftd4::dftd4)
    set(_dftd4_target "dftd4::dftd4")
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/subprojects/dftd4/CMakeLists.txt")
    set(BUILD_SHARED_LIBS FALSE CACHE BOOL "Build shared libraries" FORCE)
    set(DFTD4_WITH_API TRUE CACHE BOOL "Build the DFT-D4 C API" FORCE)
    set(DFTD4_WITH_PYTHON FALSE CACHE BOOL "Build the DFT-D4 Python API" FORCE)
    add_subdirectory("subprojects/dftd4" "subprojects/dftd4")
    set(_dftd4_target "dftd4::dftd4")
    set(_dftd4_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/subprojects/dftd4/include")
  endif()
endif()

if(NOT _dftd4_target)
  message(FATAL_ERROR
    "Could not find dftd4 >= ${PROJECT_VERSION}; install its CMake package "
    "or build from a Python sdist containing subprojects/dftd4"
  )
endif()
if(DEFINED DFTD4_WITH_API AND NOT DFTD4_WITH_API)
  message(FATAL_ERROR "The dftd4 Python extension requires the dftd4 C API")
endif()

if(NOT _dftd4_include_dir)
  get_target_property(
    _dftd4_interface_include_dirs
    "${_dftd4_target}"
    INTERFACE_INCLUDE_DIRECTORIES
  )
  unset(_dftd4_include_dir)
  unset(_dftd4_include_dir CACHE)
  find_path(
    _dftd4_include_dir
    NAMES "dftd4.h"
    HINTS ${_dftd4_interface_include_dirs}
  )
  if(NOT _dftd4_include_dir)
    message(FATAL_ERROR "Could not locate dftd4.h from ${_dftd4_target}")
  endif()
endif()

set(_dftd4_header "${CMAKE_CURRENT_SOURCE_DIR}/include/_dftd4.h")
set(_dftd4_preprocessed_header "${CMAKE_CURRENT_BINARY_DIR}/_libdftd4.h")
set(_dftd4_cffi_source "${CMAKE_CURRENT_BINARY_DIR}/_libdftd4.c")

# CFFI cannot parse all preprocessor constructs in the public C header.
add_custom_command(
  OUTPUT "${_dftd4_preprocessed_header}"
  COMMAND
  "${CMAKE_C_COMPILER}"
  "-I${_dftd4_include_dir}"
  "-DDFTD4_CFFI"
  "-E"
  "${_dftd4_header}"
  "-o"
  "${_dftd4_preprocessed_header}"
  DEPENDS
  "${_dftd4_header}"
  "${_dftd4_include_dir}/dftd4.h"
  VERBATIM
)

add_custom_command(
  OUTPUT "${_dftd4_cffi_source}"
  COMMAND
  "${Python3_EXECUTABLE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/ffibuilder.py"
  "${_dftd4_preprocessed_header}"
  "_libdftd4"
  DEPENDS
  "${_dftd4_preprocessed_header}"
  "${CMAKE_CURRENT_SOURCE_DIR}/ffibuilder.py"
  VERBATIM
)

Python3_add_library("_libdftd4" MODULE WITH_SOABI "${_dftd4_cffi_source}")
target_include_directories("_libdftd4" PRIVATE "${_dftd4_include_dir}")
target_link_libraries("_libdftd4" PRIVATE "${_dftd4_target}")
if(CMAKE_Fortran_COMPILER_LOADED)
  set_property(TARGET "_libdftd4" PROPERTY LINKER_LANGUAGE "Fortran")
endif()

execute_process(
  COMMAND
    "${Python3_EXECUTABLE}"
    -c
    "import sysconfig, sys; print(sysconfig.get_path('platlib', vars={'base': '', 'platbase': '', 'platlibdir': sys.argv[1]}).lstrip('/'))"
    "${CMAKE_INSTALL_LIBDIR}"
  OUTPUT_VARIABLE _dftd4_python_install_root
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE _dftd4_sysconfig_status
)

if(NOT _dftd4_sysconfig_status EQUAL 0)
  message(FATAL_ERROR "Could not determine the Python platform-library directory")
endif()
set(
  DFTD4_PYTHON_INSTALL_DIR
  "${_dftd4_python_install_root}/dftd4"
  CACHE PATH
  "Directory for the dftd4 Python package"
)

install(
  TARGETS "_libdftd4"
  LIBRARY DESTINATION "${DFTD4_PYTHON_INSTALL_DIR}"
  COMPONENT "python"
  RUNTIME DESTINATION "${DFTD4_PYTHON_INSTALL_DIR}"
  COMPONENT "python"
)

set(
  _dftd4_python_sources
  "dftd4/__init__.py"
  "dftd4/ase.py"
  "dftd4/data.py"
  "dftd4/interface.py"
  "dftd4/library.py"
  "dftd4/parameters.py"
  "dftd4/pyscf.py"
  "dftd4/qcschema.py"
  "dftd4/references.json"
  "dftd4/test_ase.py"
  "dftd4/test_interface.py"
  "dftd4/test_library.py"
  "dftd4/test_parameters.py"
  "dftd4/test_pyscf.py"
  "dftd4/test_qcschema.py"
)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/dftd4/parameters.toml")
  list(APPEND _dftd4_python_sources "dftd4/parameters.toml")
endif()

install(
  FILES ${_dftd4_python_sources}
  DESTINATION "${DFTD4_PYTHON_INSTALL_DIR}"
  COMPONENT "python"
)
