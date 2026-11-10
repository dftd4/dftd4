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

cmake_minimum_required(VERSION 3.22)

project(
  "dftd4"
  LANGUAGES "Fortran"
  VERSION "4.3.1"
  DESCRIPTION "Generally Applicable Atomic-Charge Dependent London Dispersion Correction"
)

include(GNUInstallDirs)
include(CMakeDependentOption)

# Build options
option(BUILD_SHARED_LIBS "Whether the libraries built should be shared" OFF)
option(BUILD_TESTING "Build unittests" ON)

# Unprefixed options are deprecated and will be removed in the future
# For now they only provide the defaults; please use DFTD4_WITH_* instead.
option(WITH_API "Enable support for C API via iso_c_binding module" ON)
option(WITH_OpenMP "Enable support for shared memory parallelisation with OpenMP" ON)
option(WITH_ILP64 "Enable support for ILP64 BLAS/LAPACK calls" OFF)
option(WITH_API_V2 "Enable the compatibility layer for the dftd4 2.5.x API" OFF)

option(DFTD4_WITH_API "Enable support for C API via iso_c_binding module" ${WITH_API})
option(DFTD4_WITH_OpenMP "Enable support for shared memory parallelisation with OpenMP" ${WITH_OpenMP})
option(DFTD4_WITH_ILP64 "Enable support for ILP64 BLAS/LAPACK calls" ${WITH_ILP64})
option(DFTD4_WITH_API_V2 "Enable the compatibility layer for the dftd4 2.5.x API" ${WITH_API_V2})
option(DFTD4_WITH_TESTS "Build unittests of dftd4" ${BUILD_TESTING})
cmake_dependent_option(DFTD4_WITH_PYTHON "Build the Python extension module"
                       FALSE "DFTD4_WITH_API" FALSE)

set(
  DFTD4_INSTALL_FMODDIR
  "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}"
  CACHE STRING "Installation directory for Fortran modules"
)

if(NOT DEFINED "${PROJECT_NAME}-dependency-method")
  set(
    "${PROJECT_NAME}-dependency-method"
    "subproject" "cmake" "pkgconf" "fetch"
  )
endif()

# Compiler-specific configurations

# gfortran "Interface mismatch" bug (119928)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "15.0"
     AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "15.2")
    message(FATAL_ERROR
      "Unsupported GNU Fortran version ${CMAKE_Fortran_COMPILER_VERSION}. "
      "This release has a known issue (fixed in 15.2). "
      "Please use GCC >= 15.2 or < 15.0, or a different Fortran compiler.")
  endif()
endif()

if(
  "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI"
  OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NVHPC"
  OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Flang"
)
  set(
    CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -Mbackslash -Mallocatable=03"
  )
endif()

# Set build type as CMake does not provide defaults
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(
    CMAKE_BUILD_TYPE "Release"
    CACHE STRING "Build type to be used."
    FORCE
  )
  message(
    STATUS
    "Setting build type to \"${CMAKE_BUILD_TYPE}\" as none was specified."
  )
endif()

# Dependency discovery
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/config/cmake")

if(NOT TARGET "OpenMP::OpenMP_Fortran" AND DFTD4_WITH_OpenMP)
  find_package("OpenMP" REQUIRED)
endif()

if(DFTD4_WITH_ILP64)
  message(STATUS "Using LAPACK/BLAS ILP64 interface")
  set(BLA_SIZEOF_INTEGER 8)
endif()

if(NOT TARGET "LAPACK::LAPACK")
  if(BLA_VENDOR MATCHES "^Intel" OR DEFINED ENV{MKLROOT})
    enable_language(C)
  endif()
  find_package(LAPACK REQUIRED)
endif()

# Collect subprojects
if(NOT TARGET "mctc-lib::mctc-lib")
  find_package("mctc-lib")
endif()
if(NOT TARGET "mstore::mstore")
  find_package("mstore")
endif()
if(NOT TARGET "multicharge::multicharge")
  find_package("multicharge")
endif()
set(
  lib-deps
  "mctc-lib::mctc-lib"
  "multicharge::multicharge"
)

# Collect source of the project
set(srcs)
add_subdirectory("src")

# DFT-D4 library target
add_library(
  "${PROJECT_NAME}"
  "${srcs}"
)
add_library("${PROJECT_NAME}::${PROJECT_NAME}" ALIAS "${PROJECT_NAME}")
set_target_properties(
  "${PROJECT_NAME}"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
target_link_libraries(
  "${PROJECT_NAME}"
  PUBLIC
  "${lib-deps}"
  "$<$<BOOL:${DFTD4_WITH_OpenMP}>:OpenMP::OpenMP_Fortran>"
)
target_include_directories(
  "${PROJECT_NAME}"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${DFTD4_INSTALL_FMODDIR}>
)
target_include_directories(
  "${PROJECT_NAME}"
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dftd4
)

if(DFTD4_WITH_ILP64)
  target_compile_definitions(
    "${PROJECT_NAME}"
    PRIVATE IK=i8)
else()
  target_compile_definitions(
    "${PROJECT_NAME}"
    PRIVATE IK=i4)
endif()

if(DFTD4_WITH_API)
  enable_language("C")
  if(DFTD4_WITH_PYTHON)
    add_subdirectory("python")
  endif()
endif()

# Add example application
add_subdirectory("app")

install(
  TARGETS
  "${PROJECT_NAME}"
  EXPORT
  "${PROJECT_NAME}-targets"
)
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${DFTD4_INSTALL_FMODDIR}"
)
if(DFTD4_WITH_API)
  install(
    DIRECTORY
    "${CMAKE_CURRENT_SOURCE_DIR}/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )
endif()
# Package license files
install(
  FILES
  "COPYING"
  "COPYING.LESSER"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)

# Package configuration and pkg-config metadata
install(
  DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/config/cmake/"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

set(DFTD4_USE_MULTICHARGE TRUE)
set(DFTD4_USE_MCTCLIB TRUE)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config/template-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMinorVersion
)
install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

if(IS_ABSOLUTE "${DFTD4_INSTALL_FMODDIR}")
  set(pkgconfig_fmoddir "${DFTD4_INSTALL_FMODDIR}")
else()
  set(pkgconfig_fmoddir "\${prefix}/${DFTD4_INSTALL_FMODDIR}")
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config/template.pc"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
  @ONLY
)
install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
)

# add the testsuite
if(DFTD4_WITH_TESTS)
  enable_testing()
  add_subdirectory("test")
endif()
