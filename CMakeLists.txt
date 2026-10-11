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

cmake_minimum_required(VERSION 3.24)

project(
  dftd4
  LANGUAGES Fortran
  VERSION 4.3.0
  DESCRIPTION
    "Generally Applicable Atomic-Charge Dependent London Dispersion Correction")

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(DFTD4_WITH_API "Build the C API" ON)
option(DFTD4_WITH_OpenMP "Enable OpenMP" ON)
option(DFTD4_WITH_ILP64 "Use the ILP64 BLAS/LAPACK interface" OFF)
option(DFTD4_WITH_API_V2 "Build the compatibility layer for the 2.5.x API" OFF)
option(DFTD4_WITH_PYTHON "Build the Python extension" OFF)

set(DFTD4_INSTALL_FMODDIR
    "${CMAKE_INSTALL_INCLUDEDIR}/dftd4/${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}"
    CACHE STRING "Installation directory for Fortran modules")

if(DFTD4_WITH_PYTHON AND NOT DFTD4_WITH_API)
  message(FATAL_ERROR "DFTD4_WITH_PYTHON requires DFTD4_WITH_API")
endif()

# GCC 15.0/15.1 generate an incompatible Fortran interface (GCC PR 119928).
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"
   AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "15.0"
   AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "15.2")
  message(
    FATAL_ERROR
      "GNU Fortran ${CMAKE_Fortran_COMPILER_VERSION} is unsupported (GCC PR 119928); use GCC 15.2 or newer"
  )
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(PGI|NVHPC|Flang)$")
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-Mbackslash>"
                      "$<$<COMPILE_LANGUAGE:Fortran>:-Mallocatable=03>")
endif()

if(PROJECT_IS_TOP_LEVEL
   AND NOT CMAKE_BUILD_TYPE
   AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      RelWithDebInfo
      CACHE STRING "Build type" FORCE)
endif()

# FindBLAS/FindLAPACK require a C compiler when using Intel MKL.
if(DFTD4_WITH_API
   OR BLA_VENDOR MATCHES "^Intel"
   OR DEFINED ENV{MKLROOT})
  enable_language(C)
endif()

if(DFTD4_WITH_OpenMP AND NOT TARGET OpenMP::OpenMP_Fortran)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

if(DFTD4_WITH_ILP64)
  set(BLA_SIZEOF_INTEGER 8)
else()
  set(BLA_SIZEOF_INTEGER 4)
endif()
if(NOT TARGET BLAS::BLAS)
  find_package(BLAS REQUIRED)
endif()
if(NOT TARGET LAPACK::LAPACK)
  find_package(LAPACK REQUIRED)
endif()

# Prefer installed packages, then bundled sources, and finally upstream releases.
# Sources under subprojects/ are included in Meson distribution archives.
function(dftd4_find_dependency package git_ref)
  # Forward settings understood by the pinned subproject releases.
  if(package STREQUAL "mctc-lib" OR package STREQUAL "multicharge")
    set(WITH_OpenMP "${DFTD4_WITH_OpenMP}")
  endif()
  if(package STREQUAL "multicharge")
    set(WITH_ILP64 "${DFTD4_WITH_ILP64}")
    set(WITH_TESTS OFF)
  endif()

  include(FetchContent)
  if(EXISTS "${PROJECT_SOURCE_DIR}/subprojects/${package}/CMakeLists.txt")
    FetchContent_Declare(
      "${package}" SOURCE_DIR "${PROJECT_SOURCE_DIR}/subprojects/${package}"
                   FIND_PACKAGE_ARGS CONFIG)
  else()
    FetchContent_Declare(
      "${package}"
      GIT_REPOSITORY "https://github.com/grimme-lab/${package}.git"
      GIT_TAG "${git_ref}"
      FIND_PACKAGE_ARGS CONFIG)
  endif()
  FetchContent_MakeAvailable("${package}")

  # Older dependency versions may only provide an unnamespaced target.
  if(NOT TARGET "${package}::${package}" AND TARGET "${package}")
    add_library("${package}::${package}" ALIAS "${package}")
  endif()
  if(NOT TARGET "${package}::${package}")
    message(FATAL_ERROR "${package} did not provide ${package}::${package}")
  endif()
endfunction()

dftd4_find_dependency(mctc-lib v0.6.1)
dftd4_find_dependency(multicharge d4888df)

set(srcs)
add_subdirectory(src)
add_library(dftd4 ${srcs})
add_library(dftd4::dftd4 ALIAS dftd4)
set_target_properties(
  dftd4
  PROPERTIES VERSION "${PROJECT_VERSION}"
             SOVERSION "${PROJECT_VERSION_MAJOR}"
             POSITION_INDEPENDENT_CODE ON
             Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
target_include_directories(
  dftd4
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
         "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>"
         "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
         "$<INSTALL_INTERFACE:${DFTD4_INSTALL_FMODDIR}>")
target_include_directories(dftd4 PRIVATE "${PROJECT_SOURCE_DIR}/src/dftd4")
target_link_libraries(
  dftd4
  PUBLIC mctc-lib::mctc-lib multicharge::multicharge
  PRIVATE BLAS::BLAS LAPACK::LAPACK)
if(DFTD4_WITH_OpenMP)
  target_link_libraries(dftd4 PRIVATE OpenMP::OpenMP_Fortran)
endif()
if(DFTD4_WITH_ILP64)
  target_compile_definitions(dftd4 PRIVATE IK=i8)
else()
  target_compile_definitions(dftd4 PRIVATE IK=i4)
endif()

add_subdirectory(app)
if(DFTD4_WITH_PYTHON)
  add_subdirectory(python)
endif()

if(BUILD_TESTING)
  dftd4_find_dependency(mstore v0.4.0)
  add_subdirectory(test)
endif()

install(TARGETS dftd4 EXPORT dftd4-targets)
install(
  EXPORT dftd4-targets
  NAMESPACE dftd4::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/dftd4")
install(DIRECTORY "${PROJECT_BINARY_DIR}/include/"
        DESTINATION "${DFTD4_INSTALL_FMODDIR}")
if(DFTD4_WITH_API)
  install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/"
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
endif()
install(FILES COPYING COPYING.LESSER
        DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/dftd4")

set(DFTD4_USE_MCTCLIB ON)
set(DFTD4_USE_MULTICHARGE ON)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/config/template-config.cmake"
  "${PROJECT_BINARY_DIR}/dftd4-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/dftd4")
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/dftd4-config-version.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMinorVersion)
install(FILES "${PROJECT_BINARY_DIR}/dftd4-config.cmake"
              "${PROJECT_BINARY_DIR}/dftd4-config-version.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/dftd4")

if(IS_ABSOLUTE "${DFTD4_INSTALL_FMODDIR}")
  set(pkgconfig_fmoddir "${DFTD4_INSTALL_FMODDIR}")
else()
  set(pkgconfig_fmoddir "\${prefix}/${DFTD4_INSTALL_FMODDIR}")
endif()
configure_file("${PROJECT_SOURCE_DIR}/config/template.pc"
               "${PROJECT_BINARY_DIR}/dftd4.pc" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/dftd4.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
