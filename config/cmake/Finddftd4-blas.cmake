# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

if ((DFTD4_BLAS MATCHES ^Intel) OR (DEFINED ENV{MKLROOT}))
  enable_language("C")
endif()

if(WITH_ILP64)
  message(STATUS "Using LAPACK/BLAS ILP64 interface")
  set(BLA_SIZEOF_INTEGER 8)
  set(_nvpl_int "_ilp64")
else()
  set(BLA_SIZEOF_INTEGER 4)
  set(_nvpl_int "_lp64")
endif()

if(NOT DFTD4_BLAS_FOUND)
  if(DFTD4_BLAS STREQUAL "NVPL")
    find_package("nvpl_blas" REQUIRED)
    set(DFTD4_BLAS_FOUND TRUE)

    if((BLA_THREAD STREQUAL "OMP") OR (BLA_THREAD STREQUAL "ANY"))
      set(_nvpl_thread "_omp")
    else()
      set(_nvpl_thread "_seq")
    endif()

    add_library("dftd4::BLAS" INTERFACE IMPORTED GLOBAL)
    target_link_libraries("dftd4::BLAS" INTERFACE "nvpl::blas${_nvpl_int}${_nvpl_thread}")
  else()
    set(BLA_VENDOR "${DFTD4_BLAS}")
    find_package("BLAS" REQUIRED)
    set(DFTD4_BLAS_FOUND ${BLAS_FOUND})
  
    if(NOT TARGET "dftd4::BLAS")
      add_library("dftd4::BLAS" INTERFACE IMPORTED GLOBAL)
      target_link_libraries("dftd4::BLAS" INTERFACE "BLAS::BLAS")
    endif()
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(dftd4-blas DEFAULT_MSG DFTD4_BLAS_FOUND)
