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

if(NOT BLAS_FOUND)
  if("${BLA_VENDOR}" MATCHES "^Intel" OR DEFINED ENV{MKLROOT})
    # C must be enabled to use MKL
    # https://cmake.org/cmake/help/v3.14/module/FindBLAS.html#result-variables
    enable_language("C")
  endif()

  if(BLA_VENDOR STREQUAL "NVPL")
    find_package("nvpl_blas" REQUIRED)
    if(BLA_SIZEOF_INTEGER EQUAL 8)
      set(_nvpl_int "_ilp64")
    else()
      set(_nvpl_int "_lp64")
    endif()

    if((BLA_THREAD STREQUAL "OMP") OR (BLA_THREAD STREQUAL "ANY"))
      set(_nvpl_thread "_omp")
    else()
      set(_nvpl_thread "_seq")
    endif()

    if(NOT TARGET "DFTD4::BLAS")
      add_library("DFTD4::BLAS" INTERFACE IMPORTED)
      target_link_libraries("DFTD4::BLAS" INTERFACE "nvpl::blas${_nvpl_int}${_nvpl_thread}")
    endif()
  else()
    find_package("BLAS" REQUIRED)
    if(NOT TARGET "DFTD4::BLAS")
      add_library("DFTD4::BLAS" INTERFACE IMPORTED)
      if(TARGET "BLAS::BLAS")
        target_link_libraries("DFTD4::BLAS" INTERFACE "BLAS::BLAS")
      else()
        target_link_libraries("DFTD4::BLAS" INTERFACE "${BLAS_LIBRARIES}")
      endif()
    endif()
  endif()
 endif()
