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

if ((BLA_VENDOR MATCHES ^Intel) OR (DEFINED ENV{MKLROOT}))
  enable_language("C")
endif()

if(WITH_ILP64)
  set(BLA_SIZEOF_INTEGER 8)
endif()

if(NOT LAPACK_FOUND)
  find_package("LAPACK")

  if(NOT TARGET "BLAS::BLAS")
    find_package("custom-blas")
  endif()

  if(NOT TARGET "LAPACK::LAPACK")
    add_library("LAPACK::LAPACK" INTERFACE IMPORTED)
    target_link_libraries("LAPACK::LAPACK" INTERFACE "${LAPACK_LIBRARIES}" "BLAS::BLAS")
  endif()
endif()
