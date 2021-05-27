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

# Include the mctc-lib project
if(NOT TARGET mctc-lib)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mctc-lib/CMakeLists.txt)
    add_subdirectory("mctc-lib")
  else()
    set("mctc-lib-url" "https://github.com/grimme-lab/mctc-lib")
    message(STATUS "Retrieving mctc-lib from ${mctc-lib-url}")
    include(FetchContent)
    FetchContent_Declare(
      "mctc-lib"
      GIT_REPOSITORY "${mctc-lib-url}"
      GIT_TAG "HEAD"
    )
    FetchContent_MakeAvailable("mctc-lib")
  endif()
endif()

list(
  APPEND lib-deps
  "mctc-lib"
)

# Include the mstore project
if(NOT TARGET mstore)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mstore/CMakeLists.txt)
    add_subdirectory("mstore")
  else()
    set("mstore-url" "https://github.com/grimme-lab/mstore")
    message(STATUS "Retrieving mstore from ${mstore-url}")
    include(FetchContent)
    FetchContent_Declare(
      "mstore"
      GIT_REPOSITORY "${mstore-url}"
      GIT_TAG "HEAD"
    )
    FetchContent_MakeAvailable("mstore")
  endif()
endif()

# Include the multicharge project
if(NOT TARGET multicharge)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/multicharge/CMakeLists.txt)
    add_subdirectory("multicharge")
  else()
    set("multicharge-url" "https://github.com/grimme-lab/multicharge")
    message(STATUS "Retrieving multicharge from ${multicharge-url}")
    include(FetchContent)
    FetchContent_Declare(
      "multicharge"
      GIT_REPOSITORY "${multicharge-url}"
      GIT_TAG "HEAD"
    )
    FetchContent_MakeAvailable("multicharge")
  endif()
endif()

list(
  APPEND lib-deps
  "multicharge"
)

set(lib-deps "${lib-deps}" PARENT_SCOPE)
