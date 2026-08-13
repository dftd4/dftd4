! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Work partitioning for externally distributed dispersion calculations
module dftd4_partition
   use, intrinsic :: iso_fortran_env, only : int64
   use mctc_env, only : error_type, fatal_error
   implicit none
   private

   public :: work_partition, new_work_partition, serial_work_partition


   !> Cyclic partition of the symmetry-reduced atom-pair work.
   !>
   !> Parts are zero based.  Every pair `(iat, jat)`, with `jat <= iat`, is
   !> assigned to exactly one part.  The caller is responsible for summing the
   !> energy and derivative contributions returned by all parts.
   type :: work_partition
      private

      !> Zero-based index of this part
      integer :: part = 0

      !> Total number of parts
      integer :: nparts = 1
   contains
      !> Whether this is the first work partition
      procedure, public :: is_first

      !> Whether this part owns a symmetry-reduced atom pair
      procedure, public :: owns_pair
   end type work_partition

   !> Work partition representing an ordinary serial calculation
   type(work_partition), parameter :: serial_work_partition = work_partition()


contains


!> Create a work partition
subroutine new_work_partition(error, partition, part, nparts)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> New work partition
   type(work_partition), intent(out) :: partition

   !> Zero-based index of this part
   integer, intent(in) :: part

   !> Total number of parts
   integer, intent(in) :: nparts

   if (nparts <= 0 .or. part < 0 .or. part >= nparts) then
      call fatal_error(error, "Invalid dispersion work partition")
      return
   end if

   partition%part = part
   partition%nparts = nparts

end subroutine new_work_partition


!> Whether this is the first work partition
elemental function is_first(self) result(first)
   class(work_partition), intent(in) :: self
   logical :: first

   first = self%part == 0

end function is_first


!> Whether this part owns a symmetry-reduced atom pair
elemental function owns_pair(self, iat, jat) result(owned)
   class(work_partition), intent(in) :: self
   integer, intent(in) :: iat
   integer, intent(in) :: jat
   logical :: owned

   integer(int64) :: pair_index

   if (self%nparts == 1) then
      owned = .true.
      return
   end if

   ! Zero-based index in the lower-triangular atom-pair sequence:
   ! (1,1), (2,1), (2,2), (3,1), ...
   pair_index = int(iat - 1, int64)*int(iat, int64)/2_int64 + int(jat - 1, int64)
   owned = modulo(pair_index, int(self%nparts, int64)) == int(self%part, int64)

end function owns_pair


end module dftd4_partition
