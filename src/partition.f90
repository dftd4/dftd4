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
   use mctc_env, only : error_type, fatal_error, i8
   implicit none
   private

   public :: work_partition, new_work_partition, serial_work_partition
   public :: owns_pair


   !> Cyclic partition of the work of a dispersion calculation.
   !>
   !> Parts are zero based. Every unit of work is assigned to exactly one part,
   !> summing the energy and derivative contributions of all parts reproduces the
   !> complete result. An absent partition owns all of the work.
   type :: work_partition
      private

      !> Zero-based index of this part
      integer :: part = 0

      !> Total number of parts
      integer :: nparts = 1
   end type work_partition

   !> Complete work of an ordinary serial calculation, equivalent to omitting
   !> the partition entirely
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


!> Whether this part owns a symmetry-reduced atom pair
elemental function owns_pair(partition, iat, jat) result(owned)

   !> Work partition, absent selects the complete work
   type(work_partition), intent(in), optional :: partition

   !> Atom indices of the pair, with jat <= iat
   integer, intent(in) :: iat, jat

   !> Whether this part owns the pair
   logical :: owned

   integer(i8) :: pair_index

   owned = .true.
   if (.not.present(partition)) return
   if (partition%nparts == 1) return

   ! zero-based index in the lower-triangular sequence (1,1), (2,1), (2,2), ...
   pair_index = int(iat - 1, i8)*int(iat, i8)/2_i8 + int(jat - 1, i8)
   owned = modulo(pair_index, int(partition%nparts, i8)) == int(partition%part, i8)

end function owns_pair


end module dftd4_partition
