! This file is part of dftd4.
!
! Copyright (C) 2017-2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module class_results
   use iso_fortran_env, wp => real64
   implicit none

   public :: dftd_results
   private

   type :: dftd_results
      real(wp), allocatable :: energy
      real(wp), allocatable :: gradient(:,:)
      real(wp), allocatable :: stress(:,:)
      real(wp), allocatable :: lattice_gradient(:,:)
      real(wp), allocatable :: charges(:)
      real(wp), allocatable :: dipole_moment(:)
      real(wp), allocatable :: polarizibilities(:)
      real(wp), allocatable :: c6_coefficients(:,:)
      real(wp), allocatable :: hessian(:,:)
   contains
      procedure :: allocate => allocate_results
      procedure :: deallocate => deallocate_results
   end type dftd_results

contains

subroutine allocate_results(self,n)
   class(dftd_results), intent(inout) :: self
   integer, intent(in) :: n
   call self%deallocate
   self%energy = 0.0_wp
   allocate(self%gradient(3,n),         source = 0.0_wp)
   allocate(self%stress(3,3),           source = 0.0_wp)
   allocate(self%lattice_gradient(3,3), source = 0.0_wp)
   allocate(self%charges(n),            source = 0.0_wp)
   allocate(self%dipole_moment(3),      source = 0.0_wp)
   allocate(self%polarizibilities(n),   source = 0.0_wp)
   allocate(self%c6_coefficients(n,n),  source = 0.0_wp)
   allocate(self%hessian(3*n,3*n),      source = 0.0_wp)
end subroutine allocate_results

subroutine deallocate_results(self)
   class(dftd_results), intent(inout) :: self
   if (allocated(self%energy))           deallocate(self%energy)
   if (allocated(self%stress))           deallocate(self%stress)
   if (allocated(self%lattice_gradient)) deallocate(self%lattice_gradient)
   if (allocated(self%gradient))         deallocate(self%gradient)
   if (allocated(self%charges))          deallocate(self%charges)
   if (allocated(self%dipole_moment))    deallocate(self%dipole_moment)
   if (allocated(self%polarizibilities)) deallocate(self%polarizibilities)
   if (allocated(self%c6_coefficients))  deallocate(self%c6_coefficients)
   if (allocated(self%hessian))          deallocate(self%hessian)
end subroutine deallocate_results

end module class_results
