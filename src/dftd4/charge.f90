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

!> Interface to the charge model
module dftd4_charge
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use multicharge, only : mchrg_model_type, new_eeq2019_model
   implicit none
   private

   public :: get_charges


contains


!> Obtain charges from electronegativity equilibration model
subroutine get_charges(mol, qvec, dqdr, dqdL)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_charges

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: qvec(:)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain deformations
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   logical :: grad
   type(mchrg_model_type) :: model
   type(error_type), allocatable :: error
   real(wp), parameter :: cn_max = 8.0_wp, cutoff = 25.0_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: rcov(:), trans(:, :)

   grad = present(dqdr) .and. present(dqdL)

   call new_eeq2019_model(mol, model, error)
   if(allocated(error)) then
      error stop "Error occurd in the coordination number setup"
   end if

   allocate(cn(mol%nat))
   if (grad) then
      allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
   end if

   call model%ncoord%get_cn(mol, cn, dcndr, dcndL)

   call model%solve(mol, cn, dcndr, dcndL, qvec=qvec, dqdr=dqdr, dqdL=dqdL)

end subroutine get_charges


end module dftd4_charge
