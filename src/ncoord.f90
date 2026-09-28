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

module dftd4_ncoord
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count
   implicit none
   private

   public :: get_coordination_number, add_coordination_number_derivs
   public :: add_coordination_number_hessian


   !> Steepness of counting function
   real(wp), parameter :: default_kcn = 7.5_wp


contains


!> Geometric fractional coordination number, supports error function counting.
subroutine get_coordination_number(mol, trans, cutoff, rcov, en, cn, dcndr, dcndL, error)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_coordination_number

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Electronegativity
   real(wp), intent(in) :: en(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out), optional :: dcndr(:, :, :)

   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out), optional :: dcndL(:, :, :)

   !> Error on failure; results are NaN if initialization fails.
   type(error_type), allocatable, intent(out), optional :: error

   class(ncoord_type), allocatable :: ncoord
   type(error_type), allocatable :: local_error

   call new_ncoord(ncoord, mol, cn_count%dftd4, &
      & kcn=default_kcn, cutoff=cutoff, rcov=rcov, en=en, error=local_error)
   if (allocated(local_error)) then
      cn = ieee_value(0.0_wp, ieee_quiet_nan)
      if (present(dcndr)) dcndr = ieee_value(0.0_wp, ieee_quiet_nan)
      if (present(dcndL)) dcndL = ieee_value(0.0_wp, ieee_quiet_nan)
      if (present(error)) call move_alloc(local_error, error)
      return
   end if

   call ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)

end subroutine get_coordination_number


subroutine add_coordination_number_derivs(mol, trans, cutoff, rcov, en, dEdcn, gradient, sigma, error)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Electronegativity
   real(wp), intent(in) :: en(:)

   !> Derivative of expression with respect to the coordination number
   real(wp), intent(in) :: dEdcn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates
   real(wp), intent(inout) :: gradient(:, :)

   !> Derivative of the CN with respect to strain deformations
   real(wp), intent(inout) :: sigma(:, :)


   !> Error on failure; results are NaN if initialization fails.
   type(error_type), allocatable, intent(out), optional :: error

   class(ncoord_type), allocatable :: ncoord
   type(error_type), allocatable :: local_error

   call new_ncoord(ncoord, mol, cn_count%dftd4, &
      & kcn=default_kcn, cutoff=cutoff, rcov=rcov, en=en, error=local_error)
   if (allocated(local_error)) then
      gradient = ieee_value(0.0_wp, ieee_quiet_nan)
      sigma = ieee_value(0.0_wp, ieee_quiet_nan)
      if (present(error)) call move_alloc(local_error, error)
      return
   end if

   call ncoord%add_coordination_number_derivs(mol, trans, dEdcn, gradient, sigma)

end subroutine add_coordination_number_derivs


!> Add the second derivative of the D4 coordination number contracted with
!> the derivative of the energy w.r.t. the coordination number.
subroutine add_coordination_number_hessian(mol, trans, cutoff, rcov, en, dEdcn, hessian, error)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Electronegativity
   real(wp), intent(in) :: en(:)

   !> Derivative of expression with respect to the coordination number
   real(wp), intent(in) :: dEdcn(:)

   !> Second derivative of the energy w.r.t. the Cartesian coordinates
   real(wp), intent(inout) :: hessian(:, :)

   !> Error on failure; results are NaN if initialization fails.
   type(error_type), allocatable, intent(out), optional :: error

   class(ncoord_type), allocatable :: ncoord
   type(error_type), allocatable :: local_error

   call new_ncoord(ncoord, mol, cn_count%dftd4, &
      & kcn=default_kcn, cutoff=cutoff, rcov=rcov, en=en, error=local_error)
   if (allocated(local_error)) then
      hessian = ieee_value(0.0_wp, ieee_quiet_nan)
      if (present(error)) call move_alloc(local_error, error)
      return
   end if

   call ncoord%add_coordination_number_hessian(mol, trans, dEdcn, hessian)

end subroutine add_coordination_number_hessian


end module dftd4_ncoord
