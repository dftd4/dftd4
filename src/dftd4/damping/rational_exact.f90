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

!> Implementation of the rational (Becke--Johnson) damping function.
module dftd4_damping_rational_exact
   use dftd4_damping, only : damping_param
   use dftd4_damping_mbd_exact, only : get_exact_atm_dispersion
   use dftd4_data, only : get_r4r2_val
   use dftd4_param, only : d4_param
   use dftd4_utils, only : triple_scale
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none
   private

   public :: new_rational_exact_damping, rational_damping_exact_param


   !> Rational (Becke-Johnson) damping model
   type, extends(rational_damping_param) :: rational_damping_exact_param
   contains

      !> Evaluate ATM three-body dispersion energy expression
      procedure :: get_dispersion3

      !> Evaluate pairwise representation of non-additive dispersion energy
      procedure :: get_pairwise_dispersion3

   end type rational_damping_exact_param


contains


!> Create new rational damping model
subroutine new_rational_exact_damping(self, param)

   !> Rational damping parameters
   type(rational_damping_exact_param), intent(out) :: self

   !> Parameters
   type(d4_param), intent(in) :: param

   self%s6 = param%s6
   self%s8 = param%s8
   self%s9 = param%s9
   self%a1 = param%a1
   self%a2 = param%a2
   self%alp = param%alp

end subroutine new_rational_exact_damping


!> Evaluation of the dispersion energy expression
subroutine get_dispersion3(self, mol, trans, cutoff, r4r2, c6, dc6dcn, dc6dq, &
      & energy, dEdcn, dEdq, gradient, sigma)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_dispersion3

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in), optional :: dc6dcn(:, :)

   !> Derivative of the C6 w.r.t. the partial charges
   real(wp), intent(in), optional :: dc6dq(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)

   !> Derivative of the energy w.r.t. the partial charges
   real(wp), intent(inout), optional :: dEdq(:)

   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   call get_exact_atm_dispersion(mol, trans, cutoff, self%s9, self%a1, &
      & self%a2, self%alp, r4r2, c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, &
      & gradient, sigma)

end subroutine get_dispersion3


!> Evaluation of the dispersion energy expression
subroutine get_pairwise_dispersion3(self, mol, trans, cutoff, r4r2, c6, energy)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_pairwise_dispersion3

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:, :)

   call get_exact_atm_dispersion(mol, trans, cutoff, self%s9, self%a1, &
   & self%a2, self%alp, r4r2, c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, &
   & gradient, sigma)

end subroutine get_pairwise_dispersion3


end module dftd4_damping_rational_exact
