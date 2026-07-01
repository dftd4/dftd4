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

!> Generic interface to define damping functions for the DFT-D4 model
module dftd4_damping
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none

   public :: damping_param, dispersion_interface


   type, abstract :: damping_param
   contains
      generic :: get_dispersion2 => get_dispersion2_impl, get_dispersion2_compat
      procedure(dispersion_interface), deferred :: get_dispersion2_impl
      procedure :: get_dispersion2_compat
      generic :: get_dispersion3 => get_dispersion3_impl, get_dispersion3_compat
      procedure(dispersion_interface), deferred :: get_dispersion3_impl
      procedure :: get_dispersion3_compat
      generic :: get_pairwise_dispersion2 => get_pairwise_dispersion2_impl, get_pairwise_dispersion2_compat
      procedure(pairwise_dispersion_interface), deferred :: get_pairwise_dispersion2_impl
      procedure :: get_pairwise_dispersion2_compat
      generic :: get_pairwise_dispersion3 => get_pairwise_dispersion3_impl, get_pairwise_dispersion3_compat
      procedure(pairwise_dispersion_interface), deferred :: get_pairwise_dispersion3_impl
      procedure :: get_pairwise_dispersion3_compat
   end type damping_param


   abstract interface
      !> Evaluation of the dispersion energy expression
      subroutine dispersion_interface(self, mol, trans, cutoff, width, r4r2, &
            & c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient, sigma)
         import :: structure_type, damping_param, wp

         !> Damping parameters
         class(damping_param), intent(in) :: self

         !> Molecular structure data
         class(structure_type), intent(in) :: mol

         !> Lattice points
         real(wp), intent(in) :: trans(:, :)

         !> Real space cutoff
         real(wp), intent(in) :: cutoff

         !> Width of smooth cutoff
         real(wp), intent(in) :: width

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
      end subroutine dispersion_interface

      !> Evaluation of the pairwise representation of the dispersion energy
      subroutine pairwise_dispersion_interface(self, mol, trans, cutoff, width, r4r2, c6, energy)
         import :: structure_type, damping_param, wp

         !> Damping parameters
         class(damping_param), intent(in) :: self

         !> Molecular structure data
         class(structure_type), intent(in) :: mol

         !> Lattice points
         real(wp), intent(in) :: trans(:, :)

         !> Real space cutoff
         real(wp), intent(in) :: cutoff

         !> Width of smooth cutoff
         real(wp), intent(in) :: width

         !> Expectation values for r4 over r2 operator
         real(wp), intent(in) :: r4r2(:)

         !> C6 coefficients for all atom pairs.
         real(wp), intent(in) :: c6(:, :)

         !> Pairwise representation of the dispersion energy
         real(wp), intent(inout) :: energy(:, :)
      end subroutine pairwise_dispersion_interface
   end interface

contains

!> Evaluation of the dispersion energy expression
subroutine get_dispersion2_compat(self, mol, trans, cutoff, r4r2, &
      & c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient, sigma)

   !> Damping parameters
   class(damping_param), intent(in) :: self

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

   call self%get_dispersion2(mol, trans, cutoff, 0.0_wp, r4r2, c6, &
      & dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient, sigma)
end subroutine get_dispersion2_compat

!> Evaluation of the dispersion energy expression
subroutine get_dispersion3_compat(self, mol, trans, cutoff, r4r2, &
      & c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient, sigma)

   !> Damping parameters
   class(damping_param), intent(in) :: self

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

   call self%get_dispersion3(mol, trans, cutoff, 0.0_wp, r4r2, c6, &
      & dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient, sigma)
end subroutine get_dispersion3_compat

!> Evaluation of the pairwise representation of the dispersion energy
subroutine get_pairwise_dispersion2_compat(self, mol, trans, cutoff, r4r2, c6, energy)

   !> Damping parameters
   class(damping_param), intent(in) :: self

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

   !> Pairwise representation of the dispersion energy
   real(wp), intent(inout) :: energy(:, :)

   call self%get_pairwise_dispersion2(mol, trans, cutoff, 0.0_wp, r4r2, c6, energy)
end subroutine get_pairwise_dispersion2_compat

!> Evaluation of the pairwise representation of the dispersion energy
subroutine get_pairwise_dispersion3_compat(self, mol, trans, cutoff, r4r2, c6, energy)

   !> Damping parameters
   class(damping_param), intent(in) :: self

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

   !> Pairwise representation of the dispersion energy
   real(wp), intent(inout) :: energy(:, :)

   call self%get_pairwise_dispersion3(mol, trans, cutoff, 0.0_wp, r4r2, c6, energy)
end subroutine get_pairwise_dispersion3_compat

end module dftd4_damping
