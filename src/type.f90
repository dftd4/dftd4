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

!> Definition of the abstract base dispersion model for the evaluation of C6 coefficients.
module dftd4_model_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use multicharge, only : mchrg_model_type
   implicit none
   private

   public :: dispersion_model, d4_qmod


   !> Abstract base dispersion model to evaluate C6 coefficients
   type, abstract :: dispersion_model

      !> Number of atoms coupled to by pairwise parameters
      integer :: ncoup

      !> Charge scaling height
      real(wp) :: ga

      !> Charge scaling steepness
      real(wp) :: gc

      !> Effective nuclear charges
      real(wp), allocatable :: zeff(:)

      !> Chemical hardness
      real(wp), allocatable :: eta(:)

      !> Electronegativity
      real(wp), allocatable :: en(:)

      !> Covalent radii for coordination number
      real(wp), allocatable :: rcov(:)

      !> Expectation values for C8 extrapolation
      real(wp), allocatable :: r4r2(:)

      !> Number of reference systems
      integer, allocatable :: ref(:)

      !> Number of Gaussian weights for each reference
      integer, allocatable :: ngw(:, :)

      !> Reference coordination numbers
      real(wp), allocatable :: cn(:, :)

      !> Reference partial charges
      real(wp), allocatable :: q(:, :)

      !> Reference dynamic polarizabilities
      real(wp), allocatable :: aiw(:, :, :)

      !> Reference C6 coefficients
      real(wp), allocatable :: c6(:, :, :, :)

      !> Multicharge model
      class(mchrg_model_type), allocatable :: mchrg 

   contains

      !> Generate weights for all reference systems
      procedure(weight_references), deferred :: weight_references

      !> Evaluate C6 coefficient
      procedure(get_atomic_c6), deferred :: get_atomic_c6

      !> Evaluate atomic polarizabilities
      procedure(get_polarizabilities), deferred :: get_polarizabilities

   end type dispersion_model

   abstract interface

      !> Calculate the weights of the reference system and the derivatives w.r.t.
      !> coordination number for later use.
      subroutine weight_references(self, mol, cn, q, gwvec, gwdcn, gwdq)
         import dispersion_model, structure_type, wp
         !> Instance of the dispersion model
         class(dispersion_model), intent(in) :: self
         !> Molecular structure data
         class(structure_type), intent(in) :: mol
         !> Coordination number of every atom: [nat]
         real(wp), intent(in) :: cn(:)
         !> Partial charge of every atom: [nat]
         real(wp), intent(in) :: q(:)
         !> weighting for the atomic reference systems: [nref, nat, ncoup]
         real(wp), intent(out) :: gwvec(:, :, :)
         !> derivative of the weighting function w.r.t. the coordination number: [nref, nat, ncoup]
         real(wp), intent(out), optional :: gwdcn(:, :, :)
         !> derivative of the weighting function w.r.t. the charge scaling: [nref, nat, ncoup]
         real(wp), intent(out), optional :: gwdq(:, :, :)
      end subroutine 

      !> Calculate atomic dispersion coefficients and their derivatives w.r.t.
      !> the coordination numbers and atomic partial charges.
      subroutine get_atomic_c6(self, mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)
         import dispersion_model, structure_type, wp
         !> Instance of the dispersion model
         class(dispersion_model), intent(in) :: self
         !> Molecular structure data
         class(structure_type), intent(in) :: mol
         !> Weighting function for the atomic reference systems
         real(wp), intent(in) :: gwvec(:, :, :)
         !> Derivative of the weighting function w.r.t. the coordination number
         real(wp), intent(in), optional :: gwdcn(:, :, :)
         !> Derivative of the weighting function w.r.t. the partial charge
         real(wp), intent(in), optional :: gwdq(:, :, :)
         !> C6 coefficients for all atom pairs.
         real(wp), intent(out) :: c6(:, :)
         !> Derivative of the C6 w.r.t. the coordination number
         real(wp), intent(out), optional :: dc6dcn(:, :)
         !> Derivative of the C6 w.r.t. the partial charge
         real(wp), intent(out), optional :: dc6dq(:, :)
      end subroutine get_atomic_c6

      !> Calculate atomic polarizabilities and their derivatives w.r.t.
      !> the coordination numbers and atomic partial charges.
      subroutine get_polarizabilities(self, mol, gwvec, gwdcn, gwdq, alpha, dadcn, dadq)
         import dispersion_model, structure_type, wp
         !> Instance of the dispersion model
         class(dispersion_model), intent(in) :: self
         !> Molecular structure data
         class(structure_type), intent(in) :: mol
         !> Weighting function for the atomic reference systems
         real(wp), intent(in) :: gwvec(:, :, :)
         !> Derivative of the weighting function w.r.t. the coordination number
         real(wp), intent(in), optional :: gwdcn(:, :, :)
         !> Derivative of the weighting function w.r.t. the partial charge
         real(wp), intent(in), optional :: gwdq(:, :, :)
         !> Static polarizabilities for all atoms.
         real(wp), intent(out) :: alpha(:)
         !> Derivative of the polarizibility w.r.t. the coordination number
         real(wp), intent(out), optional :: dadcn(:)
         !> Derivative of the polarizibility w.r.t. the partial charge
         real(wp), intent(out), optional :: dadq(:)
      end subroutine

   end interface


   !> Possible reference charges for D4
   type :: enum_qmod

      !> Electronegativity equilibration charges
      integer :: eeq = 1

      !> GFN2-xTB Mulliken partial charges
      integer :: gfn2 = 2

      !> Bond-Capcity Electronegativity equilibration charges
      integer :: eeqbc = 3

   end type enum_qmod

   !> Actual enumerator for D4 reference charges
   type(enum_qmod), parameter :: d4_qmod = enum_qmod()
   !DEC$ ATTRIBUTES DLLEXPORT :: d4_qmod

end module dftd4_model_type
