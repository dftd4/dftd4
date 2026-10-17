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

      !> Evaluate C6 coefficients and their second CN/charge derivatives
      procedure :: get_atomic_c6_hessian

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
      end subroutine weight_references

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
         !> Derivative of the polarizability w.r.t. the coordination number
         real(wp), intent(out), optional :: dadcn(:)
         !> Derivative of the polarizability w.r.t. the partial charge
         real(wp), intent(out), optional :: dadq(:)
      end subroutine get_polarizabilities

   end interface


   !> Possible reference charges for D4
   type :: enum_qmod

      !> Electronegativity equilibration charges
      integer :: eeq = 1

      !> GFN2-xTB Mulliken partial charges
      integer :: gfn2 = 2

      !> Bond-Capacity Electronegativity equilibration charges
      integer :: eeqbc = 3

   end type enum_qmod

   !> Actual enumerator for D4 reference charges
   type(enum_qmod), parameter :: d4_qmod = enum_qmod()
   !DEC$ ATTRIBUTES DLLEXPORT :: d4_qmod

contains


!> Calculate C6 coefficients and their first and second CN/charge derivatives.
subroutine get_atomic_c6_hessian(self, mol, gwvec, gwdcn, gwd2cn, c6, &
      & dc6dcn, d2c6dcn2, d2c6dcnij, gwdq, gwd2q, gwdcnq, dc6dq, &
      & d2c6dq2, d2c6dqij, d2c6dcnq, d2c6dcnqij)

   !> Instance of the dispersion model
   class(dispersion_model), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Reference weights and their CN derivatives
   real(wp), intent(in) :: gwvec(:, :, :), gwdcn(:, :, :), gwd2cn(:, :, :)

   !> C6 coefficients and their CN derivatives
   real(wp), intent(out) :: c6(:, :), dc6dcn(:, :)
   real(wp), intent(out) :: d2c6dcn2(:, :), d2c6dcnij(:, :)

   !> Optional charge response: supply all eight arrays together.
   !> For pair (i,j), same-end derivatives act on i; mixed derivatives act on
   !> CN_i and q_j. Diagonal pairs still describe two independent endpoints.
   real(wp), intent(in), optional :: gwdq(:, :, :), gwd2q(:, :, :), gwdcnq(:, :, :)
   real(wp), intent(out), optional :: dc6dq(:, :), d2c6dq2(:, :), d2c6dqij(:, :)
   real(wp), intent(out), optional :: d2c6dcnq(:, :), d2c6dcnqij(:, :)

   logical :: charge
   integer :: icoup, jcoup
   real(wp) :: dq(2), dqq(2), dqqij, dcq(2), dcqij(2)
   integer :: ipair, npair, iat, jat, izp, jzp, iref, jref
   real(wp) :: refc6, c6ij, dc6dcni, dc6dcnj, d2c6dcni, d2c6dcnj, d2c6mix

   charge = present(gwdq) .and. present(gwd2q) .and. present(gwdcnq) &
      & .and. present(dc6dq) .and. present(d2c6dq2) .and. present(d2c6dqij) &
      & .and. present(d2c6dcnq) .and. present(d2c6dcnqij)

   ! Flatten the lower triangle so OpenMP can schedule O(N^2) pair tasks.
   npair = mol%nat*(mol%nat + 1)/2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(c6, dc6dcn, d2c6dcn2, d2c6dcnij, mol, self, gwvec, gwdcn, gwd2cn, npair) &
   !$omp shared(charge, gwdq, gwd2q, gwdcnq, dc6dq, d2c6dq2, d2c6dqij, d2c6dcnq, d2c6dcnqij) &
   !$omp private(icoup, jcoup, dq, dqq, dqqij, dcq, dcqij) &
   !$omp private(ipair, iat, jat, izp, jzp, iref, jref, refc6, c6ij, dc6dcni, dc6dcnj, &
   !$omp& d2c6dcni, d2c6dcnj, d2c6mix)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp) + 1.0_wp) - 1.0_wp))
      if (iat*(iat + 1)/2 < ipair) iat = iat + 1
      jat = ipair - iat*(iat - 1)/2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      ! D4 weights are atom-local; D4S weights depend on the partner element.
      icoup = merge(jat, 1, self%ncoup > 1)
      jcoup = merge(iat, 1, self%ncoup > 1)
      dq = 0.0_wp
      dqq = 0.0_wp
      dqqij = 0.0_wp
      dcq = 0.0_wp
      dcqij = 0.0_wp
      c6ij = 0.0_wp
      dc6dcni = 0.0_wp
      dc6dcnj = 0.0_wp
      d2c6dcni = 0.0_wp
      d2c6dcnj = 0.0_wp
      d2c6mix = 0.0_wp
      do iref = 1, self%ref(izp)
         do jref = 1, self%ref(jzp)
            refc6 = self%c6(iref, jref, izp, jzp)
            c6ij = c6ij + gwvec(iref, iat, icoup)*gwvec(jref, jat, jcoup)*refc6
            dc6dcni = dc6dcni + gwdcn(iref, iat, icoup)*gwvec(jref, jat, jcoup)*refc6
            dc6dcnj = dc6dcnj + gwvec(iref, iat, icoup)*gwdcn(jref, jat, jcoup)*refc6
            d2c6dcni = d2c6dcni + gwd2cn(iref, iat, icoup)*gwvec(jref, jat, jcoup)*refc6
            d2c6dcnj = d2c6dcnj + gwvec(iref, iat, icoup)*gwd2cn(jref, jat, jcoup)*refc6
            d2c6mix = d2c6mix + gwdcn(iref, iat, icoup)*gwdcn(jref, jat, jcoup)*refc6
            if (charge) then
               dq(1) = dq(1) + gwdq(iref, iat, icoup)*gwvec(jref, jat, jcoup)*refc6
               dq(2) = dq(2) + gwvec(iref, iat, icoup)*gwdq(jref, jat, jcoup)*refc6
               dqq(1) = dqq(1) + gwd2q(iref, iat, icoup)*gwvec(jref, jat, jcoup)*refc6
               dqq(2) = dqq(2) + gwvec(iref, iat, icoup)*gwd2q(jref, jat, jcoup)*refc6
               dqqij = dqqij + gwdq(iref, iat, icoup)*gwdq(jref, jat, jcoup)*refc6
               dcq(1) = dcq(1) + gwdcnq(iref, iat, icoup)*gwvec(jref, jat, jcoup)*refc6
               dcq(2) = dcq(2) + gwvec(iref, iat, icoup)*gwdcnq(jref, jat, jcoup)*refc6
               dcqij(1) = dcqij(1) + gwdcn(iref, iat, icoup)*gwdq(jref, jat, jcoup)*refc6
               dcqij(2) = dcqij(2) + gwdq(iref, iat, icoup)*gwdcn(jref, jat, jcoup)*refc6
            end if
         end do
      end do
      c6(iat, jat) = c6ij
      c6(jat, iat) = c6ij
      dc6dcn(iat, jat) = dc6dcni
      dc6dcn(jat, iat) = dc6dcnj
      d2c6dcn2(iat, jat) = d2c6dcni
      d2c6dcn2(jat, iat) = d2c6dcnj
      d2c6dcnij(iat, jat) = d2c6mix
      d2c6dcnij(jat, iat) = d2c6mix
      if (charge) then
         dc6dq(iat, jat) = dq(1)
         dc6dq(jat, iat) = dq(2)
         d2c6dq2(iat, jat) = dqq(1)
         d2c6dq2(jat, iat) = dqq(2)
         d2c6dqij(iat, jat) = dqqij
         d2c6dqij(jat, iat) = dqqij
         d2c6dcnq(iat, jat) = dcq(1)
         d2c6dcnq(jat, iat) = dcq(2)
         d2c6dcnqij(iat, jat) = dcqij(1)
         d2c6dcnqij(jat, iat) = dcqij(2)
      end if
   end do

end subroutine get_atomic_c6_hessian


end module dftd4_model_type
