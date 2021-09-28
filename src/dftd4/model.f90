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

!> Definition of the D4 dispersion model for the evaluation of C6 coefficients.
module dftd4_model
   use ieee_arithmetic, only : ieee_is_nan
   use dftd4_data, only : get_covalent_rad, get_r4r2_val, get_effective_charge, &
      get_electronegativity, get_hardness
   use dftd4_reference
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   implicit none
   private

   public :: d4_model, new_d4_model, d4_ref


   !> Base D4 dispersion model to evaluate C6 coefficients
   type :: d4_model

      !> Charge scaling height
      real(wp) :: ga

      !> Charge scaling steepness
      real(wp) :: gc

      !> Weighting factor for CN interpolation
      real(wp) :: wf

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

      !> Reference dynamic polarizibilities
      real(wp), allocatable :: aiw(:, :, :)

      !> Reference C6 coefficients
      real(wp), allocatable :: c6(:, :, :, :)

   contains

      !> Generate weights for all reference systems
      procedure :: weight_references

      !> Evaluate C6 coefficient
      procedure :: get_atomic_c6

      !> Evaluate atomic polarizibilities
      procedure :: get_polarizibilities

   end type d4_model


   !> Default maximum charge scaling height for partial charge extrapolation
   real(wp), parameter :: ga_default = 3.0_wp

   !> Default charge scaling steepness for partial charge extrapolation
   real(wp), parameter :: gc_default = 2.0_wp

   !> Default weighting factor for coordination number interpolation
   real(wp), parameter :: wf_default = 6.0_wp


   !> Possible reference charges for D4
   type :: enum_ref

      !> Electronegativity equilibration charges
      integer :: eeq = 1

      !> GFN2-xTB Mulliken partial charges
      integer :: gfn2 = 2

   end type enum_ref

   !> Actual enumerator for D4 reference charges
   type(enum_ref), parameter :: d4_ref = enum_ref()


contains


!> Create new dispersion model from molecular structure input
subroutine new_d4_model(self, mol, ga, gc, wf, ref)

   !> Instance of the dispersion model
   type(d4_model), intent(out) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Charge scaling height
   real(wp), intent(in), optional :: ga

   !> Charge scaling steepness
   real(wp), intent(in), optional :: gc

   !> Weighting factor for coordination number interpolation
   real(wp), intent(in), optional :: wf

   !> Reference charge selection
   integer, intent(in), optional :: ref

   integer :: isp, izp, iref, jsp, jzp, jref
   integer :: mref, ref_charge
   real(wp) :: aiw(23), c6
   real(wp), parameter :: thopi = 3.0_wp/pi

   if (present(ref)) then
      ref_charge = ref
   else
      ref_charge = d4_ref%eeq
   end if

   if (present(ga)) then
      self%ga = ga
   else
      self%ga = ga_default
   end if

   if (present(gc)) then
      self%gc = gc
   else
      self%gc = gc_default
   end if

   if (present(wf)) then
      self%wf = wf
   else
      self%wf = wf_default
   end if

   allocate(self%rcov(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      self%rcov(isp) = get_covalent_rad(izp)
   end do

   allocate(self%en(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      self%en(isp) = get_electronegativity(izp)
   end do

   allocate(self%zeff(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      self%zeff(isp) = get_effective_charge(izp)
   end do

   allocate(self%eta(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      self%eta(isp) = get_hardness(izp)
   end do

   allocate(self%r4r2(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      self%r4r2(isp) = get_r4r2_val(izp)
   end do

   allocate(self%ref(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      self%ref(isp) = get_nref(izp)
   end do

   mref = maxval(self%ref)
   allocate(self%cn(mref, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      call set_refcn(self%cn(:, isp), izp)
   end do

   allocate(self%q(mref, mol%nid))
   select case(ref_charge)
   case default
      do isp = 1, mol%nid
         izp = mol%num(isp)
         call set_refq_eeq(self%q(:, isp), izp)
      end do
   case(d4_ref%gfn2)
      do isp = 1, mol%nid
         izp = mol%num(isp)
         call set_refq_gfn2(self%q(:, isp), izp)
      end do
   end select

   allocate(self%ngw(mref, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      call set_refgw(self%ngw(:, isp), izp)
   end do

   allocate(self%aiw(23, mref, mol%nid))
   select case(ref_charge)
   case default
      do isp = 1, mol%nid
         izp = mol%num(isp)
         call set_refalpha_eeq(self%aiw(:, :, isp), self%ga, self%gc, izp)
      end do
   case(d4_ref%gfn2)
      do isp = 1, mol%nid
         izp = mol%num(isp)
         call set_refalpha_gfn2(self%aiw(:, :, isp), self%ga, self%gc, izp)
      end do
   end select

   allocate(self%c6(mref, mref, mol%nid, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, isp
         jzp = mol%num(jsp)
         do iref = 1, self%ref(isp)
            do jref = 1, self%ref(jsp)
               aiw(:) = self%aiw(:, iref, isp) * self%aiw(:, jref, jsp)
               c6 = thopi * trapzd(aiw)
               self%c6(jref, iref, jsp, isp) = c6
               self%c6(iref, jref, isp, jsp) = c6
            end do
         end do
      end do
   end do

end subroutine new_d4_model


!> Calculate the weights of the reference system and the derivatives w.r.t.
!> coordination number for later use.
subroutine weight_references(self, mol, cn, q, gwvec, gwdcn, gwdq)

   !> Instance of the dispersion model
   class(d4_model), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Coordination number of every atom
   real(wp), intent(in) :: cn(:)

   !> Partial charge of every atom
   real(wp), intent(in) :: q(:)

   !> weighting for the atomic reference systems
   real(wp), intent(out) :: gwvec(:, :)

   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(out), optional :: gwdcn(:, :)

   !> derivative of the weighting function w.r.t. the charge scaling
   real(wp), intent(out), optional :: gwdq(:, :)

   integer :: iat, izp, iref, igw
   real(wp) :: norm, dnorm, gw, expw, expd, gwk, dgwk, wf, zi, gi

   if (present(gwdcn) .and. present(gwdq)) then
      gwvec(:, :) = 0.0_wp
      gwdcn(:, :) = 0.0_wp
      gwdq(:, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(gwvec, gwdcn, gwdq, mol, self, cn, q) private(iat, izp, iref, &
      !$omp& igw, norm, dnorm, gw, expw, expd, gwk, dgwk, wf, zi, gi)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         zi = self%zeff(izp)
         gi = self%eta(izp) * self%gc
         norm = 0.0_wp
         dnorm = 0.0_wp
         do iref = 1, self%ref(izp)
            do igw = 1, self%ngw(iref, izp)
               wf = igw * self%wf
               gw = weight_cn(wf, cn(iat), self%cn(iref, izp))
               norm = norm + gw
               dnorm = dnorm + 2*wf * (self%cn(iref, izp) - cn(iat)) * gw
            end do
         end do
         norm = 1.0_wp / norm
         do iref = 1, self%ref(izp)
            expw = 0.0_wp
            expd = 0.0_wp
            do igw = 1, self%ngw(iref, izp)
               wf = igw * self%wf
               gw = weight_cn(wf, cn(iat), self%cn(iref, izp))
               expw = expw + gw
               expd = expd + 2*wf * (self%cn(iref, izp) - cn(iat)) * gw
            end do
            gwk = expw * norm
            if (is_exceptional(gwk)) then
               if (maxval(self%cn(:self%ref(izp), izp)) == self%cn(iref, izp)) then
                  gwk = 1.0_wp
               else
                  gwk = 0.0_wp
               end if
            end if
            gwvec(iref, iat) = gwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
            gwdq(iref, iat) = gwk * dzeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)

            dgwk = norm * (expd - expw * dnorm * norm)
            if (is_exceptional(dgwk)) then
               dgwk = 0.0_wp
            end if
            gwdcn(iref, iat) = dgwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
         end do
      end do

   else

      gwvec(:, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(gwvec, mol, self, cn, q) &
      !$omp private(iat, izp, iref, igw, norm, gw, expw, gwk, wf, zi, gi)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         zi = self%zeff(izp)
         gi = self%eta(izp) * self%gc
         norm = 0.0_wp
         do iref = 1, self%ref(izp)
            do igw = 1, self%ngw(iref, izp)
               wf = igw * self%wf
               gw = weight_cn(wf, cn(iat), self%cn(iref, izp))
               norm = norm + gw
            end do
         end do
         norm = 1.0_wp / norm
         do iref = 1, self%ref(izp)
            expw = 0.0_wp
            do igw = 1, self%ngw(iref, izp)
               wf = igw * self%wf
               expw = expw + weight_cn(wf, cn(iat), self%cn(iref, izp))
            end do
            gwk = expw * norm
            if (is_exceptional(gwk)) then
               if (maxval(self%cn(:self%ref(izp), izp)) == self%cn(iref, izp)) then
                  gwk = 1.0_wp
               else
                  gwk = 0.0_wp
               end if
            end if
            gwvec(iref, iat) = gwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
         end do
      end do
   end if

end subroutine weight_references


!> Check whether we are dealing with an exceptional value, NaN or Inf
elemental function is_exceptional(val)
   real(wp), intent(in) :: val
   logical :: is_exceptional
   is_exceptional = ieee_is_nan(val) .or. abs(val) > huge(val)
end function is_exceptional


!> Calculate atomic dispersion coefficients and their derivatives w.r.t.
!> the coordination numbers and atomic partial charges.
subroutine get_atomic_c6(self, mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   !> Instance of the dispersion model
   class(d4_model), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :)

   !> Derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in), optional :: gwdcn(:, :)

   !> Derivative of the weighting function w.r.t. the partial charge
   real(wp), intent(in), optional :: gwdq(:, :)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(out) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(out), optional :: dc6dcn(:, :)

   !> Derivative of the C6 w.r.t. the partial charge
   real(wp), intent(out), optional :: dc6dq(:, :)

   integer :: iat, jat, izp, jzp, iref, jref
   real(wp) :: refc6, dc6, dc6dcni, dc6dcnj, dc6dqi, dc6dqj

   if (present(gwdcn).and.present(dc6dcn) &
      & .and.present(gwdq).and.present(dc6dq)) then
      c6(:, :) = 0.0_wp
      dc6dcn(:, :) = 0.0_wp
      dc6dq(:, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(c6, dc6dcn, dc6dq, mol, self, gwvec, gwdcn, gwdq) &
      !$omp private(iat, jat, izp, jzp, iref, jref, refc6, dc6, dc6dqi, dc6dqj, &
      !$omp& dc6dcni, dc6dcnj)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            dc6 = 0.0_wp
            dc6dcni = 0.0_wp
            dc6dcnj = 0.0_wp
            dc6dqi = 0.0_wp
            dc6dqj = 0.0_wp
            do iref = 1, self%ref(izp)
               do jref = 1, self%ref(jzp)
                  refc6 = self%c6(iref, jref, izp, jzp)
                  dc6 = dc6 + gwvec(iref, iat) * gwvec(jref, jat) * refc6
                  dc6dcni = dc6dcni + gwdcn(iref, iat) * gwvec(jref, jat) * refc6
                  dc6dcnj = dc6dcnj + gwvec(iref, iat) * gwdcn(jref, jat) * refc6
                  dc6dqi = dc6dqi + gwdq(iref, iat) * gwvec(jref, jat) * refc6
                  dc6dqj = dc6dqj + gwvec(iref, iat) * gwdq(jref, jat) * refc6
               end do
            end do
            c6(iat, jat) = dc6
            c6(jat, iat) = dc6
            dc6dcn(iat, jat) = dc6dcni
            dc6dcn(jat, iat) = dc6dcnj
            dc6dq(iat, jat) = dc6dqi
            dc6dq(jat, iat) = dc6dqj
         end do
      end do

   else

      c6(:, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(c6, mol, self, gwvec) &
      !$omp private(iat, jat, izp, jzp, iref, jref, refc6, dc6)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            dc6 = 0.0_wp
            do iref = 1, self%ref(izp)
               do jref = 1, self%ref(jzp)
                  refc6 = self%c6(iref, jref, izp, jzp)
                  dc6 = dc6 + gwvec(iref, iat) * gwvec(jref, jat) * refc6
               end do
            end do
            c6(iat, jat) = dc6
            c6(jat, iat) = dc6
         end do
      end do
   end if

end subroutine get_atomic_c6


!> Calculate atomic polarizibilities and their derivatives w.r.t.
!> the coordination numbers and atomic partial charges.
subroutine get_polarizibilities(self, mol, gwvec, gwdcn, gwdq, alpha, dadcn, dadq)

   !> Instance of the dispersion model
   class(d4_model), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :)

   !> Derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in), optional :: gwdcn(:, :)

   !> Derivative of the weighting function w.r.t. the partial charge
   real(wp), intent(in), optional :: gwdq(:, :)

   !> Static polarizibilities for all atoms.
   real(wp), intent(out) :: alpha(:)

   !> Derivative of the polarizibility w.r.t. the coordination number
   real(wp), intent(out), optional :: dadcn(:)

   !> Derivative of the polarizibility w.r.t. the partial charge
   real(wp), intent(out), optional :: dadq(:)

   integer :: iat, izp, iref
   real(wp) :: refa, da, dadcni, dadqi

   if (present(gwdcn).and.present(dadcn) &
      & .and.present(gwdq).and.present(dadq)) then
      alpha(:) = 0.0_wp
      dadcn(:) = 0.0_wp
      dadq(:) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(alpha, dadcn, dadq, mol, self, gwvec, gwdcn, gwdq) &
      !$omp private(iat, izp, iref, refa, da, dadqi, dadcni)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         da = 0.0_wp
         dadcni = 0.0_wp
         dadqi = 0.0_wp
         do iref = 1, self%ref(izp)
            refa = self%aiw(1, iref, izp)
            da = da + gwvec(iref, iat) * refa
            dadcni = dadcni + gwdcn(iref, iat) * refa
            dadqi = dadqi + gwdq(iref, iat) * refa
         end do
         alpha(iat) = da
         dadcn(iat) = dadcni
         dadq(iat) = dadqi
      end do

   else

      alpha(:) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(alpha, mol, self, gwvec) private(iat, izp, iref, refa, da)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         da = 0.0_wp
         do iref = 1, self%ref(izp)
            da = da + gwvec(iref, iat) * self%aiw(1, iref, izp)
         end do
         alpha(iat) = da
      end do
   end if

end subroutine get_polarizibilities


elemental function weight_cn(wf,cn,cnref) result(cngw)
   real(wp),intent(in) :: wf, cn, cnref
   real(wp) :: cngw
   intrinsic :: exp
   cngw = exp ( -wf * ( cn - cnref )**2 )
end function weight_cn

!> charge scaling function
elemental function zeta(a, c, qref, qmod)
   real(wp), intent(in) :: a
   real(wp), intent(in) :: c
   real(wp), intent(in) :: qref
   real(wp), intent(in) :: qmod
   real(wp) :: zeta

   intrinsic :: exp

   if (qmod < 0.0_wp) then
      zeta = exp( a )
   else
      zeta = exp( a * ( 1.0_wp - exp( c * ( 1.0_wp - qref/qmod ) ) ) )
   endif

end function zeta

!> derivative of charge scaling function w.r.t. charge
elemental function dzeta(a, c, qref, qmod)
   real(wp), intent(in) :: a
   real(wp), intent(in) :: c
   real(wp), intent(in) :: qref
   real(wp), intent(in) :: qmod
   real(wp) :: dzeta

   intrinsic :: exp

   if (qmod < 0.0_wp) then
      dzeta = 0.0_wp
   else
      dzeta = - a * c * exp( c * ( 1.0_wp - qref/qmod ) ) &
         & * zeta(a,c,qref,qmod) * qref / ( qmod**2 )
   endif

end function dzeta

!> numerical Casimir--Polder integration
pure function trapzd(pol)
   real(wp), intent(in) :: pol(23)
   real(wp) :: trapzd

   real(wp), parameter :: freq(23) = [ &
      & 0.000001_wp, 0.050000_wp, 0.100000_wp, &
      & 0.200000_wp, 0.300000_wp, 0.400000_wp, &
      & 0.500000_wp, 0.600000_wp, 0.700000_wp, &
      & 0.800000_wp, 0.900000_wp, 1.000000_wp, &
      & 1.200000_wp, 1.400000_wp, 1.600000_wp, &
      & 1.800000_wp, 2.000000_wp, 2.500000_wp, &
      & 3.000000_wp, 4.000000_wp, 5.000000_wp, &
      & 7.500000_wp, 10.00000_wp]
   real(wp), parameter :: weights(23) = 0.5_wp * [ &
      &  ( freq (2) - freq (1) ),  &
      &  ( freq (2) - freq (1) ) + ( freq (3) - freq (2) ),  &
      &  ( freq (3) - freq (2) ) + ( freq (4) - freq (3) ),  &
      &  ( freq (4) - freq (3) ) + ( freq (5) - freq (4) ),  &
      &  ( freq (5) - freq (4) ) + ( freq (6) - freq (5) ),  &
      &  ( freq (6) - freq (5) ) + ( freq (7) - freq (6) ),  &
      &  ( freq (7) - freq (6) ) + ( freq (8) - freq (7) ),  &
      &  ( freq (8) - freq (7) ) + ( freq (9) - freq (8) ),  &
      &  ( freq (9) - freq (8) ) + ( freq(10) - freq (9) ),  &
      &  ( freq(10) - freq (9) ) + ( freq(11) - freq(10) ),  &
      &  ( freq(11) - freq(10) ) + ( freq(12) - freq(11) ),  &
      &  ( freq(12) - freq(11) ) + ( freq(13) - freq(12) ),  &
      &  ( freq(13) - freq(12) ) + ( freq(14) - freq(13) ),  &
      &  ( freq(14) - freq(13) ) + ( freq(15) - freq(14) ),  &
      &  ( freq(15) - freq(14) ) + ( freq(16) - freq(15) ),  &
      &  ( freq(16) - freq(15) ) + ( freq(17) - freq(16) ),  &
      &  ( freq(17) - freq(16) ) + ( freq(18) - freq(17) ),  &
      &  ( freq(18) - freq(17) ) + ( freq(19) - freq(18) ),  &
      &  ( freq(19) - freq(18) ) + ( freq(20) - freq(19) ),  &
      &  ( freq(20) - freq(19) ) + ( freq(21) - freq(20) ),  &
      &  ( freq(21) - freq(20) ) + ( freq(22) - freq(21) ),  &
      &  ( freq(22) - freq(21) ) + ( freq(23) - freq(22) ),  &
      &  ( freq(23) - freq(22) ) ]

   trapzd = sum(pol*weights)

end function trapzd


end module dftd4_model
