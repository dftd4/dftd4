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

!> Definition of the D4S dispersion model for the evaluation of C6 coefficients.
module dftd4_model_d4s
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
   use ieee_arithmetic, only : ieee_is_nan
   use dftd4_model_type, only : dispersion_model, d4_qmod
   use dftd4_data, only : get_covalent_rad, get_r4r2_val, get_wfpair_val, &
      & get_effective_charge, get_electronegativity, get_hardness
   use dftd4_reference
   use dftd4_model_utils
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use multicharge, only : new_eeq2019_model, new_eeqbc2025_model
   implicit none
   private

   public :: d4s_model, new_d4s_model


   !> D4S dispersion model to evaluate C6 coefficients
   type, extends(dispersion_model) :: d4s_model

      !> Weighting factors for CN interpolation
      real(wp), allocatable :: wf(:, :)

   contains

      !> Generate weights for all reference systems
      procedure :: weight_references

      !> Evaluate C6 coefficient
      procedure :: get_atomic_c6

      !> Evaluate atomic polarizabilities
      procedure :: get_polarizabilities

   end type d4s_model


   !> Default maximum charge scaling height for partial charge extrapolation
   real(wp), parameter :: ga_default = 3.0_wp

   !> Default charge scaling steepness for partial charge extrapolation
   real(wp), parameter :: gc_default = 2.0_wp

contains


!> Create new D4S dispersion model from molecular structure input
subroutine new_d4s_model(error, d4, mol, ga, gc, qmod)
   !DEC$ ATTRIBUTES DLLEXPORT :: new_d4_model

   !> Instance of the dispersion model
   type(d4s_model), intent(out) :: d4

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Charge scaling height
   real(wp), intent(in), optional :: ga

   !> Charge scaling steepness
   real(wp), intent(in), optional :: gc

   !> Charge model selection
   integer, intent(in), optional :: qmod

   integer :: isp, izp, iref, jsp, jzp, jref
   integer :: mref, tmp_qmod
   real(wp) :: aiw(23), c6
   real(wp), parameter :: thopi = 3.0_wp/pi

   ! check for unsupported elements (104 (Rf) - 111 (Rg))
   do isp = 1, mol%nid
      if (mol%num(isp) > 103 .and. mol%num(isp) < 112) then
         call fatal_error(error, "Structure contains unsupported element '"//trim(mol%sym(isp))//"'")
         return
      end if
   end do

   d4%ncoup = mol%nat

   if (present(ga)) then
      d4%ga = ga
   else
      d4%ga = ga_default
   end if

   if (present(gc)) then
      d4%gc = gc
   else
      d4%gc = gc_default
   end if

   allocate(d4%wf(mol%nid, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, mol%nid
         jzp = mol%num(jsp)
         d4%wf(isp, jsp) = get_wfpair_val(izp, jzp)
      end do 
   end do

   allocate(d4%rcov(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      d4%rcov(isp) = get_covalent_rad(izp)
   end do

   allocate(d4%en(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      d4%en(isp) = get_electronegativity(izp)
   end do

   allocate(d4%zeff(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      d4%zeff(isp) = get_effective_charge(izp)
   end do

   allocate(d4%eta(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      d4%eta(isp) = get_hardness(izp)
   end do

   allocate(d4%r4r2(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      d4%r4r2(isp) = get_r4r2_val(izp)
   end do

   allocate(d4%ref(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      d4%ref(isp) = get_nref(izp)
   end do

   mref = maxval(d4%ref)
   allocate(d4%cn(mref, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      call set_refcn(d4%cn(:, isp), izp)
   end do

   if (present(qmod)) then
      tmp_qmod = qmod
   else
      tmp_qmod = d4_qmod%eeq
   end if
   
   allocate(d4%q(mref, mol%nid))
   allocate(d4%aiw(23, mref, mol%nid))
   select case(tmp_qmod)
   case default
      call fatal_error(error, "Unsupported option for charge model.")
      return
   case(d4_qmod%eeq)
      do isp = 1, mol%nid
         izp = mol%num(isp)
         call set_refq_eeq(d4%q(:, isp), izp)
         call set_refalpha_eeq(d4%aiw(:, :, isp), d4%ga, d4%gc, izp)
      end do
      ! Setup EEQ model
      call new_eeq2019_model(mol, d4%mchrg, error)
      if(allocated(error)) return
   case(d4_qmod%eeqbc)
      do isp = 1, mol%nid
         izp = mol%num(isp)
         call set_refq_eeqbc(d4%q(:, isp), izp)
         call set_refalpha_eeqbc(d4%aiw(:, :, isp), d4%ga, d4%gc, izp)
      end do
      ! Setup EEQBC model
      call new_eeqbc2025_model(mol, d4%mchrg, error)  
      if(allocated(error)) return
   case(d4_qmod%gfn2)
      do isp = 1, mol%nid
         izp = mol%num(isp)
         call set_refq_gfn2(d4%q(:, isp), izp)
         call set_refalpha_gfn2(d4%aiw(:, :, isp), d4%ga, d4%gc, izp)
      end do
   end select

   allocate(d4%ngw(mref, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      call set_refgw(d4%ngw(:, isp), izp)
   end do

   allocate(d4%c6(mref, mref, mol%nid, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, isp
         jzp = mol%num(jsp)
         do iref = 1, d4%ref(isp)
            do jref = 1, d4%ref(jsp)
               aiw(:) = d4%aiw(:, iref, isp) * d4%aiw(:, jref, jsp)
               c6 = thopi * trapzd(aiw)
               d4%c6(jref, iref, jsp, isp) = c6
               d4%c6(iref, jref, isp, jsp) = c6
            end do
         end do
      end do
   end do

end subroutine new_d4s_model


!> Calculate the weights of the reference system and the derivatives w.r.t.
!> coordination number for later use.
subroutine weight_references(self, mol, cn, q, gwvec, gwdcn, gwdq)
   !DEC$ ATTRIBUTES DLLEXPORT :: weight_references

   !> Instance of the dispersion model
   class(d4s_model), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Coordination number of every atom
   real(wp), intent(in) :: cn(:)

   !> Partial charge of every atom
   real(wp), intent(in) :: q(:)

   !> Pairwise weighting for the atomic reference systems
   real(wp), intent(out) :: gwvec(:, :, :)

   !> derivative of the pairwise weighting function w.r.t. the coordination number
   real(wp), intent(out), optional :: gwdcn(:, :, :)

   !> derivative of the pairwise weighting function w.r.t. the charge scaling
   real(wp), intent(out), optional :: gwdq(:, :, :)

   integer :: iat, izp, iref, igw, jat, jzp
   real(wp) :: norm, dnorm, gw, expw, expd, gwk, dgwk, wf, zi, gi, maxcn
   
   if (present(gwdcn) .and. present(gwdq)) then
      gwvec(:, :, :) = 0.0_wp
      gwdcn(:, :, :) = 0.0_wp
      gwdq(:, :, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(gwvec, gwdcn, gwdq, mol, self, cn, q) &
      !$omp private(iat, izp, iref, igw, zi, gi, jat, jzp) &
      !$omp private(norm, dnorm, gw, expw, expd, gwk, dgwk, wf, maxcn)      
      do iat = 1, mol%nat
         izp = mol%id(iat)
         zi = self%zeff(izp)
         gi = self%eta(izp) * self%gc

         do jat = 1, mol%nat
            jzp = mol%id(jat)

            norm = 0.0_wp
            dnorm = 0.0_wp
            do iref = 1, self%ref(izp)
               do igw = 1, self%ngw(iref, izp)
                  wf = igw * self%wf(izp, jzp)
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
                  wf = igw * self%wf(izp, jzp)
                  gw = weight_cn(wf, cn(iat), self%cn(iref, izp))
                  expw = expw + gw
                  expd = expd + 2*wf * (self%cn(iref, izp) - cn(iat)) * gw
               end do
               gwk = expw * norm
               if (is_exceptional(gwk)) then
                  maxcn = maxval(self%cn(:self%ref(izp), izp))
                  if (abs(maxcn - self%cn(iref, izp)) < 1e-12_wp) then
                     gwk = 1.0_wp
                  else
                     gwk = 0.0_wp
                  end if
               end if
               gwvec(iref, iat, jat) = gwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
               gwdq(iref, iat, jat) = gwk * dzeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
               
               dgwk = norm * (expd - expw * dnorm * norm)
               if (is_exceptional(dgwk)) then
                  dgwk = 0.0_wp
               end if
               gwdcn(iref, iat, jat) = dgwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
            end do

         end do 
      end do
   else

      gwvec(:, :, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(gwvec, mol, self, cn, q) &
      !$omp private(iat, izp, iref, igw, zi, gi, jat, jzp) &
      !$omp private(norm, gw, expw, gwk, wf, maxcn)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         zi = self%zeff(izp)
         gi = self%eta(izp) * self%gc

         do jat = 1, mol%nat
            jzp = mol%id(jat)

            norm = 0.0_wp
            do iref = 1, self%ref(izp)
               do igw = 1, self%ngw(iref, izp)
                  wf = igw * self%wf(izp, jzp)
                  norm = norm + weight_cn(wf, cn(iat), self%cn(iref, izp))
               end do
            end do
            norm = 1.0_wp / norm
            
            do iref = 1, self%ref(izp)
               expw = 0.0_wp
               do igw = 1, self%ngw(iref, izp)
                  wf = igw * self%wf(izp, jzp)
                  expw = expw + weight_cn(wf, cn(iat), self%cn(iref, izp))
               end do
               gwk = expw * norm
               if (is_exceptional(gwk)) then
                  maxcn = maxval(self%cn(:self%ref(izp), izp))
                  if (abs(maxcn - self%cn(iref, izp)) < 1e-12_wp) then
                     gwk = 1.0_wp
                  else
                     gwk = 0.0_wp
                  end if
               end if
               gwvec(iref, iat, jat) = gwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
            end do

         end do 
      end do
   end if

end subroutine weight_references


!> Calculate atomic dispersion coefficients and their derivatives w.r.t.
!> the coordination numbers and atomic partial charges.
subroutine get_atomic_c6(self, mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_atomic_c6

   !> Instance of the dispersion model
   class(d4s_model), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Pairwise weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :, :)

   !> Derivative of the pairwise weighting function w.r.t. the coordination number
   real(wp), intent(in), optional :: gwdcn(:, :, :)

   !> Derivative of the pairwise weighting function w.r.t. the partial charge
   real(wp), intent(in), optional :: gwdq(:, :, :)

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
                  dc6 = dc6 + gwvec(iref, iat, jat) * gwvec(jref, jat, iat) * refc6
                  dc6dcni = dc6dcni + gwdcn(iref, iat, jat) * gwvec(jref, jat, iat) * refc6
                  dc6dcnj = dc6dcnj + gwvec(iref, iat, jat) * gwdcn(jref, jat, iat) * refc6
                  dc6dqi = dc6dqi + gwdq(iref, iat, jat) * gwvec(jref, jat, iat) * refc6
                  dc6dqj = dc6dqj + gwvec(iref, iat, jat) * gwdq(jref, jat, iat) * refc6
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
                  dc6 = dc6 + gwvec(iref, iat, jat) * gwvec(jref, jat, iat) * refc6
               end do
            end do
            c6(iat, jat) = dc6
            c6(jat, iat) = dc6
         end do
      end do
   end if

end subroutine get_atomic_c6


!> Calculate atomic polarizabilities and their derivatives w.r.t.
!> the coordination numbers and atomic partial charges.
subroutine get_polarizabilities(self, mol, gwvec, gwdcn, gwdq, alpha, dadcn, dadq)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_polarizabilities

   !> Instance of the dispersion model
   class(d4s_model), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Pairwise weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :, :)

   !> Derivative of the pairwise weighting function w.r.t. the coordination number
   real(wp), intent(in), optional :: gwdcn(:, :, :)

   !> Derivative of the pairwise weighting function w.r.t. the partial charge
   real(wp), intent(in), optional :: gwdq(:, :, :)

   !> Static polarizabilities for all atoms.
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
            da = da + gwvec(iref, iat, iat) * refa
            dadcni = dadcni + gwdcn(iref, iat, iat) * refa
            dadqi = dadqi + gwdq(iref, iat, iat) * refa
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
            da = da + gwvec(iref, iat, iat) * self%aiw(1, iref, izp)
         end do
         alpha(iat) = da
      end do
   end if

end subroutine get_polarizabilities


end module dftd4_model_d4s
