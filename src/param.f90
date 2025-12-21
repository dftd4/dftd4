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

module dftd4_param
   use mctc_env, only : wp
   use dftd4_damping, only : damping_param
   use dftd4_damping_rational, only : rational_damping_param
   use dftd4_utils, only : lowercase
   implicit none
   private

   public :: functional_group
   public :: get_rational_damping, get_functionals, get_functional_id
   public :: p_r2scan_3c


   enum, bind(C)
      enumerator :: p_invalid, &
         & p_hf, p_blyp, p_bpbe, p_bp, p_bpw, p_lb94, p_mpwlyp, p_mpwpw, &
         & p_olyp, p_opbe, p_pbe, p_rpbe, p_revpbe, p_pw86pbe, &
         & p_rpw86pbe, p_pw91, p_pwp, p_xlyp, p_b97, p_tpss, p_revtpss, &
         & p_scan, p_rscan, p_r2scan, p_b1lyp, p_b3lyp, p_bhlyp, p_b1p, &
         & p_b3p, p_b1pw, p_b3pw, p_o3lyp, p_revpbe0, p_revpbe38, &
         & p_pbe0, p_pwp1, p_pw1pw, p_mpw1pw, p_mpw1lyp, p_pw6b95, &
         & p_tpssh, p_tpss0, p_x3lyp, p_m06l, p_m06, p_b97d, &
         & p_wb97, p_wb97x_2008, p_b97m, p_wb97m, p_camb3lyp, p_lcblyp, &
         & p_lh07tsvwn, p_lh07ssvwn, p_lh12ctssirpw92, p_lh12ctssifpw92, &
         & p_lh14tcalpbe, p_lh20t, p_b2plyp, p_b2gpplyp, p_mpw2plyp, p_pwpb95, &
         & p_dsdblyp, p_dsdpbe, p_dsdpbeb95, p_dsdpbep86, p_dsdsvwn, &
         & p_dodblyp, p_dodpbe, p_dodpbeb95, p_dodpbep86, p_dodsvwn, &
         & p_pbe0_2, p_pbe0_dh, p_hsesol, p_dftb_3ob, p_dftb_mio, p_dftb_ob2, &
         & p_dftb_matsci, p_dftb_pbc, p_b1b95, p_pbesol, p_hse06, p_mpwb1k, &
         & p_hse03, p_revtpssh, p_mn12sx, p_glyp, p_mpw1b95, &
         & p_revpbe0dh, p_revtpss0, p_revdsdpbep86, p_revdsdpbe, &
         & p_revdsdblyp, p_revdodpbep86, p_am05, p_hse12, p_hse12s, &
         & p_r2scanh, p_r2scan0, p_r2scan50, p_r2scan_3c, p_camqtp01, &
         & p_lcwpbe, p_lcwpbeh, p_wb97x_rev, p_wb97m_rev, &
         & p_wb97x_3c, p_wr2scan, p_r2scan0_dh, p_r2scan_cidh, &
         & p_r2scan_qidh, p_r2scan0_2, p_pr2scan50, p_pr2scan69, &
         & p_kpr2scan50, p_wpr2scan50, p_wb97x, p_last
   end enum
   integer, parameter :: df_enum = kind(p_invalid)


   !> Group different spellings/names of functionals
   type functional_group
      character(len=:), allocatable :: names(:)
   end type functional_group

   !> Retrieve rational damping parameters from functional name or ID
   interface get_rational_damping
      module procedure :: get_rational_damping_name
      module procedure :: get_rational_damping_id
   end interface get_rational_damping

contains


!> Create a new group of functional names
function new_funcgroup(input_names) result(group)

   !> List of spellings/names of the functional
   character(len=*), intent(in) :: input_names(:)

   !> Functional with possibly different spellings
   type(functional_group) :: group
   
   integer :: n, i, max_len
   n = size(input_names)

   ! Determine the length of the longest name
   max_len = 0
   do i = 1, n
      max_len = max(max_len, len_trim(input_names(i)))
   end do

   ! Allocate based on the longest name's length
   allocate(character(len=max_len) :: group%names(n))
   do i = 1, n
      group%names(i) = trim(input_names(i))
   end do
end function new_funcgroup


!> Collect all supported functionals
subroutine get_functionals(funcs)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_functionals

   !> Collection of functionals with possibly different spellings/names
   type(functional_group), allocatable, intent(out) :: funcs(:)

   allocate(funcs(p_last - 1))

   funcs(p_hf) = new_funcgroup([character(len=20) :: 'hf'])
   funcs(p_am05) = new_funcgroup([character(len=20) :: 'am05'])
   funcs(p_blyp) = new_funcgroup([character(len=20) :: 'b-lyp', 'blyp'])
   funcs(p_bpbe) = new_funcgroup([character(len=20) :: 'bpbe'])
   funcs(p_bp) = new_funcgroup([character(len=20) :: 'b-p', 'bp86', 'bp', 'b-p86'])
   funcs(p_bpw) = new_funcgroup([character(len=20) :: 'bpw', 'b-pw'])
   funcs(p_lb94) = new_funcgroup([character(len=20) :: 'lb94'])
   funcs(p_mpwlyp) = new_funcgroup([character(len=20) :: 'mpwlyp', 'mpw-lyp'])
   funcs(p_mpwpw) = new_funcgroup([character(len=20) :: 'mpwpw', 'mpw-pw', 'mpwpw91'])
   funcs(p_olyp) = new_funcgroup([character(len=20) :: 'o-lyp', 'olyp'])
   funcs(p_opbe) = new_funcgroup([character(len=20) :: 'opbe'])
   funcs(p_pbe) = new_funcgroup([character(len=20) :: 'pbe'])
   funcs(p_rpbe) = new_funcgroup([character(len=20) :: 'rpbe'])
   funcs(p_revpbe) = new_funcgroup([character(len=20) :: 'revpbe'])
   funcs(p_pw86pbe) = new_funcgroup([character(len=20) :: 'pw86pbe'])
   funcs(p_rpw86pbe) = new_funcgroup([character(len=20) :: 'rpw86pbe'])
   funcs(p_pw91) = new_funcgroup([character(len=20) :: 'pw91'])
   funcs(p_pwp) = new_funcgroup([character(len=20) :: 'pwp', 'pw-p', 'pw91p86'])
   funcs(p_xlyp) = new_funcgroup([character(len=20) :: 'x-lyp', 'xlyp'])
   funcs(p_b97) = new_funcgroup([character(len=20) :: 'b97'])
   funcs(p_tpss) = new_funcgroup([character(len=20) :: 'tpss'])
   funcs(p_revtpss) = new_funcgroup([character(len=20) :: 'revtpss'])
   funcs(p_scan) = new_funcgroup([character(len=20) :: 'scan'])
   funcs(p_rscan) = new_funcgroup([character(len=20) :: 'rscan'])
   funcs(p_r2scan) = new_funcgroup([character(len=20) :: 'r2scan', 'r²scan'])
   funcs(p_r2scanh) = new_funcgroup([character(len=20) :: 'r2scanh', 'r²scanh'])
   funcs(p_r2scan0) = new_funcgroup([character(len=20) :: 'r2scan0', 'r²scan0'])
   funcs(p_r2scan50) = new_funcgroup([character(len=20) :: 'r2scan50', 'r²scan50'])
   funcs(p_r2scan_3c) = new_funcgroup([character(len=20) :: 'r2scan-3c', &
      & 'r²scan-3c', 'r2scan_3c', 'r²scan_3c', 'r2scan3c'])
   funcs(p_wr2scan) = new_funcgroup([character(len=20) :: 'wr2scan', 'wr²scan'])
   funcs(p_r2scan0_dh) = new_funcgroup([character(len=20) :: 'r2scan0-dh', &
      & 'r²scan0-dh', 'r2scan0dh', 'r²scan0dh'])
   funcs(p_r2scan_cidh) = new_funcgroup([character(len=20) :: 'r2scan-cidh', &
      & 'r²scan-cidh', 'r2scancidh', 'r²scancidh'])
   funcs(p_r2scan_qidh) = new_funcgroup([character(len=20) :: 'r2scan-qidh', &  
      & 'r²scan-qidh', 'r2scanqidh', 'r²scanqidh'])
    funcs(p_r2scan0_2) = new_funcgroup([character(len=20) :: 'r2scan0-2', &
      & 'r²scan0-2', 'r2scan02', 'r²scan02'])
   funcs(p_pr2scan50) = new_funcgroup([character(len=20) :: 'pr2scan50', &
      & 'pr²scan50', 'pr2scan50', 'pr²scan50'])
   funcs(p_pr2scan69) = new_funcgroup([character(len=20) :: 'pr2scan69', &
      & 'pr²scan69', 'pr2scan69', 'pr²scan69'])
   funcs(p_kpr2scan50) = new_funcgroup([character(len=20) :: 'kpr2scan50', & 
      & 'kpr²scan50', 'kpr2scan50', 'kpr²scan50'])
   funcs(p_wpr2scan50) = new_funcgroup([character(len=20) :: 'wpr2scan50', &
      & 'wpr²scan50', 'wpr2scan50', 'wpr²scan50'])
   funcs(p_b1lyp) = new_funcgroup([character(len=20) :: 'b1lyp', 'b1-lyp'])
   funcs(p_b3lyp) = new_funcgroup([character(len=20) :: 'b3-lyp', 'b3lyp'])
   funcs(p_bhlyp) = new_funcgroup([character(len=20) :: 'bh-lyp', 'bhlyp'])
   funcs(p_b1p) = new_funcgroup([character(len=20) :: 'b1p', 'b1-p', 'b1p86'])
   funcs(p_b3p) = new_funcgroup([character(len=20) :: 'b3p', 'b3-p', 'b3p86'])
   funcs(p_b1pw) = new_funcgroup([character(len=20) :: 'b1pw', 'b1-pw', 'b1pw91'])
   funcs(p_b3pw) = new_funcgroup([character(len=20) :: 'b3pw', 'b3-pw', 'b3pw91'])
   funcs(p_o3lyp) = new_funcgroup([character(len=20) :: 'o3-lyp', 'o3lyp'])
   funcs(p_revpbe0) = new_funcgroup([character(len=20) :: 'revpbe0'])
   funcs(p_revpbe38) = new_funcgroup([character(len=20) :: 'revpbe38'])
   funcs(p_pbe0) = new_funcgroup([character(len=20) :: 'pbe0'])
   funcs(p_pwp1) = new_funcgroup([character(len=20) :: 'pwp1'])
   funcs(p_pw1pw) = new_funcgroup([character(len=20) :: 'pw1pw', 'pw1-pw'])
   funcs(p_mpw1pw) = new_funcgroup([character(len=20) :: 'mpw1pw', 'mpw1-pw', 'mpw1pw91'])
   funcs(p_mpw1lyp) = new_funcgroup([character(len=20) :: 'mpw1lyp', 'mpw1-lyp'])
   funcs(p_pw6b95) = new_funcgroup([character(len=20) :: 'pw6b95'])
   funcs(p_tpssh) = new_funcgroup([character(len=20) :: 'tpssh'])
   funcs(p_tpss0) = new_funcgroup([character(len=20) :: 'tpss0'])
   funcs(p_x3lyp) = new_funcgroup([character(len=20) :: 'x3-lyp', 'x3lyp'])
   funcs(p_m06) = new_funcgroup([character(len=20) :: 'm06'])
   funcs(p_m06l) = new_funcgroup([character(len=20) :: 'm06l'])
   funcs(p_mn12sx) = new_funcgroup([character(len=20) :: 'mn12sx', 'mn12-sx'])
   funcs(p_b97d) = new_funcgroup([character(len=20) :: 'b97d'])
   funcs(p_lh07tsvwn) = new_funcgroup([character(len=20) :: 'lh07tsvwn', 'lh07t-svwn'])
   funcs(p_lh07ssvwn) = new_funcgroup([character(len=20) :: 'lh07ssvwn', 'lh07s-svwn'])
   funcs(p_lh12ctssirpw92) = new_funcgroup([character(len=20) :: 'lh12ctssirpw92', 'lh12ct-ssirpw92'])
   funcs(p_lh12ctssifpw92) = new_funcgroup([character(len=20) :: 'lh12ctssifpw92', 'lh12ct-ssifpw92'])
   funcs(p_lh14tcalpbe) = new_funcgroup([character(len=20) :: 'lh14tcalpbe', 'lh14t-calpbe'])
   funcs(p_lh20t) = new_funcgroup([character(len=20) :: 'lh20t'])
   funcs(p_b2plyp) = new_funcgroup([character(len=20) :: 'b2plyp', 'b2-plyp'])
   funcs(p_b2gpplyp) = new_funcgroup([character(len=20) :: 'b2gpplyp', 'b2gp-plyp'])
   funcs(p_mpw2plyp) = new_funcgroup([character(len=20) :: 'mpw2plyp'])
   funcs(p_pwpb95) = new_funcgroup([character(len=20) :: 'pwpb95'])
   funcs(p_dsdblyp) = new_funcgroup([character(len=20) :: 'dsdblyp', 'dsd-blyp'])
   funcs(p_dsdpbe) = new_funcgroup([character(len=20) :: 'dsdpbe', 'dsd-pbe'])
   funcs(p_dsdpbeb95) = new_funcgroup([character(len=20) :: 'dsdpbeb95', 'dsd-pbeb95'])
   funcs(p_dsdpbep86) = new_funcgroup([character(len=20) :: 'dsdpbep86', 'dsd-pbep86'])
   funcs(p_dsdsvwn) = new_funcgroup([character(len=20) :: 'dsdsvwn', 'dsd-svwn'])
   funcs(p_dodblyp) = new_funcgroup([character(len=20) :: 'dodblyp', 'dod-blyp'])
   funcs(p_dodpbe) = new_funcgroup([character(len=20) :: 'dodpbe', 'dod-pbe'])
   funcs(p_dodpbeb95) = new_funcgroup([character(len=20) :: 'dodpbeb95', 'dod-pbeb95'])
   funcs(p_dodpbep86) = new_funcgroup([character(len=20) :: 'dodpbep86', 'dod-pbep86'])
   funcs(p_dodsvwn) = new_funcgroup([character(len=20) :: 'dodsvwn', 'dod-svwn'])
   funcs(p_pbe0_2) = new_funcgroup([character(len=20) :: 'pbe02', 'pbe0-2'])
   funcs(p_pbe0_dh) = new_funcgroup([character(len=20) :: 'pbe0dh', 'pbe0-dh'])
   funcs(p_dftb_3ob) = new_funcgroup([character(len=20) :: 'dftb3', 'dftb(3ob)'])
   funcs(p_dftb_mio) = new_funcgroup([character(len=20) :: 'dftb(mio)'])
   funcs(p_dftb_pbc) = new_funcgroup([character(len=20) :: 'dftb(pbc)'])
   funcs(p_dftb_matsci) = new_funcgroup([character(len=20) :: 'dftb(matsci)'])
   funcs(p_dftb_ob2) = new_funcgroup([character(len=20) :: 'lc-dftb', 'dftb(ob2)'])
   funcs(p_b1b95) = new_funcgroup([character(len=20) :: 'b1b95'])
   funcs(p_pbesol) = new_funcgroup([character(len=20) :: 'pbesol'])
   funcs(p_mpwb1k) = new_funcgroup([character(len=20) :: 'mpwb1k'])
   funcs(p_mpw1b95) = new_funcgroup([character(len=20) :: 'mpw1b95'])
   funcs(p_hse03) = new_funcgroup([character(len=20) :: 'hse03'])
   funcs(p_hse06) = new_funcgroup([character(len=20) :: 'hse06'])
   funcs(p_hse12) = new_funcgroup([character(len=20) :: 'hse12'])
   funcs(p_hse12s) = new_funcgroup([character(len=20) :: 'hse12s'])
   funcs(p_hsesol) = new_funcgroup([character(len=20) :: 'hsesol'])
   funcs(p_revtpssh) = new_funcgroup([character(len=20) :: 'revtpssh'])
   funcs(p_glyp) = new_funcgroup([character(len=20) :: 'glyp', 'g-lyp'])
   funcs(p_revpbe0dh) = new_funcgroup([character(len=20) :: 'revpbe0dh', 'revpbe0-dh'])
   funcs(p_revtpss0) = new_funcgroup([character(len=20) :: 'revtpss0'])
   funcs(p_revdsdpbep86) = new_funcgroup([character(len=20) :: 'revdsd-pbep86', 'revdsdpbep86'])
   funcs(p_revdsdpbe) = new_funcgroup([character(len=20) :: 'revdsd-pbe', 'revdsd-pbepbe', 'revdsdpbe', 'revdsdpbepbe'])
   funcs(p_revdsdblyp) = new_funcgroup([character(len=20) :: 'revdsd-blyp', 'revdsdblyp'])
   funcs(p_revdodpbep86) = new_funcgroup([character(len=20) :: 'revdod-pbep86', 'revdodpbep86'])
   funcs(p_b97m) = new_funcgroup([character(len=20) :: 'b97m'])
   funcs(p_wb97m) = new_funcgroup([character(len=20) :: 'wb97m', 'ωb97m', 'omegab97m'])
   funcs(p_wb97m_rev) = new_funcgroup([character(len=20) :: 'wb97m-rev', &
      & 'ωb97m-rev', 'omegab97m-rev', 'wb97m_rev', 'ωb97m_rev', 'omegab97m_rev'])
   funcs(p_wb97) = new_funcgroup([character(len=20) :: 'wb97', 'ωb97', 'omegab97'])
   funcs(p_wb97x_2008) = new_funcgroup([character(len=20) :: 'wb97x_2008', &
      & 'ωb97x_2008', 'omegab97x_2008', 'wb97x-2008', 'ωb97x-2008', &
      & 'omegab97x-2008'])
   funcs(p_wb97x) = new_funcgroup([character(len=20) :: 'wb97x', 'ωb97x', &
      & 'omegab97x'])
   funcs(p_wb97x_rev) = new_funcgroup([character(len=20) :: 'wb97x-rev', &
      & 'ωb97x-rev', 'omegab97x-rev', 'wb97x_rev', 'ωb97x_rev', 'omegab97x_rev'])
   funcs(p_wb97x_3c) = new_funcgroup([character(len=20) :: 'wb97x-3c', &
      & 'ωb97x-3c', 'omegab97x-3c', 'wb97x_3c', 'ωb97x_3c', 'omegab97x_3c'])
   funcs(p_camb3lyp) = new_funcgroup([character(len=20) :: 'cam-b3lyp', 'camb3lyp'])
   funcs(p_camqtp01) = new_funcgroup([character(len=20) :: 'cam-qtp01', &
      & 'camqtp01', 'camqtp(01)', 'cam-qtp(01)'])
   funcs(p_lcblyp) = new_funcgroup([character(len=20) :: 'lc-blyp', 'lcblyp'])
   funcs(p_lcwpbe) = new_funcgroup([character(len=20) :: 'lc-wpbe', &
      & 'lcwpbe', 'lc-ωpbe', 'lcωpbe', 'lc-omegapbe', 'lcomegapbe'])
   funcs(p_lcwpbeh) = new_funcgroup([character(len=20) :: 'lc-wpbeh', &
      & 'lcwpbeh', 'lc-ωpbeh', 'lcωpbeh', 'lc-omegapbeh', 'lcomegapbeh'])

end subroutine get_functionals


!> Retrieve rational damping parameters from functional name
subroutine get_rational_damping_name(functional, param, s9)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_rational_damping_name

   !> Functional name for which to retrieve the damping parameters
   character(len=*), intent(in) :: functional

   !> Damping parameters for the functional
   class(damping_param), allocatable, intent(out) :: param

   !> Scaling factor for the three-body term
   real(wp), intent(in), optional :: s9

   character(len=:), allocatable :: fname
   integer :: is, id

   is = index(functional, '/')
   if (is == 0) is = len_trim(functional) + 1
   fname = lowercase(functional(:is-1))

   id = get_functional_id(fname)

   call get_rational_damping_id(id, param, s9=s9)

end subroutine get_rational_damping_name


!> Retrieve rational damping parameters from functional ID
subroutine get_rational_damping_id(id, param, s9)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_rational_damping_id

   !> Functional ID for which to retrieve the damping parameters
   integer, intent(in) :: id

   !> Damping parameters for the functional
   class(damping_param), allocatable, intent(out) :: param

   !> Scaling factor for the three-body term
   real(wp), intent(in), optional :: s9

   logical :: mbd

   mbd = .true.
   if (present(s9)) mbd = abs(s9) > epsilon(s9)

   if (mbd) then
      call get_d4eeq_bjatm_parameter(id, param, s9)
      if (.not.allocated(param)) then
         call get_d4eeq_bj_parameter(id, param, s9)
      end if
   else
      call get_d4eeq_bj_parameter(id, param, s9)
      if (.not.allocated(param)) then
         call get_d4eeq_bjatm_parameter(id, param, s9)
      end if
   end if

end subroutine get_rational_damping_id


subroutine get_d4eeq_bj_parameter(dfnum, param, s9)
   integer(df_enum), intent(in) :: dfnum
   class(damping_param), allocatable, intent(out) :: param
   real(wp), intent(in), optional :: s9
   select case(dfnum)
   case(p_dftb_3ob)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=0.4727337_wp, a1=0.5467502_wp, a2=4.4955068_wp)
   case(p_dftb_matsci)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=2.7711819_wp, a1=0.4681712_wp, a2=5.2918629_wp)
   case(p_dftb_mio)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=1.1948145_wp, a1=0.6074567_wp, a2=4.9336133_wp)
   case(p_dftb_ob2)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=2.7611320_wp, a1=0.6037249_wp, a2=5.3900004_wp)
   case(p_dftb_pbc)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=1.7303734_wp, a1=0.5546548_wp, a2=4.7973454_wp)
   end select

contains

   pure function dftd_param(s6, s8, a1, a2, alp) result(par)
      real(wp), intent(in) :: s8, a1, a2
      real(wp), intent(in), optional :: s6, alp
      type(rational_damping_param) :: par
      real(wp) :: s6_, alp_, s9_

      s6_ = 1.0_wp
      if (present(s6)) s6_ = s6
      s9_ = 0.0_wp
      if (present(s9)) s9_ = s9
      alp_ = 16.0_wp
      if (present(alp)) alp_ = alp

      par = rational_damping_param(&
         & s6=s6_, &
         & s8=s8, a1=a1, a2=a2, &
         & s9=s9_, &
         & alp=alp_)
   end function dftd_param

end subroutine get_d4eeq_bj_parameter

subroutine get_d4eeq_bjatm_parameter(dfnum, param, s9)
   integer(df_enum), intent(in) :: dfnum
   class(damping_param), allocatable, intent(out) :: param
   real(wp), intent(in), optional :: s9
   select case(dfnum)
   case(p_b1b95)
      param = dftd_param ( & ! (SAW190107)
         &  s6=1.0000_wp, s8=1.27701162_wp, a1=0.40554715_wp, a2=4.63323074_wp )
      !  Fitset: MD= 0.22852 MAD= 0.35189 RMSD= 0.46982
   case(p_b1lyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.98553711_wp, a1=0.39309040_wp, a2=4.55465145_wp )
      !  Fitset: MD= -0.04797 MAD= 0.25597 RMSD= 0.38778
   case(p_b1p)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.36115015_wp, a1=0.48665293_wp, a2=5.05219572_wp )
      !  Fitset: MD= -0.01406 MAD= 0.27441 RMSD= 0.47328
   case(p_b1pw)
      param = dftd_param ( & ! (SAW190107)
         &  s6=1.0000_wp, s8=3.02227550_wp, a1=0.47396846_wp, a2=4.49845309_wp )
      !  Fitset: MD= 0.10485 MAD= 0.32175 RMSD= 0.48508
   case(p_b2gpplyp)
      param = dftd_param ( & ! (SAW190107)
         &  s6=0.5600_wp, s8=0.94633372_wp, a1=0.42907301_wp, a2=5.18802602_wp )
      !  Fitset: MD= -0.05248 MAD= 0.18110 RMSD= 0.27365
   case(p_b2plyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.6400_wp, s8=1.16888646_wp, a1=0.44154604_wp, a2=4.73114642_wp )
      !  Fitset: MD= -0.03761 MAD= 0.18247 RMSD= 0.27109
   case(p_b3lyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.02929367_wp, a1=0.40868035_wp, a2=4.53807137_wp )
      !  Fitset: MD= -0.05892 MAD= 0.26117 RMSD= 0.40531
   case(p_b3p)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.08822155_wp, a1=0.47324238_wp, a2=4.98682134_wp )
      !  Fitset: MD= -0.02970 MAD= 0.26962 RMSD= 0.46761
   case(p_b3pw)
      param = dftd_param ( & ! (SAW190107)
         &  s6=1.0000_wp, s8=2.88364295_wp, a1=0.46990860_wp, a2=4.51641422_wp )
      !  Fitset: MD= 0.06643 MAD= 0.29151 RMSD= 0.45541
   case(p_b97)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.87854260_wp, a1=0.29319126_wp, a2=4.51647719_wp )
      !  Fitset: MD= -0.13017 MAD= 0.24778 RMSD= 0.36116
   case(p_bhlyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.65281646_wp, a1=0.27263660_wp, a2=5.48634586_wp )
      !  Fitset: MD= -0.15832 MAD= 0.34132 RMSD= 0.57342
   case(p_blyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.34076671_wp, a1=0.44488865_wp, a2=4.09330090_wp )
      !  Fitset: MD= 0.04801 MAD= 0.28161 RMSD= 0.38321
   case(p_bpbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.64405246_wp, a1=0.52905620_wp, a2=4.11311891_wp )
      !  Fitset: MD= 0.19316 MAD= 0.41912 RMSD= 0.60452
   case(p_bp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.35497927_wp, a1=0.43645861_wp, a2=4.92406854_wp )
      !  Fitset: MD= 0.08252 MAD= 0.32681 RMSD= 0.47063
   case(p_bpw)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.24571506_wp, a1=0.50050454_wp, a2=4.12346483_wp )
      !  Fitset: MD= 0.20607 MAD= 0.41941 RMSD= 0.59589
   case(p_camb3lyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.66041301_wp, a1=0.40267156_wp, a2=5.17432195_wp )
      !  Fitset: MD= -0.19675 MAD= 0.34901 RMSD= 0.59087
   case(p_camqtp01)
      param = dftd_param ( & ! (10.1021/acs.jctc.3c00717)
         &  s6=1.0000_wp, s8=1.156_wp, a1=0.461_wp, a2=6.375_wp )
   case(p_dodblyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.4700_wp, s8=1.31146043_wp, a1=0.43407294_wp, a2=4.27914360_wp )
      !  Fitset: MD= 0.03323 MAD= 0.13858 RMSD= 0.20861
   case(p_dodpbeb95)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.5600_wp, s8=0.01574635_wp, a1=0.43745720_wp, a2=3.69180763_wp )
      !  Fitset: MD= 0.03704 MAD= 0.13343 RMSD= 0.18278
   case(p_dodpbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.4800_wp, s8=0.92051454_wp, a1=0.43037052_wp, a2=4.38067238_wp )
      !  Fitset: MD= 0.01065 MAD= 0.13414 RMSD= 0.21424
   case(p_dodpbep86)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.4600_wp, s8=0.71405681_wp, a1=0.42408665_wp, a2=4.52884439_wp )
      !  Fitset: MD= -0.03740 MAD= 0.12467 RMSD= 0.18127
   case(p_dodsvwn)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.4200_wp, s8=0.94500207_wp, a1=0.47449026_wp, a2=5.05316093_wp )
      !  Fitset: MD= -0.07427 MAD= 0.16970 RMSD= 0.25286
   case(p_dsdblyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.5400_wp, s8=0.63018237_wp, a1=0.47591835_wp, a2=4.73713781_wp )
      !  Fitset: MD= -0.01981 MAD= 0.14823 RMSD= 0.21530
   case(p_dsdpbeb95)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.5400_wp, s8=-0.14668670_wp, a1=0.46394587_wp, a2=3.64913860_wp )
      !  Fitset: MD= 0.02996 MAD= 0.12414 RMSD= 0.16860
   case(p_dsdpbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.4500_wp, s8=0.70584116_wp, a1=0.45787085_wp, a2=4.44566742_wp )
      !  Fitset: MD= 0.00866 MAD= 0.13406 RMSD= 0.21380
   case(p_dsdpbep86)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.4700_wp, s8=0.37586675_wp, a1=0.53698768_wp, a2=5.13022435_wp )
      !  Fitset: MD= -0.05273 MAD= 0.14259 RMSD= 0.21271
   case(p_dsdsvwn)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.4100_wp, s8=0.72914436_wp, a1=0.51347412_wp, a2=5.11858541_wp )
      !  Fitset: MD= -0.08974 MAD= 0.32285 RMSD= 0.43146
   case(p_glyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=4.23798924_wp, a1=0.38426465_wp, a2=4.38412863_wp )
      !  Fitset: MD= 0.63466 MAD= 0.89568 RMSD= 1.11309
   case(p_hf)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.61679827_wp, a1=0.44959224_wp, a2=3.35743605_wp )
      !  Fitset: MD= -0.02597 MAD= 0.34732 RMSD= 0.49719
   case(p_lb94)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.59538499_wp, a1=0.42088944_wp, a2=3.28193223_wp )
      !  Fitset: MD= 0.31701 MAD= 0.53196 RMSD= 0.74553
   case(p_lcblyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.60344180_wp, a1=0.45769839_wp, a2=7.86924893_wp )
      !  Fitset: MD= -0.39724 MAD= 0.72327 RMSD= 1.18218
   case(p_lcwpbe)
      param = dftd_param ( & ! (10.1021/acs.jctc.3c00717)
         &  s6=1.0000_wp, s8=1.170_wp, a1=0.378_wp, a2=4.816_wp )
   case(p_lcwpbeh)
      param = dftd_param ( & ! (10.1021/acs.jctc.3c00717)
         &  s6=1.0000_wp, s8=1.318_wp, a1=0.386_wp, a2=5.010_wp )
   case(p_lh07ssvwn)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.16675531_wp, a1=0.35965552_wp, a2=4.31947614_wp )
      !  Fitset: MD= 0.32224 MAD= 0.59006 RMSD= 0.86272
   case(p_lh07tsvwn)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.09333001_wp, a1=0.35025189_wp, a2=4.34166515_wp )
      !  Fitset: MD= 0.24243 MAD= 0.43497 RMSD= 0.61671
   case(p_lh12ctssifpw92)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.68467610_wp, a1=0.34190416_wp, a2=3.91039666_wp )
      !  Fitset: MD= 0.55106 MAD= 0.80783 RMSD= 1.11048
   case(p_lh12ctssirpw92)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.48973402_wp, a1=0.34026075_wp, a2=3.96948081_wp )
      !  Fitset: MD= 0.47785 MAD= 0.71188 RMSD= 0.98422
   case(p_lh14tcalpbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.28130770_wp, a1=0.38822021_wp, a2=4.92501211_wp )
      !  Fitset: MD= -0.02105 MAD= 0.22968 RMSD= 0.36045
   case(p_lh20t)
      param = dftd_param ( & ! (10.1021/acs.jctc.0c00498)
         & s6=1.000_wp, s8=0.113_wp, a1=0.479_wp, a2=4.635_wp )
   case(p_m06)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.16366729_wp, a1=0.53456413_wp, a2=6.06192174_wp )
      !  Fitset: MD= 0.01788 MAD= 0.24914 RMSD= 0.38604
   case(p_m06l)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.59493760_wp, a1=0.71422359_wp, a2=6.35314182_wp )
      !  Fitset: MD= 0.08395 MAD= 0.24888 RMSD= 0.34879
   case(p_mn12sx)
      param = dftd_param ( & ! (SAW211021)
         &  s6=1.0000_wp, s8=0.85964873_wp, a1=0.62662681_wp, a2=5.62088906_wp )
      !  Fitset: MD= 0.16131 MAD= 0.34142 RMSD= 0.47113
   case(p_mpw1b95)
      param = dftd_param ( & ! (SAW190107)
         &  s6=1.0000_wp, s8=0.50093024_wp, a1=0.41585097_wp, a2=4.99154869_wp )
      !  Fitset: MD= 0.00585 MAD= 0.15695 RMSD= 0.21297
   case(p_mpw1lyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.15591153_wp, a1=0.25603493_wp, a2=5.32083895_wp )
      !  Fitset: MD= -0.26979 MAD= 0.41542 RMSD= 0.60678
   case(p_mpw1pw)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.80841716_wp, a1=0.42961819_wp, a2=4.68892341_wp )
      !  Fitset: MD= -0.08840 MAD= 0.26815 RMSD= 0.45231
   case(p_mpw2plyp)
      param = dftd_param ( & ! (SAW190107)
         &  s6=0.7500_wp, s8=0.45788846_wp, a1=0.42997704_wp, a2=5.07650682_wp )
      !  Fitset: MD= -0.18921 MAD= 0.30115 RMSD= 0.44049
   case(p_mpwb1k)
      param = dftd_param ( & ! (SAW190107)
         &  s6=1.0000_wp, s8=0.57338313_wp, a1=0.44687975_wp, a2=5.21266777_wp )
      !  Fitset: MD= -0.00870 MAD= 0.17226 RMSD= 0.23614
   case(p_mpwlyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.25842942_wp, a1=0.25773894_wp, a2=5.02319542_wp )
      !  Fitset: MD= -0.24426 MAD= 0.39145 RMSD= 0.54503
   case(p_mpwpw)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.82596836_wp, a1=0.34526745_wp, a2=4.84620734_wp )
      !  Fitset: MD= -0.06278 MAD= 0.27913 RMSD= 0.43988
   case(p_o3lyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.75762508_wp, a1=0.10348980_wp, a2=6.16233282_wp )
      !  Fitset: MD= -0.19268 MAD= 0.38577 RMSD= 0.62168
   case(p_olyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.74836820_wp, a1=0.60184498_wp, a2=2.53292167_wp )
      !  Fitset: MD= 0.12352 MAD= 0.37113 RMSD= 0.58291
   case(p_opbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.06917417_wp, a1=0.68267534_wp, a2=2.22849018_wp )
      !  Fitset: MD= 0.26699 MAD= 0.55308 RMSD= 0.85023
   case(p_pbe0_2)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.5000_wp, s8=0.64299082_wp, a1=0.76542115_wp, a2=5.78578675_wp )
      !  Fitset: MD= -0.04260 MAD= 0.21186 RMSD= 0.34045
   case(p_pbe0)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.20065498_wp, a1=0.40085597_wp, a2=5.02928789_wp )
      !  Fitset: MD= -0.17892 MAD= 0.30557 RMSD= 0.51050
   case(p_pbe0_dh)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.8750_wp, s8=0.96811578_wp, a1=0.47592488_wp, a2=5.08622873_wp )
      !  Fitset: MD= -0.13857 MAD= 0.27919 RMSD= 0.47256
   case(p_pbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.95948085_wp, a1=0.38574991_wp, a2=4.80688534_wp )
      !  Fitset: MD= -0.20544 MAD= 0.33635 RMSD= 0.51168
   case(p_pbesol)
      param = dftd_param ( & ! (SAW211021)
         &  s6=1.0000_wp, s8=1.71885698_wp, a1=0.47901421_wp, a2=5.96771589_wp )
      !  Fitset: MD= -0.28899 MAD= 0.52215 RMSD= 0.93584
   case(p_am05)
      param = dftd_param ( & ! (SAW211021)
         &  s6=1.0000_wp, s8=1.71885838_wp, a1=0.47901431_wp, a2=5.96771581_wp )
      !  Fitset: MD= -0.28899 MAD= 0.52215 RMSD= 0.93584
   case(p_pw1pw)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.96850170_wp, a1=0.42427511_wp, a2=5.02060636_wp )
      !  Fitset: MD= -0.27325 MAD= 0.42206 RMSD= 0.64119
   case(p_pw6b95)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=-0.31926054_wp, a1=0.04142919_wp, a2=5.84655608_wp )
      !  Fitset: MD= -0.04767 MAD= 0.14330 RMSD= 0.18958
   case(p_pw86pbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.21362856_wp, a1=0.40510366_wp, a2=4.66737724_wp )
      !  Fitset: MD= -0.11505 MAD= 0.24691 RMSD= 0.38101
   case(p_pw91)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.77283111_wp, a1=0.39581542_wp, a2=4.93405761_wp )
      !  Fitset: MD= -0.33019 MAD= 0.48611 RMSD= 0.68110
   case(p_pwp1)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.60492565_wp, a1=0.46855837_wp, a2=5.76921413_wp )
      !  Fitset: MD= -0.35321 MAD= 0.54026 RMSD= 0.86629
   case(p_pwpb95)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.8200_wp, s8=-0.34639127_wp, a1=0.41080636_wp, a2=3.83878274_wp )
      !  Fitset: MD= 0.02143 MAD= 0.13040 RMSD= 0.17599
   case(p_pwp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.32801227_wp, a1=0.35874687_wp, a2=6.05861168_wp )
      !  Fitset: MD= -0.42482 MAD= 0.62607 RMSD= 0.91840
   case(p_revpbe0)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.57185414_wp, a1=0.38705966_wp, a2=4.11028876_wp )
      !  Fitset: MD= 0.02724 MAD= 0.21587 RMSD= 0.36040
   case(p_revpbe0dh)
      param = dftd_param ( & ! (SAW190103)
         &  s6=0.8750_wp, s8=1.24456037_wp, a1=0.36730560_wp, a2=4.71126482_wp )
      !  Fitset: MD= -0.01089 MAD= 0.20910 RMSD= 0.33564
   case(p_revpbe38)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.66597472_wp, a1=0.39476833_wp, a2=4.39026628_wp )
      !  Fitset: MD= -0.01326 MAD= 0.22598 RMSD= 0.36210
   case(p_revpbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.74676530_wp, a1=0.53634900_wp, a2=3.07261485_wp )
      !  Fitset: MD= 0.05649 MAD= 0.25212 RMSD= 0.40863
   case(p_revtpss0)
      param = dftd_param ( & ! (SAW190107)
         &  s6=1.0000_wp, s8=1.54664499_wp, a1=0.45890964_wp, a2=4.78426405_wp )
      !  Fitset: MD= -0.05298 MAD= 0.19965 RMSD= 0.32081
   case(p_revtpss)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.53089454_wp, a1=0.44880597_wp, a2=4.64042317_wp )
      !  Fitset: MD= -0.01904 MAD= 0.19568 RMSD= 0.29618
   case(p_revtpssh)
      param = dftd_param ( & ! (SAW190107)
         &  s6=1.0000_wp, s8=1.52740307_wp, a1=0.45161957_wp, a2=4.70779483_wp )
      !  Fitset: MD= -0.03731 MAD= 0.19133 RMSD= 0.29091
   case(p_rpbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.31183787_wp, a1=0.46169493_wp, a2=3.15711757_wp )
      !  Fitset: MD= -0.07156 MAD= 0.26348 RMSD= 0.38671
   case(p_rpw86pbe)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.12624034_wp, a1=0.38151218_wp, a2=4.75480472_wp )
      !  Fitset: MD= -0.12740 MAD= 0.26294 RMSD= 0.40614
   case(p_scan)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.46126056_wp, a1=0.62930855_wp, a2=6.31284039_wp )
      !  Fitset: MD= -0.13170 MAD= 0.28640 RMSD= 0.51183
   case(p_rscan)
      param = dftd_param ( & ! (10.1063/5.0041008)
         &  s6=1.0000_wp, s8=0.87728975_wp, a1=0.49116966_wp, a2=5.75859346_wp )
   case(p_r2scan)
      param = dftd_param ( & ! (10.1063/5.0041008)
         &  s6=1.0000_wp, s8=0.60187490_wp, a1=0.51559235_wp, a2=5.77342911_wp )
   case(p_r2scanh)
      param = dftd_param ( & ! (10.1063/5.0086040)
         & s6=1.0_wp, s8=0.8324_wp, a1=0.4944_wp, a2=5.9019_wp)
   case(p_r2scan0)
      param = dftd_param ( & ! (10.1063/5.0086040)
         & s6=1.0_wp, s8=0.8992_wp, a1=0.4778_wp, a2=5.8779_wp)
   case(p_r2scan50)
      param = dftd_param ( & ! (10.1063/5.0086040)
         & s6=1.0_wp, s8=1.0471_wp, a1=0.4574_wp, a2=5.8969_wp)
   case(p_r2scan_3c)
      param = dftd_param ( & ! (10.1063/5.0040021)
         & s6=1.0_wp, s8=0.00_wp, a1=0.42_wp, a2=5.65_wp)
   case(p_tpss0)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.62438102_wp, a1=0.40329022_wp, a2=4.80537871_wp )
      !  Fitset: MD= -0.09569 MAD= 0.26733 RMSD= 0.44767
   case(p_tpss)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.76596355_wp, a1=0.42822303_wp, a2=4.54257102_wp )
      !  Fitset: MD= -0.09296 MAD= 0.27505 RMSD= 0.42537
   case(p_tpssh)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.85897750_wp, a1=0.44286966_wp, a2=4.60230534_wp )
      !  Fitset: MD=  0.02238 MAD= 0.16042 RMSD= 0.33519
   case(p_b97d)
      param = dftd_param ( & ! (SAW201029)
         &  s6=1.0000_wp, s8=1.69460052_wp, a1=0.28904684_wp, a2=4.13407323_wp )
      !  Fitset: MD= -0.09858 MAD= 0.26757 RMSD= 0.42380
   case(p_wb97)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=6.55792598_wp, a1=0.76666802_wp, a2=8.36027334_wp )
      !  Fitset: MD= -0.12779 MAD= 0.36152 RMSD= 0.49991
   case(p_wb97x_2008)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=-0.07519516_wp, a1=0.45094893_wp, a2=6.78425255_wp )
      !  S22x5: MD= 0.05 MAD= 0.16 RMSD= 0.22
      !  S66x8: MD= 0.06 MAD= 0.16 RMSD= 0.21
      !  NCI10: MD= 0.08 MAD= 0.15 RMSD= 0.25
   case(p_wb97x)
      param = dftd_param ( & ! (10.1002/jcc.26411)
         &  s6=1.0000_wp, s8=0.5093_wp, a1=0.0662_wp, a2=5.4487_wp )
   case(p_wb97x_rev)
      param = dftd_param ( & ! (10.1063/5.0133026)
         &  s6=1.0000_wp, s8=0.4485_wp, a1=0.3306_wp, a2=4.279_wp )
   case(p_wb97x_3c)
      param = dftd_param ( & ! (10.1063/5.0133026)
         &  s6=1.0000_wp, s8=0.0_wp, a1=0.2464_wp, a2=4.737_wp )
   case(p_b97m)
      param = dftd_param ( & ! (10.1002/jcc.26411)
         &  s6=1.0000_wp, s8=0.6633_wp, a1=0.4288_wp, a2=3.9935_wp )
      !  S22x5: MD= 0.03 MAD= 0.12 RMSD= 0.18
      !  S66x8: MD= 0.09 MAD= 0.17 RMSD= 0.22
      !  NCI10: MD= 0.09 MAD= 0.15 RMSD= 0.32
   case(p_wb97m)
      param = dftd_param ( & ! (10.1002/jcc.26411)
         &  s6=1.0000_wp, s8=0.7761_wp, a1=0.7514_wp, a2=2.7099_wp )
      !  Fitset: MD= -0.20216 MAD= 0.34696 RMSD= 0.53641
   case(p_wb97m_rev) 
      param = dftd_param ( & ! (10.1021/acs.jctc.3c00717)
         &  s6=1.0000_wp, s8=0.842_wp, a1=0.359_wp, a2=4.668_wp )
   case(p_x3lyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.54701429_wp, a1=0.20318443_wp, a2=5.61852648_wp )
      !  Fitset: MD= -0.15607 MAD= 0.31342 RMSD= 0.49546
   case(p_xlyp)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.62972054_wp, a1=0.11268673_wp, a2=5.40786417_wp )
      !  Fitset: MD= -0.03900 MAD= 0.27562 RMSD= 0.38491
   case(p_revdsdpbep86)
      param = dftd_param ( & ! (WTMAD2)
         &  s6=0.5132_wp, s8=0.00000000_wp, a1=0.44000000_wp, a2=3.60000000_wp )
   case(p_revdsdpbe)
      param = dftd_param ( & ! (WTMAD2)
         &  s6=0.6706_wp, s8=0.00000000_wp, a1=0.40000000_wp, a2=3.60000000_wp )
   case(p_revdsdblyp)
      param = dftd_param ( & !(WTMAD2)
         &  s6=0.6141_wp, s8=0.00000000_wp, a1=0.38000000_wp, a2=3.52000000_wp )
   case(p_revdodpbep86)
      param = dftd_param ( & !(WTMAD2)
         &  s6=0.5552_wp, s8=0.00000000_wp, a1=0.44000000_wp, a2=3.60000000_wp )
   case(p_dftb_3ob)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=0.6635015_wp, a1=0.5523240_wp, a2=4.3537076_wp)
   case(p_dftb_matsci)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=3.3157614_wp, a1=0.4826330_wp, a2=5.3811976_wp)
   case(p_dftb_mio)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=1.2916225_wp, a1=0.5965326_wp, a2=4.8778602_wp)
   case(p_dftb_ob2)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=2.9692689_wp, a1=0.6068916_wp, a2=5.4476789_wp)
   case(p_dftb_pbc)
      param = dftd_param( & ! (SAW191202)
         &  s6=1.0_wp, s8=2.1667394_wp, a1=0.5646391_wp, a2=4.9576353_wp)
   case(p_hse03)
      param = dftd_param( & ! (SAW211107)
         &  s6=1.0_wp, s8=1.19812280_wp, a1=0.38662939_wp, a2=5.22925796_wp)
   case(p_hse06)
      param = dftd_param( & ! (SAW211107)
         &  s6=1.0_wp, s8=1.19528249_wp, a1=0.38663183_wp, a2=5.19133469_wp)
   case(p_hse12)
      param = dftd_param( & ! (SAW211107)
         &  s6=1.0_wp, s8=1.23500792_wp, a1=0.39226921_wp, a2=5.22036266_wp)
   case(p_hse12s)
      param = dftd_param( & ! (SAW211107)
         &  s6=1.0_wp, s8=1.23767762_wp, a1=0.39989137_wp, a2=5.34809245_wp)
   case(p_hsesol)
      param = dftd_param( & ! (SAW211107)
         &  s6=1.0_wp, s8=1.82207807_wp, a1=0.45646268_wp, a2=5.59662251_wp)
   case(p_wr2scan) ! (10.1063/5.0174988)
      param = dftd_param ( &
         & s6=1.0_wp, s8=1.0_wp, a1=0.3834_wp, a2=5.7889_wp)
   case(p_r2scan0_dh) ! (10.1063/5.0174988)
      param = dftd_param ( &
         & s6=0.9424_wp, s8=0.3856_wp, a1=0.4271_wp, a2=5.8565_wp)
   case(p_r2scan_cidh) ! (10.1063/5.0174988)
      param = dftd_param ( &
         & s6=0.8666_wp, s8=0.5336_wp, a1=0.4171_wp, a2=5.9125_wp)
   case(p_r2scan_qidh) ! (10.1063/5.0174988)
      param = dftd_param ( &
         & s6=0.7867_wp, s8=0.2955_wp, a1=0.4001_wp, a2=5.8300_wp)
   case(p_r2scan0_2) ! (10.1063/5.0174988)
      param = dftd_param ( &
         & s6=0.7386_wp, s8=0.0000_wp, a1=0.4030_wp, a2=5.5142_wp)
   case(p_pr2scan50) ! (10.1063/5.0174988)
      param = dftd_param ( &
         & s6=0.7964_wp, s8=0.3421_wp, a1=0.4663_wp, a2=5.7916_wp)
   case(p_pr2scan69) ! (10.1063/5.0174988)
      param = dftd_param ( &
         & s6=0.7167_wp, s8=0.0000_wp, a1=0.4644_wp, a2=5.2563_wp)
   case(p_kpr2scan50) ! (10.1063/5.0174988)
      param = dftd_param ( &
         & s6=0.8402_wp, s8=0.1212_wp, a1=0.4382_wp, a2=5.8232_wp)
   case(p_wpr2scan50) ! (10.1063/5.0174988)
      param = dftd_param ( &
         & s6=0.8143_wp, s8=0.3842_wp, a1=0.4135_wp, a2=5.8773_wp)
   end select

contains

   pure function dftd_param(s6, s8, a1, a2, alp) result(par)
      real(wp), intent(in) :: s8, a1, a2
      real(wp), intent(in), optional :: s6, alp
      type(rational_damping_param) :: par
      real(wp) :: s6_, alp_, s9_

      s6_ = 1.0_wp
      if (present(s6)) s6_ = s6
      s9_ = 1.0_wp
      if (present(s9)) s9_ = s9
      alp_ = 16.0_wp
      if (present(alp)) alp_ = alp

      par = rational_damping_param(&
         & s6=s6_, &
         & s8=s8, a1=a1, a2=a2, &
         & s9=s9_, &
         & alp=alp_)
   end function dftd_param

end subroutine get_d4eeq_bjatm_parameter

!> Get the unique identifier for most functionals, returns none if
!> the functional was not known at the time I implemented this mapping
pure function get_functional_id(df) result(num)
   integer(df_enum) :: num
   character(len=*), intent(in) :: df
   select case(df)
   case default
      num = p_invalid
   case('hf')
      num = p_hf
   case('am05', 'gga_x_am05:gga_c_am05')
      num = p_am05
   case('b-lyp', 'blyp', 'gga_x_b88:gga_c_lyp')
      num = p_blyp
   case('bpbe', 'gga_x_b88:gga_c_pbe')
      num = p_bpbe
   case('b-p', 'bp86', 'bp', 'b-p86', 'gga_x_b88:gga_c_p86')
      num = p_bp
   case('bpw', 'b-pw', 'gga_x_b88:gga_c_pw91')
      num = p_bpw
   case('lb94', 'gga_x_lb') ! no gga_c_lb
      num = p_lb94
   case('mpwlyp', 'mpw-lyp', 'gga_x_mpw91:gga_c_lyp')
      num = p_mpwlyp
   case('mpwpw', 'mpw-pw', 'mpwpw91', 'gga_x_mpw91:gga_c_pw91')
      num = p_mpwpw
   case('o-lyp', 'olyp', 'gga_x_optx:gga_c_lyp')
      num = p_olyp
   case('opbe', 'gga_x_optx:gga_c_pbe')
      num = p_opbe
   case('pbe', 'gga_x_pbe:gga_c_pbe')
      num = p_pbe
   case('rpbe', 'gga_x_rpbe:gga_c_pbe')
      num = p_rpbe
   case('revpbe', 'gga_x_pbe_r:gga_c_pbe')
      num = p_revpbe
   case('pbesol', 'gga_x_pbe_sol:gga_c_pbe_sol')
      num = p_pbesol
   case('pw86pbe', 'gga_x_pw86:gga_c_pbe')
      num = p_pw86pbe
   case('rpw86pbe', 'gga_x_rpw86:gga_c_pbe')
      num = p_rpw86pbe
   case('pw91', 'gga_x_pw91:gga_c_pw91')
      num = p_pw91
   case('pwp', 'pw-p', 'pw91p86', 'gga_x_pw91:gga_c_p86')
      num = p_pwp
   case('x-lyp', 'xlyp', 'gga_xc_xlyp')
      num = p_xlyp
   case('b97', 'hyb_gga_xc_b97')
      num = p_b97
   case('b97d', 'gga_xc_b97_d')
      num = p_b97d
   case('tpss', 'mgga_c_tpss:mgga_x_tpss')
      num = p_tpss
   case('revtpss', 'mgga_c_revtpss:mgga_x_revtpss')
      num = p_revtpss
   case('scan', 'mgga_x_scan:mgga_c_scan')
      num = p_scan
   case('rscan', 'mgga_x_rscan:mgga_c_rscan')
      num = p_rscan
   case('r2scan', 'r²scan', 'mgga_x_r2scan:mgga_c_r2scan')
      num = p_r2scan
   case('r2scanh', 'r²scanh', 'hyb_mgga_xc_r2scanh')
      num = p_r2scanh
   case('r2scan0', 'r²scan0', 'hyb_mgga_xc_r2scan0')
      num = p_r2scan0
   case('r2scan50', 'r²scan50', 'hyb_mgga_xc_r2scan50')
      num = p_r2scan50
   case('r2scan-3c', 'r²scan-3c', 'r2scan_3c', 'r²scan_3c', 'r2scan3c')
      num = p_r2scan_3c
   case('b1lyp', 'b1-lyp', 'hyb_gga_xc_b1lyp')
      num = p_b1lyp
   case('b3-lyp', 'b3lyp', 'hyb_gga_xc_b3lyp', 'hyb_gga_xc_b3lyp3', 'hyb_gga_xc_b3lyp5')
      num = p_b3lyp
   case('bh-lyp', 'bhlyp', 'hyb_gga_xc_bhandh', 'hyb_gga_xc_bhandhlyp')
      num = p_bhlyp
   case('b1p', 'b1-p', 'b1p86') ! 0.75 b88 + 0.25 hf; p86 (nonloc) + pw81 (loc)
      num = p_b1p
   case('b3p', 'b3-p', 'b3p86', 'hyb_gga_xc_b3p86', 'hyb_gga_xc_b3p86_nwchem')
      num = p_b3p
   case('b1pw', 'b1-pw', 'b1pw91', 'hyb_gga_xc_b1pw91')
      num = p_b1pw
   case('b3pw', 'b3-pw', 'b3pw91', 'hyb_gga_xc_b3pw91')
      num = p_b3pw
   case('o3-lyp', 'o3lyp', 'hyb_gga_xc_o3lyp')
      num = p_o3lyp
   case('revpbe0') ! no libxc
      num = p_revpbe0
   case('revpbe38') ! no libxc
      num = p_revpbe38
   case('pbe0', 'hyb_gga_xc_pbeh')
      num = p_pbe0
   case('pwp1') ! no libxc
      num = p_pwp1
   case('pw1pw', 'pw1-pw') ! no libxc
      num = p_pw1pw
   case('mpw1pw', 'mpw1-pw', 'mpw1pw91', 'hyb_gga_xc_mpw1pw')
      num = p_mpw1pw
   case('mpw1lyp', 'mpw1-lyp', 'hyb_gga_xc_mpw1lyp')
      num = p_mpw1lyp
   case('pw6b95', 'hyb_mgga_xc_pw6b95')
      num = p_pw6b95
   case('tpssh', 'hyb_mgga_xc_tpssh')
      num = p_tpssh
   case('tpss0', 'hyb_mgga_xc_tpss0')
      num = p_tpss0
   case('x3-lyp', 'x3lyp', 'hyb_gga_xc_x3lyp')
      num = p_x3lyp
   case('m06', 'mgga_x_m06:mgga_c_m06')
      num = p_m06
   case('m06l', 'mgga_x_m06_l:mgga_c_m06_l')
      num = p_m06l
   case('mn12sx', 'mn12-sx', 'mgga_c_mn12_sx:mgga_c_mn12_sx')
      num = p_mn12sx
   case('cam-b3lyp', 'camb3lyp', 'hyb_gga_xc_cam_b3lyp')
      num = p_camb3lyp
   case('cam-qtp01', 'camqtp01', 'camqtp(01)', 'cam-qtp(01)', &
      & 'hyb_gga_xc_cam_qtp_01')
      num = p_camqtp01
   case('lc-blyp', 'lcblyp', 'hyb_gga_xc_lc_blyp')
      num = p_lcblyp
   case('lc-wpbe', 'lcwpbe', 'lc-ωpbe', 'lcωpbe', 'lc-omegapbe', 'lcomegapbe', &
      & 'hyb_gga_xc_lc_wpbe', 'hyb_gga_xc_lc_wpbe08_whs', &
      & 'hyb_gga_xc_lc_wpbe_whs', 'hyb_gga_xc_lrc_wpbe')
      num = p_lcwpbe
   case('lc-wpbeh', 'lcwpbeh', 'lc-ωpbeh', 'lcωpbeh', 'lc-omegapbeh', &
      & 'lcomegapbeh', 'hyb_gga_xc_lc_wpbeh_whs', 'hyb_gga_xc_lrc_wpbeh')
      num = p_lcwpbeh
   case('lh07tsvwn', 'lh07t-svwn') ! no libxc
      num = p_lh07tsvwn
   case('lh07ssvwn', 'lh07s-svwn') ! no libxc
      num = p_lh07ssvwn
   case('lh12ctssirpw92', 'lh12ct-ssirpw92') ! no libxc
      num = p_lh12ctssirpw92
   case('lh12ctssifpw92', 'lh12ct-ssifpw92') ! no libxc
      num = p_lh12ctssifpw92
   case('lh14tcalpbe', 'lh14t-calpbe') ! no libxc
      num = p_lh14tcalpbe
   case('lh20t') ! no libxc
      num = p_lh20t
   case('b2plyp', 'b2-plyp', 'xc_hyb_gga_xc_b2plyp') ! only in code
      num = p_b2plyp
   case('b2gpplyp', 'b2gp-plyp', 'xc_hyb_gga_xc_b2gpplyp') ! only in code
      num = p_b2gpplyp
   case('mpw2plyp') ! no libxc
      num = p_mpw2plyp
   case('pwpb95') ! no libxc
      num = p_pwpb95
   case('dsdblyp', 'dsd-blyp') ! no libxc
      num = p_dsdblyp
   case('dsdpbe', 'dsd-pbe') ! no libxc
      num = p_dsdpbe
   case('dsdpbeb95', 'dsd-pbeb95') ! no libxc
      num = p_dsdpbeb95
   case('dsdpbep86', 'dsd-pbep86') ! no libxc
      num = p_dsdpbep86
   case('dsdsvwn', 'dsd-svwn') ! no libxc
      num = p_dsdsvwn
   case('dodblyp', 'dod-blyp') ! no libxc
      num = p_dodblyp
   case('dodpbe', 'dod-pbe') ! no libxc
      num = p_dodpbe
   case('dodpbeb95', 'dod-pbeb95') ! no libxc
      num = p_dodpbeb95
   case('dodpbep86', 'dod-pbep86') ! no libxc
      num = p_dodpbep86
   case('dodsvwn', 'dod-svwn') ! no libxc
      num = p_dodsvwn
   case('pbe02', 'pbe0-2') ! no libxc
      num = p_pbe0_2
   case('pbe0dh', 'pbe0-dh') ! no libxc
      num = p_pbe0_dh
   case('dftb3', 'dftb(3ob)') ! no libxc
      num = p_dftb_3ob
   case('dftb(mio)') ! no libxc
      num = p_dftb_mio
   case('dftb(pbc)') ! no libxc
      num = p_dftb_pbc
   case('dftb(matsci)') ! no libxc
      num = p_dftb_matsci
   case('lc-dftb', 'dftb(ob2)') ! no libxc
      num = p_dftb_ob2
   case('b1b95', 'hyb_mgga_xc_b88b95')
      num = p_b1b95
   case('mpwb1k', 'hyb_mgga_xc_mpwb1k')
      num = p_mpwb1k
   case('mpw1b95', 'hyb_mgga_xc_mpw1b95')
      num = p_mpw1b95
   case('hse03', 'hyb_gga_xc_hse03')
      num = p_hse03
   case('hse06', 'hyb_gga_xc_hse06')
      num = p_hse06
   case('hse12', 'hyb_gga_xc_hse12')
      num = p_hse12
   case('hse12s', 'hyb_gga_xc_hse12s')
      num = p_hse12s
   case('hsesol', 'hyb_gga_xc_hse_sol')
      num = p_hsesol
   case('glyp', 'g-lyp', 'gga_x_g96:gga_c_lyp')
      num = p_glyp
   case('revpbe0dh', 'revpbe0-dh') ! no libxc
      num = p_revpbe0dh
   case('revtpssh', 'hyb_mgga_xc_revtpssh')
      num = p_revtpssh
    case('revtpss0') ! no libxc
      num = p_revtpss0
   case('revdsd-pbep86', 'revdsdpbep86') ! no libxc
      num = p_revdsdpbep86
   case('revdsd-pbe', 'revdsd-pbepbe', 'revdsdpbe', 'revdsdpbepbe') ! no libxc
      num = p_revdsdpbe
   case('revdsd-blyp', 'revdsdblyp') ! no libxc
      num = p_revdsdblyp
   case('revdod-pbep86', 'revdodpbep86') ! no libxc
      num = p_revdodpbep86
   case('b97m', 'mgga_xc_b97m_v')
      num = p_b97m
   case('wb97m', 'ωb97m', 'omegab97m', 'hyb_mgga_xc_wb97m_v')
      num = p_wb97m
   case('wb97m-rev', 'ωb97m-rev', 'omegab97m-rev', 'wb97m_rev', 'ωb97m_rev', &
      & 'omegab97m_rev') ! D4 re-parametrization
      num = p_wb97m_rev
   case('wb97', 'ωb97', 'omegab97', 'hyb_gga_xc_wb97')
      num = p_wb97
   case('wb97x-2008', 'ωb97x-2008', 'omegab97x-2008', 'hyb_gga_xc_wb97x', &
      & 'wb97x_2008', 'ωb97x_2008', 'omegab97x_2008')
      num = p_wb97x_2008
   case('wb97x', 'ωb97x', 'omegab97x', 'hyb_gga_xc_wb97x_v')
      num = p_wb97x
   case('wb97x-rev', 'ωb97x-rev', 'omegab97x-rev', 'wb97x_rev', 'ωb97x_rev', &
      & 'omegab97x_rev') ! D4 re-parametrization
      num = p_wb97x_rev
   case('wb97x-3c', 'ωb97x-3c', 'omegab97x-3c', 'wb97x_3c', 'ωb97x_3c', &
      & 'omegab97x_3c') ! no libxc
      num = p_wb97x_3c
   case('wr2scan', 'wr²scan') ! no libxc
      num = p_wr2scan
   case('r2scan0-dh', 'r²scan0-dh', 'r2scan0dh', 'r²scan0dh') ! no libxc
      num = p_r2scan0_dh
   case('r2scan-cidh', 'r²scan-cidh', 'r2scancidh', 'r²scancidh') ! no libxc
      num = p_r2scan_cidh
   case('r2scan-qidh', 'r²scan-qidh', 'r2scanqidh', 'r²scanqidh') ! no libxc
      num = p_r2scan_qidh
   case('r2scan0-2', 'r²scan0-2', 'r2scan02', 'r²scan02') ! no libxc
      num = p_r2scan0_2
   case('pr2scan50', 'pr²scan50') ! no libxc
      num = p_pr2scan50
   case('pr2scan69', 'pr²scan69') ! no libxc
      num = p_pr2scan69
   case('kpr2scan50', 'kpr²scan50') ! no libxc
      num = p_kpr2scan50
   case('wpr2scan50', 'wpr²scan50') ! no libxc
      num = p_wpr2scan50
   end select
end function get_functional_id


end module dftd4_param
