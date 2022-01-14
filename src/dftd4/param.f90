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
   implicit none
   private

   public :: get_rational_damping


   enum, bind(C)
      enumerator :: p_invalid, &
         & p_hf, p_blyp, p_bpbe, p_bp, p_bpw, p_lb94, p_mpwlyp, p_mpwpw, &
         & p_olyp, p_opbe, p_pbe, p_rpbe, p_revpbe, p_pw86pbe, &
         & p_rpw86pbe, p_pw91, p_pwp, p_xlyp, p_b97, p_tpss, p_revtpss, &
         & p_scan, p_rscan, p_r2scan, p_b1lyp, p_b3lyp, p_bhlyp, p_b1p, &
         & p_b3p, p_b1pw, p_b3pw, p_o3lyp, p_revpbe0, p_revpbe38, &
         & p_pbe0, p_pwp1, p_pw1pw, p_mpw1pw, p_mpw1lyp, p_pw6b95, &
         & p_tpssh, p_tpss0, p_x3lyp, p_m06l, p_m06, p_m062x, p_b97d, &
         & p_wb97, p_wb97x, p_b97m, p_wb97m, p_camb3lyp, p_lcblyp, &
         & p_lh07tsvwn, p_lh07ssvwn, p_lh12ctssirpw92, p_lh12ctssifpw92, &
         & p_lh14tcalpbe, p_lh20t, &
         & p_b2plyp, p_b2gpplyp, p_mpw2plyp, p_pwpb95, &
         & p_dsdblyp, p_dsdpbe, p_dsdpbeb95, p_dsdpbep86, p_dsdsvwn, &
         & p_dodblyp, p_dodpbe, p_dodpbeb95, p_dodpbep86, p_dodsvwn, &
         & p_pbe0_2, p_pbe0_dh, p_hf3c, p_hf3cv, p_pbeh3c, p_b973c, &
         & p_hsesol, p_pwgga, p_dftb_3ob, p_dftb_mio, p_dftb_ob2, &
         & p_dftb_matsci, p_dftb_pbc, p_hcth120, p_ptpss, p_lcwpbe, &
         & p_bmk, p_b1b95, p_pwb6k, p_otpss, p_ssb, p_revssb, &
         & p_pbesol, p_hse06, p_pbexalpha, p_pbehpbe, p_hcth407, &
         & p_n12, p_pkzb, p_thcth, p_m11l, p_mn15l, p_mpwb1k, &
         & p_mpw1kcis, p_mpwkcis1k, p_pbeh1pbe, p_pbe1kcis, p_b97_1, &
         & p_b97_2, p_b98, p_hiss, p_hse03, p_revtpssh, p_tpss1kcis, &
         & p_m05, p_m052x, p_m08hx, p_lcwhpbe, p_mn12l, p_tauhcthhyb, &
         & p_sogga11x, p_n12sx, p_mn12sx, p_mn15, p_glyp, p_bop, &
         & p_mpw1b95, p_revpbe0dh, p_revtpss0, p_revdsdpbep86, p_revdsdpbe, &
         & p_revdsdblyp, p_revdodpbep86, p_am05, p_hse12, p_hse12s, &
         & p_r2scanh, p_r2scan0, p_r2scan50
   end enum
   integer, parameter :: df_enum = kind(p_invalid)

contains

subroutine get_rational_damping(functional, param, s9)
   character(len=*), intent(in) :: functional
   class(damping_param), allocatable, intent(out) :: param
   real(wp), intent(in), optional :: s9

   character(len=:), allocatable :: fname
   integer :: is, id
   logical :: mbd

   mbd = merge(s9 /= 0.0_wp, .true., present(s9))

   is = index(functional, '/')
   if (is == 0) is = len_trim(functional) + 1
   fname = lowercase(functional(:is-1))

   id = get_functional_id(fname)

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

end subroutine get_rational_damping

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

   pure function dftd_param(s6, s8, a1, a2, alp) result(param)
      real(wp), intent(in) :: s8, a1, a2
      real(wp), intent(in), optional :: s6, alp
      type(rational_damping_param) :: param
      param = rational_damping_param(&
         & s6=merge(s6, 1.0_wp, present(s6)), &
         & s8=s8, a1=a1, a2=a2, &
         & s9=merge(s9, 0.0_wp, present(s9)), &
         & alp=merge(alp, 16.0_wp, present(alp)))
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
      param = dftd_param (s6=1.0_wp, s8=0.8324_wp, a1=0.4944_wp, a2=5.9019_wp)
   case(p_r2scan0)
      param = dftd_param (s6=1.0_wp, s8=0.8992_wp, a1=0.4778_wp, a2=5.8779_wp)
   case(p_r2scan50)
      param = dftd_param (s6=1.0_wp, s8=1.0471_wp, a1=0.4574_wp, a2=5.8969_wp)
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
   case(p_wb97x)
      param = dftd_param ( & ! (SAW190103)
         &  s6=1.0000_wp, s8=-0.07519516_wp, a1=0.45094893_wp, a2=6.78425255_wp )
      !  S22x5: MD= 0.05 MAD= 0.16 RMSD= 0.22
      !  S66x8: MD= 0.06 MAD= 0.16 RMSD= 0.21
      !  NCI10: MD= 0.08 MAD= 0.15 RMSD= 0.25
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
   end select

contains

   pure function dftd_param(s6, s8, a1, a2, alp) result(param)
      real(wp), intent(in) :: s8, a1, a2
      real(wp), intent(in), optional :: s6, alp
      type(rational_damping_param) :: param
      param = rational_damping_param(&
         & s6=merge(s6, 1.0_wp, present(s6)), &
         & s8=s8, a1=a1, a2=a2, &
         & s9=merge(s9, 1.0_wp, present(s9)), &
         & alp=merge(alp, 16.0_wp, present(alp)))
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
   case('b-lyp', 'blyp')
      num = p_blyp
   case('bpbe')
      num = p_bpbe
   case('b-p', 'bp86', 'bp', 'b-p86')
      num = p_bp
   case('bpw', 'b-pw')
      num = p_bpw
   case('lb94')
      num = p_lb94
   case('mpwlyp', 'mpw-lyp')
      num = p_mpwlyp
   case('mpwpw', 'mpw-pw', 'mpwpw91')
      num = p_mpwpw
   case('o-lyp', 'olyp')
      num = p_olyp
   case('opbe')
      num = p_opbe
   case('pbe')
      num = p_pbe
   case('rpbe')
      num = p_rpbe
   case('revpbe')
      num = p_revpbe
   case('pw86pbe')
      num = p_pw86pbe
   case('rpw86pbe')
      num = p_rpw86pbe
   case('pw91')
      num = p_pw91
   case('pwp', 'pw-p', 'pw91p86')
      num = p_pwp
   case('x-lyp', 'xlyp')
      num = p_xlyp
   case('b97')
      num = p_b97
   case('tpss')
      num = p_tpss
   case('revtpss')
      num = p_revtpss
   case('scan')
      num = p_scan
   case('rscan')
      num = p_rscan
   case('r2scan', 'r²scan')
      num = p_r2scan
   case('r2scanh', 'r²scanh')
      num = p_r2scanh
   case('r2scan0', 'r²scan0')
      num = p_r2scan0
   case('r2scan50', 'r²scan50')
      num = p_r2scan50
   case('b1lyp', 'b1-lyp')
      num = p_b1lyp
   case('b3-lyp', 'b3lyp')
      num = p_b3lyp
   case('bh-lyp', 'bhlyp')
      num = p_bhlyp
   case('b1p', 'b1-p', 'b1p86')
      num = p_b1p
   case('b3p', 'b3-p', 'b3p86')
      num = p_b3p
   case('b1pw', 'b1-pw', 'b1pw91')
      num = p_b1pw
   case('b3pw', 'b3-pw', 'b3pw91')
      num = p_b3pw
   case('o3-lyp', 'o3lyp')
      num = p_o3lyp
   case('revpbe0')
      num = p_revpbe0
   case('revpbe38')
      num = p_revpbe38
   case('pbe0')
      num = p_pbe0
   case('pwp1')
      num = p_pwp1
   case('pw1pw', 'pw1-pw')
      num = p_pw1pw
   case('mpw1pw', 'mpw1-pw', 'mpw1pw91')
      num = p_mpw1pw
   case('mpw1lyp', 'mpw1-lyp')
      num = p_mpw1lyp
   case('pw6b95')
      num = p_pw6b95
   case('tpssh')
      num = p_tpssh
   case('tpss0')
      num = p_tpss0
   case('x3-lyp', 'x3lyp')
      num = p_x3lyp
   case('m06l')
      num = p_m06l
   case('m06')
      num = p_m06
   case('m06-2x', 'm062x')
      num = p_m062x
   case('wb97', 'ωb97', 'omegab97')
      num = p_wb97
   case('wb97x', 'ωb97x', 'omegab97x')
      num = p_wb97x
   case('cam-b3lyp')
      num = p_camb3lyp
   case('lc-blyp')
      num = p_lcblyp
   case('lh07tsvwn', 'lh07t-svwn')
      num = p_lh07tsvwn
   case('lh07ssvwn', 'lh07s-svwn')
      num = p_lh07ssvwn
   case('lh12ctssirpw92', 'lh12ct-ssirpw92')
      num = p_lh12ctssirpw92
   case('lh12ctssifpw92', 'lh12ct-ssifpw92')
      num = p_lh12ctssifpw92
   case('lh14tcalpbe', 'lh14t-calpbe')
      num = p_lh14tcalpbe
   case('lh20t')
      num = p_lh20t
   case('b2plyp', 'b2-plyp')
      num = p_b2plyp
   case('b2gpplyp', 'b2gp-plyp')
      num = p_b2gpplyp
   case('mpw2plyp')
      num = p_mpw2plyp
   case('pwpb95')
      num = p_pwpb95
   case('dsdblyp', 'dsd-blyp')
      num = p_dsdblyp
   case('dsdpbe', 'dsd-pbe')
      num = p_dsdpbe
   case('dsdpbeb95', 'dsd-pbeb95')
      num = p_dsdpbeb95
   case('dsdpbep86', 'dsd-pbep86')
      num = p_dsdpbep86
   case('dsdsvwn', 'dsd-svwn')
      num = p_dsdsvwn
   case('dodblyp', 'dod-blyp')
      num = p_dodblyp
   case('dodpbe', 'dod-pbe')
      num = p_dodpbe
   case('dodpbeb95', 'dod-pbeb95')
      num = p_dodpbeb95
   case('dodpbep86', 'dod-pbep86')
      num = p_dodpbep86
   case('dodsvwn', 'dod-svwn')
      num = p_dodsvwn
   case('pbe02', 'pbe0-2')
      num = p_pbe0_2
   case('pbe0dh', 'pbe0-dh')
      num = p_pbe0_dh
   case('hf-3c', 'hf3c')
      num = p_hf3c
   case('hf-3cv', 'hf3cv')
      num = p_hf3cv
   case('pbeh3c', 'pbeh-3c')
      num = p_pbeh3c
   case('b973c', 'b97-3c')
      num = p_b973c
   case('pwgga')
      num = p_pwgga
   case('dftb3', 'dftb(3ob)')
      num = p_dftb_3ob
   case('dftb(mio)')
      num = p_dftb_mio
   case('dftb(pbc)')
      num = p_dftb_pbc
   case('dftb(matsci)')
      num = p_dftb_matsci
   case('lc-dftb', 'dftb(ob2)')
      num = p_dftb_ob2
   case('hcth120')
      num = p_hcth120
   case('ptpss')
      num = p_ptpss
   case('lc-wpbe', 'lcwpbe')
      num = p_lcwpbe
   case('bmk')
      num = p_bmk
   case('b1b95')
      num = p_b1b95
   case('bwb6k')
      num = p_pwb6k
   case('otpss')
      num = p_otpss
   case('ssb')
      num = p_ssb
   case('revssb')
      num = p_revssb
   case('pbesol')
      num = p_pbesol
   case('pbexalpha')
      num = p_pbexalpha
   case('pbehpbe')
      num = p_pbehpbe
   case('hcth407')
      num = p_hcth407
   case('n12')
      num = p_n12
   case('pkzb')
      num = p_pkzb
   case('thcth', 'tauhctc')
      num = p_thcth
   case('m11l')
      num = p_m11l
   case('mn15l')
      num = p_mn15l
   case('mpwb1k')
      num = p_mpwb1k
   case('mpw1kcis')
      num = p_mpw1kcis
   case('mpwkcis1k')
      num = p_mpwkcis1k
   case('pbeh1pbe')
      num = p_pbeh1pbe
   case('pbe1kcis')
      num = p_pbe1kcis
   case('b97-1')
      num = p_b97_1
   case('b97-2')
      num = p_b97_2
   case('b98')
      num = p_b98
   case('hiss')
      num = p_hiss
   case('hse03')
      num = p_hse03
   case('hse06')
      num = p_hse06
   case('hse12')
      num = p_hse12
   case('hse12s')
      num = p_hse12s
   case('hsesol')
      num = p_hsesol
   case('revtpssh')
      num = p_revtpssh
   case('tpss1kcis')
      num = p_tpss1kcis
   case('m05')
      num = p_m05
   case('m052x', 'm05-2x')
      num = p_m052x
   case('m08hx', 'm08-hx')
      num = p_m08hx
   case('lcwhpbe', 'lc-whpbe')
      num = p_lcwhpbe
   case('mn12l')
      num = p_mn12l
   case('tauhcthhyb')
      num = p_tauhcthhyb
   case('sogga11x')
      num = p_sogga11x
   case('n12sx')
      num = p_n12sx
   case('mn12sx')
      num = p_mn12sx
   case('mn15')
      num = p_mn15
   case('glyp', 'g-lyp')
      num = p_glyp
   case('revpbe0dh', 'revpbe0-dh')
      num = p_revpbe0dh
   case('revtpss0')
      num = p_revtpss0
   case('revdsd-pbep86', 'revdsdpbep86')
      num = p_revdsdpbep86
   case('revdsd-pbe', 'revdsd-pbepbe', 'revdsdpbe', 'revdsdpbepbe')
      num = p_revdsdpbe
   case('revdsd-blyp', 'revdsdblyp')
      num = p_revdsdblyp
   case('revdod-pbep86', 'revdodpbep86')
      num = p_revdodpbep86
   case('b97m')
      num = p_b97m
   case('wb97m', 'ωb97m', 'omegab97m')
      num = p_wb97m
   case('am05')
      num = p_am05
   end select
end function get_functional_id

!> Convert string to lower case
pure function lowercase(str) result(lcstr)
   character(len=*), intent(in)  :: str
   character(len=len_trim(str)) :: lcstr
   integer :: ilen, ioffset, iquote, i, iav, iqc

   ilen=len_trim(str)
   ioffset=iachar('A')-iachar('a')
   iquote=0
   lcstr=str
   do i=1, ilen
      iav=iachar(str(i:i))
      if(iquote==0 .and. (iav==34 .or.iav==39)) then
         iquote=1
         iqc=iav
        cycle
      endif
      if(iquote==1 .and. iav==iqc) then
         iquote=0
         cycle
      endif
      if (iquote==1) cycle
      if(iav >= iachar('A') .and. iav <= iachar('Z')) then
         lcstr(i:i)=achar(iav-ioffset)
      else
         lcstr(i:i)=str(i:i)
      endif
   enddo

end function lowercase


end module dftd4_param
