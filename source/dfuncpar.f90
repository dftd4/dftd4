! This file is part of dftd4.
!
! Copyright (C) 2017-2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> saves all the damping parameters obtained by fitting the DFT-D4
!  dispersion energies to CCSD(T) interaction energies obtained
!  on different dissociation curves
module dfuncpar
   use iso_fortran_env, only : wp => real64
   use class_param, only : dftd_parameter
   implicit none
   !> unique idenifiers for all functionals known to me -- SAW
   enum, bind(C)
      enumerator :: &
         &  p_df_none, &
         &  p_df_hf, &
         &  p_df_blyp, &
         &  p_df_bpbe, &
         &  p_df_bp, &
         &  p_df_bpw, &
         &  p_df_lb94, &
         &  p_df_mpwlyp, &
         &  p_df_mpwpw, &
         &  p_df_olyp, &
         &  p_df_opbe, &
         &  p_df_pbe, &
         &  p_df_rpbe, &
         &  p_df_revpbe, &
         &  p_df_pw86pbe, &
         &  p_df_rpw86pbe, &
         &  p_df_pw91, &
         &  p_df_pwp, &
         &  p_df_xlyp, &
         &  p_df_b97, &
         &  p_df_tpss, &
         &  p_df_revtpss, &
         &  p_df_scan, &
         &  p_df_b1lyp, &
         &  p_df_b3lyp, &
         &  p_df_bhlyp, &
         &  p_df_b1p, &
         &  p_df_b3p, &
         &  p_df_b1pw, &
         &  p_df_b3pw, &
         &  p_df_o3lyp, &
         &  p_df_revpbe0, &
         &  p_df_revpbe38, &
         &  p_df_pbe0, &
         &  p_df_pwp1, &
         &  p_df_pw1pw, &
         &  p_df_mpw1pw, &
         &  p_df_mpw1lyp, &
         &  p_df_pw6b95, &
         &  p_df_tpssh, &
         &  p_df_tpss0, &
         &  p_df_x3lyp, &
         &  p_df_m06l, &
         &  p_df_m06, &
         &  p_df_m062x, &
         &  p_df_wb97, &
         &  p_df_wb97x, &
         &  p_df_camb3lyp, &
         &  p_df_lcblyp, &
         &  p_df_lh07tsvwn, &
         &  p_df_lh07ssvwn, &
         &  p_df_lh12ctssirpw92, &
         &  p_df_lh12ctssifpw92, &
         &  p_df_lh14tcalpbe, &
         &  p_df_b2plyp, &
         &  p_df_b2gpplyp, &
         &  p_df_mpw2plyp, &
         &  p_df_pwpb95, &
         &  p_df_dsdblyp, &
         &  p_df_dsdpbe, &
         &  p_df_dsdpbeb95, &
         &  p_df_dsdpbep86, &
         &  p_df_dsdsvwn, &
         &  p_df_dodblyp, &
         &  p_df_dodpbe, &
         &  p_df_dodpbeb95, &
         &  p_df_dodpbep86, &
         &  p_df_dodsvwn, &
         &  p_df_pbe0_2, &
         &  p_df_pbe0_dh, &
         &  p_df_hf3c, &
         &  p_df_hf3cv, &
         &  p_df_pbeh3c, &
         &  p_df_b973c, &
         &  p_df_hsesol, &
         &  p_df_pwgga, &
         &  p_df_dftb_3ob, &
         &  p_df_dftb_mio, &
         &  p_df_dftb_ob2, &
         &  p_df_dftb_matsci, &
         &  p_df_dftb_pbc, &
         &  p_df_hcth120, &
         &  p_df_ptpss, &
         &  p_df_lcwpbe, &
         &  p_df_bmk, &
         &  p_df_b1b95, &
         &  p_df_pwb6k, &
         &  p_df_otpss, &
         &  p_df_ssb, &
         &  p_df_revssb, &
         &  p_df_pbesol, &
         &  p_df_hse06, &
         &  p_df_pbexalpha, &
         &  p_df_pbehpbe, &
         &  p_df_hcth407, &
         &  p_df_n12, &
         &  p_df_pkzb, &
         &  p_df_thcth, &
         &  p_df_m11l, &
         &  p_df_mn15l, &
         &  p_df_mpwb1k, &
         &  p_df_mpw1kcis, &
         &  p_df_mpwkcis1k, &
         &  p_df_pbeh1pbe, &
         &  p_df_pbe1kcis, &
         &  p_df_b97_1, &
         &  p_df_b97_2, &
         &  p_df_b98, &
         &  p_df_hiss, &
         &  p_df_hse03, &
         &  p_df_revtpssh, &
         &  p_df_tpss1kcis, &
         &  p_df_m05, &
         &  p_df_m052x, &
         &  p_df_m08hx, &
         &  p_df_lcwhpbe, &
         &  p_df_mn12l, &
         &  p_df_tauhcthhyb, &
         &  p_df_sogga11x, &
         &  p_df_n12sx, &
         &  p_df_mn12sx, &
         &  p_df_mn15, &
         &  p_df_glyp, &
         &  p_df_bop, &
         &  p_df_mpw1b95, &
         &  p_df_revpbe0dh, &
         &  p_df_revtpss0, &
         &  p_df_revdsdpbep86, &
         &  p_df_revdsdpbe, &
         &  p_df_revdsdblyp, &
         &  p_df_revdodpbep86
   end enum
   integer, parameter :: df_enum = kind(p_df_none)
   enum, bind(C)
      enumerator :: &
         &  p_bas_default, &
         &  p_bas_631gd, &
         &  p_bas_mixed, &
         &  p_bas_sv, &
         &  p_bas_minis
   end enum
   integer, parameter :: bas_enum = kind(p_bas_default)

contains

! ========================================================================
!> DFT-D4(EEQ)/def2-QZVP fitted on NCIBLIND10, S22x5, S66x8
subroutine get_d4eeqbj_2019_parameter(dfnum,bsnum,param)
   implicit none
   integer(df_enum),intent(in) :: dfnum
   integer(bas_enum),intent(in) :: bsnum
   type(dftd_parameter),allocatable :: param
   select case(dfnum)
   case default; return
   case(p_df_dftb_3ob); param = dftd_parameter( & ! (SAW191202)
   &  s6=1.0_wp, s8=0.4727337_wp, a1=0.5467502_wp, a2=4.4955068_wp, s9=0.0_wp )
   case(p_df_dftb_matsci); param = dftd_parameter( & ! (SAW191202)
   &  s6=1.0_wp, s8=2.7711819_wp, a1=0.4681712_wp, a2=5.2918629_wp, s9=0.0_wp )
   case(p_df_dftb_mio); param = dftd_parameter( & ! (SAW191202)
   &  s6=1.0_wp, s8=1.1948145_wp, a1=0.6074567_wp, a2=4.9336133_wp, s9=0.0_wp )
   case(p_df_dftb_ob2); param = dftd_parameter( & ! (SAW191202)
   &  s6=1.0_wp, s8=2.7611320_wp, a1=0.6037249_wp, a2=5.3900004_wp, s9=0.0_wp )
   case(p_df_dftb_pbc); param = dftd_parameter( & ! (SAW191202)
   &  s6=1.0_wp, s8=1.7303734_wp, a1=0.5546548_wp, a2=4.7973454_wp, s9=0.0_wp )
   end select
! ------------------------------------------------------------------------
end subroutine get_d4eeqbj_2019_parameter
! ========================================================================

! ========================================================================
!> DFT-D4(EEQ)-ATM/def2-QZVP fitted on NCIBLIND10, S22x5, S66x8
subroutine get_d4eeqbjatm_2019_parameter(dfnum,bsnum,param)
   implicit none
   integer(df_enum),intent(in) :: dfnum
   integer(bas_enum),intent(in) :: bsnum
   type(dftd_parameter),allocatable :: param
   select case(dfnum)
   case default; return
! ------------------------------------------------------------------------
   case(p_df_b1b95); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=1.27701162_wp, a1=0.40554715_wp, a2=4.63323074_wp )
!  Fitset: MD= 0.22852 MAD= 0.35189 RMSD= 0.46982
   case(p_df_b1lyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.98553711_wp, a1=0.39309040_wp, a2=4.55465145_wp )
!  Fitset: MD= -0.04797 MAD= 0.25597 RMSD= 0.38778
   case(p_df_b1p); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=3.36115015_wp, a1=0.48665293_wp, a2=5.05219572_wp )
!  Fitset: MD= -0.01406 MAD= 0.27441 RMSD= 0.47328
   case(p_df_b1pw); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=3.02227550_wp, a1=0.47396846_wp, a2=4.49845309_wp )
!  Fitset: MD= 0.10485 MAD= 0.32175 RMSD= 0.48508
   case(p_df_b2gpplyp); param = dftd_parameter ( & ! (SAW190107)
   &  s6=0.5600_wp, s8=0.94633372_wp, a1=0.42907301_wp, a2=5.18802602_wp )
!  Fitset: MD= -0.05248 MAD= 0.18110 RMSD= 0.27365
   case(p_df_b2plyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.6400_wp, s8=1.16888646_wp, a1=0.44154604_wp, a2=4.73114642_wp )
!  Fitset: MD= -0.03761 MAD= 0.18247 RMSD= 0.27109
   case(p_df_b3lyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.02929367_wp, a1=0.40868035_wp, a2=4.53807137_wp )
!  Fitset: MD= -0.05892 MAD= 0.26117 RMSD= 0.40531
   case(p_df_b3p); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=3.08822155_wp, a1=0.47324238_wp, a2=4.98682134_wp )
!  Fitset: MD= -0.02970 MAD= 0.26962 RMSD= 0.46761
   case(p_df_b3pw); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=2.88364295_wp, a1=0.46990860_wp, a2=4.51641422_wp )
!  Fitset: MD= 0.06643 MAD= 0.29151 RMSD= 0.45541
   case(p_df_b97); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.87854260_wp, a1=0.29319126_wp, a2=4.51647719_wp )
!  Fitset: MD= -0.13017 MAD= 0.24778 RMSD= 0.36116
   case(p_df_bhlyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.65281646_wp, a1=0.27263660_wp, a2=5.48634586_wp )
!  Fitset: MD= -0.15832 MAD= 0.34132 RMSD= 0.57342
   case(p_df_blyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.34076671_wp, a1=0.44488865_wp, a2=4.09330090_wp )
!  Fitset: MD= 0.04801 MAD= 0.28161 RMSD= 0.38321
   case(p_df_bpbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=3.64405246_wp, a1=0.52905620_wp, a2=4.11311891_wp )
!  Fitset: MD= 0.19316 MAD= 0.41912 RMSD= 0.60452
   case(p_df_bp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=3.35497927_wp, a1=0.43645861_wp, a2=4.92406854_wp )
!  Fitset: MD= 0.08252 MAD= 0.32681 RMSD= 0.47063
   case(p_df_bpw); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=3.24571506_wp, a1=0.50050454_wp, a2=4.12346483_wp )
!  Fitset: MD= 0.20607 MAD= 0.41941 RMSD= 0.59589
   case(p_df_camb3lyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.66041301_wp, a1=0.40267156_wp, a2=5.17432195_wp )
!  Fitset: MD= -0.19675 MAD= 0.34901 RMSD= 0.59087
   case(p_df_dodblyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.4700_wp, s8=1.31146043_wp, a1=0.43407294_wp, a2=4.27914360_wp )
!  Fitset: MD= 0.03323 MAD= 0.13858 RMSD= 0.20861
   case(p_df_dodpbeb95); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.5600_wp, s8=0.01574635_wp, a1=0.43745720_wp, a2=3.69180763_wp )
!  Fitset: MD= 0.03704 MAD= 0.13343 RMSD= 0.18278
   case(p_df_dodpbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.4800_wp, s8=0.92051454_wp, a1=0.43037052_wp, a2=4.38067238_wp )
!  Fitset: MD= 0.01065 MAD= 0.13414 RMSD= 0.21424
   case(p_df_dodpbep86); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.4600_wp, s8=0.71405681_wp, a1=0.42408665_wp, a2=4.52884439_wp )
!  Fitset: MD= -0.03740 MAD= 0.12467 RMSD= 0.18127
   case(p_df_dodsvwn); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.4200_wp, s8=0.94500207_wp, a1=0.47449026_wp, a2=5.05316093_wp )
!  Fitset: MD= -0.07427 MAD= 0.16970 RMSD= 0.25286
   case(p_df_dsdblyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.5400_wp, s8=0.63018237_wp, a1=0.47591835_wp, a2=4.73713781_wp )
!  Fitset: MD= -0.01981 MAD= 0.14823 RMSD= 0.21530
   case(p_df_dsdpbeb95); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.5400_wp, s8=-0.14668670_wp, a1=0.46394587_wp, a2=3.64913860_wp )
!  Fitset: MD= 0.02996 MAD= 0.12414 RMSD= 0.16860
   case(p_df_dsdpbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.4500_wp, s8=0.70584116_wp, a1=0.45787085_wp, a2=4.44566742_wp )
!  Fitset: MD= 0.00866 MAD= 0.13406 RMSD= 0.21380
   case(p_df_dsdpbep86); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.4700_wp, s8=0.37586675_wp, a1=0.53698768_wp, a2=5.13022435_wp )
!  Fitset: MD= -0.05273 MAD= 0.14259 RMSD= 0.21271
   case(p_df_dsdsvwn); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.4100_wp, s8=0.72914436_wp, a1=0.51347412_wp, a2=5.11858541_wp )
!  Fitset: MD= -0.08974 MAD= 0.32285 RMSD= 0.43146
   case(p_df_glyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=4.23798924_wp, a1=0.38426465_wp, a2=4.38412863_wp )
!  Fitset: MD= 0.63466 MAD= 0.89568 RMSD= 1.11309
   case(p_df_hf); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.61679827_wp, a1=0.44959224_wp, a2=3.35743605_wp )
!  Fitset: MD= -0.02597 MAD= 0.34732 RMSD= 0.49719
   case(p_df_lb94); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.59538499_wp, a1=0.42088944_wp, a2=3.28193223_wp )
!  Fitset: MD= 0.31701 MAD= 0.53196 RMSD= 0.74553
   case(p_df_lcblyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.60344180_wp, a1=0.45769839_wp, a2=7.86924893_wp )
!  Fitset: MD= -0.39724 MAD= 0.72327 RMSD= 1.18218
   case(p_df_lh07ssvwn); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=3.16675531_wp, a1=0.35965552_wp, a2=4.31947614_wp )
!  Fitset: MD= 0.32224 MAD= 0.59006 RMSD= 0.86272
   case(p_df_lh07tsvwn); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.09333001_wp, a1=0.35025189_wp, a2=4.34166515_wp )
!  Fitset: MD= 0.24243 MAD= 0.43497 RMSD= 0.61671
   case(p_df_lh12ctssifpw92); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.68467610_wp, a1=0.34190416_wp, a2=3.91039666_wp )
!  Fitset: MD= 0.55106 MAD= 0.80783 RMSD= 1.11048
   case(p_df_lh12ctssirpw92); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.48973402_wp, a1=0.34026075_wp, a2=3.96948081_wp )
!  Fitset: MD= 0.47785 MAD= 0.71188 RMSD= 0.98422
   case(p_df_lh14tcalpbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.28130770_wp, a1=0.38822021_wp, a2=4.92501211_wp )
!  Fitset: MD= -0.02105 MAD= 0.22968 RMSD= 0.36045
   case(p_df_m06); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.16366729_wp, a1=0.53456413_wp, a2=6.06192174_wp )
!  Fitset: MD= 0.01788 MAD= 0.24914 RMSD= 0.38604
   case(p_df_m06l); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.59493760_wp, a1=0.71422359_wp, a2=6.35314182_wp )
!  Fitset: MD= 0.08395 MAD= 0.24888 RMSD= 0.34879
   case(p_df_mpw1b95); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=0.50093024_wp, a1=0.41585097_wp, a2=4.99154869_wp )
!  Fitset: MD= 0.00585 MAD= 0.15695 RMSD= 0.21297
   case(p_df_mpw1lyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.15591153_wp, a1=0.25603493_wp, a2=5.32083895_wp )
!  Fitset: MD= -0.26979 MAD= 0.41542 RMSD= 0.60678
   case(p_df_mpw1pw); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.80841716_wp, a1=0.42961819_wp, a2=4.68892341_wp )
!  Fitset: MD= -0.08840 MAD= 0.26815 RMSD= 0.45231
   case(p_df_mpw2plyp); param = dftd_parameter ( & ! (SAW190107)
   &  s6=0.7500_wp, s8=0.45788846_wp, a1=0.42997704_wp, a2=5.07650682_wp )
!  Fitset: MD= -0.18921 MAD= 0.30115 RMSD= 0.44049
   case(p_df_mpwb1k); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=0.57338313_wp, a1=0.44687975_wp, a2=5.21266777_wp )
!  Fitset: MD= -0.00870 MAD= 0.17226 RMSD= 0.23614
   case(p_df_mpwlyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.25842942_wp, a1=0.25773894_wp, a2=5.02319542_wp )
!  Fitset: MD= -0.24426 MAD= 0.39145 RMSD= 0.54503
   case(p_df_mpwpw); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.82596836_wp, a1=0.34526745_wp, a2=4.84620734_wp )
!  Fitset: MD= -0.06278 MAD= 0.27913 RMSD= 0.43988
   case(p_df_o3lyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.75762508_wp, a1=0.10348980_wp, a2=6.16233282_wp )
!  Fitset: MD= -0.19268 MAD= 0.38577 RMSD= 0.62168
   case(p_df_olyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.74836820_wp, a1=0.60184498_wp, a2=2.53292167_wp )
!  Fitset: MD= 0.12352 MAD= 0.37113 RMSD= 0.58291
   case(p_df_opbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=3.06917417_wp, a1=0.68267534_wp, a2=2.22849018_wp )
!  Fitset: MD= 0.26699 MAD= 0.55308 RMSD= 0.85023
   case(p_df_pbe0_2); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.5000_wp, s8=0.64299082_wp, a1=0.76542115_wp, a2=5.78578675_wp )
!  Fitset: MD= -0.04260 MAD= 0.21186 RMSD= 0.34045
   case(p_df_pbe0); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.20065498_wp, a1=0.40085597_wp, a2=5.02928789_wp )
!  Fitset: MD= -0.17892 MAD= 0.30557 RMSD= 0.51050
   case(p_df_pbe0_dh); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.8750_wp, s8=0.96811578_wp, a1=0.47592488_wp, a2=5.08622873_wp )
!  Fitset: MD= -0.13857 MAD= 0.27919 RMSD= 0.47256
   case(p_df_pbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.95948085_wp, a1=0.38574991_wp, a2=4.80688534_wp )
!  Fitset: MD= -0.20544 MAD= 0.33635 RMSD= 0.51168
   case(p_df_pw1pw); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.96850170_wp, a1=0.42427511_wp, a2=5.02060636_wp )
!  Fitset: MD= -0.27325 MAD= 0.42206 RMSD= 0.64119
   case(p_df_pw6b95); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=-0.31926054_wp, a1=0.04142919_wp, a2=5.84655608_wp )
!  Fitset: MD= -0.04767 MAD= 0.14330 RMSD= 0.18958
   case(p_df_pw86pbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.21362856_wp, a1=0.40510366_wp, a2=4.66737724_wp )
!  Fitset: MD= -0.11505 MAD= 0.24691 RMSD= 0.38101
   case(p_df_pw91); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.77283111_wp, a1=0.39581542_wp, a2=4.93405761_wp )
!  Fitset: MD= -0.33019 MAD= 0.48611 RMSD= 0.68110
   case(p_df_pwp1); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.60492565_wp, a1=0.46855837_wp, a2=5.76921413_wp )
!  Fitset: MD= -0.35321 MAD= 0.54026 RMSD= 0.86629
   case(p_df_pwpb95); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.8200_wp, s8=-0.34639127_wp, a1=0.41080636_wp, a2=3.83878274_wp )
!  Fitset: MD= 0.02143 MAD= 0.13040 RMSD= 0.17599
   case(p_df_pwp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.32801227_wp, a1=0.35874687_wp, a2=6.05861168_wp )
!  Fitset: MD= -0.42482 MAD= 0.62607 RMSD= 0.91840
   case(p_df_revpbe0); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.57185414_wp, a1=0.38705966_wp, a2=4.11028876_wp )
!  Fitset: MD= 0.02724 MAD= 0.21587 RMSD= 0.36040
   case(p_df_revpbe0dh); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.8750_wp, s8=1.24456037_wp, a1=0.36730560_wp, a2=4.71126482_wp )
!  Fitset: MD= -0.01089 MAD= 0.20910 RMSD= 0.33564
   case(p_df_revpbe38); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.66597472_wp, a1=0.39476833_wp, a2=4.39026628_wp )
!  Fitset: MD= -0.01326 MAD= 0.22598 RMSD= 0.36210
   case(p_df_revpbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.74676530_wp, a1=0.53634900_wp, a2=3.07261485_wp )
!  Fitset: MD= 0.05649 MAD= 0.25212 RMSD= 0.40863
   case(p_df_revtpss0); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=1.54664499_wp, a1=0.45890964_wp, a2=4.78426405_wp )
!  Fitset: MD= -0.05298 MAD= 0.19965 RMSD= 0.32081
   case(p_df_revtpss); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.53089454_wp, a1=0.44880597_wp, a2=4.64042317_wp )
!  Fitset: MD= -0.01904 MAD= 0.19568 RMSD= 0.29618
   case(p_df_revtpssh); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=1.52740307_wp, a1=0.45161957_wp, a2=4.70779483_wp )
!  Fitset: MD= -0.03731 MAD= 0.19133 RMSD= 0.29091
   case(p_df_rpbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.31183787_wp, a1=0.46169493_wp, a2=3.15711757_wp )
!  Fitset: MD= -0.07156 MAD= 0.26348 RMSD= 0.38671
   case(p_df_rpw86pbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.12624034_wp, a1=0.38151218_wp, a2=4.75480472_wp )
!  Fitset: MD= -0.12740 MAD= 0.26294 RMSD= 0.40614
   case(p_df_scan); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.46126056_wp, a1=0.62930855_wp, a2=6.31284039_wp )
!  Fitset: MD= -0.13170 MAD= 0.28640 RMSD= 0.51183
   case(p_df_tpss0); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.62438102_wp, a1=0.40329022_wp, a2=4.80537871_wp )
!  Fitset: MD= -0.09569 MAD= 0.26733 RMSD= 0.44767
   case(p_df_tpss); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.76596355_wp, a1=0.42822303_wp, a2=4.54257102_wp )
!  Fitset: MD= -0.09296 MAD= 0.27505 RMSD= 0.42537
   case(p_df_tpssh); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.85897750_wp, a1=0.44286966_wp, a2=4.60230534_wp )
!  Fitset: MD= -0.09858 MAD= 0.26757 RMSD= 0.42380
   case(p_df_wb97); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=6.55792598_wp, a1=0.76666802_wp, a2=8.36027334_wp )
!  Fitset: MD= -0.12779 MAD= 0.36152 RMSD= 0.49991
   case(p_df_wb97x); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=-0.07519516_wp, a1=0.45094893_wp, a2=6.78425255_wp )
!  Fitset: MD= -0.20216 MAD= 0.34696 RMSD= 0.53641
   case(p_df_x3lyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.54701429_wp, a1=0.20318443_wp, a2=5.61852648_wp )
!  Fitset: MD= -0.15607 MAD= 0.31342 RMSD= 0.49546
   case(p_df_xlyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.62972054_wp, a1=0.11268673_wp, a2=5.40786417_wp )
!  Fitset: MD= -0.03900 MAD= 0.27562 RMSD= 0.38491
   case(p_df_revdsdpbep86); param = dftd_parameter ( & ! (WTMAD2)
   &  s6=0.5132_wp, s8=0.00000000_wp, a1=0.44000000_wp, a2=3.60000000_wp )
   case(p_df_revdsdpbe); param = dftd_parameter ( & ! (WTMAD2)
   &  s6=0.6706_wp, s8=0.00000000_wp, a1=0.40000000_wp, a2=3.60000000_wp )
   case(p_df_revdsdblyp); param = dftd_parameter ( & !(WTMAD2)
   &  s6=0.6141_wp, s8=0.00000000_wp, a1=0.38000000_wp, a2=3.52000000_wp )
   case(p_df_revdodpbep86); param = dftd_parameter ( & !(WTMAD2)
   &  s6=0.5552_wp, s8=0.00000000_wp, a1=0.44000000_wp, a2=3.60000000_wp )
   case(p_df_dftb_3ob); param = dftd_parameter( & ! (SAW191202)
   &  s6=1.0_wp, s8=0.6635015_wp, a1=0.5523240_wp, a2=4.3537076_wp, s9=1.0_wp )
   case(p_df_dftb_matsci); param = dftd_parameter( & ! (SAW191202)
   &  s6=1.0_wp, s8=3.3157614_wp, a1=0.4826330_wp, a2=5.3811976_wp, s9=1.0_wp )
   case(p_df_dftb_mio); param = dftd_parameter( & ! (SAW191202)
   &  s6=1.0_wp, s8=1.2916225_wp, a1=0.5965326_wp, a2=4.8778602_wp, s9=1.0_wp )
   case(p_df_dftb_ob2); param = dftd_parameter( & ! (SAW191202)
   &  s6=1.0_wp, s8=2.9692689_wp, a1=0.6068916_wp, a2=5.4476789_wp, s9=1.0_wp )
   case(p_df_dftb_pbc); param = dftd_parameter( & ! (SAW191202)
   &  s6=1.0_wp, s8=2.1667394_wp, a1=0.5646391_wp, a2=4.9576353_wp, s9=1.0_wp )
   end select
end subroutine get_d4eeqbjatm_2019_parameter
! ========================================================================

! ========================================================================
!> DFT-D4(EEQ)-MBD/def2-QZVP fitted on NCIBLIND10, S22x5, S66x8
subroutine get_d4eeqbjmbd_2019_parameter(dfnum,bsnum,param)
   implicit none
   integer(df_enum),intent(in) :: dfnum
   integer(bas_enum),intent(in) :: bsnum
   type(dftd_parameter), allocatable :: param
   select case(dfnum)
   case default; return
! ------------------------------------------------------------------------
   case(p_df_b1b95); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=1.19549420_wp, a1=0.39241474_wp, a2=4.60397611_wp )
!  Fitset: MD= 0.21329 MAD= 0.33289 RMSD= 0.44693
   case(p_df_b1lyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.94609514_wp, a1=0.38643351_wp, a2=4.54135968_wp )
!  Fitset: MD= -0.06493 MAD= 0.25607 RMSD= 0.39776
   case(p_df_b1p); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=3.38693011_wp, a1=0.48478615_wp, a2=5.04361224_wp )
!  Fitset: MD= -0.02348 MAD= 0.27543 RMSD= 0.48014
   case(p_df_b1pw); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=2.98402204_wp, a1=0.46862950_wp, a2=4.48637849_wp )
!  Fitset: MD= 0.09181 MAD= 0.31824 RMSD= 0.48100
   case(p_df_b2gpplyp); param = dftd_parameter ( & ! (SAW190107)
   &  s6=0.5600_wp, s8=1.00494214_wp, a1=0.42447353_wp, a2=5.19461329_wp )
!  Fitset: MD= -0.06339 MAD= 0.18624 RMSD= 0.28316
   case(p_df_b2plyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.6400_wp, s8=1.15117773_wp, a1=0.42666167_wp, a2=4.73635790_wp )
!  Fitset: MD= -0.05031 MAD= 0.18506 RMSD= 0.28010
   case(p_df_b3lyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.00246246_wp, a1=0.40276191_wp, a2=4.52778320_wp )
!  Fitset: MD= -0.07554 MAD= 0.26205 RMSD= 0.41586
   case(p_df_b3p); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=3.14456298_wp, a1=0.47187947_wp, a2=4.98624258_wp )
!  Fitset: MD= -0.04085 MAD= 0.27080 RMSD= 0.47542
   case(p_df_b3pw); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=2.85656268_wp, a1=0.46491801_wp, a2=4.50601452_wp )
!  Fitset: MD= 0.05302 MAD= 0.28885 RMSD= 0.45409
   case(p_df_b97); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.81171211_wp, a1=0.28461283_wp, a2=4.48691468_wp )
!  Fitset: MD= -0.15506 MAD= 0.26269 RMSD= 0.38008
   case(p_df_bhlyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.68082973_wp, a1=0.26835837_wp, a2=5.48847218_wp )
!  Fitset: MD= -0.17297 MAD= 0.34819 RMSD= 0.58920
   case(p_df_blyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.33971306_wp, a1=0.44733688_wp, a2=4.06583931_wp )
!  Fitset: MD= 0.02464 MAD= 0.27084 RMSD= 0.37647
   case(p_df_bpbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=3.65322996_wp, a1=0.49933501_wp, a2=4.24294852_wp )
!  Fitset: MD= 0.18169 MAD= 0.41723 RMSD= 0.60190
   case(p_df_bp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=3.33728176_wp, a1=0.43220330_wp, a2=4.91443061_wp )
!  Fitset: MD= 0.07108 MAD= 0.32331 RMSD= 0.46959
   case(p_df_bpw); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=3.23137432_wp, a1=0.49955226_wp, a2=4.10411084_wp )
!  Fitset: MD= 0.18835 MAD= 0.41014 RMSD= 0.58149
   case(p_df_camb3lyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.74407961_wp, a1=0.40137870_wp, a2=5.18731225_wp )
!  Fitset: MD= -0.20897 MAD= 0.35917 RMSD= 0.60661
   case(p_df_dodblyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.4700_wp, s8=1.17809956_wp, a1=0.40252428_wp, a2=4.25096555_wp )
!  Fitset: MD= 0.01949 MAD= 0.13138 RMSD= 0.19558
   case(p_df_dodpbeb95); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.5400_wp, s8=-0.15702803_wp, a1=0.30629389_wp, a2=3.69170956_wp )
!  Fitset: MD= -0.00076 MAD= 0.09552 RMSD= 0.12939
   case(p_df_dodpbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.4800_wp, s8=0.83908332_wp, a1=0.40655901_wp, a2=4.33601239_wp )
!  Fitset: MD= -0.00580 MAD= 0.12659 RMSD= 0.20494
   case(p_df_dodpbep86); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.4600_wp, s8=0.68309910_wp, a1=0.40600975_wp, a2=4.50011772_wp )
!  Fitset: MD= -0.05518 MAD= 0.12509 RMSD= 0.18597
   case(p_df_dodsvwn); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.4200_wp, s8=1.01890345_wp, a1=0.46167459_wp, a2=5.11121382_wp )
!  Fitset: MD= -0.08458 MAD= 0.17447 RMSD= 0.26263
   case(p_df_dsdblyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.5400_wp, s8=0.65438817_wp, a1=0.46549574_wp, a2=4.73449899_wp )
!  Fitset: MD= -0.03244 MAD= 0.15073 RMSD= 0.22004
   case(p_df_dsdpbeb95); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.5400_wp, s8=-0.24336862_wp, a1=0.32697409_wp, a2=3.69767540_wp )
!  Fitset: MD= -0.00830 MAD= 0.09184 RMSD= 0.12546
   case(p_df_dsdpbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.4500_wp, s8=0.66116783_wp, a1=0.43565915_wp, a2=4.41110670_wp )
!  Fitset: MD= -0.00698 MAD= 0.12881 RMSD= 0.20652
   case(p_df_dsdpbep86); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.4700_wp, s8=0.51157821_wp, a1=0.53889789_wp, a2=5.18645943_wp )
!  Fitset: MD= -0.06077 MAD= 0.14727 RMSD= 0.21987
   case(p_df_dsdsvwn); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.4100_wp, s8=0.90084457_wp, a1=0.51106529_wp, a2=5.22490148_wp )
!  Fitset: MD= -0.08917 MAD= 0.32295 RMSD= 0.43170
   case(p_df_glyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=3.83861584_wp, a1=0.36343954_wp, a2=4.32875183_wp )
!  Fitset: MD= 0.63864 MAD= 0.89437 RMSD= 1.11594
   case(p_df_hf); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.46001146_wp, a1=0.43186901_wp, a2=3.34116014_wp )
!  Fitset: MD= -0.05320 MAD= 0.34504 RMSD= 0.50301
   case(p_df_lb94); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.36461524_wp, a1=0.41518379_wp, a2=3.19365471_wp )
!  Fitset: MD= 0.28794 MAD= 0.49917 RMSD= 0.70755
   case(p_df_lcblyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.40109962_wp, a1=0.47867438_wp, a2=8.01038424_wp )
!  Fitset: MD= -0.39664 MAD= 0.72492 RMSD= 1.18496
   case(p_df_lh07ssvwn); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.92498406_wp, a1=0.34173988_wp, a2=4.28404951_wp )
!  Fitset: MD= 0.31275 MAD= 0.58894 RMSD= 0.85534
   case(p_df_lh07tsvwn); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.95389300_wp, a1=0.33511515_wp, a2=4.31853958_wp )
!  Fitset: MD= 0.22401 MAD= 0.42077 RMSD= 0.59849
   case(p_df_lh12ctssifpw92); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.41356607_wp, a1=0.31391316_wp, a2=3.88935769_wp )
!  Fitset: MD= 0.53322 MAD= 0.78975 RMSD= 1.08954
   case(p_df_lh12ctssirpw92); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.24917162_wp, a1=0.31446575_wp, a2=3.95070925_wp )
!  Fitset: MD= 0.45858 MAD= 0.69302 RMSD= 0.96196
   case(p_df_lh14tcalpbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.27677253_wp, a1=0.38128670_wp, a2=4.91698883_wp )
!  Fitset: MD= -0.03475 MAD= 0.22645 RMSD= 0.36554
   case(p_df_m06); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.22948274_wp, a1=0.52927285_wp, a2=6.06516782_wp )
!  Fitset: MD= 0.01376 MAD= 0.24790 RMSD= 0.38566
   case(p_df_m06l); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.40077779_wp, a1=0.69611405_wp, a2=6.29092087_wp )
!  Fitset: MD= 0.08204 MAD= 0.24719 RMSD= 0.34728
   case(p_df_mpw1b95); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=0.53791835_wp, a1=0.41016913_wp, a2=4.99284176_wp )
!  Fitset: MD= -0.00740 MAD= 0.15464 RMSD= 0.20903
   case(p_df_mpw1lyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.19986100_wp, a1=0.25502469_wp, a2=5.32301304_wp )
!  Fitset: MD= -0.28762 MAD= 0.43446 RMSD= 0.62777
   case(p_df_mpw1pw); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.80656973_wp, a1=0.42456967_wp, a2=4.68132317_wp )
!  Fitset: MD= -0.10273 MAD= 0.27455 RMSD= 0.46383
   case(p_df_mpw2plyp); param = dftd_parameter ( & ! (SAW190107)
   &  s6=0.7500_wp, s8=0.61161179_wp, a1=0.43748316_wp, a2=5.12540364_wp )
!  Fitset: MD= -0.20058 MAD= 0.31617 RMSD= 0.45768
   case(p_df_mpwb1k); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=0.62221146_wp, a1=0.44216745_wp, a2=5.21324659_wp )
!  Fitset: MD= -0.01872 MAD= 0.17358 RMSD= 0.23758
   case(p_df_mpwlyp); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=1.18243337_wp, a1=0.38968985_wp, a2=4.30835285_wp )
!  Fitset: MD= -0.27161 MAD= 0.42078 RMSD= 0.58498
   case(p_df_mpwpw); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.79674014_wp, a1=0.33870479_wp, a2=4.83442213_wp )
!  Fitset: MD= -0.07933 MAD= 0.28408 RMSD= 0.44655
   case(p_df_o3lyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.77793802_wp, a1=0.09961745_wp, a2=6.16089304_wp )
!  Fitset: MD= -0.20876 MAD= 0.39753 RMSD= 0.63678
   case(p_df_olyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.58717041_wp, a1=0.59759271_wp, a2=2.48760353_wp )
!  Fitset: MD= 0.09561 MAD= 0.37723 RMSD= 0.57321
   case(p_df_opbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=2.93544102_wp, a1=0.67903933_wp, a2=2.19810071_wp )
!  Fitset: MD= 0.23963 MAD= 0.55856 RMSD= 0.83684
   case(p_df_pbe0_2); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.5000_wp, s8=0.98834859_wp, a1=0.77911062_wp, a2=5.90389569_wp )
!  Fitset: MD= -0.04530 MAD= 0.21377 RMSD= 0.34234
   case(p_df_pbe0); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.26829475_wp, a1=0.39907098_wp, a2=5.03951304_wp )
!  Fitset: MD= -0.19236 MAD= 0.31907 RMSD= 0.52730
   case(p_df_pbe0_dh); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.8750_wp, s8=1.19306002_wp, a1=0.46106784_wp, a2=5.25210480_wp )
!  Fitset: MD= -0.15041 MAD= 0.29062 RMSD= 0.48405
   case(p_df_pbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.99924614_wp, a1=0.38142528_wp, a2=4.81839284_wp )
!  Fitset: MD= -0.22162 MAD= 0.35068 RMSD= 0.52976
   case(p_df_pw1pw); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.09759050_wp, a1=0.42759830_wp, a2=5.04559572_wp )
!  Fitset: MD= -0.28594 MAD= 0.43822 RMSD= 0.65926
   case(p_df_pw6b95); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=-0.31629935_wp, a1=0.03999357_wp, a2=5.83690254_wp )
!  Fitset: MD= -0.07192 MAD= 0.15341 RMSD= 0.20230
   case(p_df_pw86pbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.22842987_wp, a1=0.39998824_wp, a2=4.66739111_wp )
!  Fitset: MD= -0.13164 MAD= 0.25352 RMSD= 0.39362
   case(p_df_pw91); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.81406882_wp, a1=0.34094706_wp, a2=5.18568823_wp )
!  Fitset: MD= -0.34117 MAD= 0.49844 RMSD= 0.69615
   case(p_df_pwp1); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.95936222_wp, a1=0.48552982_wp, a2=5.84956411_wp )
!  Fitset: MD= -0.35812 MAD= 0.55157 RMSD= 0.87759
   case(p_df_pwpb95); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.8200_wp, s8=-0.46453780_wp, a1=0.29884136_wp, a2=3.87641255_wp )
!  Fitset: MD= -0.01702 MAD= 0.10667 RMSD= 0.14619
   case(p_df_pwp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.66056055_wp, a1=0.37768052_wp, a2=6.14787138_wp )
!  Fitset: MD= -0.42973 MAD= 0.63830 RMSD= 0.93118
   case(p_df_revpbe0); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.47198256_wp, a1=0.37471756_wp, a2=4.08904369_wp )
!  Fitset: MD= 0.00549 MAD= 0.21492 RMSD= 0.35464
   case(p_df_revpbe0dh); param = dftd_parameter ( & ! (SAW190103)
   &  s6=0.8750_wp, s8=1.22494188_wp, a1=0.35904781_wp, a2=4.70216012_wp )
!  Fitset: MD= -0.02829 MAD= 0.21126 RMSD= 0.33843
   case(p_df_revpbe38); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.60423529_wp, a1=0.38938475_wp, a2=4.35557832_wp )
!  Fitset: MD= -0.03169 MAD= 0.22643 RMSD= 0.36450
   case(p_df_revpbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.62543693_wp, a1=0.54031831_wp, a2=2.97965648_wp )
!  Fitset: MD= 0.02784 MAD= 0.24691 RMSD= 0.39242
   case(p_df_revtpss0); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=1.55321888_wp, a1=0.45355319_wp, a2=4.77588598_wp )
!  Fitset: MD= -0.06532 MAD= 0.20465 RMSD= 0.32888
   case(p_df_revtpss); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.51858035_wp, a1=0.44243222_wp, a2=4.62881620_wp )
!  Fitset: MD= -0.03308 MAD= 0.19729 RMSD= 0.29587
   case(p_df_revtpssh); param = dftd_parameter ( & ! (SAW190107)
   &  s6=1.0000_wp, s8=1.52542064_wp, a1=0.44570207_wp, a2=4.69883717_wp )
!  Fitset: MD= -0.05069 MAD= 0.19530 RMSD= 0.29596
   case(p_df_rpbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.11793696_wp, a1=0.44632488_wp, a2=3.08890917_wp )
!  Fitset: MD= -0.10098 MAD= 0.27012 RMSD= 0.39113
   case(p_df_rpw86pbe); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.13795871_wp, a1=0.37636536_wp, a2=4.75236384_wp )
!  Fitset: MD= -0.14431 MAD= 0.26983 RMSD= 0.41992
   case(p_df_scan); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.75408315_wp, a1=0.63571334_wp, a2=6.35690748_wp )
!  Fitset: MD= -0.13418 MAD= 0.28800 RMSD= 0.51482
   case(p_df_tpss0); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.66752698_wp, a1=0.40074746_wp, a2=4.80927196_wp )
!  Fitset: MD= -0.11178 MAD= 0.27563 RMSD= 0.45972
   case(p_df_tpss); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.91130849_wp, a1=0.43332851_wp, a2=4.56986797_wp )
!  Fitset: MD= -0.11454 MAD= 0.28494 RMSD= 0.43743
   case(p_df_tpssh); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.88783525_wp, a1=0.43968167_wp, a2=4.60342700_wp )
!  Fitset: MD= -0.11391 MAD= 0.27497 RMSD= 0.43539
   case(p_df_wb97); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=7.11022468_wp, a1=0.76423345_wp, a2=8.44559334_wp )
!  Fitset: MD= -0.12818 MAD= 0.36163 RMSD= 0.50023
   case(p_df_wb97x); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=0.38815338_wp, a1=0.47448629_wp, a2=6.91367384_wp )
!  Fitset: MD= -0.20303 MAD= 0.34880 RMSD= 0.54190
   case(p_df_x3lyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.55067492_wp, a1=0.19818545_wp, a2=5.61262748_wp )
!  Fitset: MD= -0.17267 MAD= 0.32089 RMSD= 0.51125
   case(p_df_xlyp); param = dftd_parameter ( & ! (SAW190103)
   &  s6=1.0000_wp, s8=1.51577878_wp, a1=0.10026585_wp, a2=5.37506460_wp )
!  Fitset: MD= -0.06152 MAD= 0.27318 RMSD= 0.38360
   end select
end subroutine get_d4eeqbjmbd_2019_parameter
! ========================================================================

!> get the unique identifier for most functionals, returns none if
!  the functional was not known at the time I implemented this mapping
function get_dfnum(df) result(num)
   integer(df_enum) :: num
   character(len=*),intent(in) :: df
   select case(df)
   case default;                             num = p_df_none
   case('hf');                               num = p_df_hf
   case('b-lyp','blyp');                     num = p_df_blyp
   case('bpbe');                             num = p_df_bpbe
   case('b-p','bp86','bp','b-p86');          num = p_df_bp
   case('bpw','b-pw');                       num = p_df_bpw
   case('lb94');                             num = p_df_lb94
   case('mpwlyp','mpw-lyp');                 num = p_df_mpwlyp
   case('mpwpw','mpw-pw','mpwpw91');         num = p_df_mpwpw
   case('o-lyp','olyp');                     num = p_df_olyp
   case('opbe');                             num = p_df_opbe
   case('pbe');                              num = p_df_pbe
   case('rpbe');                             num = p_df_rpbe
   case('revpbe');                           num = p_df_revpbe
   case('pw86pbe');                          num = p_df_pw86pbe
   case('rpw86pbe');                         num = p_df_rpw86pbe
   case('pw91');                             num = p_df_pw91
   case('pwp','pw-p','pw91p86');             num = p_df_pwp
   case('x-lyp','xlyp');                     num = p_df_xlyp
   case('b97');                              num = p_df_b97
   case('tpss');                             num = p_df_tpss
   case('revtpss');                          num = p_df_revtpss
   case('scan');                             num = p_df_scan
   case('b1lyp','b1-lyp');                   num = p_df_b1lyp
   case('b3-lyp','b3lyp');                   num = p_df_b3lyp
   case('bh-lyp','bhlyp');                   num = p_df_bhlyp
   case('b1p','b1-p','b1p86');               num = p_df_b1p
   case('b3p','b3-p','b3p86');               num = p_df_b3p
   case('b1pw','b1-pw','b1pw91');            num = p_df_b1pw
   case('b3pw','b3-pw','b3pw91');            num = p_df_b3pw
   case('o3-lyp','o3lyp');                   num = p_df_o3lyp
   case('revpbe0');                          num = p_df_revpbe0
   case('revpbe38');                         num = p_df_revpbe38
   case('pbe0');                             num = p_df_pbe0
   case('pwp1');                             num = p_df_pwp1
   case('pw1pw','pw1-pw');                   num = p_df_pw1pw
   case('mpw1pw','mpw1-pw','mpw1pw91');      num = p_df_mpw1pw
   case('mpw1lyp','mpw1-lyp');               num = p_df_mpw1lyp
   case('pw6b95');                           num = p_df_pw6b95
   case('tpssh');                            num = p_df_tpssh
   case('tpss0');                            num = p_df_tpss0
   case('x3-lyp','x3lyp');                   num = p_df_x3lyp
   case('m06l');                             num = p_df_m06l
   case('m06');                              num = p_df_m06
   case('m06-2x','m062x');                   num = p_df_m062x
   case('wb97','ωb97','omegab97');           num = p_df_wb97
   case('wb97x','ωb97x','omegab97x');        num = p_df_wb97x
   case('cam-b3lyp');                        num = p_df_camb3lyp
   case('lc-blyp');                          num = p_df_lcblyp
   case('lh07tsvwn','lh07t-svwn');           num = p_df_lh07tsvwn
   case('lh07ssvwn','lh07s-svwn');           num = p_df_lh07ssvwn
   case('lh12ctssirpw92','lh12ct-ssirpw92'); num = p_df_lh12ctssirpw92
   case('lh12ctssifpw92','lh12ct-ssifpw92'); num = p_df_lh12ctssifpw92
   case('lh14tcalpbe','lh14t-calpbe');       num = p_df_lh14tcalpbe
   case('b2plyp','b2-plyp');                 num = p_df_b2plyp
   case('b2gpplyp','b2gp-plyp');             num = p_df_b2gpplyp
   case('mpw2plyp');                         num = p_df_mpw2plyp
   case('pwpb95');                           num = p_df_pwpb95
   case('dsdblyp','dsd-blyp');               num = p_df_dsdblyp
   case('dsdpbe','dsd-pbe');                 num = p_df_dsdpbe
   case('dsdpbeb95','dsd-pbeb95');           num = p_df_dsdpbeb95
   case('dsdpbep86','dsd-pbep86');           num = p_df_dsdpbep86
   case('dsdsvwn','dsd-svwn');               num = p_df_dsdsvwn
   case('dodblyp','dod-blyp');               num = p_df_dodblyp
   case('dodpbe','dod-pbe');                 num = p_df_dodpbe
   case('dodpbeb95','dod-pbeb95');           num = p_df_dodpbeb95
   case('dodpbep86','dod-pbep86');           num = p_df_dodpbep86
   case('dodsvwn','dod-svwn');               num = p_df_dodsvwn
   case('pbe02','pbe0-2');                   num = p_df_pbe0_2
   case('pbe0dh','pbe0-dh');                 num = p_df_pbe0_dh
   case('hf-3c','hf3c');                     num = p_df_hf3c
   case('hf-3cv','hf3cv');                   num = p_df_hf3cv
   case('pbeh3c','pbeh-3c');                 num = p_df_pbeh3c
   case('b973c','b97-3c');                   num = p_df_b973c
   case('hsesol');                           num = p_df_hsesol
   case('pwgga');                            num = p_df_pwgga
   case('dftb3','dftb(3ob)');                num = p_df_dftb_3ob
   case('dftb(mio)');                        num = p_df_dftb_mio
   case('dftb(pbc)');                        num = p_df_dftb_pbc
   case('dftb(matsci)');                     num = p_df_dftb_matsci
   case('lc-dftb','dftb(ob2)');              num = p_df_dftb_ob2
   case('hcth120');                          num = p_df_hcth120
   case('ptpss');                            num = p_df_ptpss
   case('lc-wpbe','lcwpbe');                 num = p_df_lcwpbe
   case('bmk');                              num = p_df_bmk
   case('b1b95');                            num = p_df_b1b95
   case('bwb6k');                            num = p_df_pwb6k
   case('otpss');                            num = p_df_otpss
   case('ssb');                              num = p_df_ssb
   case('revssb');                           num = p_df_revssb
   case('pbesol');                           num = p_df_pbesol
   case('hse06');                            num = p_df_hse06
   case('pbexalpha');                        num = p_df_pbexalpha
   case('pbehpbe');                          num = p_df_pbehpbe
   case('hcth407');                          num = p_df_hcth407
   case('n12');                              num = p_df_n12
   case('pkzb');                             num = p_df_pkzb
   case('thcth','tauhctc');                  num = p_df_thcth
   case('m11l');                             num = p_df_m11l
   case('mn15l');                            num = p_df_mn15l
   case('mpwb1k');                           num = p_df_mpwb1k
   case('mpw1kcis');                         num = p_df_mpw1kcis
   case('mpwkcis1k');                        num = p_df_mpwkcis1k
   case('pbeh1pbe');                         num = p_df_pbeh1pbe
   case('pbe1kcis');                         num = p_df_pbe1kcis
   case('b97-1');                            num = p_df_b97_1
   case('b97-2');                            num = p_df_b97_2
   case('b98');                              num = p_df_b98
   case('hiss');                             num = p_df_hiss
   case('hse03');                            num = p_df_hse03
   case('revtpssh');                         num = p_df_revtpssh
   case('tpss1kcis');                        num = p_df_tpss1kcis
   case('m05');                              num = p_df_m05
   case('m052x','m05-2x');                   num = p_df_m052x
   case('m08hx','m08-hx');                   num = p_df_m08hx
   case('lcwhpbe','lc-whpbe');               num = p_df_lcwhpbe
   case('mn12l');                            num = p_df_mn12l
   case('tauhcthhyb');                       num = p_df_tauhcthhyb
   case('sogga11x');                         num = p_df_sogga11x
   case('n12sx');                            num = p_df_n12sx
   case('mn12sx');                           num = p_df_mn12sx
   case('mn15');                             num = p_df_mn15
   case('glyp','g-lyp');                     num = p_df_glyp
   case('revpbe0dh','revpbe0-dh');           num = p_df_revpbe0dh
   case('revtpss0');                         num = p_df_revtpss0
   case('revdsd-pbep86', 'revdsdpbep86');                    num = p_df_revdsdpbep86
   case('revdsd-pbe', 'revdsd-pbepbe', 'revdsdpbe', 'revdsdpbepbe');      num = p_df_revdsdpbe
   case('revdsd-blyp', 'revdsdblyp');                      num = p_df_revdsdblyp
   case('revdod-pbep86', 'revdodpbep86');                    num = p_df_revdodpbep86
   end select
end function get_dfnum

subroutine d4par(inp,dparam,lmbd,env)
   use iso_fortran_env, only : wp => real64
   use mctc_environment
   use dftd4, only : p_mbd_none,p_mbd_rpalike,p_mbd_exact_atm,p_mbd_approx_atm
   implicit none
   type(mctc_logger),intent(inout) :: env
   character(len=*),intent(in)   :: inp
   type(dftd_parameter),intent(out) :: dparam
   type(dftd_parameter),allocatable :: param
   integer,intent(in) :: lmbd
   character(len=20) :: dfa,bas
   character(len=80) :: ftmp,homedir
   character(len=99) :: dtmp
   logical :: ex
   integer :: i
   integer(df_enum) :: dfnum
!  init
   ex = .false.
   i = index(inp,'/')
   if (i.eq.0) then
      dfa = lowercase(inp)
      bas = ''
   else
      dfa = lowercase(inp(:i-1))
      bas = lowercase(inp(i+1:))
   endif
!  Get DFA parameter
   dfnum = get_dfnum(dfa)
   if (dfnum.gt.0) then
      if (lmbd.eq.p_mbd_approx_atm) then
         call get_d4eeqbjatm_2019_parameter(dfnum,p_bas_default,param)
      else if (lmbd.eq.p_mbd_none) then
         call get_d4eeqbj_2019_parameter(dfnum,p_bas_default,param)
      else if (lmbd.eq.p_mbd_rpalike) then
         call get_d4eeqbjmbd_2019_parameter(dfnum,p_bas_default,param)
      endif
      if (.not.allocated(param)) then
         call env%error(5,"No parameters for '"//trim(inp)//"' available")
      endif
   else
      call env%error(5,"Functional '"//trim(inp)//"' unknown")
   endif
   if (env%sane) dparam = param
end subroutine d4par

!> convert string to lower case
function lowercase(str) result(lcstr)
   implicit none
   character(len=*),intent(in)  :: str
   character(len=len_trim(str)) :: lcstr
   integer :: ilen,ioffset,iquote,i,iav,iqc

   ilen=len_trim(str)
   ioffset=iachar('A')-iachar('a')
   iquote=0
   lcstr=str
   do i=1,ilen
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

end module dfuncpar
