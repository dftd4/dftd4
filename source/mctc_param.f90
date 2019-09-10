! This file is part of dftd4.
!
! Copyright (C) 2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

module mctc_param
   use iso_fortran_env, wp => real64
   use mctc_econv

   implicit none

   public :: atomic_mass
   public :: covalent_radius
   public :: transition_metal_row
   public :: maingroup

   private

   integer, private, parameter :: max_elem = 118

!  covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
!  188-197), values for metals decreased by 10 %
   real(wp),public,parameter :: covalent_radius_2009(max_elem) = aatoau * [ &
   & 0.32_wp,0.46_wp, & ! H,He
   & 1.20_wp,0.94_wp,0.77_wp,0.75_wp,0.71_wp,0.63_wp,0.64_wp,0.67_wp, & ! Li-Ne
   & 1.40_wp,1.25_wp,1.13_wp,1.04_wp,1.10_wp,1.02_wp,0.99_wp,0.96_wp, & ! Na-Ar
   & 1.76_wp,1.54_wp, & ! K,Ca
   &                 1.33_wp,1.22_wp,1.21_wp,1.10_wp,1.07_wp, & ! Sc-
   &                 1.04_wp,1.00_wp,0.99_wp,1.01_wp,1.09_wp, & ! -Zn
   &                 1.12_wp,1.09_wp,1.15_wp,1.10_wp,1.14_wp,1.17_wp, & ! Ga-Kr
   & 1.89_wp,1.67_wp, & ! Rb,Sr
   &                 1.47_wp,1.39_wp,1.32_wp,1.24_wp,1.15_wp, & ! Y-
   &                 1.13_wp,1.13_wp,1.08_wp,1.15_wp,1.23_wp, & ! -Cd
   &                 1.28_wp,1.26_wp,1.26_wp,1.23_wp,1.32_wp,1.31_wp, & ! In-Xe
   & 2.09_wp,1.76_wp, & ! Cs,Ba
   &         1.62_wp,1.47_wp,1.58_wp,1.57_wp,1.56_wp,1.55_wp,1.51_wp, & ! La-Eu
   &         1.52_wp,1.51_wp,1.50_wp,1.49_wp,1.49_wp,1.48_wp,1.53_wp, & ! Gd-Yb
   &                 1.46_wp,1.37_wp,1.31_wp,1.23_wp,1.18_wp, & ! Lu-
   &                 1.16_wp,1.11_wp,1.12_wp,1.13_wp,1.32_wp, & ! -Hg
   &                 1.30_wp,1.30_wp,1.36_wp,1.31_wp,1.38_wp,1.42_wp, & ! Tl-Rn
   & 2.01_wp,1.81_wp, & ! Fr,Ra
   &         1.67_wp,1.58_wp,1.52_wp,1.53_wp,1.54_wp,1.55_wp,1.49_wp, & ! Ac-Am
   &         1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
   &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
   &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
   &                 1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp ] ! Nh-Og

!  based on "Atomic Radii of the Elements," M. Mantina, R. Valero,
!  C. J. Cramer, and D. G. Truhlar,
!  in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
!  edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
!  corrected Nov. 17, 2010 for the 92nd edition.
   real(wp),public,parameter :: covalent_radius_2010(94) = [ &
   & 0.32_wp,0.37_wp, & ! H,He
   & 1.30_wp,0.99_wp,0.84_wp,0.75_wp,0.71_wp,0.64_wp,0.60_wp,0.62_wp, & ! Li-Ne
   & 1.60_wp,1.40_wp,1.24_wp,1.14_wp,1.09_wp,1.04_wp,1.00_wp,1.01_wp, & ! Na-Ar
   & 2.00_wp,1.74_wp, & ! K,Ca
   &                 1.59_wp,1.48_wp,1.44_wp,1.30_wp,1.29_wp, & ! Sc-
   &                 1.24_wp,1.18_wp,1.17_wp,1.22_wp,1.20_wp, & ! -Zn
   &                 1.23_wp,1.20_wp,1.20_wp,1.18_wp,1.17_wp,1.16_wp, & ! Ga-Kr
   & 2.15_wp,1.90_wp, & ! Rb,Sr
   &                 1.76_wp,1.64_wp,1.56_wp,1.46_wp,1.38_wp, & ! Y-
   &                 1.36_wp,1.34_wp,1.30_wp,1.36_wp,1.40_wp, & ! -Cd
   &                 1.42_wp,1.40_wp,1.40_wp,1.37_wp,1.36_wp,1.36_wp, & ! In-Xe
   & 2.38_wp,2.06_wp, & ! Cs,Ba
   &         1.94_wp,1.84_wp,1.90_wp,1.88_wp,1.86_wp,1.85_wp,1.83_wp, & ! La-Eu
   &         1.82_wp,1.81_wp,1.80_wp,1.79_wp,1.77_wp,1.77_wp,1.78_wp, & ! Gd-Yb
   &                 1.74_wp,1.64_wp,1.58_wp,1.50_wp,1.41_wp, & ! Lu-
   &                 1.36_wp,1.32_wp,1.30_wp,1.30_wp,1.32_wp, & ! -Hg
   &                 1.44_wp,1.45_wp,1.50_wp,1.42_wp,1.48_wp,1.46_wp, & ! Tl-Rn
   & 2.42_wp,2.11_wp, & ! Fr,Ra
   &         2.01_wp,1.90_wp,1.84_wp,1.83_wp,1.80_wp,1.80_wp ] ! Ac-Pu

   real(wp),public,parameter :: atomic_mass_nist(max_elem) = amutoau * [ &
      &  1.00794075_wp,  4.00260193_wp,  6.94003660_wp,  9.01218307_wp,&
      & 10.81102805_wp, 12.01073590_wp, 14.00670321_wp, 15.99940492_wp,&
      & 18.99840316_wp, 20.18004638_wp, 22.98976928_wp, 24.30505162_wp,&
      & 26.98153853_wp, 28.08549871_wp, 30.97376200_wp, 32.06478741_wp,&
      & 35.45293758_wp, 39.94779856_wp, 39.09830091_wp, 40.07802251_wp,&
      & 44.95590828_wp, 47.86674496_wp, 50.94146504_wp, 51.99613176_wp,&
      & 54.93804391_wp, 55.84514443_wp, 58.93319429_wp, 58.69334711_wp,&
      & 63.54603995_wp, 65.37778253_wp, 69.72306607_wp, 72.62755016_wp,&
      & 74.92159457_wp, 78.95938856_wp, 79.90352778_wp, 83.79800000_wp,&
      & 85.46766360_wp, 87.61664447_wp, 88.90584030_wp, 91.22364160_wp,&
      & 92.90637300_wp, 95.95978854_wp, 97.90721240_wp,101.06494014_wp,&
      &102.90549800_wp,106.41532751_wp,107.86814963_wp,112.41155782_wp,&
      &114.81808663_wp,118.71011259_wp,121.75978367_wp,127.60312648_wp,&
      &126.90447190_wp,131.29276145_wp,132.90545196_wp,137.32689163_wp,&
      &138.90546887_wp,140.11573074_wp,140.90765760_wp,144.24159603_wp,&
      &144.91275590_wp,150.36635571_wp,151.96437813_wp,157.25213065_wp,&
      &158.92535470_wp,162.49947282_wp,164.93032880_wp,167.25908265_wp,&
      &168.93421790_wp,173.05415017_wp,174.96681496_wp,178.48497872_wp,&
      &180.94787564_wp,183.84177755_wp,186.20670455_wp,190.22485963_wp,&
      &192.21605165_wp,195.08445686_wp,196.96656879_wp,200.59916703_wp,&
      &204.38341284_wp,207.21690806_wp,208.98039910_wp,208.98243080_wp,&
      &209.98714790_wp,222.01757820_wp,223.01973600_wp,226.02541030_wp,&
      &227.02775230_wp,232.03805580_wp,231.03588420_wp,238.02891046_wp,&
      &237.04817360_wp,244.06420530_wp,243.06138130_wp,247.07035410_wp,&
      &247.07030730_wp,251.07958860_wp,252.08298000_wp,257.09510610_wp,&
      &258.09843150_wp,259.10103000_wp,262.10961000_wp,267.12179000_wp,&
      &269.12791000_wp,271.13393000_wp,270.13336000_wp,276.14846000_wp,&
      &276.15159000_wp,280.16131000_wp,282.16912000_wp,284.17416000_wp,&
      &284.17873000_wp,289.19042000_wp,288.19274000_wp,293.20449000_wp,&
      &292.20746000_wp,294.21392000_wp]

!> pauling EN's
   real(wp),public,parameter :: pauling_en(max_elem) = [ &
   & 2.20_wp,3.00_wp, & ! H,He
   & 0.98_wp,1.57_wp,2.04_wp,2.55_wp,3.04_wp,3.44_wp,3.98_wp,4.50_wp, & ! Li-Ne
   & 0.93_wp,1.31_wp,1.61_wp,1.90_wp,2.19_wp,2.58_wp,3.16_wp,3.50_wp, & ! Na-Ar
   & 0.82_wp,1.00_wp, & ! K,Ca
   &                 1.36_wp,1.54_wp,1.63_wp,1.66_wp,1.55_wp, & ! Sc-
   &                 1.83_wp,1.88_wp,1.91_wp,1.90_wp,1.65_wp, & ! -Zn
   &                 1.81_wp,2.01_wp,2.18_wp,2.55_wp,2.96_wp,3.00_wp, & ! Ga-Kr
   & 0.82_wp,0.95_wp, & ! Rb,Sr
   &                 1.22_wp,1.33_wp,1.60_wp,2.16_wp,1.90_wp, & ! Y-
   &                 2.20_wp,2.28_wp,2.20_wp,1.93_wp,1.69_wp, & ! -Cd
   &                 1.78_wp,1.96_wp,2.05_wp,2.10_wp,2.66_wp,2.60_wp, & ! In-Xe
   & 0.79_wp,0.89_wp, & ! Cs,Ba
   &         1.10_wp,1.12_wp,1.13_wp,1.14_wp,1.15_wp,1.17_wp,1.18_wp, & ! La-Eu
   &         1.20_wp,1.21_wp,1.22_wp,1.23_wp,1.24_wp,1.25_wp,1.26_wp, & ! Gd-Yb
   &                 1.27_wp,1.30_wp,1.50_wp,2.36_wp,1.90_wp, & ! Lu-
   &                 2.20_wp,2.20_wp,2.28_wp,2.54_wp,2.00_wp, & ! -Hg
   &                 1.62_wp,2.33_wp,2.02_wp,2.00_wp,2.20_wp,2.20_wp, & ! Tl-Rn
   ! only dummies below
   & 1.50_wp,1.50_wp, & ! Fr,Ra
   &         1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Ac-Am
   &         1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Cm-No
   &                 1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Rf-
   &                 1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Rf-Cn
   &                 1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp ] ! Nh-Og

!  Semiempirical Evaluation of the GlobalHardness of the Atoms of 103
!  Elements of the Periodic Table Using the Most Probable Radii as
!  their Size Descriptors DULAL C. GHOSH, NAZMUL ISLAM 2009 in
!  Wiley InterScience (www.interscience.wiley.com).
!  DOI 10.1002/qua.22202
!  values in the paper multiplied by two because
!  (ii:ii)=(IP-EA)=d^2 E/dN^2 but the hardness
!  definition they use is 1/2d^2 E/dN^2 (in Eh)
   real(wp),public, parameter :: chemical_hardness(max_elem) = [ &
  &0.47259288_Wp,0.92203391_Wp,0.17452888_Wp,0.25700733_Wp,0.33949086_wp,0.42195412_wp, & ! H-C
  &0.50438193_Wp,0.58691863_Wp,0.66931351_Wp,0.75191607_Wp,0.17964105_wp,0.22157276_wp, & ! N-Mg
  &0.26348578_Wp,0.30539645_Wp,0.34734014_Wp,0.38924725_Wp,0.43115670_wp,0.47308269_wp, & ! Al-Ar
  &0.17105469_Wp,0.20276244_Wp,0.21007322_Wp,0.21739647_Wp,0.22471039_wp,0.23201501_wp, & ! Ca-Cr
  &0.23933969_Wp,0.24665638_Wp,0.25398255_Wp,0.26128863_Wp,0.26859476_wp,0.27592565_wp, & ! Mn-Zn
  &0.30762999_Wp,0.33931580_Wp,0.37235985_Wp,0.40273549_Wp,0.43445776_wp,0.46611708_wp, & ! Ga-Kr
  &0.15585079_Wp,0.18649324_Wp,0.19356210_Wp,0.20063311_Wp,0.20770522_wp,0.21477254_wp, & ! Rb-Mo
  &0.22184614_Wp,0.22891872_Wp,0.23598621_Wp,0.24305612_Wp,0.25013018_wp,0.25719937_wp, & ! Tc-Cd
  &0.28784780_Wp,0.31848673_Wp,0.34912431_Wp,0.37976593_Wp,0.41040808_wp,0.44105777_wp, & ! In-Xe
  &0.05019332_Wp,0.06762570_Wp,0.08504445_Wp,0.10247736_Wp,0.11991105_wp,0.13732772_wp, & ! Cs-Nd
  &0.15476297_Wp,0.17218265_Wp,0.18961288_Wp,0.20704760_Wp,0.22446752_wp,0.24189645_wp, & ! Pm-Dy
  &0.25932503_Wp,0.27676094_Wp,0.29418231_Wp,0.31159587_Wp,0.32902274_wp,0.34592298_wp, & ! Ho-Hf
  &0.36388048_Wp,0.38130586_Wp,0.39877476_Wp,0.41614298_Wp,0.43364510_wp,0.45104014_wp, & ! Ta-Pt
  &0.46848986_Wp,0.48584550_Wp,0.12526730_Wp,0.14268677_Wp,0.16011615_wp,0.17755889_wp, & ! Au-Po
  &0.19497557_Wp,0.21240778_Wp,0.07263525_Wp,0.09422158_Wp,0.09920295_wp,0.10418621_wp, & ! At-Th
  &0.14235633_Wp,0.16394294_Wp,0.18551941_Wp,0.22370139_Wp,0.00000000_wp,0.00000000_wp, & ! Pa-Cm
  &0.00000000_Wp,0.00000000_Wp,0.00000000_Wp,0.00000000_Wp,0.00000000_wp,0.00000000_wp, & ! Bk-No
  &0.00000000_Wp,0.00000000_Wp,0.00000000_Wp,0.00000000_Wp,0.00000000_wp,0.00000000_wp, & ! Rf-Mt
  &0.00000000_Wp,0.00000000_Wp,0.00000000_Wp,0.00000000_Wp,0.00000000_wp,0.00000000_wp, & ! Ds-Mc
  &0.00000000_Wp,0.00000000_Wp,0.00000000_Wp,0.00000000_Wp ] ! Lv-Og

   real(wp),public,parameter :: def2ecp_nuclear_charges(max_elem) = (/ &
   &   1,                                                 2,  & ! H-He
   &   3, 4,                               5, 6, 7, 8, 9,10,  & ! Li-Ne
   &  11,12,                              13,14,15,16,17,18,  & ! Na-Ar
   &  19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,  & ! K-Kr
   &   9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,  & ! Rb-Xe
   &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Cs-Lu
   &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26, & ! Hf-Rn
   !  just copy & paste from above
   &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Fr-Lr
   &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 /) ! Rf-Og

!  PBE0/def2-QZVP atomic values calculated by S. Grimme in Gaussian (2010)
!  rare gases recalculated by J. Mewes with PBE0/aug-cc-pVQZ in Dirac (2018)
!  He: 3.4698 -> 3.5544, Ne: 3.1036 -> 3.7943, Ar: 5.6004 -> 5.6638,
!  Kr: 6.1971 -> 6.2312, Xe: 7.5152 -> 8.8367
!  not replaced but recalculated (PBE0/cc-pVQZ) were
!   H: 8.0589 ->10.9359, Li:29.0974 ->39.7226, Be:14.8517 ->17.7460
!  also new super heavies Cn,Nh,Fl,Lv,Og
   real(wp),public,parameter :: r4_over_r2(max_elem) = (/  &
   &  8.0589_wp, 3.4698_wp, & ! H,He
   & 29.0974_wp,14.8517_wp,11.8799_wp, 7.8715_wp, 5.5588_wp, 4.7566_wp, 3.8025_wp, 3.1036_wp, & ! Li-Ne
   & 26.1552_wp,17.2304_wp,17.7210_wp,12.7442_wp, 9.5361_wp, 8.1652_wp, 6.7463_wp, 5.6004_wp, & ! Na-Ar
   & 29.2012_wp,22.3934_wp, & ! K,Ca
   &         19.0598_wp,16.8590_wp,15.4023_wp,12.5589_wp,13.4788_wp, & ! Sc-
   &         12.2309_wp,11.2809_wp,10.5569_wp,10.1428_wp, 9.4907_wp, & ! -Zn
   &                 13.4606_wp,10.8544_wp, 8.9386_wp, 8.1350_wp, 7.1251_wp, 6.1971_wp, & ! Ga-Kr
   & 30.0162_wp,24.4103_wp, & ! Rb,Sr
   &         20.3537_wp,17.4780_wp,13.5528_wp,11.8451_wp,11.0355_wp, & ! Y-
   &         10.1997_wp, 9.5414_wp, 9.0061_wp, 8.6417_wp, 8.9975_wp, & ! -Cd
   &                 14.0834_wp,11.8333_wp,10.0179_wp, 9.3844_wp, 8.4110_wp, 7.5152_wp, & ! In-Xe
   & 32.7622_wp,27.5708_wp, & ! Cs,Ba
   &         23.1671_wp,21.6003_wp,20.9615_wp,20.4562_wp,20.1010_wp,19.7475_wp,19.4828_wp, & ! La-Eu
   &         15.6013_wp,19.2362_wp,17.4717_wp,17.8321_wp,17.4237_wp,17.1954_wp,17.1631_wp, & ! Gd-Yb
   &         14.5716_wp,15.8758_wp,13.8989_wp,12.4834_wp,11.4421_wp, & ! Lu-
   &         10.2671_wp, 8.3549_wp, 7.8496_wp, 7.3278_wp, 7.4820_wp, & ! -Hg
   &                 13.5124_wp,11.6554_wp,10.0959_wp, 9.7340_wp, 8.8584_wp, 8.0125_wp, & ! Tl-Rn
   & 29.8135_wp,26.3157_wp, & ! Fr,Ra
   &         19.1885_wp,15.8542_wp,16.1305_wp,15.6161_wp,15.1226_wp,16.1576_wp, 0.0000_wp, & ! Ac-Am
   &          0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, & ! Cm-No
   &          0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, & ! Lr-
   &          0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 5.4929_wp, & ! -Cn
   &                  6.7286_wp, 6.5144_wp,10.9169_wp,10.3600_wp, 9.4723_wp, 8.6641_wp /) ! Nh-Og
   integer,private :: idum
   real(wp),public,parameter :: sqrt_z_r4_over_r2(max_elem) = &
   &  sqrt(0.5_wp*(r4_over_r2*(/(sqrt(real(idum,wp)),idum=1,max_elem)/)))

   !> D3 covalent radii used to construct the coordination number
   real(wp),public,parameter :: covalent_radius_d3(max_elem) = &
      & 4.0_wp / 3.0_wp * covalent_radius_2009

contains

pure elemental function covalent_radius(iatom) result(rcov)
   implicit none
   integer,intent(in) :: iatom
   real(wp) :: rcov

   if (iatom > 0 .and. iatom <= max_elem) then
      rcov = covalent_radius_2009(iatom)
   else
      rcov = 0.0_wp
   endif
end function covalent_radius

pure elemental function atomic_mass(iatom) result(mass)
   implicit none
   integer,intent(in) :: iatom
   real(wp) :: mass

   if (iatom > 0 .and. iatom <= max_elem) then
      mass = atomic_mass_nist(iatom)
   else
      mass = 0.0_wp
   endif
end function atomic_mass

logical function maingroup(i)
   integer, intent(in) :: i
   logical main_group(107)
   data main_group /&
      &  2*.true.,&                              ! H  - He
      &  8*.true.,&                              ! Li - Ne
      &  8*.true.,&                              ! Na - Ar
      &  2*.true., 9*.false., 7*.true.,&         ! K  - Kr
      &  2*.true., 9*.false., 7*.true.,&         ! Rb - Xe
      &  2*.true.,23*.false., 7*.true.,&         ! Cs - Rn
      & 21*.true.                               / ! Fr - Tv

   maingroup = main_group(i)

end function maingroup

pure elemental function transition_metal_row(element) result(tmmetal)
   integer, intent(in) :: element
   integer :: tmmetal

   select case(element)
   case default
      tmmetal = 0
   case(21:29)
      !if(i.gt.20.and.i.lt.30) j=1
      !if(i.gt.38.and.i.lt.48) j=2
      !if(i.gt.56.and.i.lt.80) j=3
      tmmetal = 1
   case(39:47)
      tmmetal = 2
   case(57:79)
      tmmetal = 3
   end select

end function transition_metal_row

end module mctc_param
