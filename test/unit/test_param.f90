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

module test_param
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use dftd4_param
   use dftd4
   implicit none
   private

   public :: collect_param

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   type(realspace_cutoff), parameter :: cutoff = &
      & realspace_cutoff(cn=30_wp, disp2=60.0_wp, disp3=15.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_param(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("rational-damping", test_rational_damping) &
      & ]

end subroutine collect_param


subroutine test_dftd4_gen(error, mol, param, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Expected dispersion energy
   real(wp), intent(in) :: ref

   type(d4_model) :: d4
   real(wp) :: energy

   call new_d4_model(d4, mol)
   call get_dispersion(mol, d4, param, cutoff, energy)

   call check(error, energy, ref, thr=thr)
   if (allocated(error)) then
      print'(es21.14)',energy
   end if

end subroutine test_dftd4_gen


subroutine test_rational_damping(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: func(*) = [character(len=32)::&
      & 'hf', 'b-lyp', 'bpbe', 'b-p', 'bpw', 'lb94', 'mpwlyp', 'mpwpw', &
      & 'o-lyp', 'opbe', 'pbe', 'rpbe', 'revpbe', 'pw86pbe', 'rpw86pbe', &
      & 'pw91', 'pwp', 'x-lyp', 'b97', 'tpss', 'revtpss', 'scan', 'rscan', &
      & 'r2scan', 'b1lyp', 'b3-lyp', 'bh-lyp', 'b1p', 'b3p', 'b1pw', 'b3pw', &
      & 'o3-lyp', 'revpbe0', 'revpbe38', 'pbe0', 'pwp1', 'pw1pw', 'mpw1pw', &
      & 'mpw1lyp', 'pw6b95', 'tpssh', 'tpss0', 'x3-lyp', 'm06l', 'm06', &
      & 'wb97', 'wb97x', 'cam-b3lyp', 'lc-blyp', 'lh07tsvwn', 'lh07ssvwn', &
      & 'lh12ctssirpw92', 'lh12ctssifpw92', 'lh14tcalpbe', 'lh20t', &
      & 'b2plyp', 'b2gpplyp', 'mpw2plyp', 'pwpb95', 'dsdblyp', 'dsdpbe', &
      & 'dsdpbeb95', 'dsdpbep86', 'dsdsvwn', 'dodblyp', 'dodpbe', 'dodpbeb95', &
      & 'dodpbep86', 'dodsvwn', 'pbe02', 'pbe0dh', 'dftb(3ob)', 'dftb(mio)', &
      & 'dftb(pbc)', 'dftb(matsci)', 'dftb(ob2)', 'b1b95', 'glyp', 'revpbe0dh', &
      & 'revtpss0', 'revdsd-pbep86', 'revdsd-pbe', 'revdsd-blyp', &
      & 'revdod-pbep86', 'b97m', 'wb97m', 'pbesol', 'am05', 'mn12sx', &
      & 'hse03', 'hse06', 'hse12', 'hse12s', 'hsesol', 'r2scanh', 'r2scan0', &
      & 'r2scan50']
   real(wp), parameter :: ref(*) = [&
      &-2.82477943738524E-1_wp,-1.83931932418458E-1_wp,-1.67883732655799E-1_wp, &
      &-1.19260344932822E-1_wp,-1.73553644083274E-1_wp,-4.64187011491173E-1_wp, &
      &-1.33350632679087E-1_wp,-1.29120541441992E-1_wp,-4.23162979118617E-1_wp, &
      &-4.19510902600366E-1_wp,-8.78277690565286E-2_wp,-2.88574076971878E-1_wp, &
      &-2.59731167367680E-1_wp,-9.87868660020315E-2_wp,-9.80588492098119E-2_wp, &
      &-7.25055148336794E-2_wp,-3.38536837308831E-2_wp,-1.99101155252712E-1_wp, &
      &-1.51912676673518E-1_wp,-1.17214547943120E-1_wp,-9.43596374590356E-2_wp, &
      &-1.84970292487715E-2_wp,-3.25233368429052E-2_wp,-2.79219863681671E-2_wp, &
      &-1.41170942861054E-1_wp,-1.36036601060694E-1_wp,-9.86753156979172E-2_wp, &
      &-8.99492843498988E-2_wp,-9.44139526338724E-2_wp,-1.35742320819228E-1_wp, &
      &-1.32216888108482E-1_wp,-1.14737687058121E-1_wp,-1.86535362608147E-1_wp, &
      &-1.46887479247351E-1_wp,-7.64961642342119E-2_wp,-3.19815918817886E-2_wp, &
      &-6.61290132842028E-2_wp,-1.05230601191541E-1_wp,-1.02669100564886E-1_wp, &
      &-7.17939223473286E-2_wp,-1.08455785302610E-1_wp,-1.01080258586891E-1_wp, &
      &-1.12546243554503E-1_wp,-1.23576092834163E-2_wp,-2.00031673564059E-2_wp, &
      &-8.37504030553571E-3_wp,-1.58515233166313E-2_wp,-7.77668958299319E-2_wp, &
      &-1.29932494731440E-2_wp,-2.07804102792114E-1_wp,-2.62377727009727E-1_wp, &
      &-3.39614227696319E-1_wp,-3.74698168802945E-1_wp,-8.84516858200918E-2_wp, &
      &-5.21275385809976E-2_wp,-6.24504150834387E-2_wp,-3.94277426207638E-2_wp, &
      &-4.14986339439292E-2_wp,-6.11949139991988E-2_wp,-3.85836296996462E-2_wp, &
      &-4.77329014912804E-2_wp,-4.36022696940550E-2_wp,-1.90151395280036E-2_wp, &
      &-2.31910927374136E-2_wp,-8.51189317714713E-2_wp,-6.56921633268103E-2_wp, &
      &-5.96170926790428E-2_wp,-5.17434306954188E-2_wp,-3.15736134106408E-2_wp, &
      &-8.31688594037642E-3_wp,-4.85454966544324E-2_wp,-6.17262149162204E-2_wp, &
      &-4.55955216246633E-2_wp,-5.85633742996495E-2_wp,-7.13638871727705E-2_wp, &
      &-4.28867617727652E-2_wp,-1.03264577526582E-1_wp,-2.71073701507995E-1_wp, &
      &-1.02702213711187E-1_wp,-8.21356035564771E-2_wp,-5.65853022691401E-2_wp, &
      &-8.56296238144821E-2_wp,-8.91479952796072E-2_wp,-6.13524150208560E-2_wp, &
      &-1.24018193150396E-1_wp,-1.05459213596423E-1_wp,-3.62056545113500E-2_wp, &
      &-3.62056561428035E-2_wp,-2.40058303330640E-2_wp,-6.96544056298022E-2_wp, &
      &-7.14846530839187E-2_wp,-6.94907084631183E-2_wp,-6.20062460023045E-2_wp, &
      &-5.03611953456781E-2_wp,-2.92623319798973E-2_wp,-3.16651213971901E-2_wp, &
      &-3.45136901370079E-2_wp]
   class(damping_param), allocatable :: param
   type(structure_type) :: mol
   integer :: ii

   call get_structure(mol, "UPU23", "0a")
   do ii = 1, size(func)
      call get_rational_damping(trim(func(ii)), param, s9=1.0_wp)
      call check(error, allocated(param))
      if (allocated(error)) exit
      call test_dftd4_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_rational_damping

end module test_param
