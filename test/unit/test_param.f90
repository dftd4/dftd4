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
      &-2.82477942248756E-1_wp,-1.83931931447376E-1_wp,-1.67883731806544E-1_wp, &
      &-1.19260344286630E-1_wp,-1.73553643190541E-1_wp,-4.64187008858948E-1_wp, &
      &-1.33350631907425E-1_wp,-1.29120540723649E-1_wp,-4.23162977083251E-1_wp, &
      &-4.19510900741786E-1_wp,-8.78277686099131E-2_wp,-2.88574075471549E-1_wp, &
      &-2.59731166101835E-1_wp,-9.87868654968457E-2_wp,-9.80588487018907E-2_wp, &
      &-7.25055144755650E-2_wp,-3.38536835705987E-2_wp,-1.99101153940959E-1_wp, &
      &-1.51912675830051E-1_wp,-1.17214547334950E-1_wp,-9.43596369853342E-2_wp, &
      &-1.84970291718116E-2_wp,-3.25233366958794E-2_wp,-2.79219862478753E-2_wp, &
      &-1.41170942097673E-1_wp,-1.36036600333643E-1_wp,-9.86753151279879E-2_wp, &
      &-8.99492838831500E-2_wp,-9.44139521417514E-2_wp,-1.35742320110483E-1_wp, &
      &-1.32216887418510E-1_wp,-1.14737686310415E-1_wp,-1.86535361600070E-1_wp, &
      &-1.46887478462845E-1_wp,-7.64961638449203E-2_wp,-3.19815917382690E-2_wp, &
      &-6.61290129604851E-2_wp,-1.05230600646517E-1_wp,-1.02669099979520E-1_wp, &
      &-7.17939219769430E-2_wp,-1.08455784745289E-1_wp,-1.01080258057282E-1_wp, &
      &-1.12546242873133E-1_wp,-1.23576092386848E-2_wp,-2.00031672761610E-2_wp, &
      &-8.37504026938291E-3_wp,-1.58515232511485E-2_wp,-7.77668954242432E-2_wp, &
      &-1.29932494116969E-2_wp,-2.07804101611958E-1_wp,-2.62377725485470E-1_wp, &
      &-3.39614225692246E-1_wp,-3.74698166577280E-1_wp,-8.84516853613503E-2_wp, &
      &-5.21275383612525E-2_wp,-6.24504147646583E-2_wp,-3.94277424205698E-2_wp, &
      &-4.14986337497841E-2_wp,-6.11949137677711E-2_wp,-3.85836295170037E-2_wp, &
      &-4.77329012547881E-2_wp,-4.36022695279866E-2_wp,-1.90151394463137E-2_wp, &
      &-2.31910926279841E-2_wp,-8.51189313158795E-2_wp,-6.56921629858802E-2_wp, &
      &-5.96170924208417E-2_wp,-5.17434304314999E-2_wp,-3.15736132524860E-2_wp, &
      &-8.31688591017636E-3_wp,-4.85454964254385E-2_wp,-6.17262146525898E-2_wp, &
      &-4.55955214285284E-2_wp,-5.85633740272026E-2_wp,-7.13638868029425E-2_wp, &
      &-4.28867615756140E-2_wp,-1.03264576995975E-1_wp,-2.71073699938782E-1_wp, &
      &-1.02702213161577E-1_wp,-8.21356031484763E-2_wp,-5.65853020266326E-2_wp, &
      &-8.56296234267241E-2_wp,-8.91479948664164E-2_wp,-6.13524147568788E-2_wp, &
      &-1.24018192547142E-1_wp,-1.05459213220861E-1_wp,-3.62056543364222E-2_wp, &
      &-3.62056559678757E-2_wp,-2.40058302367699E-2_wp,-6.96544052724497E-2_wp, &
      &-7.14846527169669E-2_wp,-6.94907081073225E-2_wp,-6.20062456876537E-2_wp, &
      &-5.03611950945288E-2_wp,-2.92623318488036E-2_wp,-3.16651212523792E-2_wp, &
      &-3.45136899744303E-2_wp]
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
