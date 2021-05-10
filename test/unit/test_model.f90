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

module test_model
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, &
      & test_failed
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use dftd4_charge, only : get_charges
   use dftd4_cutoff, only : get_lattice_points
   use dftd4_data, only : get_covalent_rad
   use dftd4_ncoord, only : get_coordination_number
   use dftd4_model
   implicit none
   private

   public :: collect_model

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_model(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("gw-mb01", test_gw_mb01), &
      & new_unittest("gw-mb02", test_gw_mb02), &
      & new_unittest("gw-mb03", test_gw_mb03), &
      & new_unittest("dgw-mb04", test_dgw_mb04), &
      & new_unittest("dgw-mb05", test_dgw_mb05), &
      & new_unittest("dgw-mb06", test_dgw_mb06), &
      & new_unittest("gw-gfn2", test_gw_mb07), &
      & new_unittest("dgw-gfn2", test_dgw_mb08) &
      & ]

end subroutine collect_model


subroutine test_gw_gen(error, mol, ref, with_cn, with_q)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Reference Gaussian weights
   real(wp), intent(in) :: ref(:, :)

   logical, intent(in) :: with_cn
   logical, intent(in) :: with_q

   type(d4_model) :: d4
   real(wp), allocatable :: cn(:), q(:), gwvec(:, :)
   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :)

   call new_d4_model(d4, mol)

   allocate(cn(mol%nat), q(mol%nat), gwvec(maxval(d4%ref), mol%nat))
   cn(:) = 0.0_wp
   q(:) = 0.0_wp

   if (with_cn) then
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call get_coordination_number(mol, lattr, cutoff, d4%rcov, d4%en, cn)
   end if
   if (with_q) then
      call get_charges(mol, q)
   end if

   call d4%weight_references(mol, cn, q, gwvec)

   if (any(abs(gwvec - ref) > thr2)) then
      call test_failed(error, "Gaussian weights do not match")
      where(abs(gwvec) < thr) gwvec = 0.0_wp
      print'(3(es20.13,"_wp,"), " &")', gwvec
   end if

end subroutine test_gw_gen


subroutine test_dgw_gen(error, mol, with_cn, with_q)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   logical, intent(in) :: with_cn
   logical, intent(in) :: with_q

   integer :: iat, mref
   type(d4_model) :: d4
   real(wp), allocatable :: cn(:), q(:), gwvec(:, :), gwdcn(:, :), gwdq(:, :)
   real(wp), allocatable :: gwr(:, :), gwl(:, :), numdcn(:, :), numdq(:, :)
   real(wp), parameter :: cutoff = 30.0_wp, lattr(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp

   call new_d4_model(d4, mol)

   mref = maxval(d4%ref)
   allocate(cn(mol%nat), q(mol%nat), gwvec(mref, mol%nat), &
      & gwdcn(mref, mol%nat), gwdq(mref, mol%nat), &
      & gwr(mref, mol%nat), gwl(mref, mol%nat), &
      & numdcn(mref, mol%nat), numdq(mref, mol%nat))
   cn(:) = 0.0_wp
   q(:) = 0.0_wp

   if (with_cn) then
      call get_coordination_number(mol, lattr, cutoff, d4%rcov, d4%en, cn)
   end if
   if (with_q) then
      call get_charges(mol, q)
   end if

   if (with_cn) then
      do iat = 1, mol%nat
         cn(iat) = cn(iat) + step
         call d4%weight_references(mol, cn, q, gwr)
         cn(iat) = cn(iat) - 2*step
         call d4%weight_references(mol, cn, q, gwl)
         cn(iat) = cn(iat) + step
         gwdcn(:, :) = 0.5_wp*(gwr - gwl)/step
         numdcn(:, iat) = gwdcn(:, iat)
         gwdcn(:, iat) = 0.0_wp
         if (any(abs(gwdcn) > thr)) then
            call test_failed(error, "Unexpected non-zero gradient element found")
            exit
         end if
      end do
      if (allocated(error)) return
   end if

   if (with_q) then
      do iat = 1, mol%nat
         q(iat) = q(iat) + step
         call d4%weight_references(mol, cn, q, gwr)
         q(iat) = q(iat) - 2*step
         call d4%weight_references(mol, cn, q, gwl)
         q(iat) = q(iat) + step
         gwdq(:, :) = 0.5_wp*(gwr - gwl)/step
         numdq(:, iat) = gwdq(:, iat)
         gwdq(:, iat) = 0.0_wp
         if (any(abs(gwdq) > thr)) then
            call test_failed(error, "Unexpected non-zero gradient element found")
            exit
         end if
      end do
      if (allocated(error)) return
   end if

   call d4%weight_references(mol, cn, q, gwvec, gwdcn, gwdq)

   if (with_cn .and. any(abs(gwdcn - numdcn) > thr2)) then
      call test_failed(error, "Gaussian weights derivatives do not match")
      print'(3es21.14)', gwdcn
      print'("---")'
      print'(3es21.14)', numdcn
      print'("---")'
      print'(3es21.14)', gwdcn - numdcn
   end if

   if (with_q .and. any(abs(gwdq - numdq) > thr2)) then
      call test_failed(error, "Gaussian weights derivatives do not match")
      print'(3es21.14)', gwdq
      print'("---")'
      print'(3es21.14)', numdq
      print'("---")'
      print'(3es21.14)', gwdq - numdq
   end if

end subroutine test_dgw_gen


subroutine test_gw_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(5, 16) = reshape([&
      & 1.2088598797363E-11_wp, 1.0375924114817E+00_wp, 2.4440673480025E-09_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 5.5306528158155E-03_wp, &
      & 9.9446934718418E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 3.1404786132959E-06_wp, 6.3956949849691E-02_wp, &
      & 4.7846310051123E-01_wp, 2.9468375847232E-01_wp, 0.0000000000000E+00_wp, &
      & 7.4311632283534E-05_wp, 9.9992568836772E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 4.4395785875153E-02_wp, &
      & 8.5383544991033E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.1950287975766E-02_wp, 9.8804971202423E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 5.6993600838495E-05_wp, 9.9994300639916E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 5.6167990420147E-07_wp, &
      & 2.9078974710698E-02_wp, 6.1859557270865E-01_wp, 1.5104290756010E-01_wp, &
      & 0.0000000000000E+00_wp, 2.4997061606460E-26_wp, 4.9771731646875E-14_wp, &
      & 1.3962125880500E-04_wp, 6.8484995443647E-01_wp, 1.6516600010230E-12_wp, &
      & 1.9880179303418E-02_wp, 9.8011982069658E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.1966391399283E-02_wp, &
      & 9.8803360860072E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.2874084623676E-05_wp, 9.7726740452908E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.4569926587417E-47_wp, 4.3446630202006E-29_wp, 2.5850639306886E-15_wp, &
      & 4.1492997927124E-06_wp, 8.6768202812026E-01_wp, 1.6508416802282E-24_wp, &
      & 4.3169770775040E-12_wp, 3.6465263662462E-04_wp, 8.8903502908917E-01_wp, &
      & 1.2720250177275E-05_wp, 8.4045163736307E-30_wp, 2.3928873165684E-16_wp, &
      & 1.4066071431515E-05_wp, 6.8494819900529E-01_wp, 1.1610740025363E-14_wp, &
      & 2.6396812522623E-41_wp, 6.3321179479930E-24_wp, 1.9575180950269E-10_wp, &
      & 1.0064486114128E+00_wp, 0.0000000000000E+00_wp], shape(ref))

   call get_structure(mol, "MB16-43", "01")
   call test_gw_gen(error, mol, ref, with_cn=.true., with_q=.false.)

end subroutine test_gw_mb01


subroutine test_gw_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(5, 16) = reshape([&
      & 8.1529927732913E-01_wp, 2.2414425943560E-03_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.0229633868523E+00_wp, &
      & 2.0291440964487E-03_wp, 4.4863014126871E-10_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.1578430860412E+00_wp, 2.0553787043636E-03_wp, &
      & 3.1900095420129E-10_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.3687710077901E+00_wp, 8.2192922249920E-03_wp, 5.9806186563013E-08_wp, &
      & 1.4414431799851E-03_wp, 0.0000000000000E+00_wp, 9.1646420559464E-01_wp, &
      & 5.3618567209686E-03_wp, 1.7524886447677E-08_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 5.5812105408862E-01_wp, 1.5344013397624E-03_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 7.5659187685708E-01_wp, 2.0800426377008E-03_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 7.7216223321749E-01_wp, &
      & 2.1228490781405E-03_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 9.9399937131066E-01_wp, 2.8186797340819E-03_wp, &
      & 7.1776130192014E-10_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 5.5220460326555E-01_wp, 1.5181356748084E-03_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.1425732736036E+00_wp, &
      & 2.0284311993372E-03_wp, 3.1484304344932E-10_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 8.1572048274912E-01_wp, 6.4067286692506E-03_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.2272138452089E+00_wp, 1.3612868750987E-02_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 8.3285820461366E-01_wp, &
      & 2.2897160672036E-03_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 7.4607914206254E-01_wp, 2.0511407458350E-03_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.0523038115105E+00_wp, 2.0872155589017E-03_wp, 4.6145584751225E-10_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp], shape(ref))

   call get_structure(mol, "MB16-43", "02")
   call test_gw_gen(error, mol, ref, with_cn=.false., with_q=.true.)

end subroutine test_gw_mb02


subroutine test_gw_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(7, 16) = reshape([&
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 3.9505293783468E-09_wp, &
      & 1.0097520141925E-02_wp, 9.3220033286238E-01_wp, 2.5639784869775E-02_wp, &
      & 0.0000000000000E+00_wp, 1.4797736763222E-11_wp, 1.9052338447722E-04_wp, &
      & 1.1116142263240E+00_wp, 2.9351828841454E-03_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 5.3735049395955E-03_wp, &
      & 8.8839283385315E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 9.7753602561687E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 9.8036040945136E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.2907115615400E-09_wp, 9.9812405004339E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 3.9707631235001E-14_wp, &
      & 3.4070545718944E-05_wp, 1.0390095281076E+00_wp, 1.4775740267284E-04_wp, &
      & 0.0000000000000E+00_wp, 1.4709668915570E-03_wp, 8.7779710280724E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 4.9950695043173E-03_wp, &
      & 8.2586382594523E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 4.7331810125278E-03_wp, 7.8936053108503E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 6.1261156722753E-04_wp, 1.1351788897547E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.4636944706824E-06_wp, 1.0081631982011E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.5989800945240E-05_wp, 1.0653935533649E+00_wp, 7.4689286200916E-05_wp, &
      & 0.0000000000000E+00_wp, 9.7209060353668E-03_wp, 5.7638167950617E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 5.8433813639092E-10_wp, &
      & 9.6868781448195E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.8479075138024E-04_wp, 9.0297478986230E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp], shape(ref))

   call get_structure(mol, "MB16-43", "03")
   call test_gw_gen(error, mol, ref, with_cn=.true., with_q=.true.)

end subroutine test_gw_mb03


subroutine test_dgw_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_dgw_gen(error, mol, with_cn=.true., with_q=.false.)

end subroutine test_dgw_mb04


subroutine test_dgw_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_dgw_gen(error, mol, with_cn=.false., with_q=.true.)

end subroutine test_dgw_mb05


subroutine test_dgw_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "06")
   call test_dgw_gen(error, mol, with_cn=.true., with_q=.true.)

end subroutine test_dgw_mb06


subroutine test_gw_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: qat(16) = [&
      &-1.57324183192355E-1_wp, 1.65228028672395E-1_wp, 3.22320366812437E-1_wp, &
      & 3.63579860576008E-2_wp, 4.85677294229959E-2_wp,-3.59193331069774E-1_wp, &
      &-1.93844259127416E-1_wp,-3.86492630088218E-1_wp, 3.10104713332147E-1_wp, &
      & 8.34803229863654E-2_wp,-3.62667644019899E-1_wp, 3.64142434058147E-1_wp, &
      & 3.34644499696670E-1_wp,-4.69889877462762E-1_wp,-1.89224201365947E-1_wp, &
      & 4.53790045287620E-1_wp]
   real(wp), parameter :: ref(7, 16) = reshape([&
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 4.0696276249573E-09_wp, &
      & 4.1017933802195E-02_wp, 1.0950334442720E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 8.2493287561538E-08_wp, &
      & 6.2778573148766E-02_wp, 7.1426028532379E-01_wp, 9.0944312586520E-07_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.7678546474024E-03_wp, &
      & 4.5688574836408E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 4.8156277971682E-03_wp, 5.0705021767586E-01_wp, 2.6271362332948E-03_wp, &
      & 3.8726296026747E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 6.2942218738658E-13_wp, 1.0411557228148E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 4.9030037605503E-05_wp, &
      & 3.4343632271501E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.0798558970910E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.3113431982261E-02_wp, 3.8190728854670E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 5.2050837196991E-08_wp, &
      & 4.0261222021330E-02_wp, 9.1473844350972E-01_wp, 6.4613328643664E-04_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.2249494480500E-09_wp, 1.0367677493893E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 4.1505840413994E-02_wp, 3.4418284274292E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.5562886172976E-03_wp, &
      & 4.2019468075114E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.6236623287431E-03_wp, 4.4567473574759E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.9601157716485E-05_wp, 5.4850578977952E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 8.5241970390889E-09_wp, 1.5557325390935E-02_wp, &
      & 9.9022169633399E-01_wp, 3.7155741362071E-02_wp, 0.0000000000000E+00_wp, &
      & 7.0489201281527E-06_wp, 3.5717430460362E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp], shape(ref))

   type(structure_type) :: mol
   type(d4_model) :: d4
   real(wp), allocatable :: cn(:), gwvec(:, :)
   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :)

   call get_structure(mol, "MB16-43", "06")
   call new_d4_model(d4, mol, ref=d4_ref%gfn2)

   allocate(cn(mol%nat), gwvec(maxval(d4%ref), mol%nat))
   cn(:) = 0.0_wp

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call get_coordination_number(mol, lattr, cutoff, d4%rcov, d4%en, cn)

   call d4%weight_references(mol, cn, qat, gwvec)

   if (any(abs(gwvec - ref) > thr2)) then
      call test_failed(error, "Gaussian weights do not match")
      where(abs(gwvec) < thr) gwvec = 0.0_wp
      print'(3(es20.13,"_wp,"), " &")', gwvec
   end if

end subroutine test_gw_mb07


subroutine test_dgw_mb08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: qat(16) = [&
      &-2.05667914412001E-1_wp,-3.99554326663093E-1_wp, 3.29243862111419E-1_wp, &
      &-3.11738256025803E-1_wp, 3.58849862618190E-2_wp, 3.21889825709581E-1_wp, &
      & 4.14746199807777E-2_wp, 2.95730046338104E-2_wp,-5.06348926564523E-1_wp, &
      & 3.43067357217139E-1_wp, 6.88375720293390E-1_wp, 7.03350634832100E-2_wp, &
      &-9.62426987152087E-2_wp,-1.32210876939567E-1_wp, 9.79003738112971E-2_wp, &
      &-3.05981814182260E-1_wp]

   type(structure_type) :: mol
   integer :: iat, mref
   type(d4_model) :: d4
   real(wp), allocatable :: cn(:), q(:), gwvec(:, :), gwdcn(:, :), gwdq(:, :)
   real(wp), allocatable :: gwr(:, :), gwl(:, :), numdcn(:, :), numdq(:, :)
   real(wp), parameter :: cutoff = 30.0_wp, lattr(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp

   call get_structure(mol, "MB16-43", "08")
   call new_d4_model(d4, mol, ref=d4_ref%gfn2)

   mref = maxval(d4%ref)
   allocate(cn(mol%nat), q(mol%nat), gwvec(mref, mol%nat), &
      & gwdcn(mref, mol%nat), gwdq(mref, mol%nat), &
      & gwr(mref, mol%nat), gwl(mref, mol%nat), &
      & numdcn(mref, mol%nat), numdq(mref, mol%nat))
   cn(:) = 0.0_wp
   q(:) = qat

   call get_coordination_number(mol, lattr, cutoff, d4%rcov, d4%en, cn)

   do iat = 1, mol%nat
      cn(iat) = cn(iat) + step
      call d4%weight_references(mol, cn, q, gwr)
      cn(iat) = cn(iat) - 2*step
      call d4%weight_references(mol, cn, q, gwl)
      cn(iat) = cn(iat) + step
      gwdcn(:, :) = 0.5_wp*(gwr - gwl)/step
      numdcn(:, iat) = gwdcn(:, iat)
      gwdcn(:, iat) = 0.0_wp
      if (any(abs(gwdcn) > thr)) then
         call test_failed(error, "Unexpected non-zero gradient element found")
         exit
      end if
   end do
   if (allocated(error)) return

   do iat = 1, mol%nat
      q(iat) = q(iat) + step
      call d4%weight_references(mol, cn, q, gwr)
      q(iat) = q(iat) - 2*step
      call d4%weight_references(mol, cn, q, gwl)
      q(iat) = q(iat) + step
      gwdq(:, :) = 0.5_wp*(gwr - gwl)/step
      numdq(:, iat) = gwdq(:, iat)
      gwdq(:, iat) = 0.0_wp
      if (any(abs(gwdq) > thr)) then
         call test_failed(error, "Unexpected non-zero gradient element found")
         exit
      end if
   end do
   if (allocated(error)) return

   call d4%weight_references(mol, cn, q, gwvec, gwdcn, gwdq)

   if (any(abs(gwdcn - numdcn) > thr2)) then
      call test_failed(error, "Gaussian weights derivatives do not match")
      print'(3es21.14)', gwdcn
      print'("---")'
      print'(3es21.14)', numdcn
      print'("---")'
      print'(3es21.14)', gwdcn - numdcn
   end if

   if (any(abs(gwdq - numdq) > thr2)) then
      call test_failed(error, "Gaussian weights derivatives do not match")
      print'(3es21.14)', gwdq
      print'("---")'
      print'(3es21.14)', numdq
      print'("---")'
      print'(3es21.14)', gwdq - numdq
   end if

end subroutine test_dgw_mb08


end module test_model
