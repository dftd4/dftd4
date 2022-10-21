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
      & 1.2088234202205E-11_wp, 1.0375924114812E+00_wp, 2.4445602531076E-09_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 5.5306460372092E-03_wp, &
      & 9.9446935396279E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 3.1404764970412E-06_wp, 6.3956931077352E-02_wp, &
      & 4.7846317290481E-01_wp, 2.9468368515225E-01_wp, 0.0000000000000E+00_wp, &
      & 7.4310920077799E-05_wp, 9.9992568907992E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 4.4395780891660E-02_wp, &
      & 8.5383545436310E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.1950286131613E-02_wp, 9.8804971386839E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 5.6992908407601E-05_wp, 9.9994300709159E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 5.6167960886507E-07_wp, &
      & 2.9078967459651E-02_wp, 6.1859560324107E-01_wp, 1.5104287585421E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 4.9771287942905E-14_wp, &
      & 1.3962072360858E-04_wp, 6.8484995485525E-01_wp, 1.6516463249478E-12_wp, &
      & 1.9880162618842E-02_wp, 9.8011983738116E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.1966389405645E-02_wp, &
      & 9.8803361059436E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.2874045548253E-05_wp, 9.7726740456727E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 4.1491138092838E-06_wp, 8.6768202830170E-01_wp, 0.0000000000000E+00_wp, &
      & 4.3169267440676E-12_wp, 3.6465051654232E-04_wp, 8.8903503100593E-01_wp, &
      & 1.2720379452861E-05_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.4066048552144E-05_wp, 6.8494819902320E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.9574616732692E-10_wp, &
      & 1.0064486114128E+00_wp, 0.0000000000000E+00_wp], shape(ref))

   call get_structure(mol, "MB16-43", "01")
   call test_gw_gen(error, mol, ref, with_cn=.true., with_q=.false.)

end subroutine test_gw_mb01


subroutine test_gw_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(5, 16) = reshape([&
      & 8.1529926850468E-01_wp, 2.2414425700956E-03_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.0229633690838E+00_wp, &
      & 2.0291440612814E-03_wp, 4.4863013350169E-10_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.1578430925721E+00_wp, 2.0553787158892E-03_wp, &
      & 3.1900095597967E-10_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.3687709975145E+00_wp, 8.2192921625604E-03_wp, 5.9806186105627E-08_wp, &
      & 1.4414431691639E-03_wp, 0.0000000000000E+00_wp, 9.1646421111918E-01_wp, &
      & 5.3618567538002E-03_wp, 1.7524886555881E-08_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 5.5812105498128E-01_wp, 1.5344013422165E-03_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 7.5659186904499E-01_wp, 2.0800426162236E-03_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 7.7216222381926E-01_wp, &
      & 2.1228490523027E-03_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 9.9399938189677E-01_wp, 2.8186797641147E-03_wp, &
      & 7.1776130956295E-10_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 5.5220460179834E-01_wp, 1.5181356707747E-03_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.1425733041093E+00_wp, &
      & 2.0284312531719E-03_wp, 3.1484305175579E-10_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 8.1572049134917E-01_wp, 6.4067287425100E-03_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.2272138342948E+00_wp, 1.3612868628584E-02_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 8.3285816475829E-01_wp, &
      & 2.2897159576321E-03_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 7.4607907232618E-01_wp, 2.0511405541139E-03_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.0523038058719E+00_wp, 2.0872155477415E-03_wp, 4.6145584504737E-10_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp], shape(ref))

   call get_structure(mol, "MB16-43", "02")
   call test_gw_gen(error, mol, ref, with_cn=.false., with_q=.true.)

end subroutine test_gw_mb02


subroutine test_gw_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(7, 16) = reshape([&
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 3.9504240794198E-09_wp, &
      & 1.0097386820563E-02_wp, 9.3220073517801E-01_wp, 2.5639475147259E-02_wp, &
      & 0.0000000000000E+00_wp, 1.4797717135069E-11_wp, 1.9052325554662E-04_wp, &
      & 1.1116141991349E+00_wp, 2.9351813145404E-03_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 5.3735027174640E-03_wp, &
      & 8.8839284731689E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 9.7753605500113E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 9.8036041832883E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.2906805265616E-09_wp, 9.9812405470369E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 3.9706733315821E-14_wp, &
      & 3.4070161345981E-05_wp, 1.0390095052930E+00_wp, 1.4775589359887E-04_wp, &
      & 0.0000000000000E+00_wp, 1.4709587472355E-03_wp, 8.7779708361772E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 4.9950676848882E-03_wp, &
      & 8.2586385264927E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 4.7331797500540E-03_wp, 7.8936056984439E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 6.1261133783377E-04_wp, 1.1351788724746E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.4636823892691E-06_wp, 1.0081631910738E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.5989624842058E-05_wp, 1.0653935240443E+00_wp, 7.4688541496264E-05_wp, &
      & 0.0000000000000E+00_wp, 9.7208858575852E-03_wp, 5.7638169663657E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 5.8433000707751E-10_wp, &
      & 9.6868781982789E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.8478831017041E-04_wp, 9.0297474706410E-01_wp, 0.0000000000000E+00_wp, &
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
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 4.0693983485343E-09_wp, &
      & 4.1016483181479E-02_wp, 1.0950348404655E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 8.2493063323491E-08_wp, &
      & 6.2778502874922E-02_wp, 7.1426034737421E-01_wp, 9.0944083534497E-07_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.7678536584132E-03_wp, &
      & 4.5688574935307E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 4.8156277656763E-03_wp, 5.0705021724268E-01_wp, 2.6271362451520E-03_wp, &
      & 3.8726296078651E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 6.2939729329860E-13_wp, 1.0411557228148E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 4.9029610641780E-05_wp, &
      & 3.4343632275771E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.0798558970910E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.3113423827832E-02_wp, 3.8190728936215E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 5.2050674938049E-08_wp, &
      & 4.0261160521426E-02_wp, 9.1473850489949E-01_wp, 6.4613435497592E-04_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.2249110023524E-09_wp, 1.0367677493894E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 4.1505834538896E-02_wp, 3.4418284333043E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.5562876550469E-03_wp, &
      & 4.2019468171339E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.6236585470333E-03_wp, 4.4567473952930E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.9600876599811E-05_wp, 5.4850578980763E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 8.5241700290816E-09_wp, 1.5557301075470E-02_wp, &
      & 9.9022177224599E-01_wp, 3.7155688096309E-02_wp, 0.0000000000000E+00_wp, &
      & 7.0488497407158E-06_wp, 3.5717430467401E-01_wp, 0.0000000000000E+00_wp, &
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
