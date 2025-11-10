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
   use dftd4_cutoff, only : get_lattice_points
   use dftd4_data, only : get_covalent_rad
   use dftd4_model, only : dispersion_model, new_dispersion_model, d4_qmod
   use dftd4_model_d4, only : d4_model, new_d4_model
   use dftd4_model_d4s, only : d4s_model, new_d4s_model
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, &
      & test_failed
   use mctc_io_structure, only : new, structure_type
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use mstore, only : get_structure
   use multicharge, only : get_charges
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
      & new_unittest("gw-D4-mb01", test_gw_d4_mb01), &
      & new_unittest("gw-D4S-mb01", test_gw_d4s_mb01), &
      & new_unittest("gw-D4-mb02", test_gw_d4_mb02), &
      & new_unittest("gw-D4-EEQBC-mb02", test_gw_d4_eeqbc_mb02), &
      & new_unittest("gw-D4-mb03", test_gw_d4_mb03), &
      & new_unittest("dgw_D4-mb04", test_dgw_d4_mb04), &
      & new_unittest("dgw_D4S-mb04", test_dgw_d4s_mb04), &
      & new_unittest("dgw_D4-mb05", test_dgw_d4_mb05), &
      & new_unittest("dgw_D4S-mb05", test_dgw_d4s_mb05), &
      & new_unittest("dgw_D4S-EEQBC-mb05", test_dgw_d4s_eeqbc_mb05), &
      & new_unittest("dgw_D4-mb06", test_dgw_d4_mb06), &
      & new_unittest("dgw_D4S-mb06", test_dgw_d4s_mb06), &
      & new_unittest("gw-D4-gfn2", test_gw_d4_mb07), &
      & new_unittest("dgw-D4-gfn2", test_dgw_d4_mb08), &
      & new_unittest("dgw-D4S-gfn2", test_dgw_d4s_mb08), &
      & new_unittest("pol-D4-mb09", test_pol_d4_mb09), &
      & new_unittest("pol-D4s-mb09", test_pol_d4s_mb09), &
      & new_unittest("dpol-D4-mb10", test_dpol_d4_mb10), &
      & new_unittest("dpol-D4S-mb10", test_dpol_d4s_mb10), &
      & new_unittest("model-D4-error", test_d4_model_error, should_fail=.true.), &
      & new_unittest("model-D4S-error", test_d4s_model_error, should_fail=.true.), &
      & new_unittest("model-wrapper", test_model_wrapper), &
      & new_unittest("model-wrapper-fail", test_model_wrapper_fail, should_fail=.true.) &
      & ]

end subroutine collect_model


subroutine test_gw_gen(error, mol, d4, ref, with_cn, with_q, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Reference Gaussian weights
   real(wp), intent(in) :: ref(:, :, :)

   !> Calculate coordination number
   logical, intent(in) :: with_cn

   !> Calculate atomic charges
   logical, intent(in) :: with_q

   !> Atomic charges
   real(wp), optional, intent(in) :: qat(:)

   real(wp), allocatable :: cn(:), q(:), gwvec(:, :, :)
   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :)
   class(ncoord_type), allocatable :: ncoord

   allocate(cn(mol%nat), q(mol%nat), gwvec(maxval(d4%ref), mol%nat, d4%ncoup))
   cn(:) = 0.0_wp
   q(:) = 0.0_wp

   if (with_cn) then
      call new_ncoord(ncoord, mol, cn_count%dftd4, &
         & cutoff=cutoff, rcov=d4%rcov, en=d4%en, error=error)
      if (allocated(error)) return
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call ncoord%get_coordination_number(mol, lattr, cn)
   end if
   if (with_q) then
      if(present(qat)) then
         q(:) = qat
      else
         call get_charges(d4%mchrg, mol, error, q)
         if (allocated(error)) return
      end if
   end if

   call d4%weight_references(mol, cn, q, gwvec)
   if (any(abs(gwvec - ref) > thr2)) then
      call test_failed(error, "Gaussian weights do not match")
      where(abs(gwvec) < thr) gwvec = 0.0_wp
      print'(3(es20.13,"_wp,"), " &")', gwvec
   end if

end subroutine test_gw_gen


subroutine test_dgw_gen(error, mol, d4, with_cn, with_q, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Calculate coordination number
   logical, intent(in) :: with_cn

   !> Calculate atomic charges
   logical, intent(in) :: with_q

   !> Atomic charges
   real(wp), optional, intent(in) :: qat(:)

   integer :: iat, mref, ncoup
   real(wp), allocatable :: cn(:), q(:), gwvec(:, :, :), gwdcn(:, :, :), gwdq(:, :, :)
   real(wp), allocatable :: gwr(:, :, :), gwl(:, :, :), numdcn(:, :, :), numdq(:, :, :)
   real(wp), parameter :: cutoff = 30.0_wp, lattr(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   class(ncoord_type), allocatable :: ncoord

   mref = maxval(d4%ref)
   ncoup = d4%ncoup
   allocate(cn(mol%nat), q(mol%nat), gwvec(mref, mol%nat, ncoup), &
      & gwdcn(mref, mol%nat, ncoup), gwdq(mref, mol%nat, ncoup), &
      & gwr(mref, mol%nat, ncoup), gwl(mref, mol%nat, ncoup), &
      & numdcn(mref, mol%nat, ncoup), numdq(mref, mol%nat, ncoup))
   cn(:) = 0.0_wp
   q(:) = 0.0_wp

   if (with_cn) then
      call new_ncoord(ncoord, mol, cn_count%dftd4, &
         & cutoff=cutoff, rcov=d4%rcov, en=d4%en, error=error)
      if (allocated(error)) return
      call ncoord%get_coordination_number(mol, lattr, cn)
   end if
   if (with_q) then
      if(present(qat)) then
         q(:) = qat
      else
         call get_charges(d4%mchrg, mol, error, q)
         if (allocated(error)) return
      end if
   end if

   if (with_cn) then
      do iat = 1, mol%nat
         cn(iat) = cn(iat) + step
         call d4%weight_references(mol, cn, q, gwr)
         cn(iat) = cn(iat) - 2*step
         call d4%weight_references(mol, cn, q, gwl)
         cn(iat) = cn(iat) + step
         gwdcn(:, :, :) = 0.5_wp*(gwr - gwl)/step
         numdcn(:, iat, :) = gwdcn(:, iat, :)
         gwdcn(:, iat, :) = 0.0_wp 
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
         gwdq(:, :, :) = 0.5_wp*(gwr - gwl)/step
         numdq(:, iat, :) = gwdq(:, iat, :)
         gwdq(:, iat, :) = 0.0_wp
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


subroutine test_pol_gen(error, mol, d4, ref, with_cn, with_q, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Reference polarizabilities
   real(wp), intent(in) :: ref(:)

   !> Calculate coordination number
   logical, intent(in) :: with_cn

   !> Calculate atomic charges
   logical, intent(in) :: with_q

   !> Atomic charges
   real(wp), optional, intent(in) :: qat(:)

   real(wp), allocatable :: cn(:), q(:), alpha(:), gwvec(:, :, :)
   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :)
   class(ncoord_type), allocatable :: ncoord

   allocate(cn(mol%nat), q(mol%nat), alpha(mol%nat), &
      & gwvec(maxval(d4%ref), mol%nat, d4%ncoup))
   cn(:) = 0.0_wp
   q(:) = 0.0_wp

   if (with_cn) then
      call new_ncoord(ncoord, mol, cn_count%dftd4, &
         & cutoff=cutoff, rcov=d4%rcov, en=d4%en, error=error)
      if (allocated(error)) return
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call ncoord%get_coordination_number(mol, lattr, cn)
   end if

   if (with_q) then
      if(present(qat)) then
         q(:) = qat
      else
         call get_charges(d4%mchrg, mol, error, q)
         if (allocated(error)) return
      end if
   end if

   call d4%weight_references(mol, cn, q, gwvec)

   call d4%get_polarizabilities(mol, gwvec, alpha=alpha)

   if (any(abs(alpha - ref) > thr2)) then
      call test_failed(error, "Polarizabilities do not match")
      where(abs(alpha) < thr) alpha = 0.0_wp
      print'(3(es20.13,"_wp,"), " &")', alpha
   end if

end subroutine test_pol_gen


subroutine test_dpol_gen(error, mol, d4, with_cn, with_q, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Calculate coordination number
   logical, intent(in) :: with_cn

   !> Calculate atomic charges
   logical, intent(in) :: with_q

   !> Atomic charges
   real(wp), optional, intent(in) :: qat(:)

   integer :: iat, mref, ncoup
   real(wp), allocatable :: cn(:), q(:), gwvec(:, :, :), gwdcn(:, :, :), gwdq(:, :, :)
   real(wp), allocatable :: alpha(:), alphadcn(:), alphadq(:)
   real(wp), allocatable :: alphar(:), alphal(:), numdcn(:), numdq(:)
   real(wp), parameter :: cutoff = 30.0_wp, lattr(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   class(ncoord_type), allocatable :: ncoord

   mref = maxval(d4%ref)
   ncoup = d4%ncoup
   allocate(cn(mol%nat), q(mol%nat), gwvec(mref, mol%nat, ncoup), &
      & gwdcn(mref, mol%nat, ncoup), gwdq(mref, mol%nat, ncoup), &
      & alpha(mol%nat), alphadcn(mol%nat), alphadq(mol%nat), &
      & alphar(mol%nat), alphal(mol%nat), numdcn(mol%nat), numdq(mol%nat))
   cn(:) = 0.0_wp
   q(:) = 0.0_wp

   if (with_cn) then
      call new_ncoord(ncoord, mol, cn_count%dftd4, &
         & cutoff=cutoff, rcov=d4%rcov, en=d4%en, error=error)
      if (allocated(error)) return
      call ncoord%get_coordination_number(mol, lattr, cn)
   end if

   if (with_q) then
      if(present(qat)) then
         q(:) = qat
      else
         call get_charges(d4%mchrg, mol, error, q)
         if (allocated(error)) return
      end if
   end if

   if (with_cn) then
      do iat = 1, mol%nat
         cn(iat) = cn(iat) + step
         call d4%weight_references(mol, cn, q, gwvec)
         call d4%get_polarizabilities(mol, gwvec, alpha=alphar)
         cn(iat) = cn(iat) - 2*step
         call d4%weight_references(mol, cn, q, gwvec)
         call d4%get_polarizabilities(mol, gwvec, alpha=alphal)
         cn(iat) = cn(iat) + step
         numdcn(iat) = 0.5_wp*(alphar(iat) - alphal(iat))/step
      end do
   end if

   if (with_q) then
      do iat = 1, mol%nat
         q(iat) = q(iat) + step
         call d4%weight_references(mol, cn, q, gwvec)
         call d4%get_polarizabilities(mol, gwvec, alpha=alphar)
         q(iat) = q(iat) - 2*step
         call d4%weight_references(mol, cn, q, gwvec)
         call d4%get_polarizabilities(mol, gwvec, alpha=alphal)
         q(iat) = q(iat) + step
         numdq(iat) = 0.5_wp*(alphar(iat) - alphal(iat))/step
      end do
   end if

   call d4%weight_references(mol, cn, q, gwvec, gwdcn, gwdq)
   call d4%get_polarizabilities(mol, gwvec, gwdcn=gwdcn, gwdq=gwdq, &
      & alpha=alpha, dadcn=alphadcn, dadq=alphadq)


   if (with_cn .and. any(abs(alphadcn - numdcn) > thr2)) then
     call test_failed(error, "Gaussian weights derivatives do not match")
     print'(3es21.14)', alphadcn
     print'("---")'
     print'(3es21.14)', numdcn
     print'("---")'
     print'(3es21.14)', alphadcn - numdcn
   end if

   if (with_q .and. any(abs(alphadq - numdq) > thr2)) then
      call test_failed(error, "Gaussian weights derivatives do not match")
      print'(3es21.14)', alphadq
      print'("---")'
      print'(3es21.14)', numdq
      print'("---")'
      print'(3es21.14)', alphadq - numdq
   end if

end subroutine test_dpol_gen


subroutine test_gw_d4_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   real(wp), parameter :: ref(5, 16, 1) = reshape([&
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
      & 1.0064486114128E+00_wp, 0.0000000000000E+00_wp], &
      & [5, 16, 1])

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_gw_gen(error, mol, d4, ref, with_cn=.true., with_q=.false.)

end subroutine test_gw_d4_mb01

subroutine test_gw_d4s_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   real(wp), parameter :: ref(5, 16, 16) = reshape([&
      & 2.3131592204844E-01_wp, 5.4031280728148E-01_wp, 2.6137589903151E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 5.1809010783127E-01_wp, &
      & 4.8190989216873E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 8.4636699292911E-02_wp, 2.8916075663289E-01_wp, &
      & 1.2624087141442E-01_wp, 4.0965002951526E-01_wp, 0.0000000000000E+00_wp, &
      & 3.3462728435151E-01_wp, 6.6537271564849E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 6.1876633063175E-01_wp, &
      & 3.4063351416268E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 5.5278384311783E-01_wp, 4.4721615688217E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 3.2486444154312E-01_wp, 6.7513555845688E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 7.5533304743754E-02_wp, &
      & 2.8014049695601E-01_wp, 1.4158516457382E-01_wp, 4.0913000753423E-01_wp, &
      & 0.0000000000000E+00_wp, 7.5382669995544E-03_wp, 5.6004612639252E-02_wp, &
      & 2.5449767768026E-01_wp, 3.7693011028422E-01_wp, 8.5938011138013E-02_wp, &
      & 5.7538171652020E-01_wp, 4.2461828347980E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 5.5284405139333E-01_wp, &
      & 4.4715594860667E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 4.9408119795623E-01_wp, 4.9442431983631E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.0184875449037E-03_wp, 1.1949848707631E-02_wp, 7.6197788846927E-02_wp, &
      & 2.6107757645099E-01_wp, 5.2942068197745E-01_wp, 1.6601987730686E-02_wp, &
      & 8.4375771909016E-02_wp, 2.4112060689684E-01_wp, 3.7174716102808E-01_wp, &
      & 1.8680908647673E-01_wp, 4.6453010724907E-03_wp, 4.0443789086525E-02_wp, &
      & 2.3923007812319E-01_wp, 4.1892702798575E-01_wp, 6.3157269210764E-02_wp, &
      & 1.3398605819820E-02_wp, 7.0902968576691E-02_wp, 2.5998769707888E-01_wp, &
      & 6.6186984945544E-01_wp, 0.0000000000000E+00_wp, 5.5310440465207E-02_wp, &
      & 8.7883262647393E-01_wp, 1.0298872576581E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 4.1485451554553E-02_wp, 9.5851454844545E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 6.7032432969941E-05_wp, 1.0329252259896E-01_wp, 3.9006942107954E-01_wp, &
      & 3.6487976107888E-01_wp, 0.0000000000000E+00_wp, 2.9607006073691E-03_wp, &
      & 9.9703929939263E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 7.5207498914890E-02_wp, 8.2630508485236E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 6.6472988033016E-02_wp, 9.3352701196698E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.5170748238388E-03_wp, &
      & 9.9748292517616E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.1199100233663E-05_wp, 6.1735949348506E-02_wp, &
      & 5.2547467287298E-01_wp, 2.3534849766078E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.1545009180054E-08_wp, 4.8156324470167E-03_wp, &
      & 6.8119098700326E-01_wp, 9.8640706987545E-08_wp, 9.0609306620312E-02_wp, &
      & 9.0939069337969E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 6.6527647187693E-02_wp, 9.3347235281231E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 5.0319522106180E-03_wp, 9.7236235992780E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.0833572554318E-11_wp, 3.0735372456843E-06_wp, 8.7114485461151E-03_wp, &
      & 8.5918441677693E-01_wp, 1.0676321498807E-09_wp, 4.8559648922814E-05_wp, &
      & 4.5302763715688E-02_wp, 8.3289153450766E-01_wp, 1.2364799599275E-02_wp, &
      & 0.0000000000000E+00_wp, 4.8553772734271E-10_wp, 1.2356257410943E-03_wp, &
      & 6.8399234763729E-01_wp, 5.1995420254894E-09_wp, 1.2027896985957E-11_wp, &
      & 5.7304401546169E-07_wp, 2.4468959606775E-03_wp, 1.0040005130011E+00_wp, &
      & 0.0000000000000E+00_wp, 8.8802202146486E-02_wp, 8.0146979806795E-01_wp, &
      & 1.4628116593891E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.2198992510853E-02_wp, 9.7780100748915E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 3.1404764970412E-06_wp, &
      & 6.3956931077352E-02_wp, 4.7846317290481E-01_wp, 2.9468368515225E-01_wp, &
      & 0.0000000000000E+00_wp, 9.4397261326925E-04_wp, 9.9905602738673E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 4.4395780891660E-02_wp, 8.5383545436310E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 3.9032141252572E-02_wp, &
      & 9.6096785874743E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 7.7726046635812E-04_wp, 9.9922273953364E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 5.6167960886507E-07_wp, 2.9078967459651E-02_wp, 6.1859560324107E-01_wp, &
      & 1.5104287585421E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 8.6218571603764E-13_wp, 3.1641253888564E-04_wp, 6.8471161857811E-01_wp, &
      & 2.0919541563050E-11_wp, 5.6611833802765E-02_wp, 9.4338816619724E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 3.9070615176354E-02_wp, 9.6092938482365E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.2874045548253E-05_wp, &
      & 9.7726740456727E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.3913444097198E-09_wp, 6.3616927572147E-04_wp, 8.6706543647995E-01_wp, &
      & 0.0000000000000E+00_wp, 1.9173135281474E-07_wp, 8.9543641410655E-03_wp, &
      & 8.7950766885962E-01_wp, 1.2071886445160E-03_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 3.9522919456175E-05_wp, 6.8492827950137E-01_wp, &
      & 2.3393071766261E-13_wp, 0.0000000000000E+00_wp, 1.3854678522329E-13_wp, &
      & 4.1660837709814E-06_wp, 1.0064444444547E+00_wp, 0.0000000000000E+00_wp, &
      & 5.5310440465207E-02_wp, 8.7883262647393E-01_wp, 1.0298872576581E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 4.1485451554553E-02_wp, &
      & 9.5851454844545E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 6.7032432969941E-05_wp, 1.0329252259896E-01_wp, &
      & 3.9006942107954E-01_wp, 3.6487976107888E-01_wp, 0.0000000000000E+00_wp, &
      & 2.9607006073691E-03_wp, 9.9703929939263E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 7.5207498914890E-02_wp, &
      & 8.2630508485236E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 6.6472988033016E-02_wp, 9.3352701196698E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.5170748238388E-03_wp, 9.9748292517616E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.1199100233663E-05_wp, &
      & 6.1735949348506E-02_wp, 5.2547467287298E-01_wp, 2.3534849766078E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.1545009180054E-08_wp, &
      & 4.8156324470167E-03_wp, 6.8119098700326E-01_wp, 9.8640706987545E-08_wp, &
      & 9.0609306620312E-02_wp, 9.0939069337969E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 6.6527647187693E-02_wp, &
      & 9.3347235281231E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 5.0319522106180E-03_wp, 9.7236235992780E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.0833572554318E-11_wp, 3.0735372456843E-06_wp, &
      & 8.7114485461151E-03_wp, 8.5918441677693E-01_wp, 1.0676321498807E-09_wp, &
      & 4.8559648922814E-05_wp, 4.5302763715688E-02_wp, 8.3289153450766E-01_wp, &
      & 1.2364799599275E-02_wp, 0.0000000000000E+00_wp, 4.8553772734271E-10_wp, &
      & 1.2356257410943E-03_wp, 6.8399234763729E-01_wp, 5.1995420254894E-09_wp, &
      & 1.2027896985957E-11_wp, 5.7304401546169E-07_wp, 2.4468959606775E-03_wp, &
      & 1.0040005130011E+00_wp, 0.0000000000000E+00_wp, 6.0951918468791E-02_wp, &
      & 8.6529295743687E-01_wp, 1.1079756841654E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.3336298784774E-02_wp, 9.8666370121523E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 3.1404764970412E-06_wp, 6.3956931077352E-02_wp, 4.7846317290481E-01_wp, &
      & 2.9468368515225E-01_wp, 0.0000000000000E+00_wp, 3.7174949263315E-04_wp, &
      & 9.9962825050737E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 4.4395780891660E-02_wp, 8.5383545436310E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.5290254001583E-02_wp, 9.7470974599842E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.9822583135455E-04_wp, &
      & 9.9970177416865E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 5.6167960886507E-07_wp, 2.9078967459651E-02_wp, &
      & 6.1859560324107E-01_wp, 1.5104287585421E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 4.9771287942905E-14_wp, 1.3962072360858E-04_wp, &
      & 6.8484995485525E-01_wp, 1.6516463249478E-12_wp, 3.8575470395886E-02_wp, &
      & 9.6142452960411E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.5318540756121E-02_wp, 9.7468145924388E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.2874045548253E-05_wp, 9.7726740456727E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.6273317494726E-10_wp, 2.3773118630915E-04_wp, &
      & 8.6745414891661E-01_wp, 0.0000000000000E+00_wp, 2.3703664758231E-08_wp, &
      & 4.8024934491719E-03_wp, 8.8424901114203E-01_wp, 4.9697715237582E-04_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.4066048552144E-05_wp, &
      & 6.8494819902320E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.3920190449891E-06_wp, 1.0064472192320E+00_wp, &
      & 0.0000000000000E+00_wp, 5.5310440465207E-02_wp, 8.7883262647393E-01_wp, &
      & 1.0298872576581E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 4.1485451554553E-02_wp, 9.5851454844545E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 6.7032432969941E-05_wp, &
      & 1.0329252259896E-01_wp, 3.9006942107954E-01_wp, 3.6487976107888E-01_wp, &
      & 0.0000000000000E+00_wp, 2.9607006073691E-03_wp, 9.9703929939263E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 7.5207498914890E-02_wp, 8.2630508485236E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 6.6472988033016E-02_wp, &
      & 9.3352701196698E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.5170748238388E-03_wp, 9.9748292517616E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.1199100233663E-05_wp, 6.1735949348506E-02_wp, 5.2547467287298E-01_wp, &
      & 2.3534849766078E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.1545009180054E-08_wp, 4.8156324470167E-03_wp, 6.8119098700326E-01_wp, &
      & 9.8640706987545E-08_wp, 9.0609306620312E-02_wp, 9.0939069337969E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 6.6527647187693E-02_wp, 9.3347235281231E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 5.0319522106180E-03_wp, &
      & 9.7236235992780E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.0833572554318E-11_wp, &
      & 3.0735372456843E-06_wp, 8.7114485461151E-03_wp, 8.5918441677693E-01_wp, &
      & 1.0676321498807E-09_wp, 4.8559648922814E-05_wp, 4.5302763715688E-02_wp, &
      & 8.3289153450766E-01_wp, 1.2364799599275E-02_wp, 0.0000000000000E+00_wp, &
      & 4.8553772734271E-10_wp, 1.2356257410943E-03_wp, 6.8399234763729E-01_wp, &
      & 5.1995420254894E-09_wp, 1.2027896985957E-11_wp, 5.7304401546169E-07_wp, &
      & 2.4468959606775E-03_wp, 1.0040005130011E+00_wp, 0.0000000000000E+00_wp, &
      & 5.5310440465207E-02_wp, 8.7883262647393E-01_wp, 1.0298872576581E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 4.1485451554553E-02_wp, &
      & 9.5851454844545E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 6.7032432969941E-05_wp, 1.0329252259896E-01_wp, &
      & 3.9006942107954E-01_wp, 3.6487976107888E-01_wp, 0.0000000000000E+00_wp, &
      & 2.9607006073691E-03_wp, 9.9703929939263E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 7.5207498914890E-02_wp, &
      & 8.2630508485236E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 6.6472988033016E-02_wp, 9.3352701196698E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.5170748238388E-03_wp, 9.9748292517616E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.1199100233663E-05_wp, &
      & 6.1735949348506E-02_wp, 5.2547467287298E-01_wp, 2.3534849766078E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.1545009180054E-08_wp, &
      & 4.8156324470167E-03_wp, 6.8119098700326E-01_wp, 9.8640706987545E-08_wp, &
      & 9.0609306620312E-02_wp, 9.0939069337969E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 6.6527647187693E-02_wp, &
      & 9.3347235281231E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 5.0319522106180E-03_wp, 9.7236235992780E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.0833572554318E-11_wp, 3.0735372456843E-06_wp, &
      & 8.7114485461151E-03_wp, 8.5918441677693E-01_wp, 1.0676321498807E-09_wp, &
      & 4.8559648922814E-05_wp, 4.5302763715688E-02_wp, 8.3289153450766E-01_wp, &
      & 1.2364799599275E-02_wp, 0.0000000000000E+00_wp, 4.8553772734271E-10_wp, &
      & 1.2356257410943E-03_wp, 6.8399234763729E-01_wp, 5.1995420254894E-09_wp, &
      & 1.2027896985957E-11_wp, 5.7304401546169E-07_wp, 2.4468959606775E-03_wp, &
      & 1.0040005130011E+00_wp, 0.0000000000000E+00_wp, 8.8802202146486E-02_wp, &
      & 8.0146979806795E-01_wp, 1.4628116593891E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.2198992510853E-02_wp, 9.7780100748915E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 3.1404764970412E-06_wp, 6.3956931077352E-02_wp, 4.7846317290481E-01_wp, &
      & 2.9468368515225E-01_wp, 0.0000000000000E+00_wp, 9.4397261326925E-04_wp, &
      & 9.9905602738673E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 4.4395780891660E-02_wp, 8.5383545436310E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 3.9032141252572E-02_wp, 9.6096785874743E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 7.7726046635812E-04_wp, &
      & 9.9922273953364E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 5.6167960886507E-07_wp, 2.9078967459651E-02_wp, &
      & 6.1859560324107E-01_wp, 1.5104287585421E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 8.6218571603764E-13_wp, 3.1641253888564E-04_wp, &
      & 6.8471161857811E-01_wp, 2.0919541563050E-11_wp, 5.6611833802765E-02_wp, &
      & 9.4338816619724E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 3.9070615176354E-02_wp, 9.6092938482365E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.2874045548253E-05_wp, 9.7726740456727E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.3913444097198E-09_wp, 6.3616927572147E-04_wp, &
      & 8.6706543647995E-01_wp, 0.0000000000000E+00_wp, 1.9173135281474E-07_wp, &
      & 8.9543641410655E-03_wp, 8.7950766885962E-01_wp, 1.2071886445160E-03_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 3.9522919456175E-05_wp, &
      & 6.8492827950137E-01_wp, 2.3393071766261E-13_wp, 0.0000000000000E+00_wp, &
      & 1.3854678522329E-13_wp, 4.1660837709814E-06_wp, 1.0064444444547E+00_wp, &
      & 0.0000000000000E+00_wp, 1.1619449180297E-01_wp, 7.4292303143376E-01_wp, &
      & 1.7688691875381E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 4.5511447319703E-02_wp, 9.5448855268030E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 9.2508392559746E-06_wp, &
      & 7.5824069699182E-02_wp, 4.4951921574350E-01_wp, 3.1887288951654E-01_wp, &
      & 0.0000000000000E+00_wp, 3.5066919289813E-03_wp, 9.9649330807102E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 4.4395780891660E-02_wp, 8.5383545436310E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 7.1925482699841E-02_wp, &
      & 9.2807451730016E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.9954258724411E-03_wp, 9.9700457412756E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.0238008000468E-06_wp, 3.8055054539114E-02_wp, 5.9074525520210E-01_wp, &
      & 1.7718982934717E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 3.4151558259365E-08_wp, 6.5639826011158E-03_wp, 6.7982280802516E-01_wp, &
      & 2.5898101368240E-07_wp, 9.7140202686718E-02_wp, 9.0285979731328E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 7.1982878187917E-02_wp, 9.2801712181208E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.1772574193346E-03_wp, &
      & 9.7612947599990E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 3.5689039339356E-13_wp, &
      & 3.8192206935584E-07_wp, 4.0712983354820E-03_wp, 8.6371381723556E-01_wp, &
      & 3.5979138804999E-11_wp, 9.7693789150509E-06_wp, 2.8536546460626E-02_wp, &
      & 8.5530890313171E-01_wp, 6.3404470879011E-03_wp, 0.0000000000000E+00_wp, &
      & 1.7396779713960E-09_wp, 1.8302935084816E-03_wp, 6.8352702270097E-01_wp, &
      & 1.6312609722908E-08_wp, 0.0000000000000E+00_wp, 7.9041698936705E-09_wp, &
      & 4.0780312005162E-04_wp, 1.0060406958243E+00_wp, 0.0000000000000E+00_wp, &
      & 5.5310440465207E-02_wp, 8.7883262647393E-01_wp, 1.0298872576581E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 4.1485451554553E-02_wp, &
      & 9.5851454844545E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 6.7032432969941E-05_wp, 1.0329252259896E-01_wp, &
      & 3.9006942107954E-01_wp, 3.6487976107888E-01_wp, 0.0000000000000E+00_wp, &
      & 2.9607006073691E-03_wp, 9.9703929939263E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 7.5207498914890E-02_wp, &
      & 8.2630508485236E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 6.6472988033016E-02_wp, 9.3352701196698E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.5170748238388E-03_wp, 9.9748292517616E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.1199100233663E-05_wp, &
      & 6.1735949348506E-02_wp, 5.2547467287298E-01_wp, 2.3534849766078E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.1545009180054E-08_wp, &
      & 4.8156324470167E-03_wp, 6.8119098700326E-01_wp, 9.8640706987545E-08_wp, &
      & 9.0609306620312E-02_wp, 9.0939069337969E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 6.6527647187693E-02_wp, &
      & 9.3347235281231E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 5.0319522106180E-03_wp, 9.7236235992780E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.0833572554318E-11_wp, 3.0735372456843E-06_wp, &
      & 8.7114485461151E-03_wp, 8.5918441677693E-01_wp, 1.0676321498807E-09_wp, &
      & 4.8559648922814E-05_wp, 4.5302763715688E-02_wp, 8.3289153450766E-01_wp, &
      & 1.2364799599275E-02_wp, 0.0000000000000E+00_wp, 4.8553772734271E-10_wp, &
      & 1.2356257410943E-03_wp, 6.8399234763729E-01_wp, 5.1995420254894E-09_wp, &
      & 1.2027896985957E-11_wp, 5.7304401546169E-07_wp, 2.4468959606775E-03_wp, &
      & 1.0040005130011E+00_wp, 0.0000000000000E+00_wp, 5.5310440465207E-02_wp, &
      & 8.7883262647393E-01_wp, 1.0298872576581E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 4.1485451554553E-02_wp, 9.5851454844545E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 6.7032432969941E-05_wp, 1.0329252259896E-01_wp, 3.9006942107954E-01_wp, &
      & 3.6487976107888E-01_wp, 0.0000000000000E+00_wp, 2.9607006073691E-03_wp, &
      & 9.9703929939263E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 7.5207498914890E-02_wp, 8.2630508485236E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 6.6472988033016E-02_wp, 9.3352701196698E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.5170748238388E-03_wp, &
      & 9.9748292517616E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.1199100233663E-05_wp, 6.1735949348506E-02_wp, &
      & 5.2547467287298E-01_wp, 2.3534849766078E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.1545009180054E-08_wp, 4.8156324470167E-03_wp, &
      & 6.8119098700326E-01_wp, 9.8640706987545E-08_wp, 9.0609306620312E-02_wp, &
      & 9.0939069337969E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 6.6527647187693E-02_wp, 9.3347235281231E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 5.0319522106180E-03_wp, 9.7236235992780E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.0833572554318E-11_wp, 3.0735372456843E-06_wp, 8.7114485461151E-03_wp, &
      & 8.5918441677693E-01_wp, 1.0676321498807E-09_wp, 4.8559648922814E-05_wp, &
      & 4.5302763715688E-02_wp, 8.3289153450766E-01_wp, 1.2364799599275E-02_wp, &
      & 0.0000000000000E+00_wp, 4.8553772734271E-10_wp, 1.2356257410943E-03_wp, &
      & 6.8399234763729E-01_wp, 5.1995420254894E-09_wp, 1.2027896985957E-11_wp, &
      & 5.7304401546169E-07_wp, 2.4468959606775E-03_wp, 1.0040005130011E+00_wp, &
      & 0.0000000000000E+00_wp, 1.7836181363501E-01_wp, 6.2456202118140E-01_wp, &
      & 2.3160354638897E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 8.7002281543997E-02_wp, 9.1299771845600E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 3.1404764970412E-06_wp, &
      & 6.3956931077352E-02_wp, 4.7846317290481E-01_wp, 2.9468368515225E-01_wp, &
      & 0.0000000000000E+00_wp, 1.1450448372045E-02_wp, 9.8854955162796E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 4.4395780891660E-02_wp, 8.5383545436310E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.2478649046186E-01_wp, &
      & 8.7521350953814E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.0110955182316E-02_wp, 9.8988904481768E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 5.6167960886507E-07_wp, 2.9078967459651E-02_wp, 6.1859560324107E-01_wp, &
      & 1.5104287585421E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.0080902930002E-08_wp, 4.6325621083875E-03_wp, 6.8133424494767E-01_wp, &
      & 8.7425251232161E-08_wp, 1.5789095742186E-01_wp, 8.4210904257814E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.2486469344048E-01_wp, 8.7513530655952E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.2874045548253E-05_wp, &
      & 9.7726740456727E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.7281235340132E-12_wp, &
      & 1.0839801911001E-06_wp, 5.9583550692028E-03_wp, 8.6187217041940E-01_wp, &
      & 1.9632749283518E-10_wp, 2.1818061591790E-05_wp, 3.6017036179768E-02_wp, &
      & 8.4547474113294E-01_wp, 8.8694849144924E-03_wp, 0.0000000000000E+00_wp, &
      & 4.1393746281189E-10_wp, 1.1763869822125E-03_wp, 6.8403870138130E-01_wp, &
      & 4.5070023665802E-09_wp, 0.0000000000000E+00_wp, 2.8208391560596E-10_wp, &
      & 1.0107828306943E-04_wp, 1.0063475070808E+00_wp, 0.0000000000000E+00_wp, &
      & 1.5052252044350E-01_wp, 6.7503900878048E-01_wp, 2.0966764053284E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.4123874345618E-01_wp, &
      & 8.5876125654382E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 3.2995570594716E-04_wp, 1.3239836014058E-01_wp, &
      & 3.3638195237351E-01_wp, 4.0116925621724E-01_wp, 0.0000000000000E+00_wp, &
      & 2.7738482922597E-02_wp, 9.7226151707740E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.2379486544581E-01_wp, &
      & 7.8289211602208E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.8811788794004E-01_wp, 8.1188211205996E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.5107682180162E-02_wp, 9.7489231781984E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.3960579256771E-04_wp, &
      & 9.0196670358265E-02_wp, 4.5739781308729E-01_wp, 2.9149367395236E-01_wp, &
      & 0.0000000000000E+00_wp, 5.8317983398498E-12_wp, 1.3182770972301E-06_wp, &
      & 1.8538244187819E-02_wp, 6.7044774462985E-01_wp, 6.6838194641483E-06_wp, &
      & 2.2640318425337E-01_wp, 7.7359681574663E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.8821108711455E-01_wp, &
      & 8.1178891288545E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.0126387421950E-02_wp, 9.6738367035521E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.9953637089135E-12_wp, 9.2422656035061E-08_wp, 2.2493337765938E-04_wp, &
      & 4.1001518442525E-02_wp, 8.2747340204614E-01_wp, 1.1163506217946E-06_wp, &
      & 1.2540172388574E-03_wp, 1.0999576750594E-01_wp, 7.3478924996110E-01_wp, &
      & 4.5973534870136E-02_wp, 1.7733574164893E-13_wp, 1.2868664921572E-07_wp, &
      & 6.8690163375128E-03_wp, 6.7958369863349E-01_wp, 7.7071301908717E-07_wp, &
      & 2.5715775485838E-07_wp, 1.6957699914291E-04_wp, 2.6115927712490E-02_wp, &
      & 9.8015581415428E-01_wp, 0.0000000000000E+00_wp, 1.5052252044350E-01_wp, &
      & 6.7503900878048E-01_wp, 2.0966764053284E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.4123874345618E-01_wp, 8.5876125654382E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 3.2995570594716E-04_wp, 1.3239836014058E-01_wp, 3.3638195237351E-01_wp, &
      & 4.0116925621724E-01_wp, 0.0000000000000E+00_wp, 2.7738482922597E-02_wp, &
      & 9.7226151707740E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.2379486544581E-01_wp, 7.8289211602208E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.8811788794004E-01_wp, 8.1188211205996E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.5107682180162E-02_wp, &
      & 9.7489231781984E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.3960579256771E-04_wp, 9.0196670358265E-02_wp, &
      & 4.5739781308729E-01_wp, 2.9149367395236E-01_wp, 0.0000000000000E+00_wp, &
      & 5.8317983398498E-12_wp, 1.3182770972301E-06_wp, 1.8538244187819E-02_wp, &
      & 6.7044774462985E-01_wp, 6.6838194641483E-06_wp, 2.2640318425337E-01_wp, &
      & 7.7359681574663E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.8821108711455E-01_wp, 8.1178891288545E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.0126387421950E-02_wp, 9.6738367035521E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.9953637089135E-12_wp, &
      & 9.2422656035061E-08_wp, 2.2493337765938E-04_wp, 4.1001518442525E-02_wp, &
      & 8.2747340204614E-01_wp, 1.1163506217946E-06_wp, 1.2540172388574E-03_wp, &
      & 1.0999576750594E-01_wp, 7.3478924996110E-01_wp, 4.5973534870136E-02_wp, &
      & 1.7733574164893E-13_wp, 1.2868664921572E-07_wp, 6.8690163375128E-03_wp, &
      & 6.7958369863349E-01_wp, 7.7071301908717E-07_wp, 2.5715775485838E-07_wp, &
      & 1.6957699914291E-04_wp, 2.6115927712490E-02_wp, 9.8015581415428E-01_wp, &
      & 0.0000000000000E+00_wp, 1.1619449180297E-01_wp, 7.4292303143376E-01_wp, &
      & 1.7688691875381E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 4.5511447319703E-02_wp, 9.5448855268030E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 9.2508392559746E-06_wp, &
      & 7.5824069699182E-02_wp, 4.4951921574350E-01_wp, 3.1887288951654E-01_wp, &
      & 0.0000000000000E+00_wp, 3.5066919289813E-03_wp, 9.9649330807102E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 4.4395780891660E-02_wp, 8.5383545436310E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 7.1925482699841E-02_wp, &
      & 9.2807451730016E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.9954258724411E-03_wp, 9.9700457412756E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.0238008000468E-06_wp, 3.8055054539114E-02_wp, 5.9074525520210E-01_wp, &
      & 1.7718982934717E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 3.4151558259365E-08_wp, 6.5639826011158E-03_wp, 6.7982280802516E-01_wp, &
      & 2.5898101368240E-07_wp, 9.7140202686718E-02_wp, 9.0285979731328E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 7.1982878187917E-02_wp, 9.2801712181208E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.1772574193346E-03_wp, &
      & 9.7612947599990E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 3.5689039339356E-13_wp, &
      & 3.8192206935584E-07_wp, 4.0712983354820E-03_wp, 8.6371381723556E-01_wp, &
      & 3.5979138804999E-11_wp, 9.7693789150509E-06_wp, 2.8536546460626E-02_wp, &
      & 8.5530890313171E-01_wp, 6.3404470879011E-03_wp, 0.0000000000000E+00_wp, &
      & 1.7396779713960E-09_wp, 1.8302935084816E-03_wp, 6.8352702270097E-01_wp, &
      & 1.6312609722908E-08_wp, 0.0000000000000E+00_wp, 7.9041698936705E-09_wp, &
      & 4.0780312005162E-04_wp, 1.0060406958243E+00_wp, 0.0000000000000E+00_wp, &
      & 2.1174800276996E-01_wp, 5.6964871898640E-01_wp, 2.5219928685690E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.4403797218710E-01_wp, &
      & 7.5596202781290E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 4.8890988223564E-04_wp, 1.4079823595826E-01_wp, &
      & 3.2239231943401E-01_wp, 4.0957756984420E-01_wp, 0.0000000000000E+00_wp, &
      & 7.5862834744969E-02_wp, 9.2413716525503E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.5224097067826E-01_wp, &
      & 7.5747542917585E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.9773509090075E-01_wp, 7.0226490909925E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 7.0607838000571E-02_wp, 9.2939216199943E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.2189803166994E-04_wp, &
      & 9.8892646451660E-02_wp, 4.3831806676090E-01_wp, 3.0633584406942E-01_wp, &
      & 0.0000000000000E+00_wp, 1.0892600979524E-09_wp, 1.9452547146752E-05_wp, &
      & 3.9425621360989E-02_wp, 6.5404367987276E-01_wp, 7.3196915833069E-05_wp, &
      & 3.3779631811082E-01_wp, 6.6220368188918E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 2.9783635509721E-01_wp, &
      & 7.0216364490279E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 9.6336440102884E-03_wp, 9.6786521862972E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 2.4358117918347E-08_wp, 2.1339725696048E-05_wp, 3.5282099660186E-03_wp, &
      & 1.0654366986150E-01_wp, 7.6039718337984E-01_wp, 9.5562006186987E-05_wp, &
      & 9.4981643353880E-03_wp, 1.7700942131591E-01_wp, 6.0914264505196E-01_wp, &
      & 9.7809437645041E-02_wp, 6.8362324068111E-11_wp, 3.0917858879577E-06_wp, &
      & 1.8135933744590E-02_wp, 6.7075660296166E-01_wp, 1.3284665444268E-05_wp, &
      & 2.8388831053016E-13_wp, 6.7329811561295E-08_wp, 9.9935157277430E-04_wp, &
      & 1.0054489358913E+00_wp, 0.0000000000000E+00_wp], &
      & [5, 16, 16])

   call get_structure(mol, "MB16-43", "01")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_gw_gen(error, mol, d4s, ref, with_cn=.true., with_q=.false.)

end subroutine test_gw_d4s_mb01


subroutine test_gw_d4_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   real(wp), parameter :: ref(5, 16, 1) = reshape([&
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
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp], &
      & [5, 16, 1])

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_gw_gen(error, mol, d4, ref, with_cn=.false., with_q=.true.)

end subroutine test_gw_d4_mb02


subroutine test_gw_d4_eeqbc_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   real(wp), parameter :: ref(5, 16, 1) = reshape([&
      & 1.2595414807089E+00_wp, 3.4627651498328E-03_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.0254439271401E+00_wp, &
      & 2.0298475069163E-03_wp, 4.5954321305318E-10_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.0048156894089E+00_wp, 1.8796029452983E-03_wp, &
      & 3.0506495005117E-10_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.1020354818356E+00_wp, 6.7928862292002E-03_wp, 5.5690274492868E-08_wp, &
      & 1.1605458622024E-03_wp, 0.0000000000000E+00_wp, 9.4187732431476E-01_wp, &
      & 5.4832647054101E-03_wp, 1.7870363924972E-08_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 6.9350750294192E-01_wp, 1.9066093885080E-03_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.0650204821510E+00_wp, 2.9279828143293E-03_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.1303254622557E+00_wp, &
      & 3.1075210134921E-03_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 9.6837394592909E-01_wp, 2.7666277695184E-03_wp, &
      & 7.1117229823468E-10_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 7.1529148033840E-01_wp, 1.9664984822049E-03_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 9.9541854542085E-01_wp, &
      & 1.8620015159966E-03_wp, 3.0220664483454E-10_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 8.9906533280080E-01_wp, 6.8945505568332E-03_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.1487780707026E+00_wp, 1.3509553061686E-02_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.0926349093356E+00_wp, &
      & 3.0039011366331E-03_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.1619247510000E+00_wp, 3.1943946238488E-03_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.0447172264156E+00_wp, 2.0679073977604E-03_wp, 4.6816982905992E-10_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp], &
      & [5, 16, 1])

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol, qmod=d4_qmod%eeqbc)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_gw_gen(error, mol, d4, ref, with_cn=.false., with_q=.true.)

end subroutine test_gw_d4_eeqbc_mb02


subroutine test_gw_d4_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   real(wp), parameter :: ref(7, 16, 1) = reshape([&
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
      & 0.0000000000000E+00_wp], &
      & [7, 16, 1])

   call get_structure(mol, "MB16-43", "03")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_gw_gen(error, mol, d4, ref, with_cn=.true., with_q=.true.)

end subroutine test_gw_d4_mb03


subroutine test_dgw_d4_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4

   call get_structure(mol, "MB16-43", "04")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_dgw_gen(error, mol, d4, with_cn=.true., with_q=.false.)

end subroutine test_dgw_d4_mb04


subroutine test_dgw_d4s_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s

   call get_structure(mol, "MB16-43", "04")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_dgw_gen(error, mol, d4s, with_cn=.true., with_q=.false.)

end subroutine test_dgw_d4s_mb04


subroutine test_dgw_d4_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4

   call get_structure(mol, "MB16-43", "05")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_dgw_gen(error, mol, d4, with_cn=.false., with_q=.true.)

end subroutine test_dgw_d4_mb05

subroutine test_dgw_d4s_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s

   call get_structure(mol, "MB16-43", "05")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_dgw_gen(error, mol, d4s, with_cn=.false., with_q=.true.)

end subroutine test_dgw_d4s_mb05


subroutine test_dgw_d4s_eeqbc_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s

   call get_structure(mol, "MB16-43", "05")
   call new_d4s_model(error, d4s, mol, qmod=d4_qmod%eeqbc)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_dgw_gen(error, mol, d4s, with_cn=.false., with_q=.true.)

end subroutine test_dgw_d4s_eeqbc_mb05


subroutine test_dgw_d4_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4

   call get_structure(mol, "MB16-43", "06")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_dgw_gen(error, mol, d4, with_cn=.true., with_q=.true.)

end subroutine test_dgw_d4_mb06

subroutine test_dgw_d4s_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s

   call get_structure(mol, "MB16-43", "06")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_dgw_gen(error, mol, d4s, with_cn=.true., with_q=.true.)

end subroutine test_dgw_d4s_mb06


subroutine test_gw_d4_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: qat(16) = [&
      &-1.57324183192355E-1_wp, 1.65228028672395E-1_wp, 3.22320366812437E-1_wp, &
      & 3.63579860576008E-2_wp, 4.85677294229959E-2_wp,-3.59193331069774E-1_wp, &
      &-1.93844259127416E-1_wp,-3.86492630088218E-1_wp, 3.10104713332147E-1_wp, &
      & 8.34803229863654E-2_wp,-3.62667644019899E-1_wp, 3.64142434058147E-1_wp, &
      & 3.34644499696670E-1_wp,-4.69889877462762E-1_wp,-1.89224201365947E-1_wp, &
      & 4.53790045287620E-1_wp]
   real(wp), parameter :: ref(7, 16, 1) = reshape([&
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
      & 0.0000000000000E+00_wp], &
      & [7, 16, 1])

   type(structure_type) :: mol
   type(d4_model) :: d4

   call get_structure(mol, "MB16-43", "06")
   call new_d4_model(error, d4, mol, qmod=d4_qmod%gfn2)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_gw_gen(error, mol, d4, ref, with_cn=.true., with_q=.true., qat=qat)

end subroutine test_gw_d4_mb07


subroutine test_dgw_d4_mb08(error)

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
   type(d4_model) :: d4

   call get_structure(mol, "MB16-43", "08")
   call new_d4_model(error, d4, mol, qmod=d4_qmod%gfn2)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_dgw_gen(error, mol, d4, with_cn=.true., with_q=.true., qat=qat)

end subroutine test_dgw_d4_mb08

subroutine test_dgw_d4s_mb08(error)

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
   type(d4s_model) :: d4s

   call get_structure(mol, "MB16-43", "08")
   call new_d4s_model(error, d4s, mol, qmod=d4_qmod%gfn2)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_dgw_gen(error, mol, d4s, with_cn=.true., with_q=.true., qat=qat)

end subroutine test_dgw_d4s_mb08


subroutine test_pol_d4_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   real(wp), parameter :: ref(16) = [&
      & 2.7207807469636E+00_wp, 2.7207697669106E+00_wp, 2.7340295240471E+00_wp, &
      & 2.7373201241818E+00_wp, 2.0051462957754E+01_wp, 2.7348043044018E+00_wp, &
      & 6.9923962318719E+00_wp, 9.2345547677329E+00_wp, 2.7338175995796E+00_wp, &
      & 2.7347959316428E+00_wp, 2.3018401095949E+01_wp, 2.7347926176379E+00_wp, &
      & 1.5079046616252E+01_wp, 3.5315073387745E+00_wp, 2.7340548152913E+00_wp, &
      & 9.2442543133947E+00_wp]

   call get_structure(mol, "MB16-43", "09")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_pol_gen(error, mol, d4, ref, with_cn=.true., with_q=.false.)

end subroutine test_pol_d4_mb09

subroutine test_pol_d4s_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   real(wp), parameter :: ref(16) = [&
      & 2.0796638857659E+00_wp, 2.0317560362159E+00_wp, 2.1362479991008E+00_wp, &
      & 2.4840420351444E+00_wp, 3.5144617579624E+01_wp, 2.1798399760895E+00_wp, &
      & 8.7285804463012E+00_wp, 1.0754720091309E+01_wp, 2.1068749047750E+00_wp, &
      & 2.2167964014782E+00_wp, 2.4223007635619E+01_wp, 2.1626725547918E+00_wp, &
      & 1.5342832501948E+01_wp, 4.0831999012216E+00_wp, 2.2319948533110E+00_wp, &
      & 1.0723200212288E+01_wp]

   call get_structure(mol, "MB16-43", "09")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_pol_gen(error, mol, d4s, ref, with_cn=.true., with_q=.true.)

end subroutine test_pol_d4s_mb09

subroutine test_dpol_d4_mb10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4

   call get_structure(mol, "MB16-43", "10")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_dpol_gen(error, mol, d4, with_cn=.true., with_q=.true.)

end subroutine test_dpol_d4_mb10

subroutine test_dpol_d4s_mb10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s

   call get_structure(mol, "MB16-43", "10")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_dpol_gen(error, mol, d4s, with_cn=.true., with_q=.true.)

end subroutine test_dpol_d4s_mb10


subroutine test_d4_model_error(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4

   integer, parameter :: nat = 2
   character(len=*), parameter :: sym(nat) = [character(len=4) :: "Rf", "Db"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp],&
      & [3, nat])

   call new(mol, sym, xyz)
   call new_d4_model(error, d4, mol)
  
end subroutine test_d4_model_error

subroutine test_d4s_model_error(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s

   integer, parameter :: nat = 2
   character(len=*), parameter :: sym(nat) = [character(len=4) :: "H", "Db"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp],&
      & [3, nat])

   call new(mol, sym, xyz)
   call new_d4s_model(error, d4s, mol)
  
end subroutine test_d4s_model_error

subroutine test_model_wrapper(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_model), allocatable :: disp

   integer, parameter :: nat = 2
   character(len=*), parameter :: sym(nat) = [character(len=4) :: "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp],&
      & [3, nat])

   call new(mol, sym, xyz)

   call new_dispersion_model(error, disp, mol, "d4")
   if (allocated(error)) then
      call test_failed(error, "D4 model could not be created")
      return
   end if

   ! check if the model is of type d4_model
   select type (disp)
      type is (d4_model)
      class default
         call test_failed(error, "Expected d4_model but got something else")
         return
   end select

   deallocate(disp)
   call new_dispersion_model(error, disp, mol, "D4S")
   if (allocated(error)) then
      call test_failed(error, "D4S model could not be created")
      return
   end if

   ! check if the model is of type d4s_model
   select type (disp)
      type is (d4s_model)
      class default
         call test_failed(error, "Expected d4s_model but got something else")
         return
   end select

end subroutine test_model_wrapper

subroutine test_model_wrapper_fail(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_model), allocatable :: disp

   integer, parameter :: nat = 2
   character(len=*), parameter :: sym(nat) = [character(len=4) :: "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp],&
      & [3, nat])

   call new(mol, sym, xyz)

   call new_dispersion_model(error, disp, mol, "wrong")

   if (.not. allocated(error)) then
      call test_failed(error, "Expected an error for key 'wrong'")
      return
   end if

   if (allocated(disp)) then
      call test_failed(error, "Model should not be allocated with invalid key")
      return
   end if

end subroutine test_model_wrapper_fail

end module test_model
