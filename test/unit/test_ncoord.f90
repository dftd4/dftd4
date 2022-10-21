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

module test_ncoord
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, &
      & test_failed
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use dftd4_cutoff, only : get_lattice_points
   use dftd4_data, only : get_covalent_rad, get_electronegativity
   use dftd4_ncoord
   implicit none
   private

   public :: collect_ncoord

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_ncoord(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("cn-mb01", test_cn_mb01), &
      & new_unittest("cn-mb02", test_cn_mb02), &
      & new_unittest("cn-mb03", test_cn_mb03), &
      & new_unittest("cn-acetic", test_cn_acetic), &
      & new_unittest("dcndr-mb04", test_dcndr_mb04), &
      & new_unittest("dcndr-mb05", test_dcndr_mb05), &
      & new_unittest("dcndr-ammonia", test_dcndr_ammonia), &
      & new_unittest("dcndL-mb06", test_dcndL_mb06), &
      & new_unittest("dcndL-mb07", test_dcndL_mb07), &
      & new_unittest("dcndL-antracene", test_dcndL_anthracene) &
      & ]

end subroutine collect_ncoord


subroutine test_cn_gen(error, mol, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Reference CNs
   real(wp), intent(in) :: ref(:)

   real(wp), allocatable :: cn(:), rcov(:), en(:)
   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :)

   allocate(rcov(mol%nid), en(mol%nid), cn(mol%nat))
   rcov(:) = get_covalent_rad(mol%num)
   en(:) = get_electronegativity(mol%num)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call get_coordination_number(mol, lattr, cutoff, rcov, en, cn)

   if (any(abs(cn - ref) > thr)) then
      call test_failed(error, "Coordination numbers do not match")
      print'(*(6x,"&",3(es20.14e1,"_wp",:,","),1x,"&",/))', cn
   end if

end subroutine test_cn_gen


subroutine test_numgrad(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   integer :: iat, ic, mat
   real(wp), allocatable :: cn(:), rcov(:), en(:), cnr(:), cnl(:)
   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: numdr(:, :, :)
   real(wp), allocatable :: lattr(:, :)
   real(wp), parameter :: cutoff = 20.0_wp
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(rcov(mol%nid), en(mol%nid), cn(mol%nat), cnr(mol%nat), cnl(mol%nat), &
      & dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & numdr(3, mol%nat, mol%nat))
   rcov(:) = get_covalent_rad(mol%num)
   en(:) = get_electronegativity(mol%num)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   if (any(mol%periodic)) then
      mat = min(mol%nat, 5)
   else
      mat = mol%nat
   end if

   do iat = 1, mat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_coordination_number(mol, lattr, cutoff, rcov, en, cnr)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_coordination_number(mol, lattr, cutoff, rcov, en, cnl)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numdr(ic, iat, :) = 0.5_wp*(cnr - cnl)/step
      end do
   end do

   call get_coordination_number(mol, lattr, cutoff, rcov, en, cn, dcndr, dcndL)

   if (any(abs(dcndr(:, :mat, :) - numdr(:, :mat, :)) > thr2)) then
      call test_failed(error, "Derivative of coordination number does not match")
   end if

end subroutine test_numgrad


subroutine test_numsigma(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   integer :: ic, jc
   real(wp) :: eps(3, 3)
   real(wp), allocatable :: cn(:), rcov(:), en(:), cnr(:), cnl(:), xyz(:, :)
   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: numdL(:, :, :)
   real(wp), allocatable :: lattr(:, :), trans(:, :)
   real(wp), parameter :: cutoff = 20.0_wp
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(rcov(mol%nid), en(mol%nid), cn(mol%nat), cnr(mol%nat), cnl(mol%nat), &
      & dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), xyz(3, mol%nat), &
      & numdL(3, 3, mol%nat))
   rcov(:) = get_covalent_rad(mol%num)
   en(:) = get_electronegativity(mol%num)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   trans = lattr
   do ic = 1, 3
      do jc = 1, 3
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call get_coordination_number(mol, lattr, cutoff, rcov, en, cnr)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call get_coordination_number(mol, lattr, cutoff, rcov, en, cnl)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         lattr(:, :) = trans
         numdL(jc, ic, :) = 0.5_wp*(cnr - cnl)/step
      end do
   end do

   call get_coordination_number(mol, lattr, cutoff, rcov, en, cn, dcndr, dcndL)

   if (any(abs(dcndL - numdL) > thr2)) then
      call test_failed(error, "Derivative of coordination number does not match")
   end if

end subroutine test_numsigma


subroutine test_cn_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      & 3.07349677110402E+0_wp, 9.31461605116103E-1_wp, 1.43709439375839E+0_wp, &
      & 1.33309431581960E+0_wp, 7.20743527030337E-1_wp, 8.59659004770982E-1_wp, &
      & 1.35782158177921E+0_wp, 1.53940006996025E+0_wp, 3.19400368195259E+0_wp, &
      & 8.12162111631342E-1_wp, 8.59533443784854E-1_wp, 1.53347108155587E+0_wp, &
      & 4.23314989525721E+0_wp, 3.03048504567396E+0_wp, 3.45229319488306E+0_wp, &
      & 4.28478289652264E+0_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb01


subroutine test_cn_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      & 9.20259141516190E-1_wp, 3.29216939906043E+0_wp, 3.51944438412931E+0_wp, &
      & 2.25877973040028E+0_wp, 4.46999073626179E+0_wp, 8.18916367808423E-1_wp, &
      & 9.28914937407466E-1_wp, 9.30833050893587E-1_wp, 4.60708718003244E+0_wp, &
      & 8.18343168300509E-1_wp, 3.70959638795740E+0_wp, 2.87405845608016E+0_wp, &
      & 1.24015900552686E+0_wp, 9.11079070954527E-1_wp, 1.57258868344791E+0_wp, &
      & 1.67284525339418E+0_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb02


subroutine test_cn_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      & 3.70329575672631E+0_wp, 2.11476582851627E+0_wp, 9.23682826697708E-1_wp, &
      & 3.88999767055963E+0_wp, 6.17490081567969E+0_wp, 4.10595506858888E+0_wp, &
      & 4.22916534938705E+0_wp, 1.04287687415719E+0_wp, 9.23686985143960E-1_wp, &
      & 9.24487848931390E-1_wp, 1.22927636800717E+0_wp, 2.59853001985457E+0_wp, &
      & 4.30015470969650E+0_wp, 8.29081650895103E-1_wp, 2.65239010637793E+0_wp, &
      & 1.19840431336618E+0_wp]

   call get_structure(mol, "MB16-43", "03")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb03


subroutine test_cn_acetic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(32) = [&
      & 7.97759158536926E-1_wp, 7.97780729307072E-1_wp, 7.97759174018738E-1_wp, &
      & 7.97777897618605E-1_wp, 9.24187491895463E-1_wp, 9.24182909365931E-1_wp, &
      & 9.24173453029854E-1_wp, 9.23551144552296E-1_wp, 9.23548004418282E-1_wp, &
      & 9.23546232108523E-1_wp, 9.23543089852042E-1_wp, 9.23659744859844E-1_wp, &
      & 9.23675500925197E-1_wp, 9.23670884769485E-1_wp, 9.23665145186263E-1_wp, &
      & 9.24185332194868E-1_wp, 2.68674571323560E+0_wp, 2.68674659839631E+0_wp, &
      & 2.68674144410079E+0_wp, 2.68674612746527E+0_wp, 3.74955132990418E+0_wp, &
      & 3.74956146465094E+0_wp, 3.74954089184621E+0_wp, 3.74954806564687E+0_wp, &
      & 1.64848046237829E+0_wp, 1.64850720288453E+0_wp, 1.64848691881326E+0_wp, &
      & 1.64850679513430E+0_wp, 8.60766216408348E-1_wp, 8.60767776978061E-1_wp, &
      & 8.60769276517625E-1_wp, 8.60775180678160E-1_wp]

   call get_structure(mol, "X23", "acetic")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_acetic


subroutine test_dcndr_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol)

end subroutine test_dcndr_mb04


subroutine test_dcndr_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numgrad(error, mol)

end subroutine test_dcndr_mb05


subroutine test_dcndr_ammonia(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "ammonia")
   call test_numgrad(error, mol)

end subroutine test_dcndr_ammonia


subroutine test_dcndL_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "06")
   call test_numsigma(error, mol)

end subroutine test_dcndL_mb06


subroutine test_dcndL_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "07")
   call test_numsigma(error, mol)

end subroutine test_dcndL_mb07


subroutine test_dcndL_anthracene(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "anthracene")
   call test_numsigma(error, mol)

end subroutine test_dcndL_anthracene


end module test_ncoord
