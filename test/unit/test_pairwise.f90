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

module test_pairwise
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use dftd4
   implicit none
   private

   public :: collect_pairwise

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   type(realspace_cutoff), parameter :: cutoff = &
      & realspace_cutoff(cn=30_wp, disp2=60.0_wp, disp3=15.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_pairwise(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("PBE-D4", test_pbed4_mb01), &
      & new_unittest("B97-D4", test_b97d4_mb02), &
      & new_unittest("TPSS-D4", test_tpssd4_ammonia) &
      & ]

end subroutine collect_pairwise


subroutine test_dftd4_pairwise(error, mol, param)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   type(d4_model) :: d4
   real(wp) :: energy
   real(wp), allocatable :: energy2(:, :), energy3(:, :)

   allocate(energy2(mol%nat, mol%nat), energy3(mol%nat, mol%nat))
   call new_d4_model(d4, mol)
   call get_dispersion(mol, d4, param, cutoff, energy)
   call get_pairwise_dispersion(mol, d4, param, cutoff, energy2, energy3)

   call check(error, energy, sum(energy2) + sum(energy3), thr=thr)
   if (allocated(error)) then
      print*,energy, sum(energy2) + sum(energy3)
   end if

end subroutine test_dftd4_pairwise


subroutine test_pbed4_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "MB16-43", "01")
   call test_dftd4_pairwise(error, mol, param)

end subroutine test_pbed4_mb01


subroutine test_b97d4_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.69460052_wp, a1 = 0.28904684_wp, a2 = 4.13407323_wp)

   call get_structure(mol, "MB16-43", "02")
   call test_dftd4_pairwise(error, mol, param)

end subroutine test_b97d4_mb02


subroutine test_tpssd4_ammonia(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.76596355_wp, a1 = 0.42822303_wp, a2 = 4.54257102_wp )

   call get_structure(mol, "X23", "ammonia")
   call test_dftd4_pairwise(error, mol, param)

end subroutine test_tpssd4_ammonia


end module test_pairwise
