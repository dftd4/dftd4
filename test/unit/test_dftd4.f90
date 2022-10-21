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

module test_dftd4
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use dftd4
   implicit none
   private

   public :: collect_dftd4

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_dftd4(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("PBE-D4", test_pbed4_mb01), &
      & new_unittest("B97-D4", test_b97d4_mb02), &
      & new_unittest("TPSS-D4", test_tpssd4_mb03), &
      & new_unittest("PWPB95-D4", test_pwpb95d4_mb04), &
      & new_unittest("B2PLYP-D4", test_b2plypd4_mb05), &
      & new_unittest("PW6B95-D4", test_pw6b95d4_mb06), &
      & new_unittest("OLYP-D4", test_olypd4_mb07), &
      & new_unittest("PBE0-D4", test_pbe0d4_mb08), &
      & new_unittest("RPBE-D4-ATM", test_rpbed4atm_mb09), &
      & new_unittest("B2GPPLYP-D4-ATM", test_b2gpplypd4atm_mb10), &
      & new_unittest("LH14t-calPBE-D4-ATM", test_lh14tcalpbed4atm_mb11), &
      & new_unittest("B1B95-D4-ATM", test_b1b95d4atm_mb12), &
      & new_unittest("M06L-D4-ATM", test_m06ld4atm_mb13), &
      & new_unittest("TPSSh-D4-ATM", test_tpsshd4atm_mb14), &
      & new_unittest("HF-D4-ATM", test_hfd4atm_mb15), &
      & new_unittest("CAM-B3LYP-D4-ATM", test_camb3lypd4atm_mb16) &
      & ]

end subroutine collect_dftd4


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
   call get_dispersion(mol, d4, param, realspace_cutoff(), energy)

   call check(error, energy, ref, thr=thr)
   if (allocated(error)) then
      print*,energy
   end if

end subroutine test_dftd4_gen


subroutine test_numgrad(error, mol, param)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   integer :: iat, ic
   type(d4_model) :: d4
   real(wp) :: energy, er, el, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   call new_d4_model(d4, mol)

   do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_dispersion(mol, d4, param, realspace_cutoff(), er)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_dispersion(mol, d4, param, realspace_cutoff(), el)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(er - el)/step
      end do
   end do

   call get_dispersion(mol, d4, param, realspace_cutoff(), energy, gradient, sigma)

   if (any(abs(gradient - numgrad) > thr2)) then
      call test_failed(error, "Gradient of dispersion energy does not match")
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_numgrad


subroutine test_numsigma(error, mol, param)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   integer :: ic, jc
   type(d4_model) :: d4
   real(wp) :: energy, er, el, sigma(3, 3), eps(3, 3), numsigma(3, 3)
   real(wp), allocatable :: gradient(:, :), xyz(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), xyz(3, mol%nat))
   call new_d4_model(d4, mol)

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   do ic = 1, 3
      do jc = 1, 3
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call get_dispersion(mol, d4, param, realspace_cutoff(), er)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call get_dispersion(mol, d4, param, realspace_cutoff(), el)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         numsigma(jc, ic) = 0.5_wp*(er - el)/step
      end do
   end do

   call get_dispersion(mol, d4, param, realspace_cutoff(), energy, gradient, sigma)

   if (any(abs(sigma - numsigma) > thr2)) then
      call test_failed(error, "Strain derivatives do not match")
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_numsigma


subroutine test_pbed4_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "MB16-43", "01")
   call test_dftd4_gen(error, mol, param, -1.8578752883363366E-002_wp)

end subroutine test_pbed4_mb01


subroutine test_b97d4_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.69460052_wp, a1 = 0.28904684_wp, a2 = 4.13407323_wp)

   call get_structure(mol, "MB16-43", "02")
   call test_dftd4_gen(error, mol, param, -8.9181168937810723E-002_wp)

end subroutine test_b97d4_mb02


subroutine test_tpssd4_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.91130849_wp, a1 = 0.43332851_wp, a2 = 4.56986797_wp)

   call get_structure(mol, "MB16-43", "03")
   call test_dftd4_gen(error, mol, param, -2.4695638764787930E-002_wp)

end subroutine test_tpssd4_mb03


subroutine test_pwpb95d4_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 0.82_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = -0.34639127_wp, a1 = 0.41080636_wp, a2 = 3.83878274_wp)

   call get_structure(mol, "MB16-43", "04")
   call test_dftd4_gen(error, mol, param, -9.5128100471706181E-003_wp)

end subroutine test_pwpb95d4_mb04


subroutine test_b2plypd4_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 0.64_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.15117773_wp, a1 = 0.42666167_wp, a2 = 4.73635790_wp)

   call get_structure(mol, "MB16-43", "05")
   call test_numgrad(error, mol, param)

end subroutine test_b2plypd4_mb05


subroutine test_pw6b95d4_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = -0.31629935_wp, a1 = 0.03999357_wp, a2 = 5.83690254_wp)

   call get_structure(mol, "MB16-43", "06")
   call test_numgrad(error, mol, param)

end subroutine test_pw6b95d4_mb06


subroutine test_olypd4_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 2.74836820_wp, a1 = 0.60184498_wp, a2 = 2.53292167_wp)

   call get_structure(mol, "MB16-43", "07")
   call test_numgrad(error, mol, param)

end subroutine test_olypd4_mb07


subroutine test_pbe0d4_mb08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.20065498_wp, a1 = 0.40085597_wp, a2 = 5.02928789_wp)

   call get_structure(mol, "MB16-43", "08")
   call test_numsigma(error, mol, param)

end subroutine test_pbe0d4_mb08


subroutine test_rpbed4atm_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.31183787_wp, a1 = 0.46169493_wp, a2 = 3.15711757_wp)

   call get_structure(mol, "MB16-43", "09")
   call test_dftd4_gen(error, mol, param, -4.5140422485299259E-002_wp)

end subroutine test_rpbed4atm_mb09


subroutine test_b2gpplypd4atm_mb10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 0.56_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 0.94633372_wp, a1 = 0.42907301_wp, a2 = 5.18802602_wp)

   call get_structure(mol, "MB16-43", "10")
   call test_dftd4_gen(error, mol, param, -9.6812427202205668E-003_wp)

end subroutine test_b2gpplypd4atm_mb10


subroutine test_lh14tcalpbed4atm_mb11(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.27677253_wp, a1 = 0.38128670_wp, a2 = 4.91698883_wp)

   call get_structure(mol, "MB16-43", "11")
   call test_dftd4_gen(error, mol, param, -1.7460015867914524E-002_wp)

end subroutine test_lh14tcalpbed4atm_mb11


subroutine test_b1b95d4atm_mb12(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.27701162_wp, a1 = 0.40554715_wp, a2 = 4.63323074_wp)

   call get_structure(mol, "MB16-43", "12")
   call test_dftd4_gen(error, mol, param, -2.5712178361964221E-002_wp)

end subroutine test_b1b95d4atm_mb12


subroutine test_m06ld4atm_mb13(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 0.59493760_wp, a1 = 0.71422359_wp, a2 = 6.35314182_wp)

   call get_structure(mol, "MB16-43", "13")
   call test_numgrad(error, mol, param)

end subroutine test_m06ld4atm_mb13


subroutine test_tpsshd4atm_mb14(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.85897750_wp, a1 = 0.44286966_wp, a2 = 4.60230534_wp)

   call get_structure(mol, "MB16-43", "14")
   call test_numgrad(error, mol, param)

end subroutine test_tpsshd4atm_mb14


subroutine test_hfd4atm_mb15(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.61679827_wp, a1 = 0.44959224_wp, a2 = 3.35743605_wp)

   call get_structure(mol, "MB16-43", "15")
   call test_numgrad(error, mol, param)

end subroutine test_hfd4atm_mb15


subroutine test_camb3lypd4atm_mb16(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.74407961_wp, a1 = 0.40137870_wp, a2 = 5.18731225_wp)

   call get_structure(mol, "MB16-43", "16")
   call test_numsigma(error, mol, param)

end subroutine test_camb3lypd4atm_mb16


end module test_dftd4
