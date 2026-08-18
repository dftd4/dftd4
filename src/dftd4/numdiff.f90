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

!> Numerical differentiation of DFT-D4 model
module dftd4_numdiff
   use, intrinsic :: iso_fortran_env, only : error_unit
   use dftd4_cutoff, only : realspace_cutoff
   use dftd4_damping, only : damping_param
   use dftd4_damping_rational, only : rational_damping_param
   use dftd4_disp, only : get_dispersion, get_dispersion2, get_dispersion3_hessian
   use dftd4_model, only : dispersion_model, d4_model, d4s_model
   use dftd4_partition, only : work_partition
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   implicit none
   private

   public :: get_dispersion_hessian


contains


!> Evaluate Hessian matrix using an analytical ATM contribution when available
subroutine get_dispersion_hessian(mol, disp, param, cutoff, hessian, partition)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_dispersion_hessian

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Dispersion Hessian
   real(wp), intent(out) :: hessian(:, :, :, :)

   !> Work partition of the interaction loops, defaults to the complete work
   type(work_partition), intent(in), optional :: partition

   integer :: iat, ix, jat, jx, ii, jj, ndim
   logical :: analytical_atm
   real(wp), parameter :: step = 1.0e-4_wp
   type(structure_type) :: displ
   type(error_type), allocatable :: error
   real(wp) :: el, er
   real(wp), allocatable :: gl(:, :), gr(:, :), sl(:, :), sr(:, :), hessian3(:, :)

   analytical_atm = .false.
   select type(param)
   class is (rational_damping_param)
      select type(disp)
      type is (d4_model)
         analytical_atm = .true.
      type is (d4s_model)
         analytical_atm = .true.
      end select
   end select

   hessian(:, :, :, :) = 0.0_wp
   !$omp parallel default(none) &
   !$omp private(iat, ix, displ, er, el, gr, gl, sr, sl) &
   !$omp shared(mol, disp, param, cutoff, hessian, partition, analytical_atm)
   displ = mol
   allocate(gl(3, mol%nat), gr(3, mol%nat), sl(3, 3), sr(3, 3))
   !$omp do schedule(dynamic) collapse(2)
   do iat = 1, mol%nat
      do ix = 1, 3
         displ%xyz(ix, iat) = mol%xyz(ix, iat) + step
         if (analytical_atm) then
            call get_dispersion2(displ, disp, param, cutoff, el, gl, sl, partition)
         else
            call get_dispersion(displ, disp, param, cutoff, el, gl, sl, partition)
         end if

         displ%xyz(ix, iat) = mol%xyz(ix, iat) - step
         if (analytical_atm) then
            call get_dispersion2(displ, disp, param, cutoff, er, gr, sr, partition)
         else
            call get_dispersion(displ, disp, param, cutoff, er, gr, sr, partition)
         end if

         displ%xyz(ix, iat) = mol%xyz(ix, iat)
         hessian(:, :, ix, iat) = (gl - gr) / (2 * step)
      end do
   end do
   !$omp end parallel

   if (analytical_atm) then
      ndim = 3*mol%nat
      allocate(hessian3(ndim, ndim))
      call get_dispersion3_hessian(error, mol, disp, param, cutoff, hessian3, partition)
      if (allocated(error)) then
         write(error_unit, '("[Error]:", 1x, a)') error%message
         error stop
      end if

      !$omp parallel do collapse(4) schedule(static) default(none) &
      !$omp shared(mol, hessian, hessian3) private(jat, jx, iat, ix, ii, jj)
      do jat = 1, mol%nat
         do jx = 1, 3
            do iat = 1, mol%nat
               do ix = 1, 3
                  jj = 3*(jat - 1) + jx
                  ii = 3*(iat - 1) + ix
                  hessian(ix, iat, jx, jat) = hessian(ix, iat, jx, jat) + hessian3(ii, jj)
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end if

end subroutine get_dispersion_hessian

end module dftd4_numdiff
