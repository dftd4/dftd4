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

!> Analytical Hessian dispatch and numerical reference for the DFT-D4 model
module dftd4_numdiff
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use dftd4_cutoff, only : realspace_cutoff
   use dftd4_damping, only : damping_param
   use dftd4_damping_rational, only : rational_damping_param
   use dftd4_disp, only : get_dispersion, get_dispersion2_hessian, get_dispersion3_hessian
   use dftd4_model, only : dispersion_model, d4_model, d4s_model
   use dftd4_partition, only : work_partition
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
   implicit none
   private

   public :: get_dispersion_hessian, get_dispersion_hessian_numdiff


contains


!> Dispersion Hessian: analytic D4/D4S rational damping, numerical fallback otherwise.
subroutine get_dispersion_hessian(mol, disp, param, cutoff, hessian, partition, error)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_dispersion_hessian
   class(structure_type), intent(in) :: mol
   class(dispersion_model), intent(in) :: disp
   class(damping_param), intent(in) :: param
   type(realspace_cutoff), intent(in) :: cutoff
   real(wp), intent(out) :: hessian(:, :, :, :)
   type(work_partition), intent(in), optional :: partition
   type(error_type), allocatable, intent(out), optional :: error

   logical :: analytic
   integer :: ndim
   real(wp), allocatable :: hess2(:, :), hess3(:, :)
   type(error_type), allocatable :: local_error

   analytic = .false.
   select type(param)
   class is (rational_damping_param)
      select type(disp)
      type is (d4_model)
         analytic = .true.
      type is (d4s_model)
         analytic = .true.
      class default
         continue
      end select
   class default
      continue
   end select
   analytic = analytic .and. (all(mol%periodic) .or. .not. any(mol%periodic))
   if (.not. analytic) then
      call get_dispersion_hessian_numdiff(mol, disp, param, cutoff, hessian, partition, error)
      return
   end if

   ndim = 3*mol%nat
   allocate(hess2(ndim, ndim), hess3(ndim, ndim))
   call get_dispersion2_hessian(local_error, mol, disp, param, cutoff, hess2, partition)
   if (.not. allocated(local_error)) then
      call get_dispersion3_hessian(local_error, mol, disp, param, cutoff, hess3, partition)
   end if
   if (allocated(local_error)) then
      hessian = ieee_value(0.0_wp, ieee_quiet_nan)
      if (present(error)) call move_alloc(local_error, error)
      return
   end if
   hess2(:, :) = hess2 + hess3
   hessian(:, :, :, :) = reshape(hess2, shape(hessian))
end subroutine get_dispersion_hessian


!> Independent central difference of the full dispersion gradient.
subroutine get_dispersion_hessian_numdiff(mol, disp, param, cutoff, hessian, partition, error, step)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_dispersion_hessian_numdiff

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

   !> Error on failure; the entire Hessian is set to NaN on failure.
   type(error_type), allocatable, intent(out), optional :: error

   !> Displacement in bohr (default 1e-4)
   real(wp), intent(in), optional :: step

   integer :: iat, ix
   real(wp) :: delta
   type(structure_type) :: displ
   type(error_type), allocatable :: local_error, first_error
   real(wp) :: el, er
   real(wp), allocatable :: gl(:, :), gr(:, :), sl(:, :), sr(:, :)

   delta = 1.0e-4_wp
   if (present(step)) delta = step
   if (delta <= 0.0_wp) then
      hessian = ieee_value(0.0_wp, ieee_quiet_nan)
      if (present(error)) call fatal_error(error, "Hessian finite-difference step must be positive")
      return
   end if

   hessian(:, :, :, :) = 0.0_wp
   !$omp parallel default(none) &
   !$omp private(iat, ix, displ, er, el, gr, gl, sr, sl, local_error) &
   !$omp shared(mol, disp, param, cutoff, hessian, partition, delta, first_error)
   displ = mol
   allocate(gl(3, mol%nat), gr(3, mol%nat), sl(3, 3), sr(3, 3))
   !$omp do schedule(dynamic) collapse(2)
   do iat = 1, mol%nat
      do ix = 1, 3
         ! A failed worker skips its remaining displacements, retaining its error.
         if (allocated(local_error)) cycle
         displ%xyz(ix, iat) = mol%xyz(ix, iat) + delta
         call get_dispersion(displ, disp, param, cutoff, el, gl, sl, partition, local_error)
         if (allocated(local_error)) then
            displ%xyz(ix, iat) = mol%xyz(ix, iat)
            cycle
         end if

         displ%xyz(ix, iat) = mol%xyz(ix, iat) - delta
         call get_dispersion(displ, disp, param, cutoff, er, gr, sr, partition, local_error)

         displ%xyz(ix, iat) = mol%xyz(ix, iat)
         if (allocated(local_error)) cycle
         hessian(:, :, ix, iat) = (gl - gr) / (2 * delta)
      end do
   end do
   !$omp end do
   if (allocated(local_error)) then
      !$omp critical(dftd4_hessian_error)
      if (.not. allocated(first_error)) first_error = local_error
      !$omp end critical(dftd4_hessian_error)
   end if
   !$omp end parallel

   if (allocated(first_error)) then
      hessian = ieee_value(0.0_wp, ieee_quiet_nan)
      if (present(error)) call move_alloc(first_error, error)
      return
   end if

end subroutine get_dispersion_hessian_numdiff

end module dftd4_numdiff
