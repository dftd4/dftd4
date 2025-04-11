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

!> This is a compatibility module for dftd4 2.5.0 reproducing enough of the old API
!> to compile the interface with Vasp.
module dftd4_compat
   use mctc_env, only : wp
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use dftd4, only: structure_type, new, d4_model, new_d4_model, rational_damping_param, &
      & damping_param, get_rational_damping, get_dispersion, realspace_cutoff
   implicit none
   private

   public :: dftd_options, dftd_parameter, dftd_results, molecule, mctc_logger, ws_cell
   public :: dlat_to_cell, dlat_to_dvol, dlat_to_rlat, d4par, d4_calculation

   type :: dftd_options
      integer :: lmbd = 3
      integer :: refq = 5
      real(wp) :: wf = 6.0_wp
      real(wp) :: g_a
      real(wp) :: g_c
      logical :: lmolpol
      logical :: lenergy
      logical :: lgradient
      logical :: lhessian
      integer :: print_level
   end type

   type :: dftd_parameter
      real(wp) :: s6
      real(wp) :: s8
      real(wp) :: s9 = 1.0_wp
      real(wp) :: a1
      real(wp) :: a2
      real(wp) :: alp = 16.0_wp
   end type

   type :: dftd_results
      real(wp), allocatable :: energy
      real(wp), allocatable :: gradient(:, :)
      real(wp), allocatable :: lattice_gradient(:, :)
   end type

   type :: ws_cell
   end type

   type :: molecule
      integer, allocatable :: at(:)
      real(wp), allocatable :: xyz(:, :)
      real(wp), allocatable :: lattice(:, :)
      real(wp), allocatable :: chrg
      integer, allocatable :: npbc
      logical, allocatable :: pbc
      real(wp), allocatable :: volume
      real(wp), allocatable :: cellpar(:)
      real(wp), allocatable :: rec_lat(:, :)
      type(ws_cell) :: wsc
   contains
      procedure :: allocate
      procedure :: wrap_back
      procedure :: calculate_distances
   end type

   type :: mctc_logger
      logical :: sane = .true.
   end type

contains

subroutine allocate(self, n, l)
   class(molecule), intent(inout) :: self
   integer, intent(in) :: n
   logical, intent(in) :: l
end subroutine allocate

subroutine wrap_back(self)
   class(molecule), intent(inout) :: self
end subroutine wrap_back

subroutine calculate_distances(self)
   class(molecule), intent(inout) :: self
end subroutine calculate_distances

subroutine d4_calculation(io, env, options, mol_, param_, res)
   integer, intent(in) :: io
   type(mctc_logger), intent(inout) :: env
   type(dftd_options), intent(in) :: options
   type(molecule), intent(in) :: mol_
   type(dftd_parameter), intent(in) :: param_
   type(dftd_results), intent(out) :: res

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(rational_damping_param) :: param
   real(wp), allocatable :: energy, gradient(:, :), sigma(:, :)

   call new(mol, mol_%at, mol_%xyz, lattice=mol_%lattice)
   call new_d4_model(d4, mol, ga=options%g_a, gc=options%g_c, wf=options%wf)
   param = rational_damping_param(&
      & s6=param_%s6, &
      & s8=param_%s8, &
      & s9=param_%s9, &
      & a1=param_%a1, &
      & a2=param_%a2, &
      & alp=param_%alp)

   allocate(energy, gradient(3, mol%nat), sigma(3, 3))
   call get_dispersion(mol, d4, param, realspace_cutoff(), energy, gradient, sigma)
   call move_alloc(energy, res%energy)
   call move_alloc(gradient, res%gradient)
   res%lattice_gradient = matmul(sigma, transpose(matinv_3x3(mol%lattice)))

end subroutine d4_calculation

subroutine d4par(fname, param_, lmbd, env)
   character(len=*), intent(in) :: fname
   type(dftd_parameter), intent(out) :: param_
   integer, intent(in) :: lmbd
   type(mctc_logger) :: env

   class(damping_param), allocatable :: param

   call get_rational_damping(fname, param, merge(1.0_wp, 0.0_wp, lmbd == 3))
   if (allocated(param)) then
      select type(param)
      type is(rational_damping_param)
         env%sane = .true.
         param_%s6 = param%s6
         param_%s8 = param%s8
         param_%s9 = param%s9
         param_%a1 = param%a1
         param_%a2 = param%a2
         param_%alp = param%alp
      class default
         env%sane = .false.
      end select
   else
      env%sane = .false.
   end if
end subroutine d4par

subroutine dlat_to_cell(lattice, cellpar)
   real(wp), intent(in) :: lattice(3, 3)
   real(wp), intent(out), optional :: cellpar(6)
end subroutine dlat_to_cell

subroutine dlat_to_rlat(lattice, reclatt)
   real(wp), intent(in) :: lattice(3, 3)
   real(wp), intent(out), optional :: reclatt(3, 3)
end subroutine dlat_to_rlat

function dlat_to_dvol(lattice) result(vol)
   real(wp), intent(in) :: lattice(3, 3)
   real(wp) :: vol
   vol = matdet_3x3(lattice)
end function dlat_to_dvol

end module dftd4_compat


subroutine generate_wsc(mol, wsc)
   use dftd4_compat, only : molecule, ws_cell
   type(molecule), intent(inout) :: mol
   type(ws_cell), intent(inout) :: wsc
end subroutine generate_wsc

module class_set
   use dftd4_compat, only : dftd_options
end module class_set

module class_param
   use dftd4_compat, only : dftd_parameter
end module class_param

module class_molecule
   use dftd4_compat, only : molecule
end module class_molecule

module class_results
   use dftd4_compat, only : dftd_results
end module class_results

module class_wsc
   use dftd4_compat, only : ws_cell
end module class_wsc

module mctc_environment
   use dftd4_compat, only : mctc_logger
end module mctc_environment

module dispersion_calculator
   use dftd4_compat, only : d4_calculation
end module dispersion_calculator

module pbc_tools
   use dftd4_compat, only: dlat_to_cell, dlat_to_dvol, dlat_to_rlat
end module pbc_tools

module dfuncpar
   use dftd4_compat, only : d4par
end module dfuncpar
