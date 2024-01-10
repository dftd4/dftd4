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

!> Entry point for running single point calculations with dftd4
module dftd4_driver
   use, intrinsic :: iso_fortran_env, only : output_unit, input_unit
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type, read_structure, filetype
   use dftd4, only : get_dispersion, d4_model, new_d4_model, &
      realspace_cutoff, get_lattice_points, get_coordination_number, &
      damping_param, rational_damping_param, get_rational_damping, &
      get_properties, get_pairwise_dispersion, get_dispersion_hessian
   use dftd4_charge, only : get_charges
   use dftd4_output
   use dftd4_utils
   use dftd4_cli, only : cli_config, param_config, run_config
   use dftd4_help, only : header
   use dftd4_param, only : functional_group, get_functionals, &
      get_functional_id, p_r2scan_3c
   implicit none
   private

   public :: main

contains


!> Main entry point for the driver
subroutine main(config, error)

   !> Configuration for this driver
   class(cli_config), intent(in) :: config

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   select type(config)
   class default
      call fatal_error(error, "Unknown runtime selected")
   type is(run_config)
      call run_main(config, error)
   type is(param_config)
      call run_param(config, error)
   end select
end subroutine main


!> Entry point for the single point driver
subroutine run_main(config, error)

   !> Configuration for this driver
   type(run_config), intent(in) :: config

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   character(len=:), allocatable :: filename
   character(len=:), allocatable :: functional
   class(damping_param), allocatable :: param
   type(d4_model) :: d4
   real(wp) :: charge
   real(wp), allocatable :: energy, gradient(:, :), sigma(:, :), hessian(:, :, :, :)
   real(wp), allocatable :: pair_disp2(:, :), pair_disp3(:, :)
   real(wp), allocatable :: cn(:), q(:), c6(:, :), alpha(:)
   real(wp), allocatable :: s9
   real(wp) :: ga, gc
   integer :: stat, unit, is, id
   logical :: exist

   if (config%verbosity > 1) then
      call header(output_unit)
   end if

   if (config%input == "-") then
      if (.not.allocated(config%input_format)) then
         call read_structure(mol, input_unit, filetype%xyz, error)
      else
         call read_structure(mol, input_unit, config%input_format, error)
      end if
   else
      call read_structure(mol, config%input, error, config%input_format)
   end if
   if (allocated(error)) return

   if (allocated(config%charge)) then
      mol%charge = config%charge
   else
      filename = join(dirname(config%input), ".CHRG")
      if (exists(filename)) then
         open(file=filename, newunit=unit)
         read(unit, *, iostat=stat) charge
         if (stat == 0) then
            mol%charge = charge
            if (config%verbosity > 0) write(output_unit, '(a)') &
               "[Info] Molecular charge read from '"//filename//"'"
         else
            if (config%verbosity > 0) write(output_unit, '(a)') &
               "[Warn] Could not read molecular charge read from '"//filename//"'"
         end if
         close(unit)
      end if
   end if

   if (config%wrap) then
      call wrap_to_central_cell(mol%xyz, mol%lattice, mol%periodic)
   end if

   ga = config%ga
   gc = config%gc
   if (config%mbdscale) s9 = config%inp%s9
   if (config%rational) then
      if (config%has_param) then
         param = config%inp
      else
         is = index(config%method, '/')
         if (is == 0) is = len_trim(config%method) + 1
         functional = lowercase(config%method(:is-1))
         id = get_functional_id(functional)

         ! special case: r2SCAN-3c (modifies s9, ga, gc)
         if (id == p_r2scan_3c) then
            if (.not.config%mbdscale) then
               s9 = 2.0_wp
            end if
            if (.not.config%zeta) then
               ga = 2.0_wp
               gc = 1.0_wp
            end if
         end if

         call get_rational_damping(functional, param, s9)
         if (.not.allocated(param)) then
            call fatal_error(error, "No parameters for '"//config%method//"' available")
            return
         end if
      end if
   end if

   if (allocated(param) .and. config%verbosity > 0) then
      call ascii_damping_param(output_unit, param, config%method)
   end if

   if (allocated(param)) then
      energy = 0.0_wp
      if (config%grad) then
         allocate(gradient(3, mol%nat), sigma(3, 3))
      end if
      if (config%hessian) then
         allocate(hessian(3, mol%nat, 3, mol%nat))
      end if
   end if

   call new_d4_model(error, d4, mol, ga=ga, gc=gc, wf=config%wf)
   if (allocated(error)) return

   if (config%properties) then
      if (config%verbosity > 1) then
         call ascii_atomic_radii(output_unit, mol, d4)
         call ascii_atomic_references(output_unit, mol, d4)
      end if
      allocate(cn(mol%nat), q(mol%nat), c6(mol%nat, mol%nat), alpha(mol%nat))
      call get_properties(mol, d4, realspace_cutoff(), cn, q, c6, alpha)

      if (config%verbosity > 0) then
         call ascii_system_properties(output_unit, mol, d4, cn, q, c6)
      end if
   end if

   if (allocated(param)) then
      call get_dispersion(mol, d4, param, realspace_cutoff(), energy, gradient, &
         & sigma)
      if (config%pair_resolved) then
         allocate(pair_disp2(mol%nat, mol%nat), pair_disp3(mol%nat, mol%nat))
         call get_pairwise_dispersion(mol, d4, param, realspace_cutoff(), pair_disp2, &
            & pair_disp3)
      end if
      if (config%hessian) then
         call get_dispersion_hessian(mol, d4, param, realspace_cutoff(), hessian)
      end if
      if (config%verbosity > 0) then
         call ascii_results(output_unit, mol, energy, gradient, sigma)
         if (config%pair_resolved) then
            call ascii_pairwise(output_unit, mol, pair_disp2, pair_disp3)
         end if
      end if
      if (config%tmer) then
         if (config%verbosity > 0) then
            write(output_unit, '(a)') "[Info] Dispersion energy written to .EDISP"
         end if
         open(file=".EDISP", newunit=unit)
         write(unit, '(f24.14)') energy
         close(unit)
      end if
      if (config%grad) then
         open(file=config%grad_output, newunit=unit)
         call tagged_result(unit, energy, gradient, sigma, hessian)
         close(unit)
         if (config%verbosity > 0) then
            write(output_unit, '(a)') &
               & "[Info] Dispersion results written to '"//config%grad_output//"'"
         end if

         inquire(file="gradient", exist=exist)
         if (exist) then
            call turbomole_gradient(mol, "gradient", energy, gradient, stat)
            if (config%verbosity > 0) then
               if (stat == 0) then
                  write(output_unit, '(a)') &
                     & "[Info] Dispersion gradient added to Turbomole gradient file"
               else
                  write(output_unit, '(a)') &
                     & "[Warn] Could not add to Turbomole gradient file"
               end if
            end if
         end if
         inquire(file="gradlatt", exist=exist)
         if (exist) then
            call turbomole_gradlatt(mol, "gradlatt", energy, sigma, stat)
            if (config%verbosity > 0) then
               if (stat == 0) then
                  write(output_unit, '(a)') &
                     & "[Info] Dispersion virial added to Turbomole gradlatt file"
               else
                  write(output_unit, '(a)') &
                     & "[Warn] Could not add to Turbomole gradlatt file"
               end if
            end if
         end if
      end if

   end if

   if (config%json) then
      open(file=config%json_output, newunit=unit)
      call json_results(unit, "  ", energy=energy, gradient=gradient, sigma=sigma, &
         & hessian=hessian, &
         & cn=cn, q=q, c6=c6, alpha=alpha, &
         & pairwise_energy2=pair_disp2, pairwise_energy3=pair_disp3)
      close(unit)
      if (config%verbosity > 0) then
         write(output_unit, '(a)') &
            & "[Info] JSON dump of results written to '"//config%json_output//"'"
      end if
   end if

end subroutine run_main


subroutine run_param(config, error)

   !> Configuration for this driver
   type(param_config), intent(in) :: config

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (config%list) then
      block
         type(functional_group), allocatable :: funcs(:)
         character(len=:), allocatable :: temp_names(:)
         integer, parameter :: MAX_LEN = 20
    
         integer :: i, j, nfuncs
         integer :: size_j, size_jp1
 
         call get_functionals(funcs)
         nfuncs = size(funcs)

         ! Bubble sort based on the first name in each group of funcs
         do i = 1, nfuncs - 1
            do j = 1, nfuncs - i
               if (funcs(j)%names(1) > funcs(j+1)%names(1)) then
                  size_j = size(funcs(j)%names)
                  size_jp1 = size(funcs(j+1)%names)

                  allocate(character(len=MAX_LEN) :: temp_names(size_j))
                  temp_names = funcs(j)%names

                  ! De- and reallocate before swap
                  if (allocated(funcs(j)%names)) deallocate(funcs(j)%names)
                  allocate(character(len=MAX_LEN) :: funcs(j)%names(size_jp1))
                  funcs(j)%names = funcs(j+1)%names

                  if (allocated(funcs(j+1)%names)) deallocate(funcs(j+1)%names)
                  allocate(character(len=MAX_LEN) :: funcs(j+1)%names(size_j))
                  funcs(j+1)%names = temp_names

                  deallocate(temp_names)
               end if
            end do
         end do
         
         write(output_unit, '(a)') "List of available functionals:"
         
         do i = 1, nfuncs
            associate(names => funcs(i)%names)
               do j = 1, size(names)
                  if (len_trim(names(j)) > 0) then
                     write(output_unit, '(a)', advance='no') trim(names(j)) // " "
                     
                     ! new line if last in list
                     if (size(names) == j) then
                        write(output_unit, *)
                     end if
                  end if
               end do
            end associate
         end do

      end block
   end if

end subroutine run_param


!> Construct path by joining strings with os file separator
function join(a1, a2) result(path)
   use mctc_env_system, only : is_windows
   character(len=*), intent(in) :: a1, a2
   character(len=:), allocatable :: path
   character :: filesep

   if (is_windows()) then
      filesep = '\'
   else
      filesep = '/'
   end if

   path = a1 // filesep // a2
end function join


!> test if pathname already exists
function exists(filename)
    character(len=*), intent(in) :: filename
    logical :: exists
    inquire(file=filename, exist=exists)
end function exists


!> Extract dirname from path
function dirname(filename)
   character(len=*), intent(in) :: filename
   character(len=:), allocatable :: dirname

   dirname = filename(1:scan(filename, "/\", back=.true.))
   if (len_trim(dirname) == 0) dirname = "."
end function dirname


end module dftd4_driver
