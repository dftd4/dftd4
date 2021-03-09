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

program main
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use mctc_env, only : error_type, fatal_error, get_argument, wp
   use mctc_io, only : structure_type, read_structure, filetype, get_filetype
   use dftd4, only : get_dispersion, d4_model, new_d4_model, &
      realspace_cutoff, get_lattice_points, get_coordination_number, &
      damping_param, rational_damping_param, get_rational_damping, &
      get_dftd4_version, get_pairwise_dispersion
   use dftd4_charge, only : get_charges
   use dftd4_output
   use dftd4_utils
   implicit none
   character(len=*), parameter :: prog_name = "dftd4"

   type :: d4_config
      logical :: json = .false.
      character(len=:), allocatable :: json_output
      logical :: wrap = .true.
      logical :: tmer = .true.
      logical :: properties = .false.
      logical :: mbdscale = .false.
      logical :: grad = .false.
      character(len=:), allocatable :: grad_output
      logical :: rational = .false.
      logical :: has_param = .false.
      integer :: verbosity = 2
      real(wp), allocatable :: charge
      real(wp) :: ga = 3.0_wp
      real(wp) :: gc = 2.0_wp
      real(wp) :: wf = 6.0_wp
      logical :: pair_resolved = .false.
   end type d4_config
   type(d4_config) :: config

   character(len=:), allocatable :: input
   integer, allocatable :: input_format
   type(structure_type) :: mol
   class(damping_param), allocatable :: param
   type(rational_damping_param) :: inp
   type(d4_model) :: d4
   type(error_type), allocatable :: error
   real(wp) :: energy, sigma(3, 3), charge
   real(wp), allocatable :: gradient(:, :), pair_disp2(:, :), pair_disp3(:, :)
   character(len=:), allocatable :: method
   real(wp), allocatable :: s9
   integer :: stat, unit
   logical :: exist

   call get_arguments(input, input_format, config, method, inp, error)
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   if (config%verbosity > 1) then
      call header(output_unit)
   end if

   if (input == "-") then
      if (.not.allocated(input_format)) input_format = filetype%xyz
      call read_structure(mol, input_unit, input_format, error)
   else
      call read_structure(mol, input, error, input_format)
   end if
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if
   if (allocated(config%charge)) then
      mol%charge = config%charge
   else
      inquire(file='.CHRG', exist=exist)
      if (exist) then
         open(file='.CHRG', newunit=unit)
         read(unit, *, iostat=stat) charge
         if (stat == 0) then
            mol%charge = charge
            if (config%verbosity > 0) write(output_unit, '(a)') &
               "[Info] Molecular charge read from .CHRG"
         else
            if (config%verbosity > 0) write(output_unit, '(a)') &
               "[Warn] Could not read molecular charge read from .CHRG"
         end if
         close(unit)
      end if
   end if
   if (config%wrap) then
      call wrap_to_central_cell(mol%xyz, mol%lattice, mol%periodic)
   end if

   if (config%mbdscale) s9 = inp%s9
   if (config%rational) then
      if (config%has_param) then
         param = inp
      else
         call get_rational_damping(method, param, s9)
         if (.not.allocated(param)) then
            write(error_unit, '("[Error]", 1x, a)') &
               "No parameters for '"//method//"' available"
            error stop
         end if
      end if
   end if

   if (allocated(param) .and. config%verbosity > 0) then
      call ascii_damping_param(output_unit, param, method)
   end if

   if (config%grad) then
      allocate(gradient(3, mol%nat))
   end if

   call new_d4_model(d4, mol, ga=config%ga, gc=config%gc, wf=config%wf)

   if (config%properties) then
      call property_calc(output_unit, mol, d4, config%verbosity)
   end if

   if (allocated(param)) then
      call get_dispersion(mol, d4, param, realspace_cutoff(), energy, gradient, &
         & sigma)
      if (config%pair_resolved) then
         allocate(pair_disp2(mol%nat, mol%nat), pair_disp3(mol%nat, mol%nat))
         call get_pairwise_dispersion(mol, d4, param, realspace_cutoff(), pair_disp2, &
            & pair_disp3)
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
         call tagged_result(unit, energy, gradient, sigma)
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

      if (config%json) then
         open(file=config%json_output, newunit=unit)
         if (config%grad) then
            call json_results(unit, "  ", energy, gradient, sigma, &
               & pairwise_energy2=pair_disp2, pairwise_energy3=pair_disp3)
         else
            call json_results(unit, "  ", energy, &
               & pairwise_energy2=pair_disp2, pairwise_energy3=pair_disp3)
         end if
         close(unit)
         if (config%verbosity > 0) then
            write(output_unit, '(a)') &
               & "[Info] JSON dump of results written to '"//config%json_output//"'"
         end if
      end if
   end if


contains


subroutine property_calc(unit, mol, disp, verbosity)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d4_model), intent(in) :: disp

   !> Printout verbosity
   integer, intent(in) :: verbosity

   integer :: mref
   real(wp), allocatable :: cn(:), q(:), gwvec(:, :), c6(:, :), lattr(:, :)

   if (verbosity > 1) then
      call ascii_atomic_radii(unit, mol, disp)
      call ascii_atomic_references(unit, mol, disp)
   end if

   mref = maxval(disp%ref)
   allocate(cn(mol%nat), q(mol%nat), gwvec(mref, mol%nat), c6(mol%nat, mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, 30.0_wp, lattr)
   call get_coordination_number(mol, lattr, 30.0_wp, disp%rcov, disp%en, cn)
   call get_charges(mol, q)
   call disp%weight_references(mol, cn, q, gwvec)
   call disp%get_atomic_c6(mol, gwvec, c6=c6)

   if (verbosity > 0) then
      call ascii_system_properties(unit, mol, disp, cn, q, c6)
   end if

end subroutine property_calc


subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options] <input>"

   write(unit, '(a)') &
      "", &
      "Generally Applicable Atomic-Charge Dependent London Dispersion Correction.", &
      "Takes an geometry input to calculate the D4 dispersion correction.", &
      "Periodic calculations are performed automatically for periodic input formats.", &
      "Specify the functional to select the correct parameters.", &
      ""

   write(unit, '(2x, a, t25, a)') &
      "-c, --charge <real>", "Set charge to molecule, overwrites .CHRG file", &
      "-i, --input <format>", "Hint for the format of the input file", &
      "-f, --func <method>", "Use damping parameters for given functional", &
      "    --param <list>", "Specify parameters for rational damping,", &
      "", "expected order is s6, s8, a1, a2 (requires four arguments)", &
      "    --mbdscale <s9>", "Use scaled ATM three-body dispersion", &
      "    --zeta <list>", "Adjust charge scaling parameters, takes two reals", &
      "", "expected order is ga, gc (default: 3.0, 2.0)", &
      "    --wfactor <real>", "Adjust weighting factor for interpolation", &
      "", "(default: 6.0)", &
      "-g, --grad", "Evaluate molecular gradient and virial", &
      "", "write results to file (default: dftd4.txt),", &
      "", "attempts to add to Turbomole gradient and gradlatt files", &
      "    --property", "Show dispersion related atomic and system properties", &
      "    --pair-resolved", "Calculate pairwise representation of dispersion energy", &
      "    --noedisp", "Disable writing of dispersion energy to .EDISP file", &
      "    --json [file]", "Dump results to JSON output (default: dftd4.json)", &
      "    --grad [file]", "Request gradient evaluation,", &
      "-v, --verbose", "Show more, can be used multiple times", &
      "-s, --silent", "Show less, use twice to supress all output", &
      "    --version", "Print program version and exit", &
      "    --citation", "Print citation information and exit", &
      "    --license", "Print license header and exit", &
      "-h, --help", "Show this help message"

   write(unit, '(a)')

end subroutine help

subroutine header(unit)
   integer, intent(in) :: unit

   write(unit,'(a)') &
      "                    ____  _____ _____     ____  _  _",&
      "      -------------|  _ \|  ___|_   _|---|  _ \| || |------------",&
      "     |             | | | | |_    | | ___ | | | | || |_           |",&
      "     |             | |_| |  _|   | ||___|| |_| |__   _|          |",&
      "     |             |____/|_|     |_|     |____/   |_|            |",&
      "     |             ===================================           |",&
      "     |            E. Caldeweyher, S. Ehlert & S. Grimme          |",&
      "     |          Mulliken Center for Theoretical Chemistry        |",&
      "     |                    University of Bonn                     |",&
      "      ----------------------------------------------------------- ",""
end subroutine header


subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_dftd4_version(string=version_string)
   write(unit, '(a, *(1x, a))') &
      & prog_name, "version", version_string

end subroutine version


subroutine citation(unit)
   integer, intent(in) :: unit

   write(unit, '(a)') &
      "Please include the appropriate citations when using DFTD4 in your work.", &
      "", &
      "Original DFTD4 idea:", &
      "Eike Caldeweyher, Christoph Bannwarth and Stefan Grimme,", &
      "J. Chem. Phys., 2017, 147, 034112.", &
      "DOI: 10.1063/1.4993215", &
      "", &
      "DFTD4 model:", &
      "Eike Caldeweyher, Sebastian Ehlert, Andreas Hansen, Hagen Neugebauer,", &
      "Sebastian Spicher, Christoph Bannwarth and Stefan Grimme,", &
      "J. Chem Phys, 2019, 150, 154122.", &
      "DOI: 10.1063/1.5090222", &
      "ChemRxiv: 10.26434/chemrxiv.7430216.v2", &
      "", &
      "Periodic DFTD4 model:", &
      "Eike Caldeweyher, Jan-Michael Mewes, Sebastian Ehlert", &
      "and Stefan Grimme, Phys. Chem. Chem. Phys., 2020, 22, 8499-8512.", &
      "DOI: 10.1039/D0CP00502A", &
      "ChemRxiv: 10.26434/chemrxiv.10299428.v1", &
      ""

end subroutine citation


subroutine license(unit)
   integer, intent(in) :: unit

   write(unit, '(a)') &
      "dftd4 is free software: you can redistribute it and/or modify it under", &
      "the terms of the Lesser GNU General Public License as published by", &
      "the Free Software Foundation, either version 3 of the License, or", &
      "(at your option) any later version.", &
      "", &
      "dftd4 is distributed in the hope that it will be useful,", &
      "but WITHOUT ANY WARRANTY; without even the implied warranty of", &
      "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the", &
      "Lesser GNU General Public License for more details.", &
      "", &
      "You should have received a copy of the Lesser GNU General Public License", &
      "along with dftd4.  If not, see <https://www.gnu.org/licenses/>."
end subroutine license


subroutine get_argument_as_real(iarg, val, error)

   !> Index of command line argument, range [0:command_argument_count()]
   integer, intent(in) :: iarg

   !> Real value
   real(wp), intent(out) :: val

   !> Error handling
   type(error_type), allocatable :: error

   integer :: stat
   character(len=:), allocatable :: arg

   call get_argument(iarg, arg)
   if (.not.allocated(arg)) then
      call fatal_error(error, "Cannot read real value, argument missing")
      return
   end if
   read(arg, *, iostat=stat) val
   if (stat /= 0) then
      call fatal_error(error, "Cannot read real value from '"//arg//"'")
      return
   end if

end subroutine get_argument_as_real


subroutine get_arguments(input, input_format, config, method, inp, error)

   !> Input file name
   character(len=:), allocatable :: input

   !> Input file format
   integer, allocatable, intent(out) :: input_format

   !> Configuation data
   type(d4_config), intent(out) :: config

   !> Method name
   character(len=:), allocatable, intent(out) :: method

   !> Damping parameters
   type(rational_damping_param), intent(out) :: inp

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   character(len=:), allocatable :: arg

   iarg = 0
   narg = command_argument_count()
   do while(iarg < narg)
      iarg = iarg + 1
      call get_argument(iarg, arg)
      select case(arg)
      case("-h", "--help")
         call help(output_unit)
         stop
      case("--version")
         call version(output_unit)
         stop
      case("--citation")
         call citation(output_unit)
         stop
      case("--license")
         call license(output_unit)
         stop
      case("-v", "--verbose")
         config%verbosity = config%verbosity + 1
      case("-s", "--silent")
         config%verbosity = config%verbosity - 1
      case default
         if (.not.allocated(input)) then
            call move_alloc(arg, input)
            cycle
         end if
         if (arg(1:1) == "-") then
            call fatal_error(error, "Unknown command option '"//arg//"'")
         else
            call fatal_error(error, "Too many positional arguments present")
         end if
         exit
      case("-i", "--input")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         input_format = get_filetype("."//arg)
      case("-c", "--charge")
         iarg = iarg + 1
         allocate(config%charge)
         call get_argument_as_real(iarg, config%charge, error)
         if (allocated(error)) exit
      case("--json")
         config%json = .true.
         config%json_output = "dftd4.json"
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%json_output)
         end if
      case("--property")
         config%properties = .true.
      case("--pair-resolved")
         config%pair_resolved = .true.
      case("--noedisp")
         config%tmer = .false.
      case("--nowrap")
         config%wrap = .false.
      case("--grad")
         config%grad = .true.
         config%grad_output = "dftd4.txt"
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%grad_output)
         end if
      case("--mbdscale")
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%s9, error)
         if (allocated(error)) exit
         config%mbdscale = .true.
      case("--zeta")
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%ga, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%gc, error)
         if (allocated(error)) exit
      case("--wfactor")
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%wf, error)
         if (allocated(error)) exit
      case("-f", "--func")
         config%rational = .true.
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
         call move_alloc(arg, method)
      case("--param")
         config%rational = .true.
         config%has_param = .true.
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%s6, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%s8, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%a1, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%a2, error)
         if (allocated(error)) exit
      end select
   end do
   if (allocated(error)) return

   if (.not.config%has_param .and. .not.allocated(method)) then
      config%properties = .true.
   end if

   if (.not.allocated(input)) then
      if (.not.allocated(error)) then
         call help(output_unit)
         error stop
      end if
   end if

end subroutine get_arguments

end program main
