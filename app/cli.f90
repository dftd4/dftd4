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

!> Definition of the command line interface
module dftd4_cli
   use, intrinsic :: iso_fortran_env, only : output_unit
   use mctc_env, only : error_type, fatal_error, get_argument, wp
   use mctc_io, only : get_filetype
   use dftd4, only : rational_damping_param, get_dftd4_version
   implicit none
   private

   public :: cli_config, run_config, get_arguments
   public :: prog_name, header

   !> The name of the program
   character(len=*), parameter :: prog_name = "dftd4"

   !> Base command line configuration
   type, abstract :: cli_config
   end type cli_config

   !> Configuration data for running single point calculations
   type, extends(cli_config) :: run_config
      character(len=:), allocatable :: input
      integer, allocatable :: input_format
      type(rational_damping_param) :: inp
      character(len=:), allocatable :: method
      logical :: json = .false.
      character(len=:), allocatable :: json_output
      logical :: wrap = .true.
      logical :: tmer = .true.
      logical :: properties = .false.
      logical :: mbdscale = .false.
      logical :: grad = .false.
      logical :: hessian = .false.
      character(len=:), allocatable :: grad_output
      logical :: rational = .false.
      logical :: has_param = .false.
      integer :: verbosity = 2
      real(wp), allocatable :: charge
      real(wp) :: ga = 3.0_wp
      real(wp) :: gc = 2.0_wp
      real(wp) :: wf = 6.0_wp
      logical :: pair_resolved = .false.
   end type run_config


contains


subroutine get_arguments(config, error)

   !> Configuation data
   class(cli_config), allocatable, intent(out) :: config

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   block
      type(run_config), allocatable :: tmp
      allocate(tmp)
      call get_run_arguments(tmp, error)
      call move_alloc(tmp, config)
   end block
end subroutine get_arguments


!> Read configuration for the single point driver
subroutine get_run_arguments(config, error)

   !> Configuation data
   type(run_config), intent(out) :: config

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   logical :: getopts
   character(len=:), allocatable :: arg

   iarg = 0
   getopts = .true.
   narg = command_argument_count()
   do while(iarg < narg)
      iarg = iarg + 1
      call get_argument(iarg, arg)
      if (.not.getopts) then
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      end if
      select case(arg)
      case("--")
         getopts = .false.
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
      case("--list-funcs", "--list-func")
         call list_funcs(output_unit)
         stop
      case("-v", "--verbose")
         config%verbosity = config%verbosity + 1
      case("-s", "--silent")
         config%verbosity = config%verbosity - 1
      case default
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
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
         config%input_format = get_filetype("."//arg)
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
      case("-g", "--grad")
         config%grad = .true.
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%grad_output)
         end if
      case("--hessian")
         config%grad = .true.
         config%hessian = .true.
      case("--mbdscale")
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%inp%s9, error)
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
         call move_alloc(arg, config%method)
      case("--param")
         config%rational = .true.
         config%has_param = .true.
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%inp%s6, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%inp%s8, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%inp%a1, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%inp%a2, error)
         if (allocated(error)) exit
      end select
   end do
   if (allocated(error)) return

   if (.not.config%has_param .and. .not.allocated(config%method)) then
      config%properties = .true.
   end if

   if (.not.allocated(config%grad_output)) then
      config%grad_output = "dftd4.txt"
   end if

   if (.not.allocated(config%input)) then
      if (.not.allocated(error)) then
         call help(output_unit)
         error stop
      end if
   end if

end subroutine get_run_arguments


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


subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options] <input>"

   write(unit, '(a)') &
      "", &
      "Generally Applicable Atomic-Charge Dependent London Dispersion Correction.", &
      "Takes an geometry input to calculate the D4 dispersion correction.", &
      "Periodic calculations are performed automatically for periodic input formats.", &
      "Reads .CHRG file (if present) from the same directory as the input.", &      
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
      "-g, --grad [file]", "Evaluate molecular gradient and virial", &
      "    --hessian", "Evaluate molecular hessian", &
      "", "write results to file (default: dftd4.txt),", &
      "", "attempts to add to Turbomole gradient and gradlatt files", &
      "    --property", "Show dispersion related atomic and system properties", &
      "    --pair-resolved", "Calculate pairwise representation of dispersion energy", &
      "    --noedisp", "Disable writing of dispersion energy to .EDISP file", &
      "    --json [file]", "Dump results to JSON output (default: dftd4.json)", &
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


subroutine list_funcs(unit)
   integer, intent(in) :: unit

   character(len=*), parameter :: func(*) = [character(len=32)::&
      & 'am05', 'b1-lyp', 'b1-p', 'b1-pw', 'b1pw', 'b1pw91', 'b1p', 'b1p86', &
      & 'b1lyp', 'b2-plyp', 'b2gp-plyp', 'b2plyp', 'b3-lyp', 'b3-p', 'b3-pw', &
      & 'b3lyp', 'b3p', 'b3p86', 'b3pw', 'b3pw91', 'b6k', 'bh-lyp', 'bhpbe', &
      & 'bmk', 'bp', 'bp86', 'bpbe', 'bpw', 'bwb6k', 'cam-b3lyp', 'dftb(3ob)', &
      & 'dftb(matsci)', 'dftb(mio)', 'dftb(ob2)', 'dftb(pbc)', 'dsd-blyp', &
      & 'dsd-pbe', 'dsd-pbep86', 'dsd-pbeb95', 'dsd-svwn', 'dodblyp', &
      & 'dod-pbe', 'dod-pbep86', 'dod-pbeb95', 'dod-svwn', 'dsdblyp', &
      & 'dsdpbe', 'dsdpbep86', 'dsdpbeb95', 'dsdsvwn', 'g-lyp', 'glyp', &
      & 'hcth120', 'hcth407', 'hcthhyb', 'hf', 'hf-3c', 'hf-3cv', 'hf3c', &
      & 'hf3cv', 'hiss', 'hse03', 'hse06', 'hse12', 'hse12s', 'hsesol', &
      & 'lb94', 'lc-blyp', 'lc-dftb', 'lc-wpbe', 'lc-whpbe', 'lh07ssvwn', &
      & 'lh07s-svwn', 'lh07tsvwn', 'lh07t-svwn', 'lh12ctssifpw92', &
      & 'lh12ct-ssifpw92', 'lh12ctssirpw92', 'lh12ct-ssirpw92', 'lh14tcalpbe', &
      & 'lh14t-calpbe', 'lh20t', 'lcwpbe', 'm05', 'm052x', 'm05-2x', 'm06', &
      & 'm06-2x', 'm06l', 'm08hx', 'm08-hx', 'm11l', 'mn12l', 'mn12sx', &
      & 'mn15', 'mn15l', 'mpw-lyp', 'mpw1-lyp', 'mpw1-pw', 'mpw1-pw91', &
      & 'mpw1-pw', 'mpw1lyp', 'mpw1pw', 'mpw1pw91', 'mpw2plyp', 'mpwkcis1k', &
      & 'mpwlyp', 'mpw-pw', 'mpwpw', 'mpwpw91', 'mpw1kcis', 'mpwb1k', 'n12', &
      & 'n12sx', 'o-lyp', 'o3-lyp', 'o3lyp', 'olyp', 'omegab97', 'omegab97m', &
      & 'opbe', 'otpss', 'pbe', 'pbe0', 'pbe02', 'pbe0-2', 'pbe0dh', &
      & 'pbe0-dh', 'pbexalpha', 'pbeh1pbe', 'pbeh-3c', 'pbep86', 'pbepbeb95', &
      & 'pbepbe', 'pbesol', 'pwp', 'pwp1', 'pw1-pw', 'pw1pw', 'pw1pw91', &
      & 'pw6b95', 'pw86pbe', 'pw91', 'pw91p86', 'pwgga', 'revdsd-blyp', &
      & 'revdsdblyp', 'revdsd-pbe', 'revdsdpbe', 'revdsd-pbepbe', &
      & 'revdsdpbepbe', 'revdsdpbep86', 'revdsdpbeb95', 'revdsd-svwn', &
      & 'revdod-pbep86', 'revdodpbep86', 'revdod-pbeb95', 'revdodpbepb95', &
      & 'revdod-svwn', 'revpbe', 'revpbe0', 'revpbe0-dh', 'revpbe38', &
      & 'revtpss', 'revtpss0', 'r2scan', 'r²scan', 'r2scan0', 'r²scan0', &
      & 'r2scan50', 'r²scan50', 'r2scanh', 'r²scanh', 'rscan', 'scan', 'ssb', &
      & 'tauhctc', 'tauhcthhyb', 'thcth', 'tpss', 'tpss0', 'tpssh', 'wb97', &
      & 'wb97m', 'ωb97', 'ωb97m', 'x-lyp', 'x3-lyp', 'x3lyp']

   character(len=1) :: current_letter = ''
   integer :: i
   
   write(unit, '(a)') &
      & 'The following functionals are available (keywords for --func):'
  
   do i = 1, size(func)
      if (func(i)(1:1) /= current_letter) then
         write(unit, '(a)') ''
         current_letter = func(i)(1:1)
      end if
        
      write(unit, '(a)') func(i)
   end do
end subroutine list_funcs

end module dftd4_cli
