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

module dftd4_help
   use dftd4, only : get_dftd4_version
   implicit none
   private

   public :: citation, license, header, prog_name, version
   public :: help_text, help_text_run, help_text_param


   !> The name of the program
   character(len=*), parameter :: prog_name = "dftd4"

   character(len=*), parameter :: nl = new_line('a')

   character(len=*), parameter :: run_options_text = &
      "-c,--charge <real>       Set charge to molecule, overwrites .CHRG file"//nl//&
      "-i,--input <format>      Hint for the format of the input file"//nl//&
      "-f,--func <method>       Use damping parameters for given functional"//nl//&
      "   --param <list>        Specify parameters for rational damping,"//nl//&
      "                         expected order is s6, s8, a1, a2 (requires four arguments)"//nl//&
      "   --mbdscale <s9>       Use scaled ATM three-body dispersion"//nl//&
      "   --zeta <list>         Adjust charge scaling parameters, takes two reals,"//nl//&
      "                         expected order is ga, gc (default: 3.0, 2.0)"//nl//&
      "   --wfactor <real>      Adjust weighting factor for interpolation"//nl//&
      "                         (default: 6.0)"//nl//&
      "-g,--grad [file]         Evaluate molecular gradient and virial,"//nl//&
      "                         write results to file (default: dftd4.txt),"//nl//&
      "                         attempts to add to Turbomole gradient and gradlatt files"//nl//&
      "   --hessian             Evaluate molecular hessian"//nl//&
      "   --property            Show dispersion related atomic and system properties"//nl//&
      "   --pair-resolved       Calculate pairwise representation of dispersion energy"//nl//&
      "   --noedisp             Disable writing of dispersion energy to .EDISP file"//nl//&
      "   --json [file]         Dump results to JSON output (default: dftd4.json)"//nl//&
      "-v,--verbose             Show more, can be used multiple times"//nl//&
      "-s,--silent              Show less, use twice to supress all output"//nl//&
      "   --version             Print program version and exit"//nl//&
      "   --citation            Print citation information and exit"//nl//&
      "   --license             Print license header and exit"//nl//&
      "-h,--help                Show this help message"

   character(len=*), parameter :: help_text_run = &
      "Usage: "//prog_name//" [run] [options] <input>"//nl//&
      ""//nl//&
      "Takes an geometry input to calculate the D4 dispersion correction."//nl//&
      "Periodic calculations are performed automatically for periodic input formats."//nl//&
      "Reads .CHRG file (if present) from the same directory as the input."//nl//&     
      "Specify the functional to select the correct parameters."//nl//&
      ""//nl//&
      run_options_text//nl//&
      ""

   character(len=*), parameter :: param_options_text = &
      "-l,--list,--funcs         List all supported functionals"

   character(len=*), parameter :: help_text_param = &
      "Usage: "//prog_name//" param [options]"//nl//&
      ""//nl//&
      "Inspect damping parameters and supported functionals"//nl//&
      ""//nl//&
      param_options_text//nl//&
      ""

   character(len=*), parameter :: help_text = &
      "Usage: "//prog_name//" [run|param] [options] ..."//nl//&
      ""//nl//&
      !
      "Generally Applicable Atomic-Charge Dependent London Dispersion Correction."//nl//&
      "Takes an geometry input to calculate the D4 dispersion correction."//nl//&
      "Periodic calculations are performed automatically for periodic input formats."//nl//&
      "Reads .CHRG file (if present) from the same directory as the input."//nl//&     
      "Specify the functional to select the correct parameters."//nl//&
      ""//nl//&
      !
      "Commands"//nl//&
      ""//nl//&
      "  run       Evaluate dispersion correction on the provided input structure."//nl//&
      "            Periodic calculations are performed automatically for periodic inputs"//nl//&
      "            If no command is specified run is selected by default."//nl//&
      ""//nl//&
      "  param     Inspect damping parameters."//nl//&
      ""//nl//&
      !
      "Options"//nl//&
      ""//nl//&
      run_options_text//nl//&
      ""

contains


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


end module dftd4_help
