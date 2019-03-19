subroutine dftd4_header(verbose)
use iso_fortran_env, istdout => output_unit
logical,intent(in) :: verbose
if (verbose) then
write(istdout,'(a)') &
   !< < < < < < < < < < < < < < < < < < > > > > > > > > > > > > > > > > > >!
   ! okay... this is ridiculous...
!  "         _______   _________ _________     _______     /\              ",&
!  "        '_   __ `.|_   ___  |  _   _  |   '_   __ `.  / /  _           ",&
!  "      ----| |--`. \-| |---\_|_/-| |-\_|-----| |--`. \/ /--| |----      ",&
!  "     |    | |   | | | '--.      | |   _____ | |   | / /___' '_   |     ",&
!  "     |    | |   | | | .--'      | |  |_____|| |   | '_____.  _'  |     ",&
!  "     |   _| |__.' /_| |_       _| |_       _| |__.' /    _| |_   |     ",&
!  "     | =|_______.''_____'====='_____'=====|_______.'===='_____'= |     ",&
   !< < < < < < < < < < < < < < < < < < > > > > > > > > > > > > > > > > > >!
   "                    ____  _____ _____     ____  _  _                   ",&
   "      -------------|  _ \|  ___|_   _|---|  _ \| || |------------      ",&
   "     |             | | | | |_    | | ___ | | | | || |_           |     ",&
   "     |             | |_| |  _|   | ||___|| |_| |__   _|          |     ",&
   "     |             |____/|_|     |_|     |____/   |_|            |     ",&
   "     |             ===================================           |     "
   !< < < < < < < < < < < < < < < < < < > > > > > > > > > > > > > > > > > >!
else
write(istdout,'(a)') &
   !< < < < < < < < < < < < < < < < < < > > > > > > > > > > > > > > > > > >!
   "      -----------------------------------------------------------      ",&
   "     |                   =====================                   |     ",&
   "     |                        D F T - D 4                        |     ",&
   "     |                   =====================                   |     "
   !< < < < < < < < < < < < < < < < < < > > > > > > > > > > > > > > > > > >!
endif
write(istdout,'(a)') &
   !< < < < < < < < < < < < < < < < < < > > > > > > > > > > > > > > > > > >!
   "     |            E. Caldeweyher, S. Ehlert & S. Grimme          |     ",&
   "     |          Mulliken Center for Theoretical Chemistry        |     ",&
   "     |                    University of Bonn                     |     ",&
   "     |                  Version 2.0 (SAW190211)                  |     ",&
   !     |  Version number by <major>.<minor>.<rev> (<author><date>) |     !
   "      -----------------------------------------------------------      ",""
   !< < < < < < < < < < < < < < < < < < > > > > > > > > > > > > > > > > > >!
end subroutine dftd4_header

subroutine eeq_header
use iso_fortran_env, istdout => output_unit
write(istdout,'(10x,a)') &
   !< < < < < < < < < < < < < > > > > > > > > > > > > >!
   " ------------------------------------------------- ",&
   "|                      E E Q                      |",&
   "|       Electronegativity Equilibrium Model       |",&
   " ------------------------------------------------- "
   !< < < < < < < < < < < < < > > > > > > > > > > > > >!
end subroutine eeq_header

subroutine dftd4_citation
use iso_fortran_env, istdout => output_unit
write(istdout,'(3x,a)') &
   "Please cite:", &
   "E. Caldeweyher, C. Bannwarth and S. Grimme, J. Chem. Phys., 2017,", &
   "147, 034112.", &
   "and",&
   "E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher,", &
   "C. Bannwarth and S. Grimme, ChemRxiv, 2018, preprint.",&
   "http://doi.org/10.26434/chemrxiv.7430216.v2",&
   "",&
   "For GFN2-xTB:", &
   "C. Bannwarth, S. Ehlert and S. Grimme., ChemRxiv, 2018, preprint.",&
   "http://doi.org/10.26434/chemrxiv.7246238.v2",&
   "",&
   "For a general overview on dispersion corrected mean-field methods",&
   "we recommend:",&
   "S. Grimme, A. Hansen, J. G. Brandenburg, C. Bannwarth, Chem. Rev. 2016,",&
   "116, 5105−5154.",&
   "",&
   "with help from (in alphabetical order)",&
   "C. Bannwarth, P. Shushkov, and S. Spicher.",&
   ""
end subroutine dftd4_citation

subroutine gpl_license
   use iso_fortran_env, istdout => output_unit
   write(istdout,'(3x,a)') &
      "Copyright (C) 2017-2019 S. Grimme",&
      "",&
      "This program is free software: you can redistribute it and/or ",&
      "modify it under the terms of the GNU General Public License as ",&
      "published by the Free Software Foundation, either version 3 of ",&
      "the License, or (at your option) any later version.",&
      "",&
      "This program is distributed in the hope that it will be useful,",&
      "but WITHOUT ANY WARRANTY; without even the implied warranty of",&
      "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the",&
      "GNU General Public License for more details.",&
      "",&
      "You should have received a copy of the GNU General Public License",&
      "along with this program.",&
      "If not, see <http://www.gnu.org/licenses/>.",&
      ""
end subroutine gpl_license

subroutine help
   use iso_fortran_env, istdout => output_unit
   write(istdout,'(a)') &
      "Usage:",&
      "dftd4 [options] <file>",&
      "",&
      "<file> is a valid Turbomole coordinate file (coordinates in Bohr) or",&
      "in xmol format (coordinates in Ångström).",&
      "",&
      "Options:",&
      "",&
      "-c, --chrg <integer>",&
      "        set charge of the System, (.CHRG can also used as input)",&
      "-f, --func <name>[/<basis>]",&
      "        calculate DFT-D4 dispersion for the given functional",&
      "    --zeta <real> <real>",&
      "        use given parameters for the charge scaling (default: 3.0 & 2.0)",&
      "    --wfactor <real>",&
      "        use <real> as weighting factor for Gaussian weighting",&
      "        (default: 6.0)",&
      "    --param <s6> <s8> <a1> <a2>",&
      "        use userdefined damping parameters",&
      "    --mbdscale <real>",&
      "        scale non-additive dispersion by <real> (default: 1.0)",&
      "    --s10 <real>",&
      "        scale C10 coeffients by <real> (default: 0.0)",&
      "-2, --two",&
      "        only use pairwise two-body dispersion energy",&
      "-3, --abc",&
      "        use Axilrod-Teller-Muto three-body dispersion energy (default)",&
      "-m, --mbd",&
      "        use RPA-like many-body dispersion energy as non-additivity",&
      "        correction (no gradient and hessian)",&
      "-g, --grad",&
      "        calculate analytical first derivates",&
      "    --hess",&
      "        calculate numerical second derivates",&
      "    --tmer",&
      "        force tmer2++ compatible .EDISP printout",&
      "    --orca",&
      "        ORCA compatibility mode (for usage by the orca binary)",&
      "    --molc6",&
      "        calculate dispersion related properties (default)",&
   !$ "-P, --parallel <integer>",&
   !$ "        use <integer> OMP threads in parallel sections",&
      "-v, --verbose",&
      "        be verbose",&
      "-s, --silent",&
      "        clutter the screen less",&
      "-h, --help",&
      "        print this message",&
      ""
end subroutine help

subroutine generic_header(iunit,string,width,offset)
implicit none
integer,intent(in) :: iunit
integer,intent(in) :: offset
integer,intent(in) :: width
character(len=*),intent(in) :: string
character(len=width) :: dum1,dum2
character(len=2*width) :: outstring
character(len=width) :: formatstr
integer :: strlen,ifront,iback
strlen = len(string)
ifront = (width - strlen)/2
iback  = width - ifront - strlen
write(dum1,*) width
write(dum2,*) offset
write(formatstr,'(i0,"x,a,",i0,"x")') ifront,iback
write(outstring,'("|",'//formatstr//',"|")') string
write(iunit,'('//dum2//'x,1x,'//dum1//'("-"),1x)')
write(iunit,'('//dum2//'x,a)') trim(outstring)
write(iunit,'('//dum2//'x,1x,'//dum1//'("-"),1x)')
end subroutine generic_header

subroutine print_pbcsum(iunit,mol)
   use iso_fortran_env, wp => real64
   use mctc_constants
   use mctc_econv
   use class_molecule
   implicit none
   integer, intent(in)  :: iunit
   type(molecule),intent(in) :: mol

   integer  :: i
   real(wp) :: conv

   conv = autoaa

   call generic_header(iunit,"Geometry Summary",49,10)
   write(iunit,'(a)')

   ! atomic coordinates
   write(iunit,'(1x,"*",1x,i0,1x,a)') mol%nat,"atoms in unit cell"
   write(iunit,'(a)')
   write(iunit,'(5x,"#",3x,"Z",3x,32x,"position/Å",8x,"charge")')
   do i = 1, mol%nat
      write(iunit,'(i6,1x,i3,1x,a2)',advance='no') i,mol%at(i),mol%sym(i)
      write(iunit,'(3f14.7)',advance='no') mol%xyz(:,i)*conv
      write(iunit,'(f14.7)') mol%z(i)
   enddo
   write(iunit,'(a)')

   ! periodicity
   write(iunit,'(1x,"*",1x,i0,a)') mol%npbc,"D periodic system"
   write(iunit,'(a)')

   if (mol%npbc > 0) then
      ! cell parameters
      write(iunit,'(1x,"*",1x,a)') "cell parameter"
      write(iunit,'(a)')
      write(iunit,'(a12,2a15,2x,3a11)') &
         "|a|/Å", "|b|/Å", "|c|/Å", "α/°", "β/°", "γ/°"
      write(iunit,'(f13.7,2f14.7,1x,3f9.3)') &
         mol%cellpar(:3)*conv,mol%cellpar(4:)*180.0_wp/pi
      write(iunit,'(a)')

      ! direct lattice (transformation abc -> xyz)
      write(iunit,'(1x,"*",1x,a)') "direct lattice/Å"
      write(iunit,'(a)')
      write(iunit,'(12x,a,3f14.7)') "a",mol%lattice(:,1)*conv
      write(iunit,'(12x,a,3f14.7)') "b",mol%lattice(:,2)*conv
      write(iunit,'(12x,a,3f14.7)') "c",mol%lattice(:,3)*conv
      write(iunit,'(a)')

      ! reciprocal lattice
      write(iunit,'(1x,"*",1x,a)') "reciprocal lattice/Å⁻¹"
      write(iunit,'(a)')
      write(iunit,'(11x,a,3f14.7)') "a*",mol%rec_lat(:,1)/conv
      write(iunit,'(11x,a,3f14.7)') "b*",mol%rec_lat(:,2)/conv
      write(iunit,'(11x,a,3f14.7)') "c*",mol%rec_lat(:,3)/conv
      write(iunit,'(a)')

      ! geometry in fractional coordinates
      write(iunit,'(1x,"*",1x,a)') "geometry in fractional coordinates"
      write(iunit,'(a)')
      write(iunit,'(5x,"#",3x,"Z",3x,20x,"fractional coordinates")')
      do i = 1, mol%nat
         write(iunit,'(i6,1x,i3,1x,a2)',advance='no') i,mol%at(i),mol%sym(i)
         write(iunit,'(3f14.7)',advance='no') mol%abc(:,i)
         write(iunit,'(a)')
      enddo
      write(iunit,'(a)')

      ! volume of unit cell
      write(iunit,'(1x,"*",1x,a,1x,"=",f14.7)') "volume of direct unit cell/Å³", &
         mol%volume*conv**3
      write(iunit,'(a)')
   endif

end subroutine print_pbcsum

