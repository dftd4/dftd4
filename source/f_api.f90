! This file is part of dftd4.
!
! Copyright (C) 2017-2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

subroutine d4_pbc_calculation_f_api &
      &   (natoms,attyp,charge,coord,lattice,dparam_in,dopt_in,edisp,grad,glat)

   use iso_fortran_env, wp => real64, istdout => output_unit

   use class_param
   use class_set
   use class_molecule
   use class_results

   use dispersion_calculator
   use pbc_tools

   implicit none

   !> number of atoms (used to determine array dimensions)
   integer,intent(in) :: natoms
   !> atom types as ordinal numbers
   integer,intent(in) :: attyp(natoms)
   !> total molecular charge
   real(wp),intent(in) :: charge
   !> cartesian coordinates in bohr
   real(wp),intent(in) :: coord(3,natoms)
   !> lattice parameters
   real(wp),intent(in) :: lattice(3,3)
   !> damping parameters
   type(dftd_parameter),intent(in) :: dparam_in
   !> calculation options
   type(dftd_options),  intent(in) :: dopt_in
   !> final dispersion energy
   real(wp),intent(out) :: edisp
   !> molecular dispersion gradient
   !  (not referenced of dopt_in.lgradient is false)
   real(wp),intent(out) :: grad(3,natoms)
   !> dispersion lattice gradient
   !  (not referenced of dopt_in.lgradient is false)
   real(wp),intent(out) :: glat(3,3)

   type(molecule) :: mol
   type(dftd_parameter) :: dparam
   type(dftd_options)   :: dopt
   type(dftd_results)   :: dresults

   write(*,*) "You're entering the DFT-D4 API for periodic systems" 

   call mol%allocate(natoms,.false.)
   mol%at = attyp
   mol%xyz = coord
   mol%chrg = charge
   mol%npbc = 3
   mol%pbc = .true.
   mol%lattice = lattice
   mol%volume = dlat_to_dvol(mol%lattice)
   call dlat_to_cell(mol%lattice,mol%cellpar)
   call dlat_to_rlat(mol%lattice,mol%rec_lat)
   call mol%wrap_back
   call mol%calculate_distances

   call generate_wsc(mol,mol%wsc)

   dparam = dparam_in
   dopt   = dopt_in

   call d4_calculation(istdout,dopt,mol,dparam,dresults)

   write(*,*)'Dispersion energy / au: ', dresults%energy

   if (allocated(dresults%energy)) &
      edisp = dresults%energy
   if (allocated(dresults%gradient)) &
      grad = dresults%gradient
   if (allocated(dresults%lattice_gradient)) &
      glat = dresults%lattice_gradient

   call mol%deallocate
   call dresults%deallocate

end subroutine d4_pbc_calculation_f_api
