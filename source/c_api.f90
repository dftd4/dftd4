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

!> external iso-c compatible interface to the DFT-D4 program
function d4_calculation_c_api &
      &   (natoms,attyp,charge,coord,file_in,dparam_in,dopt_in,edisp,grad,hess) &
      &    bind(C,name="D4_calculation") result(status)

   use iso_fortran_env, wp => real64, istdout => output_unit
   use iso_c_binding

   use mctc_environment

   use class_param
   use class_set
   use class_molecule
   use class_results

   use dispersion_calculator

   implicit none

   !> number of atoms (used to determine array dimensions)
   integer(c_int),intent(in) :: natoms
   !> atom types as ordinal numbers
   integer(c_int),intent(in) :: attyp(natoms)
   !> total molecular charge
   real(c_double),intent(in) :: charge
   !> cartesian coordinates in bohr
   real(c_double),intent(in) :: coord(3,natoms)
   !> input file name
   character(kind=c_char),intent(in) :: file_in(*)
   !> damping parameters
   type(c_dftd_parameter),intent(in) :: dparam_in
   !> calculation options
   type(c_dftd_options),  intent(in) :: dopt_in
   !> final dispersion energy
   real(c_double),intent(out) :: edisp
   !> molecular dispersion gradient
   !  (not referenced of dopt_in.lgradient is false)
   real(c_double),intent(out) :: grad(3,natoms)
   !> molecular dispersion hessian
   !  (not referenced of dopt_in.lhessian is false)
   real(c_double),intent(out) :: hess(3*natoms,3*natoms)
   !> success status of calculation
   integer(c_int) :: status

   type(mctc_logger) :: env

   type(molecule) :: mol
   type(dftd_parameter) :: dparam
   type(dftd_options)   :: dopt
   type(dftd_results)   :: dresults
   integer :: i,iunit
   character(len=:),allocatable :: outfile

   call mol%allocate(natoms,.false.)
   mol%at = attyp
   mol%xyz = coord
   mol%chrg = charge
   call mol%update

   call generate_wsc(mol,mol%wsc)

   dparam = dparam_in
   dopt   = dopt_in

   i = 0
   outfile = ''
   do
      i = i+1
      if (file_in(i).eq.c_null_char) exit
      outfile = outfile//file_in(i)
   enddo

   if (outfile.ne.'-'.and.dopt%print_level > 0) then
      open(newunit=iunit,file=outfile)
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

   call d4_calculation(iunit,env,dopt,mol,dparam,dresults)
   if (.not.env%sane) then
      call finalize
      status = 1
      return
   endif

   if (allocated(dresults%energy)) &
      edisp = dresults%energy
   if (allocated(dresults%gradient)) &
      grad = dresults%gradient
   if (allocated(dresults%hessian)) &
      hess = dresults%hessian

   call finalize
   status = 0

contains
subroutine finalize
   call mol%deallocate
   call dresults%deallocate
   if (iunit.ne.istdout) close(iunit)
end subroutine finalize
end function d4_calculation_c_api

!  extern void D4_PBC_calculation(const int& natoms, const int* attyp,
!        const double& charge, const double* coord, const double* lattice,
!        const DFTD_parameter& dparam, const DFTD_options& dopt,
!        double& energy, double* grad, double* hess);
!> external iso-c compatible interface to the DFT-D4 program
function d4_pbc_calculation_c_api &
      &   (natoms,attyp,charge,coord,lattice,pbc,file_in,dparam_in,dopt_in, &
      &    edisp,grad,glat) &
      &    bind(C,name="D4_PBC_calculation") result(status)

   use iso_fortran_env, wp => real64, istdout => output_unit
   use iso_c_binding

   use mctc_environment

   use class_param
   use class_set
   use class_molecule
   use class_results

   use dispersion_calculator
   use pbc_tools

   implicit none

   !> number of atoms (used to determine array dimensions)
   integer(c_int),intent(in) :: natoms
   !> atom types as ordinal numbers
   integer(c_int),intent(in) :: attyp(natoms)
   !> total molecular charge
   real(c_double),intent(in) :: charge
   !> cartesian coordinates in bohr
   real(c_double),intent(in) :: coord(3,natoms)
   !> lattice parameters
   real(c_double),intent(in) :: lattice(3,3)
   !> periodic dimensions
   logical(c_bool),intent(in) :: pbc(3)
   !> input file name
   character(kind=c_char),intent(in) :: file_in(*)
   !> damping parameters
   type(c_dftd_parameter),intent(in) :: dparam_in
   !> calculation options
   type(c_dftd_options),  intent(in) :: dopt_in
   !> final dispersion energy
   real(c_double),intent(out) :: edisp
   !> molecular dispersion gradient
   !  (not referenced of dopt_in.lgradient is false)
   real(c_double),intent(out) :: grad(3,natoms)
   !> dispersion lattice gradient
   !  (not referenced of dopt_in.lgradient is false)
   real(c_double),intent(out) :: glat(3,3)
   !> success status of calculation
   integer(c_int) :: status

   type(mctc_logger) :: env
   type(molecule) :: mol
   type(dftd_parameter) :: dparam
   type(dftd_results)   :: dresults
   type(dftd_options)   :: dopt
   integer :: i,iunit
   character(len=:),allocatable :: outfile

   call mol%allocate(natoms,.false.)
   mol%at = attyp
   mol%xyz = coord
   mol%chrg = charge
   mol%npbc = count(pbc)
   mol%pbc = pbc
   mol%lattice = lattice
   call mol%update

   call generate_wsc(mol,mol%wsc)

   dparam = dparam_in
   dopt   = dopt_in

   i = 0
   outfile = ''
   do
      i = i+1
      if (file_in(i).eq.c_null_char) exit
      outfile = outfile//file_in(i)
   enddo

   if (outfile.ne.'-'.and.dopt%print_level > 0) then
      open(newunit=iunit,file=outfile)
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

   call d4_calculation(iunit,env,dopt,mol,dparam,dresults)
   if (.not.env%sane) then
      call finalize
      status = 1
      return
   endif

   if (allocated(dresults%energy)) &
      edisp = dresults%energy
   if (allocated(dresults%gradient)) &
      grad = dresults%gradient
   if (allocated(dresults%lattice_gradient)) &
      glat = dresults%lattice_gradient

   call finalize
   status = 0

contains
subroutine finalize
   call mol%deallocate
   call dresults%deallocate
   if (iunit.ne.istdout) close(iunit)
end subroutine finalize
end function d4_pbc_calculation_c_api

function d4_damping_parameters_c_api(name_in,dparam_out,lmbd) &
      &  bind(C,name="D4_damping_parameters") result(status)
   use iso_c_binding
   use class_param
   use mctc_environment
   use dfuncpar
   implicit none
   !> functional name
   character(kind=c_char),intent(in) :: name_in(*)
   !> non-additive dispersion mode
   integer(c_int),intent(in) :: lmbd
   !> damping parameters
   type(c_dftd_parameter),intent(out) :: dparam_out
   !> success status of calculation
   integer(c_int) :: status

   type(mctc_logger) :: env
   type(dftd_parameter) :: dparam
   integer :: mbd_mode,i
   character(len=:),allocatable :: name

   mbd_mode = lmbd

   i = 0
   name = ''
   do
      i = i+1
      if (name_in(i).eq.c_null_char) exit
      name = name//name_in(i)
   enddo

   if (get_dfnum(name) > 0) then
      call d4par(name,dparam,mbd_mode,env)
      if (.not.env%sane) then
         status = 2
         return
      endif

      dparam_out = dparam
      if (lmbd == 0) then
         status = -2 ! not parametrized, but D4-ATM parameters returned
      else
         status = 0
      endif
   else
      status = -1 ! not found
   endif

end function d4_damping_parameters_c_api
