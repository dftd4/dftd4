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

!> official dftd4 ISO-C-API.
!
!  Every other API (Python/C++/..., except Fortran) should use this layer
!  at some point, this API gives enough flexibility downstream such that
!  the user/developer can wrap the ISO-C layer without much overhead.
!
!  We are not using C-reserved memory except for reading into Fortran-reserved
!  memory and the other way round, the only requirement is that the memory on
!  the C side has to be continous.
module c_api
   use iso_c_binding
   use iso_fortran_env, only: wp => real64

   implicit none

   !> some overloading for convience
   interface c_return
      module procedure :: c_return_double_0d
      module procedure :: c_return_double_1d
      module procedure :: c_return_double_2d
   end interface c_return

contains

!> external iso-c compatible interface to the DFT-D4 program
function d4_calculation_c_api &
      &   (natoms,attyp,charge,coord,file_in,dparam_in,dopt_in, &
      &    edisp,grad,hess,polarizibilities,c6_coefficients,charges) &
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
   real(c_double),intent(out) :: grad(3,natoms)
   !> molecular dispersion hessian
   real(c_double),intent(out) :: hess(3*natoms,3*natoms)
   !> optionally returns all atom partioned polarizibilities
   real(c_double),intent(out) :: polarizibilities(natoms)
   !> optionally returns all C6 coefficents
   real(c_double),intent(out) :: c6_coefficients(natoms,natoms)
   !> optionally returns all atomic partial charges
   real(c_double),intent(out) :: charges(natoms)
   !> success status of calculation
   integer(c_int) :: status

   type(mctc_logger) :: env

   type(molecule) :: mol
   type(dftd_parameter) :: dparam
   type(dftd_options)   :: dopt
   type(dftd_results)   :: dresults
   integer :: iunit
   character(len=:),allocatable :: outfile

   call mol%allocate(natoms,.false.)
   mol%at = attyp
   mol%xyz = coord
   mol%chrg = charge
   call mol%update

   call generate_wsc(mol,mol%wsc)

   dparam = dparam_in
   dopt   = dopt_in

   call c_string_convert(outfile, file_in)

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

   call c_return(edisp, dresults%energy)
   call c_return(grad, dresults%gradient)
   call c_return(hess, dresults%hessian)
   call c_return(charges, dresults%charges)
   call c_return(polarizibilities, dresults%polarizibilities)
   call c_return(c6_coefficients, dresults%c6_coefficients)

   call finalize
   status = 0

contains
subroutine finalize
   call mol%deallocate
   call dresults%deallocate
   if (iunit.ne.istdout) close(iunit)
end subroutine finalize
end function d4_calculation_c_api

!> external iso-c compatible interface to the DFT-D4 program
function d4_pbc_calculation_c_api &
      &   (natoms,attyp,charge,coord,lattice,pbc,file_in,dparam_in,dopt_in, &
      &    edisp,grad,glat,stress,hess,polarizibilities,c6_coefficients,charges) &
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
   real(c_double),intent(out) :: grad(3,natoms)
   !> dispersion lattice gradient
   real(c_double),intent(out) :: glat(3,3)
   !> dispersion stress tensor
   real(c_double),intent(out) :: stress(3,3)
   !> molecular dispersion hessian
   real(c_double),intent(out) :: hess(3*natoms,3*natoms)
   !> optionally returns all atom partioned polarizibilities
   real(c_double),intent(out) :: polarizibilities(natoms)
   !> optionally returns all C6 coefficents
   real(c_double),intent(out) :: c6_coefficients(natoms,natoms)
   !> optionally returns all atomic partial charges
   real(c_double),intent(out) :: charges(natoms)
   !> success status of calculation
   integer(c_int) :: status

   type(mctc_logger) :: env
   type(molecule) :: mol
   type(dftd_parameter) :: dparam
   type(dftd_results)   :: dresults
   type(dftd_options)   :: dopt
   integer :: iunit
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

   call c_string_convert(outfile, file_in)

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

   call c_return(edisp, dresults%energy)
   call c_return(grad, dresults%gradient)
   call c_return(glat, dresults%lattice_gradient)
   call c_return(stress, dresults%stress)
   call c_return(hess, dresults%hessian)
   call c_return(charges, dresults%charges)
   call c_return(polarizibilities, dresults%polarizibilities)
   call c_return(c6_coefficients, dresults%c6_coefficients)

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
   integer :: mbd_mode
   character(len=:),allocatable :: name

   mbd_mode = lmbd

   call c_string_convert(name, name_in)

   name = lowercase(name)
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

!> external iso-c compatible interface to the DFT-D4 program
function d3_calculation_c_api &
      &   (natoms,attyp,coord,file_in,dparam_in,dopt_in, &
      &    edisp,grad,hess,polarizibilities,c6_coefficients) &
      &    bind(C,name="D3_calculation") result(status)

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
   real(c_double),intent(out) :: grad(3,natoms)
   !> molecular dispersion hessian
   real(c_double),intent(out) :: hess(3*natoms,3*natoms)
   !> optionally returns all atom partioned polarizibilities
   real(c_double),intent(out) :: polarizibilities(natoms)
   !> optionally returns all C6 coefficents
   real(c_double),intent(out) :: c6_coefficients(natoms,natoms)
   !> success status of calculation
   integer(c_int) :: status

   type(mctc_logger) :: env

   type(molecule) :: mol
   type(dftd_parameter) :: dparam
   type(dftd_options)   :: dopt
   type(dftd_results)   :: dresults
   integer :: iunit
   character(len=:),allocatable :: outfile

   call mol%allocate(natoms,.false.)
   mol%at = attyp
   mol%xyz = coord
   mol%chrg = 0.0_wp
   call mol%update

   call generate_wsc(mol,mol%wsc)

   dparam = dparam_in
   dopt   = dopt_in

   call c_string_convert(outfile, file_in)

   if (outfile.ne.'-'.and.dopt%print_level > 0) then
      open(newunit=iunit,file=outfile)
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

   call d3_calculation(iunit,env,dopt,mol,dparam,dresults)
   if (.not.env%sane) then
      call finalize
      status = 1
      return
   endif

   call c_return(edisp, dresults%energy)
   call c_return(grad, dresults%gradient)
   call c_return(hess, dresults%hessian)
   call c_return(polarizibilities, dresults%polarizibilities)
   call c_return(c6_coefficients, dresults%c6_coefficients)

   call finalize
   status = 0

contains
subroutine finalize
   call mol%deallocate
   call dresults%deallocate
   if (iunit.ne.istdout) close(iunit)
end subroutine finalize
end function d3_calculation_c_api

!> external iso-c compatible interface to the DFT-D4 program
function d3_pbc_calculation_c_api &
      &   (natoms,attyp,coord,lattice,pbc,file_in,dparam_in,dopt_in, &
      &    edisp,grad,glat,stress,hess,polarizibilities,c6_coefficients) &
      &    bind(C,name="D3_PBC_calculation") result(status)

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
   real(c_double),intent(out) :: grad(3,natoms)
   !> dispersion lattice gradient
   real(c_double),intent(out) :: glat(3,3)
   !> dispersion stress tensor
   real(c_double),intent(out) :: stress(3,3)
   !> molecular dispersion hessian
   real(c_double),intent(out) :: hess(3*natoms,3*natoms)
   !> optionally returns all atom partioned polarizibilities
   real(c_double),intent(out) :: polarizibilities(natoms)
   !> optionally returns all C6 coefficents
   real(c_double),intent(out) :: c6_coefficients(natoms,natoms)
   !> success status of calculation
   integer(c_int) :: status

   type(mctc_logger) :: env
   type(molecule) :: mol
   type(dftd_parameter) :: dparam
   type(dftd_results)   :: dresults
   type(dftd_options)   :: dopt
   integer :: iunit
   character(len=:),allocatable :: outfile

   call mol%allocate(natoms,.false.)
   mol%at = attyp
   mol%xyz = coord
   mol%chrg = 0.0_wp
   mol%npbc = count(pbc)
   mol%pbc = pbc
   mol%lattice = lattice
   call mol%update

   call generate_wsc(mol,mol%wsc)

   dparam = dparam_in
   dopt   = dopt_in

   call c_string_convert(outfile, file_in)

   if (outfile.ne.'-'.and.dopt%print_level > 0) then
      open(newunit=iunit,file=outfile)
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

   call d3_calculation(iunit,env,dopt,mol,dparam,dresults)
   if (.not.env%sane) then
      call finalize
      status = 1
      return
   endif

   call c_return(edisp, dresults%energy)
   call c_return(grad, dresults%gradient)
   call c_return(glat, dresults%lattice_gradient)
   call c_return(stress, dresults%stress)
   call c_return(hess, dresults%hessian)
   call c_return(polarizibilities, dresults%polarizibilities)
   call c_return(c6_coefficients, dresults%c6_coefficients)

   call finalize
   status = 0

contains
subroutine finalize
   call mol%deallocate
   call dresults%deallocate
   if (iunit.ne.istdout) close(iunit)
end subroutine finalize
end function d3_pbc_calculation_c_api

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran value has been calculated
pure subroutine c_return_double_0d(c_array, f_array)
   real(c_double), intent(out), target :: c_array
   real(wp), allocatable, intent(in) :: f_array
   if (c_associated(c_loc(c_array)) .and. allocated(f_array)) then
      c_array = f_array
   endif
end subroutine c_return_double_0d

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran array has been populated
pure subroutine c_return_double_1d(c_array, f_array)
   real(c_double), intent(out), target :: c_array(:)
   real(wp), allocatable, intent(in) :: f_array(:)
   if (c_associated(c_loc(c_array)) .and. allocated(f_array)) then
      c_array = f_array
   endif
end subroutine c_return_double_1d

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran array has been populated (2D version)
pure subroutine c_return_double_2d(c_array, f_array)
   real(c_double), intent(out), target :: c_array(:,:)
   real(wp), allocatable, intent(in) :: f_array(:,:)
   if (c_associated(c_loc(c_array)) .and. allocated(f_array)) then
      c_array = f_array
   endif
end subroutine c_return_double_2d

!> convert vector of chars to deferred size character
subroutine c_string_convert(f_string, c_string)
   character(c_char), dimension(*), intent(in) :: c_string
   character(len=:), allocatable, intent(out) :: f_string
   integer :: i
   i = 0
   f_string = ''
   do
      i = i+1
      if (c_string(i).eq.c_null_char) exit
      f_string = f_string//c_string(i)
   enddo
end subroutine c_string_convert

end module c_api
