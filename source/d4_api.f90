!  extern void D4_calculation(const int& natoms, const int* attyp,
!        const double& charge, const double* coord,
!        const DFTD_parameter& dparam, const DFTD_options& dopt,
!        double& energy, double* grad, double* hess);
!> external iso-c compatible interface to the DFT-D4 program
subroutine d4_calculation_c_api &
      &   (natoms,attyp,charge,coord,dparam_in,dopt_in,edisp,grad,hess) &
      &    bind(C,name="D4_calculation")

   use iso_fortran_env, wp => real64, istdout => output_unit
   use iso_c_binding

   use class_param
   use class_set
   use class_molecule

   implicit none

   !> number of atoms (used to determine array dimensions)
   integer(c_int),intent(in) :: natoms
   !> atom types as ordinal numbers
   integer(c_int),intent(in) :: attyp(natoms)
   !> total molecular charge
   real(c_double),intent(in) :: charge
   !> cartesian coordinates in bohr
   real(c_double),intent(in) :: coord(3,natoms)
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

   type(molecule) :: mol
   type(dftd_parameter) :: dparam
   type(dftd_options)   :: dopt
   real(wp) :: energy
   real(wp),allocatable :: gradient(:,:)
   real(wp),allocatable :: hessian(:,:)

   call mol%allocate(natoms,.false.)
   mol%at = attyp
   mol%xyz = coord
   mol%chrg = charge

   dparam = dparam_in
   dopt   = dopt_in

   allocate(gradient(3,mol%nat), hessian(3*mol%nat,3*mol%nat))
   energy = 0.0_wp
   gradient = 0.0_wp
   hessian = 0.0_wp

   call d4_calculation(istdout,dopt,mol,dparam,energy,gradient,hessian)

   if (dopt%lenergy .or. dopt%lgradient) &
   edisp = energy
   if (dopt%lgradient) &
   grad = gradient
   if (dopt%lhessian) &
   hess = hessian

   call mol%deallocate
   deallocate(gradient,hessian)

end subroutine d4_calculation_c_api

!  extern void D4_PBC_calculation(const int& natoms, const int* attyp,
!        const double& charge, const double* coord, const double* lattice,
!        const DFTD_parameter& dparam, const DFTD_options& dopt,
!        double& energy, double* grad, double* hess);
!> external iso-c compatible interface to the DFT-D4 program
subroutine d4_pbc_calculation_c_api &
      &   (natoms,attyp,charge,coord,lattice,dparam_in,dopt_in,edisp,grad,glat) &
      &    bind(C,name="D4_PBC_calculation")

   use iso_fortran_env, wp => real64, istdout => output_unit
   use iso_c_binding

   use class_param
   use class_set
   use class_molecule

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

   type(molecule) :: mol
   type(dftd_parameter) :: dparam
   type(dftd_options)   :: dopt
   real(wp) :: energy
   real(wp) :: lattice_grad(3,3)
   real(wp),allocatable :: gradient(:,:)
   integer, parameter :: wsc_rep(3) = [1,1,1]

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

   call generate_wsc(mol,mol%wsc,wsc_rep)

   dparam = dparam_in
   dopt   = dopt_in

   allocate(gradient(3,mol%nat))
   energy = 0.0_wp
   gradient = 0.0_wp
   lattice_grad = 0.0_wp

   call d4_pbc_calculation(istdout,dopt,mol,dparam,energy,gradient,lattice_grad)

   if (dopt%lenergy .or. dopt%lgradient) &
   edisp = energy
   if (dopt%lgradient) then
      grad = gradient
      glat = lattice_grad
   endif

   call mol%deallocate
   deallocate(gradient)

end subroutine d4_pbc_calculation_c_api

subroutine d4_pbc_calculation_f_api &
      &   (natoms,attyp,charge,coord,lattice,dparam_in,dopt_in,edisp,grad,glat)

   use iso_fortran_env, wp => real64, istdout => output_unit

   use class_param
   use class_set
   use class_molecule

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
   real(wp) :: energy
   real(wp) :: lattice_grad(3,3)
   real(wp), allocatable :: gradient(:,:)
   integer, parameter :: wsc_rep(3) = [1,1,1]

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

   call generate_wsc(mol,mol%wsc,wsc_rep)

   dparam = dparam_in
   dopt   = dopt_in

   allocate(gradient(3,mol%nat))
   energy = 0.0_wp
   gradient = 0.0_wp
   lattice_grad = 0.0_wp

   call d4_pbc_calculation(istdout,dopt,mol,dparam,energy,gradient,lattice_grad)

   write(*,*)'Dispersion energy / au: ', energy

   if (dopt%lenergy .or. dopt%lgradient) &
   edisp = energy
   if (dopt%lgradient) then
      grad = gradient
      glat = lattice_grad
   endif

   call mol%deallocate
   deallocate(gradient)

end subroutine d4_pbc_calculation_f_api
