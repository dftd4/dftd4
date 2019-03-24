!  extern void D4_calculation(const int& natoms, const int* attyp,
!        const double& charge, const double* coord,
!        const DFTD_parameter& dparam, const DFTD_options& dopt,
!        double& energy, double* grad, double* hess);
subroutine d4_calculation_c_api &
      &   (natoms,attyp,charge,coord,dparam_in,dopt_in,edisp,grad,hess) &
      &    bind(C,name="D4_calculation")

   use iso_fortran_env, wp => real64, istdout => output_unit
   use iso_c_binding

   use class_param
   use class_set
   use class_molecule

   implicit none

   integer(c_int),intent(in) :: natoms
   integer(c_int),intent(in) :: attyp(natoms)
   real(c_double),intent(in) :: charge
   real(c_double),intent(in) :: coord(3,natoms)
   type(c_dftd_parameter),intent(in) :: dparam_in
   type(c_dftd_options),  intent(in) :: dopt_in
   real(c_double),intent(out) :: edisp
   real(c_double),intent(out) :: grad(3,natoms)
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
   
   dparam = convert_dftd_parameter(dparam_in)
   dopt   = convert_dftd_options  (dopt_in)

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
