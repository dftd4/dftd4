!> @brief definition of all used parameter types
!!
!! This module provides a type for the DFT-D damping parameters,
!! the charge model parameters for the EEQ calculation and
!! the SCC parameters used in GFN-xTB
module class_param
   use iso_fortran_env, only : wp => real64
   use iso_c_binding

   public :: dftd_parameter
   public :: c_dftd_parameter
   public :: convert_dftd_parameter

   public ::  scc_parameter
   public :: chrg_parameter
   private

!> @brief damping parameters for DFT-D
   type :: dftd_parameter
      sequence
      real(wp) :: s6  = -1.0_wp !< scaling of dipole-dipole dispersion
      real(wp) :: s8  = -1.0_wp !< scaling of dipole-quadrupole dispersion
      real(wp) :: s10 =  0.0_wp !< scaling of higher order dispersion
      real(wp) :: a1  = -1.0_wp !< scaling of vdW-Radius in finite damping
      real(wp) :: a2  = -1.0_wp !< constant offset off vdW-Radius in finite damping
      real(wp) :: s9  =  1.0_wp !< scaling of non-addititive dispersion
      integer  :: alp = 16      !< exponent of zero damping
      real(wp) :: beta = 1.0_wp !< range separation parameter for Fermi-damping
   end type dftd_parameter

   type,bind(C) :: c_dftd_parameter
      real(c_double) :: s6  = -1.0_wp !< scaling of dipole-dipole dispersion
      real(c_double) :: s8  = -1.0_wp !< scaling of dipole-quadrupole dispersion
      real(c_double) :: s10 =  0.0_wp !< scaling of higher order dispersion
      real(c_double) :: a1  = -1.0_wp !< scaling of vdW-Radius in finite damping
      real(c_double) :: a2  = -1.0_wp !< constant offset off vdW-Radius in finite damping
      real(c_double) :: s9  =  1.0_wp !< scaling of non-addititive dispersion
      integer(c_int) :: alp = 16      !< exponent of zero damping
      real(c_double) :: beta = 1.0_wp !< range separation parameter for Fermi-damping
   end type c_dftd_parameter

   interface convert_dftd_parameter
      module procedure :: convert_dftd_parameter_c_to_f
      module procedure :: convert_dftd_parameter_f_to_c
   end interface convert_dftd_parameter

!> @brief parameters for self consistent charge calculation in GFN-xTB
   type scc_parameter
      real(wp) :: kspd(6)
      real(wp) :: gscal
      real(wp) :: gam3l(0:3)
      real(wp) :: kcnsh(4)
      real(wp) :: kmagic(4,4)
      type(dftd_parameter) :: disp
      real(wp) :: cn_shift
      real(wp) :: cn_expo
      real(wp) :: cn_rmax
      real(wp) :: xbdamp
      real(wp) :: xbrad
      real(wp) :: ipshift
      real(wp) :: eashift
      real(wp) :: ken1
      real(wp) :: zqf
      real(wp) :: zcnf
      real(wp) :: fpol
      real(wp) :: wllscal
      real(wp) :: tfac
      real(wp) :: kcn
      real(wp) :: alphaj
      real(wp) :: lshift
      real(wp) :: lshifta
      real(wp) :: split
      real(wp) :: kenscal
      real(wp) :: g_a
      real(wp) :: g_c
      real(wp) :: wf
   end type scc_parameter

!> @brief parameters for the extended electronegativity model
!!
!! Contains fitted parameters for electronegativity, CN scaling factor,
!! chemical hardness and atom radius for a given molecule type.
   type chrg_parameter
      integer :: n
      real(wp),allocatable :: xi(:)
      real(wp),allocatable :: gam(:)
      real(wp),allocatable :: kappa(:)
      real(wp),allocatable :: alpha(:)
   contains
      procedure :: allocate => allocate_chrg
      procedure :: deallocate => deallocate_chrg
   end type chrg_parameter

contains

!> @brief constructor for charge model
subroutine allocate_chrg(self,n)
   implicit none
   class(chrg_parameter) :: self
   integer,intent(in)    :: n
   self%n = n
   allocate(self%xi(n),    source = 0.0_wp)
   allocate(self%gam(n),   source = 0.0_wp)
   allocate(self%kappa(n), source = 0.0_wp)
   allocate(self%alpha(n), source = 0.0_wp)
end subroutine allocate_chrg

!> @brief deconstructor for charge model
subroutine deallocate_chrg(self)
   implicit none
   class(chrg_parameter) :: self
   if (allocated(self%xi))    deallocate(self%xi)
   if (allocated(self%gam))   deallocate(self%gam)
   if (allocated(self%kappa)) deallocate(self%kappa)
   if (allocated(self%alpha)) deallocate(self%alpha)
end subroutine deallocate_chrg

pure elemental function convert_dftd_parameter_c_to_f &
      (c_dparam) result(f_dparam)
   implicit none
   type(c_dftd_parameter),intent(in) :: c_dparam
   type(dftd_parameter)              :: f_dparam
   f_dparam%s6   = c_dparam%s6
   f_dparam%s8   = c_dparam%s8
   f_dparam%s10  = c_dparam%s10
   f_dparam%a1   = c_dparam%a1
   f_dparam%a2   = c_dparam%a2
   f_dparam%s9   = c_dparam%s9
   f_dparam%alp  = c_dparam%alp
   f_dparam%beta = c_dparam%beta
end function convert_dftd_parameter_c_to_f

pure elemental function convert_dftd_parameter_f_to_c &
      (f_dparam) result(c_dparam)
   implicit none
   type(c_dftd_parameter)            :: c_dparam
   type(dftd_parameter),  intent(in) :: f_dparam
   c_dparam%s6   = f_dparam%s6
   c_dparam%s8   = f_dparam%s8
   c_dparam%s10  = f_dparam%s10
   c_dparam%a1   = f_dparam%a1
   c_dparam%a2   = f_dparam%a2
   c_dparam%s9   = f_dparam%s9
   c_dparam%alp  = f_dparam%alp
   c_dparam%beta = f_dparam%beta
end function convert_dftd_parameter_f_to_c

end module class_param
