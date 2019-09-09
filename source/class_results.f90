module class_results
   use iso_fortran_env, wp => real64
   implicit none

   public :: dftd_results
   private

   type :: dftd_results
      real(wp), allocatable :: energy
      real(wp), allocatable :: gradient(:,:)
      real(wp), allocatable :: stress(:)
      real(wp), allocatable :: lattice_gradient(:,:)
      real(wp), allocatable :: charges(:)
      real(wp), allocatable :: dipole_moment(:)
      real(wp), allocatable :: polarizibilities(:)
      real(wp), allocatable :: c6_coefficients(:,:)
      real(wp), allocatable :: hessian(:,:)
   contains
      procedure :: allocate => allocate_results
      procedure :: deallocate => deallocate_results
   end type dftd_results

contains

subroutine allocate_results(self,n)
   class(dftd_results), intent(inout) :: self
   integer, intent(in) :: n
   call self%deallocate
   self%energy = 0.0_wp
   allocate(self%gradient(3,n),         source = 0.0_wp)
   allocate(self%stress(6),             source = 0.0_wp)
   allocate(self%lattice_gradient(3,3), source = 0.0_wp)
   allocate(self%charges(n),            source = 0.0_wp)
   allocate(self%dipole_moment(3),      source = 0.0_wp)
   allocate(self%polarizibilities(n),   source = 0.0_wp)
   allocate(self%c6_coefficients(n,n),  source = 0.0_wp)
   allocate(self%hessian(3*n,3*n),      source = 0.0_wp)
end subroutine allocate_results

subroutine deallocate_results(self)
   class(dftd_results), intent(inout) :: self
   if (allocated(self%energy))           deallocate(self%energy)
   if (allocated(self%stress))           deallocate(self%stress)
   if (allocated(self%lattice_gradient)) deallocate(self%lattice_gradient)
   if (allocated(self%gradient))         deallocate(self%gradient)
   if (allocated(self%charges))          deallocate(self%charges)
   if (allocated(self%dipole_moment))    deallocate(self%dipole_moment)
   if (allocated(self%polarizibilities)) deallocate(self%polarizibilities)
   if (allocated(self%c6_coefficients))  deallocate(self%c6_coefficients)
   if (allocated(self%hessian))          deallocate(self%hessian)
end subroutine deallocate_results

end module class_results
