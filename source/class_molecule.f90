module class_molecule
   use iso_fortran_env, wp => real64
   implicit none

   public :: molecule
   private

   type :: molecule
      integer  :: nat
      character(len=2),allocatable :: sym(:)
      integer, allocatable :: at(:)
      real(wp),allocatable :: xyz(:,:)
      real(wp),allocatable :: abc(:)
      real(wp) :: chrg = 0.0_wp
      real(wp) :: lattice(3,3)
   contains
      procedure :: allocate => allocate_molecule
      procedure :: deallocate => deallocate_molecule
   end type

contains

subroutine allocate_molecule(self,nat,lpbc)
   implicit none
   class(molecule) :: self
   integer,intent(in) :: nat
   logical,intent(in) :: lpbc
   self%nat = nat
   allocate(self%sym(nat),   source = '  ')
   allocate(self%at(nat),    source = 0 )
   allocate(self%xyz(3,nat), source = 0.0_wp )
   if (lpbc) then
      allocate(self%abc(nat),source = 0.0_wp )
   endif
end subroutine allocate_molecule

subroutine deallocate_molecule(self)
   implicit none
   class(molecule) :: self
   if (allocated(self%sym))  deallocate(self%sym)
   if (allocated(self%at))   deallocate(self%at)
   if (allocated(self%xyz))  deallocate(self%xyz)
   if (allocated(self%abc))  deallocate(self%abc)
end subroutine deallocate_molecule

end module class_molecule
