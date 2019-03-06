!> @brief definition of the molecular structure class
!!
!! contains number of atoms, geometry, atomtypes and molecular charge
module class_molecule
   use iso_fortran_env, wp => real64
   implicit none

   public :: molecule
   private

!> @class molecule
!> @brief molecule structure class
   type :: molecule
      integer  :: nat !< number of atoms
      character(len=2),allocatable :: sym(:) !< element symbol
      integer, allocatable :: at(:)          !< ordinal number
      real(wp),allocatable :: xyz(:,:)       !< cartesian coordinates in bohr
      real(wp),allocatable :: abc(:)
      real(wp) :: chrg = 0.0_wp              !< molecular charge in e
      real(wp) :: lattice(3,3)
   contains
      procedure :: allocate => allocate_molecule
      procedure :: deallocate => deallocate_molecule
   end type

contains

!> @brief constructor for molecule class
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

!> @brief deconstructor for molecule class
subroutine deallocate_molecule(self)
   implicit none
   class(molecule) :: self
   if (allocated(self%sym))  deallocate(self%sym)
   if (allocated(self%at))   deallocate(self%at)
   if (allocated(self%xyz))  deallocate(self%xyz)
   if (allocated(self%abc))  deallocate(self%abc)
end subroutine deallocate_molecule

end module class_molecule
