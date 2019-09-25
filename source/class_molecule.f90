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

!> definition of the molecular structure class
!
!  Contains number of atoms, geometry, atomtypes and molecular charge.
!  With periodic boundary conditions the molecular structure holds
!  information about the cell and the periodic directions as well
!  a Wigner--Seitz cell as primitive unit.
module class_molecule
   use iso_fortran_env, wp => real64
   use class_wsc
   implicit none

   public :: molecule
   private

!> molecule structure class
   type :: molecule
      integer  :: n = 0            !< number of atoms
      real(wp) :: chrg = 0.0_wp    !< molecular charge in e
      logical  :: pbc(3) = .false. !< periodic dimensions
      integer  :: npbc = 0         !< periodicity of system
      character(len=2),allocatable :: sym(:) !< element symbol
      integer, allocatable :: at(:)          !< ordinal number
      real(wp),allocatable :: xyz(:,:)       !< cartesian coordinates in bohr
      real(wp),allocatable :: abc(:,:)       !< fractional coordinates
      real(wp),allocatable :: dist(:,:)      !< interatomic distances
      real(wp),allocatable :: atmass(:)      !< atomic masses in au (!)
      real(wp),allocatable :: z(:)           !< nuclear charges
      real(wp),allocatable :: cn(:)          !< coordination number
      real(wp) :: cellpar(6) = 0.0_wp        !< cell parameters
      real(wp) :: lattice(3,3) = 0.0_wp      !< direct lattice parameters
      real(wp) :: rec_lat(3,3) = 0.0_wp      !< reciprocal lattice parameters
      real(wp) :: volume = 0.0_wp            !< volume of unit cell
      type(ws_cell) :: wsc                   !< Wigner--Seitz cell
   contains
      !> constructor for molecular structure information
      procedure :: allocate => allocate_molecule
      !> (optional) deconstructor for molecular structure information
      procedure :: deallocate => deallocate_molecule
      !> debug information of internal data of the molecular structure
      procedure :: write => write_molecule
      !> update molecular structure data from xyz and lattice
      procedure :: update
      !> retrieves atomic masses from MCTC library and stores them in au
      procedure :: assign_atomic_mass
      !> calculates all interatomic distances for molecular case
      !  or minimum image distances if periodic boundary conditions are
      !  present
      procedure :: calculate_distances
      !> wraps back the coordinates into the primitive cell for
      !  periodic boundary conditions
      procedure :: wrap_back
      !> calculates molecular mass (needs atomic masses assigned)
      procedure :: molecular_mass
      !> returns the center of geometry for the non-periodic directions
      procedure :: center_of_geometry
      !> returns the center of mass for the non-periodic directions
      procedure :: center_of_mass
      !> shifts the non-periodic directions back to the center of geometry
      procedure :: shift_to_center_of_geometry
      !> shifts the non-periodic directions back to the center of mass
      procedure :: shift_to_center_of_mass
      !> calculates the moments of inertia (only for non-periodic directions)
      procedure :: moments_of_inertia
      !> aligns the molecule to its principal axes of inertia
      procedure :: align_to_principal_axes
   end type

   !> wrap Lapack and pretend it is a pure routine to call
   interface spev
      pure subroutine sspev(jobz,uplo,n,ap,w,z,ldz,work,info)
         use iso_fortran_env, only : wp => real32
         real(wp), intent(inout) :: ap(*)
         real(wp), intent(out) :: w(*)
         character, intent(in) :: uplo
         real(wp), intent(out) :: z(ldz,*)
         integer, intent(out) :: info
         character, intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(wp), intent(in) :: work(*)
      end subroutine sspev
      pure subroutine dspev(jobz,uplo,n,ap,w,z,ldz,work,info)
         use iso_fortran_env, only : wp => real64
         real(wp), intent(inout) :: ap(*)
         real(wp), intent(out) :: w(*)
         character, intent(in) :: uplo
         real(wp), intent(out) :: z(ldz,*)
         integer, intent(out) :: info
         character, intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(wp), intent(in) :: work(*)
      end subroutine dspev
   end interface spev

contains

!> constructor for molecule class
subroutine allocate_molecule(self,nat,lpbc)
   implicit none
   class(molecule) :: self !< molecule structure information
   integer,intent(in) :: nat !< number of atoms
   logical,optional,intent(in) :: lpbc !< periodic boundary conditions
   call self%deallocate
   self%n = nat
   allocate( self%sym(nat),      source = '  ')
   allocate( self%at(nat),       source = 0 )
   allocate( self%xyz(3,nat),    source = 0.0_wp )
   allocate( self%abc(3,nat),    source = 0.0_wp )
   allocate( self%dist(nat,nat), source = 0.0_wp )
   allocate( self%atmass(nat),   source = 0.0_wp )
   allocate( self%z(nat),        source = 0.0_wp )
   allocate( self%cn(nat),       source = 0.0_wp )
end subroutine allocate_molecule

!> deconstructor for molecule class
subroutine deallocate_molecule(self)
   implicit none
   class(molecule) :: self !< molecule structure information
   self%n = 0
   self%pbc = .false.
   self%chrg = 0.0_wp
   self%lattice = 0.0_wp
   if (allocated(self%sym))    deallocate(self%sym)
   if (allocated(self%at))     deallocate(self%at)
   if (allocated(self%xyz))    deallocate(self%xyz)
   if (allocated(self%abc))    deallocate(self%abc)
   if (allocated(self%dist))   deallocate(self%dist)
   if (allocated(self%atmass)) deallocate(self%atmass)
   if (allocated(self%z))      deallocate(self%z)
   if (allocated(self%cn))     deallocate(self%cn)
end subroutine deallocate_molecule

!> print information about current molecular structure to unit
subroutine write_molecule(self,iunit,comment)
   implicit none
   class(molecule),   intent(in) :: self    !< molecular structure information
   integer,           intent(in) :: iunit   !< file handle
   character(len=*),  intent(in) :: comment !< name of the variable
   character(len=*),parameter :: dfmt = '(1x,a,1x,"=",1x,g0)'

   write(iunit,'(72(">"))')
   write(iunit,'(1x,"*",1x,a)') "Writing 'molecule' class"
   write(iunit,'(  "->",1x,a)') comment
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "status of the fields"
   write(iunit,dfmt) "integer :: nat         ",self%n
   write(iunit,dfmt) "real    :: chrg        ",self%chrg
   write(iunit,dfmt) "integer :: npbc        ",self%npbc
   write(iunit,dfmt) "logical :: pbc(1)      ",self%pbc(1)
   write(iunit,dfmt) "        &  pbc(2)      ",self%pbc(2)
   write(iunit,dfmt) "        &  pbc(3)      ",self%pbc(3)
   write(iunit,dfmt) "real    :: volume      ",self%volume
   write(iunit,dfmt) "real    :: lattice(1,1)",self%lattice(1,1)
   write(iunit,dfmt) "        &  lattice(2,1)",self%lattice(2,1)
   write(iunit,dfmt) "        &  lattice(3,1)",self%lattice(3,1)
   write(iunit,dfmt) "        &  lattice(1,2)",self%lattice(1,2)
   write(iunit,dfmt) "        &  lattice(2,2)",self%lattice(2,2)
   write(iunit,dfmt) "        &  lattice(3,2)",self%lattice(3,2)
   write(iunit,dfmt) "        &  lattice(1,3)",self%lattice(1,3)
   write(iunit,dfmt) "        &  lattice(2,3)",self%lattice(2,3)
   write(iunit,dfmt) "        &  lattice(3,3)",self%lattice(3,3)
   write(iunit,dfmt) "real    :: rec_lat(1,1)",self%rec_lat(1,1)
   write(iunit,dfmt) "        &  rec_lat(2,1)",self%rec_lat(2,1)
   write(iunit,dfmt) "        &  rec_lat(3,1)",self%rec_lat(3,1)
   write(iunit,dfmt) "        &  rec_lat(1,2)",self%rec_lat(1,2)
   write(iunit,dfmt) "        &  rec_lat(2,2)",self%rec_lat(2,2)
   write(iunit,dfmt) "        &  rec_lat(3,2)",self%rec_lat(3,2)
   write(iunit,dfmt) "        &  rec_lat(1,3)",self%rec_lat(1,3)
   write(iunit,dfmt) "        &  rec_lat(2,3)",self%rec_lat(2,3)
   write(iunit,dfmt) "        &  rec_lat(3,3)",self%rec_lat(3,3)
   write(iunit,dfmt) "real    :: cellpar(1)  ",self%cellpar(1)
   write(iunit,dfmt) "        &  cellpar(2)  ",self%cellpar(2)
   write(iunit,dfmt) "        &  cellpar(3)  ",self%cellpar(3)
   write(iunit,dfmt) "        &  cellpar(4)  ",self%cellpar(4)
   write(iunit,dfmt) "        &  cellpar(5)  ",self%cellpar(5)
   write(iunit,dfmt) "        &  cellpar(6)  ",self%cellpar(6)
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "allocation status"
   write(iunit,dfmt) "allocated? sym(:)      ",allocated(self%sym)
   write(iunit,dfmt) "allocated? at(:)       ",allocated(self%at)
   write(iunit,dfmt) "allocated? xyz(:,:)    ",allocated(self%xyz)
   write(iunit,dfmt) "allocated? abc(:,:)    ",allocated(self%abc)
   write(iunit,dfmt) "allocated? dist(:,:)   ",allocated(self%dist)
   write(iunit,dfmt) "allocated? atmass(:)   ",allocated(self%atmass)
   write(iunit,dfmt) "allocated? z(:)        ",allocated(self%z)
   write(iunit,dfmt) "allocated? cn(:)       ",allocated(self%cn)
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "size of memory allocation"
   if (allocated(self%at)) then
   write(iunit,dfmt) "size(1) :: at(*)       ",size(self%at,1)
   endif
   if (allocated(self%xyz)) then
   write(iunit,dfmt) "size(1) :: xyz(*,:)    ",size(self%xyz,1)
   write(iunit,dfmt) "size(2) :: xyz(:,*)    ",size(self%xyz,2)
   endif
   if (allocated(self%abc)) then
   write(iunit,dfmt) "size(1) :: abc(*,:)    ",size(self%abc,1)
   write(iunit,dfmt) "size(2) :: abc(:,*)    ",size(self%abc,2)
   endif
   if (allocated(self%dist)) then
   write(iunit,dfmt) "size(1) :: dist(*,:)   ",size(self%dist,1)
   write(iunit,dfmt) "size(2) :: dist(:,*)   ",size(self%dist,2)
   endif
   if (allocated(self%atmass)) then
   write(iunit,dfmt) "size(1) :: atmass(*)   ",size(self%atmass,1)
   endif
   if (allocated(self%z)) then
   write(iunit,dfmt) "size(1) :: z(*)        ",size(self%z,1)
   endif
   if (allocated(self%cn)) then
   write(iunit,dfmt) "size(1) :: cn(*)       ",size(self%cn,1)
   endif
   call self%wsc%write(iunit,comment//"%wsc")
   write(iunit,'(72("<"))')
end subroutine write_molecule

subroutine update(self)
   use iso_fortran_env, wp => real64
   use pbc_tools
   implicit none
   class(molecule),intent(inout) :: self  !< molecular structure information

   if (self%npbc > 0) then
      call dlat_to_cell(self%lattice,self%cellpar)
      call dlat_to_rlat(self%lattice,self%rec_lat)
      self%volume = dlat_to_dvol(self%lattice)

      call self%wrap_back
   endif

   call self%calculate_distances

end subroutine update

!> set masses for molecular structure
subroutine assign_atomic_mass(self)
   use iso_fortran_env, wp => real64
   use mctc_param
   implicit none
   class(molecule),intent(inout) :: self !< molecular structure information

   self%atmass = atomic_mass(self%at)

end subroutine assign_atomic_mass

!> calculates all distances for molecular structures and minimum
!  image distances for periodic structures
subroutine calculate_distances(self)
   use iso_fortran_env, wp => real64
   use pbc_tools
   implicit none
   class(molecule),intent(inout) :: self !< molecular structure information
   integer :: i,j
   if (self%npbc > 0) then
      do i = 1, self%n
         do j = 1, i-1
            self%dist(j,i) = minimum_image_distance(.false.,self%abc(:,i), &
               &              self%abc(:,j),self%lattice,self%pbc)
            self%dist(i,j) = self%dist(j,i)
         enddo
         self%dist(i,i) = minimum_image_distance(.true.,self%abc(:,i), &
            &              self%abc(:,i),self%lattice,self%pbc)
      enddo
   else
      do i = 1, self%n
         do j = 1, i-1
            self%dist(j,i) = norm2(self%xyz(:,j)-self%xyz(:,i))
            self%dist(i,j) = self%dist(j,i)
         enddo
         self%dist(i,i) = 0.0_wp
      enddo
   endif
end subroutine calculate_distances

!> wrap cartesian coordinates back into cell
!
!  This automatically done when calling xyz_to_abc, so we only have
!  to perform the transformation there and back again
subroutine wrap_back(self)
   use iso_fortran_env, wp => real64
   use pbc_tools
   implicit none
   class(molecule),intent(inout) :: self !< molecular structure information
   call xyz_to_abc(self%n,self%lattice,self%xyz,self%abc,self%pbc)
   call abc_to_xyz(self%n,self%lattice,self%abc,self%xyz)
end subroutine wrap_back

!> returns the center of geometry for the non-periodic directions of
!  the system
pure function center_of_geometry(self) result(center)
   use iso_fortran_env, wp => real64
   implicit none
   class(molecule),intent(in) :: self !< molecular structure information
   real(wp) :: center(3)
   integer  :: idir
   center = 0.0_wp

   do idir = 1, 3
      if (.not.self%pbc(idir)) &
         center(idir) = sum(self%xyz(idir,:))
   enddo
   center = center/real(self%n,wp)
end function center_of_geometry

!> shifts back the non-periodic directions to the center of geometry
pure subroutine shift_to_center_of_geometry(self)
   use iso_fortran_env, wp => real64
   implicit none
   class(molecule),intent(inout) :: self !< molecular structure information
   real(wp) :: center(3)
   integer  :: iat
   center = self%center_of_geometry()
   do iat = 1, self%n
      self%xyz(:,iat) = self%xyz(:,iat) - center
   enddo
end subroutine shift_to_center_of_geometry

!> calculates the molecular mass, needs assigned molecular masses
pure function molecular_mass(self) result(molmass)
   use iso_fortran_env, wp => real64
   implicit none
   class(molecule),intent(in) :: self !< molecular structure information
   real(wp) :: molmass
   molmass = sum(self%atmass)
end function molecular_mass

!> returns the center of mass for the non-periodic directions of the system
pure function center_of_mass(self) result(center)
   use iso_fortran_env, wp => real64
   implicit none
   class(molecule),intent(in) :: self !< molecular structure information
   real(wp) :: center(3)
   integer  :: idir
   center = 0.0_wp
   do idir = 1, 3
      if (.not.self%pbc(idir)) &
         center(idir) = sum(self%atmass*self%xyz(idir,:))
   enddo
   center = center/sum(self%atmass)
end function center_of_mass

!> shifts back the non-periodic directions to the center of mass
pure subroutine shift_to_center_of_mass(self)
   use iso_fortran_env, wp => real64
   implicit none
   class(molecule),intent(inout) :: self !< molecular structure information
   real(wp) :: center(3)
   integer  :: iat
   center = self%center_of_mass()
   do iat = 1, self%n
      self%xyz(:,iat) = self%xyz(:,iat) - center
   enddo
end subroutine shift_to_center_of_mass

!> calculates the moments of inertia (only for non-periodic directions)
pure function moments_of_inertia(self) result(moments)
   use iso_fortran_env, wp => real64
   implicit none
   class(molecule),intent(in) :: self !< molecular structure information
   real(wp) :: moments(3)
   real(wp) :: center(3),atmass,t(6),work(9)
   real(wp) :: x,x2,y,y2,z,z2
   integer  :: iat,info
   ! currently not supported
   if (self%npbc > 0) then
      moments = -1.0_wp
      return
   endif
   center = self%center_of_mass()

   t = 0.0_wp

   do iat = 1, self%n
      atmass = self%atmass(iat)
      x = self%xyz(1,iat)-center(1); x2 = x**2
      y = self%xyz(2,iat)-center(2); y2 = y**2
      z = self%xyz(3,iat)-center(3); z2 = z**2
      t(1) = t(1) + atmass * (y2+z2)
      t(2) = t(2) - atmass * x*y
      t(3) = t(3) + atmass * (x2+z2)
      t(4) = t(4) - atmass * x*z
      t(5) = t(5) - atmass * y*z
      t(6) = t(6) + atmass * (x2+y2)
   enddo

   call dspev('N','U',3,t,moments,work,3,work,info)
   if (info.ne.0) moments = -1.0_wp

end function moments_of_inertia

!> aligns the molecule to its principal axes of inertia
pure subroutine align_to_principal_axes(self,break_symmetry)
   use iso_fortran_env, wp => real64
   use pbc_tools
   implicit none
   class(molecule),intent(inout) :: self !< molecular structure information
   logical, optional, intent(in)    :: break_symmetry
   real(wp) :: moments(3),det
   real(wp) :: center(3),atmass,t(6),work(9),axes(3,3)
   real(wp) :: x,x2,y,y2,z,z2
   integer  :: i,iat,info
   ! currently not supported
   if (self%npbc > 0) then
      return
   endif
   center = self%center_of_mass()
   t = 0.0_wp
   if (present(break_symmetry)) then
      if (break_symmetry) t = [(real(i,wp)*1.0e-10_wp,i=1,6)]
   endif

   do iat = 1, self%n
      atmass = self%atmass(iat)
      x = self%xyz(1,iat)-center(1); x2 = x**2
      y = self%xyz(2,iat)-center(2); y2 = y**2
      z = self%xyz(3,iat)-center(3); z2 = z**2
      t(1) = t(1) + atmass * (y2+z2)
      t(2) = t(2) - atmass * x*y
      t(3) = t(3) + atmass * (x2+z2)
      t(4) = t(4) - atmass * x*z
      t(5) = t(5) - atmass * y*z
      t(6) = t(6) + atmass * (x2+y2)
   enddo

   call dspev('V','U',3,t,moments,axes,3,work,info)
   if (info.ne.0) return

   det = mat_det_3x3(axes)
   if (det < 0) axes(:,1) = -axes(:,1)

   call coord_trafo(self%n,axes,self%xyz)

end subroutine align_to_principal_axes

end module class_molecule
