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

subroutine out_tmer(fname,energy)
   use iso_fortran_env, wp => real64, istdout => output_unit
   implicit none
   character(len=*),intent(in) :: fname
   real(wp),intent(in)  :: energy
   integer :: ich ! file handle
   write(istdout,'(a)') "dispersion energy written to file "//fname
   open(newunit=ich,file=fname)
   write(ich,'(f30.20)') energy
   close(ich)
end subroutine out_tmer

subroutine write_gradient(mol,iunit,gradient)
   use iso_fortran_env, wp => real64, istdout => output_unit
   use class_molecule
   implicit none
   type(molecule),intent(in) :: mol
   integer, intent(in)  :: iunit !< file handle
   real(wp),intent(in)  :: gradient(3,mol%n)
   integer :: i
   do i = 1, mol%n
      write(iunit,'(3D22.13)') gradient(:,i)
   enddo
end subroutine write_gradient

subroutine out_gradlatt(mol,fname,energy,gradlatt,env)
   use iso_fortran_env, wp => real64
   use mctc_systools
   use mctc_environment
   use class_molecule
   implicit none
   type(mctc_logger),intent(inout) :: env
   type(molecule),intent(in)    :: mol
   character(len=*),intent(in)  :: fname
   real(wp),intent(in)          :: energy
   real(wp),intent(in)          :: gradlatt(3,3)
   character(len=:),allocatable :: line
   integer  :: i,icycle,line_number
   integer  :: err
   integer  :: igrad ! file handle
   logical  :: exist
   real(wp) :: escf
   real(wp) :: glat(3,3)
   real(wp) :: dlat(3,3)
   icycle = 1
   i = 0
   escf = 0.0_wp

   inquire(file=fname,exist=exist)
   if (exist) then
      open(newunit=igrad,file=fname)
      read_file: do
         call getline(igrad,line,iostat=err)
         if (err.ne.0) exit read_file
         i=i+1
         if (index(line,'cycle') > 0) line_number = i
      enddo read_file
      if (line_number < 2) then
         call env%warning(9,'Illegal gradient file to add dispersion gradient')
         return
      endif

      rewind(igrad)
      skip_lines: do i = 1, line_number-1
         read(igrad,'(a)')
      enddo skip_lines
      call getline(igrad,line)
      read(line(10:17),*,iostat=err) icycle
      read(line(33:51),*,iostat=err) escf

      do i = 1, 3
         call getline(igrad,line)
         read(line,*,iostat=err) dlat(1,i),dlat(2,i),dlat(3,i)
      enddo
      if (any(abs(dlat-mol%lattice) > 1.0e-8_wp)) then
         call env%error(9,'Lattice in gradlatt does not match actual lattice')
      endif
      do i = 1, 3
         call getline(igrad,line)
         read(line,*,iostat=err) glat(1,i),glat(2,i),glat(3,i)
      enddo
      do i = 1, 3
         backspace(igrad)
         backspace(igrad)
      enddo
      backspace(igrad)
   else
      open(newunit=igrad,file=fname)
      write(igrad,'("$gradlatt")')
   endif

   write(igrad,'(2x,"cycle =",1x,i6,4x,"SCF energy =",f18.11,3x,'//&
                   '"|dE/dlatt| =",f10.6)') &
      icycle, energy+escf, norm2(gradlatt+glat)
   do i = 1, 3
      write(igrad,'(3(F20.14,2x))') mol%lattice(1,i),mol%lattice(2,i),mol%lattice(3,i)
   enddo
   do i = 1, 3
      write(igrad,'(3D22.13)') gradlatt(1,i)+glat(1,i),gradlatt(2,i)+glat(2,i),gradlatt(3,i)+glat(3,i)
   enddo
   write(igrad,'("$end")')
   close(igrad)

end subroutine out_gradlatt


subroutine out_gradient(mol,fname,energy,gradient,env)
   use iso_fortran_env, wp => real64
   use mctc_environment
   use mctc_systools
   use class_molecule
   implicit none
   type(mctc_logger),intent(inout) :: env
   type(molecule),intent(in)    :: mol
   character(len=*),intent(in)  :: fname
   real(wp),intent(in)          :: energy
   real(wp),intent(in)          :: gradient(3,mol%n)
   character(len=:),allocatable :: line
   integer  :: i,icycle,line_number
   integer  :: err
   integer  :: igrad ! file handle
   logical  :: exist
   real(wp) :: escf
   real(wp),allocatable :: gscf(:,:)
   real(wp),allocatable :: xyz (:,:)
   allocate( gscf(3,mol%n), source = 0.0_wp )
   icycle = 1
   i = 0
   escf = 0.0_wp

   inquire(file=fname,exist=exist)
   if (exist) then
      open(newunit=igrad,file=fname)
      read_file: do
         call getline(igrad,line,iostat=err)
         if (err.ne.0) exit read_file
         i=i+1
         if (index(line,'cycle') > 0) line_number = i
      enddo read_file
      if (line_number < 2) then
         call env%warning(9,'Illegal gradient file to add dispersion gradient')
         return
      endif

      rewind(igrad)
      skip_lines: do i = 1, line_number-1
         read(igrad,'(a)')
      enddo skip_lines
      call getline(igrad,line)
      read(line(10:17),*,iostat=err) icycle
      read(line(33:51),*,iostat=err) escf

      allocate(xyz(3,mol%n))
      do i = 1, mol%n
         call getline(igrad,line)
         read(line,*,iostat=err) xyz(1,i),xyz(2,i),xyz(3,i)
      enddo
      if (any(abs(xyz-mol%xyz) > 1.0e-8_wp)) then
         call env%error(9,'Geometry in gradient does not match actual geometry')
      endif
      do i = 1, mol%n
         call getline(igrad,line)
         read(line,*,iostat=err) gscf(1,i),gscf(2,i),gscf(3,i)
      enddo
      do i = 1, mol%n
         backspace(igrad)
         backspace(igrad)
      enddo
      backspace(igrad)
   else
      open(newunit=igrad,file=fname)
      write(igrad,'("$grad")')
   endif

   write(igrad,'(2x,"cycle =",1x,i6,4x,"SCF energy =",f18.11,3x,'//&
                   '"|dE/dxyz| =",f10.6)') &
      icycle, energy+escf, norm2(gradient+gscf)
   do i = 1, mol%n
      write(igrad,'(3(F20.14,2x),4x,a2)') mol%xyz(1,i),mol%xyz(2,i),mol%xyz(3,i),mol%sym(i)
   enddo
   do i = 1, mol%n
      write(igrad,'(3D22.13)') gradient(1,i)+gscf(1,i),gradient(2,i)+gscf(2,i),gradient(3,i)+gscf(3,i)
   enddo
   write(igrad,'("$end")')
   close(igrad)

end subroutine out_gradient

subroutine orca_gradient(mol,iunit,gradient)
   use iso_fortran_env, wp => real64
   use class_molecule
   implicit none
   type(molecule),intent(in) :: mol
   integer, intent(in)  :: iunit !< file handle
   real(wp),intent(in)  :: gradient(3,mol%n)
   integer :: i
   write(iunit,'("# Gdisp fmt='' %22.14lf  %22.14lf  %22.14lf ''")')
   do i = 1, mol%n
      write(iunit,'(3(1x,f22.14,1x))') gradient(1,i),gradient(2,i),gradient(3,i)
   enddo
   write(iunit,'("#")')
end subroutine orca_gradient

subroutine orca_hessian(mol,iunit,hessian)
   use iso_fortran_env, wp => real64
   use class_molecule
   implicit none
   type(molecule),intent(in) :: mol
   integer, intent(in)  :: iunit !< file handle
   real(wp),intent(in)  :: hessian(3*mol%n,3*mol%n)
   integer :: i,j
   write(iunit,'("NOTE: Printing only the lower triangle")')
   write(iunit,'("# Hdisp fmt='' %22.14lf  %22.14lf  %22.14lf ''")')
   write(iunit,'(3(1x,f22.14,1x))') ((hessian(j,i),j=1,i),i=1,3*mol%n)
   write(iunit,'("#")')
end subroutine orca_hessian
