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

!> implements readers for different periodic geometry types
!
!  currently Xmol/xyz format, Turbomole coord format and VASP POSCAR are supported
module geometry_reader
   use iso_fortran_env, wp => real64
   implicit none

   public :: read_geometry

   public :: read_coord
   public :: read_poscar
   public :: read_xmol
   public :: read_sdf
!   public :: read_pdb

   private

   integer,private,parameter :: p_str_length = 48
   integer,private,parameter :: p_arg_length = 24

contains

!> interface to read in a geometry from a file
subroutine read_geometry(fname,mol,env)
   use iso_fortran_env, wp => real64
   use mctc_environment
   use mctc_systools, only : getline
   use class_molecule
   implicit none
   type(mctc_logger),intent(inout) :: env
   character(len=*), intent(in)    :: fname !< geometry file name
   type(molecule),intent(inout) :: mol   !< molecular structure information
   character(len=:),allocatable :: line, tmp
   integer :: ifile
   integer :: err
   integer :: idum
   integer :: idot
   logical :: exist


   inquire(file=fname,exist=exist)
   if (exist) then
      open(newunit=ifile,file=fname)
   else
      call env%error(4,"Could not open '"//fname//"'")
      return
   endif

   ! find real file name for format identification
   idum = index(fname,'/',back=.true.)
   if (idum > 0) then
      tmp = fname(idum+1:)
   else
      tmp = fname
   endif
   idot = index(tmp,'.',back=.true.)
   if (idot > 0) tmp = tmp(idot+1:)

   ! get first non-empty line
   do
      call getline(ifile,line,err)
      if ((len(line) > 0).or.err /= 0) exit
   enddo
   rewind(ifile)

   ! check if first line is an integer
   read(line,*,iostat=err) idum
   if ((err.eq.0 .and. idum > 0) &
      &.or. index(tmp,'xyz') > 0) then
      call read_xmol(ifile,mol,env)
   else if ((index(line,'$') > 0) &
      &.or. index(tmp,'coord') > 0) then
      call read_coord(ifile,mol,env)
   else if (index(tmp,'POSCAR') > 0 &
      &.or. index(tmp,'CONTCAR') > 0 &
      &.or. index(tmp,'poscar') > 0 &
      &.or. index(tmp,'contcar') > 0 &
      &.or. index(tmp,'VASP') > 0 &
      &.or. index(tmp,'vasp')   > 0) then
      call read_poscar(ifile,mol,env)
   else if (index(tmp,'sdf') > 0 &
      &.or. index(tmp,'mol') > 0) then
      call read_sdf(ifile,mol,env)
!   else if (index(tmp,'pdb') > 0) then
!      call read_pdb(ifile,mol)
!      if (present(filetype)) filetype = 'pdb'
   else
      ! format not recognized => assume Turbomole
      call read_coord(ifile,mol,env)
   endif

   close(ifile)
end subroutine read_geometry

subroutine read_sdf(iunit,mol,env)
   use iso_fortran_env, wp => real64
   use mctc_econv
   use class_molecule
   use mctc_environment
   use mctc_readin, only : getline => strip_line
   implicit none
   type(mctc_logger),intent(inout) :: env
   logical, parameter :: debug = .false.
   integer,intent(in) :: iunit !< file handle
   type(molecule),intent(inout) :: mol
   integer  :: n, i, iat, ichrg
   real(wp) :: xyz(3)
   real(wp) :: conv

   character(len=:),allocatable :: line
   character(len=10) :: chdum
   integer  :: err

   conv = aatoau

   rewind(iunit)

   read(iunit,'(a)')!skip
   read(iunit,'(a)')!skip
   read(iunit,'(a)')!skip
   read(iunit,'(i3)',iostat=err) n
   if (err.ne.0) then
      call env%error(4,"Could not read number of atoms, check format!")
      return
   endif

   if (n.lt.1) then
      call env%error(4,'Found no atoms, cannot work without atoms!')
      return
   endif

   call mol%allocate(n)
   mol%npbc = 0 ! SDF is always molecular
   mol%pbc = .false.

   do i = 1, n
      call getline(iunit,line,err)
      if (err.ne.0) then
         call env%error(4,"Could not read geometry from SDF file")
         return
      endif
      if (debug) print'(">",a)',line
      read(line,*,iostat=err) chdum, xyz(1), xyz(2), xyz(3)
      if (debug) print'("->",a)',chdum
      if (debug) print'("->",3g0)',xyz

      call elem(line,iat)
      if (debug) print'("->",g0)',iat
      if (iat > 0) then
         mol%at(i) = iat
         mol%sym(i) = line(1:2)
         mol%xyz(:,i) = xyz*conv
      else
         call env%error(4,"Unknown atom kind in SDF file")
         return
      endif
   enddo

   ! search for charge modifiers
   do
      call getline(iunit,line,err)
      if (is_iostat_end(err)) exit
      if (index(line,'$$$$') > 0) exit
      if (index(line,'M  CHG') > 0) then
         read(line(7:),*,iostat=err) ichrg
         if (err.eq.0) mol%chrg = mol%chrg + ichrg
      endif
   enddo

end subroutine read_sdf

subroutine read_xmol(iunit,mol,env)
   use iso_fortran_env, wp => real64
   use mctc_econv
   use mctc_environment
   use class_molecule
   use pbc_tools
   use mctc_readin, only : getline => strip_line
   implicit none
   type(mctc_logger),intent(inout) :: env
   logical, parameter :: debug = .false.
   integer,intent(in) :: iunit !< file handle
   type(molecule),intent(inout) :: mol
   integer  :: n, iat
   real(wp) :: xyz(3)
   real(wp) :: conv

   character(len=:),allocatable :: line
   character(len=10) :: chdum
   integer  :: err

   conv = aatoau

   rewind(iunit)

   read(iunit,*,iostat=err) n
   if (err.ne.0) then
      call env%error(4,"Could not read number of atoms, check format!")
      return
   endif

   if (n.lt.1) then
      call env%error(4,'Found no atoms, cannot work without atoms!')
      return
   endif

   call mol%allocate(n)
   mol%npbc = 0 ! Xmol is always molecular (there are extensions to this...)
   mol%pbc = .false.

   ! drop next record
   read(iunit,'(a)')

   n = 0
   do while (n < mol%n)
      call getline(iunit,line,err)
      if (is_iostat_end(err)) exit
      if (err.ne.0) then
         call env%error(4,"Could not read geometry from Xmol file")
         return
      endif
      if (debug) print'(">",a)',line
      read(line,*,iostat=err) chdum, xyz(1), xyz(2), xyz(3)
      if (debug) print'("->",a)',chdum
      if (debug) print'("->",3g0)',xyz

      call elem(line,iat)
      if (debug) print'("->",g0)',iat
      if (iat > 0) then
         n = n+1
         mol%at(n) = iat
         mol%sym(n) = line(1:2)
         mol%xyz(:,n) = xyz*conv
      endif
   enddo

   if (n.ne.mol%n) then
      call env%error(4,"Atom number missmatch in Xmol file")
      return
   endif

end subroutine read_xmol

!> read geometry as POSCAR from iunit
subroutine read_poscar(iunit,mol,env)
   use iso_fortran_env, wp => real64
   use mctc_environment
   use class_molecule
   use pbc_tools
   implicit none
   type(mctc_logger),intent(inout) :: env
   logical, parameter :: debug = .false.
   integer,intent(in) :: iunit !< file handle
   type(molecule),intent(inout) :: mol
   integer :: n

   call get_atomnumber(iunit,n)

   if (n.lt.1) then
      call env%error(4,'Found no atoms, cannot work without atoms!')
      return
   endif

   call mol%allocate(n)
   mol%npbc = 3 ! VASP is always 3D
   mol%pbc = .true.

   call get_coord(iunit,mol%lattice,mol%n,mol%xyz,mol%at,mol%sym)
   call dlat_to_cell(mol%lattice,mol%cellpar)
   call dlat_to_rlat(mol%lattice,mol%rec_lat)
   mol%volume = dlat_to_dvol(mol%lattice)

   call xyz_to_abc(mol%n,mol%lattice,mol%xyz,mol%abc,mol%pbc)

contains

!> read the coordinates from POSCAR
subroutine get_coord(iunit,lattice,n,xyz,at,sym)
   use mctc_econv
   use mctc_strings
   use mctc_systools

   implicit none

   integer, intent(in)  :: n
   real(wp),intent(out) :: xyz(3,n)
   real(wp),intent(out) :: lattice(3,3)
   integer, intent(out) :: at(n)
   character(len=2),intent(out) :: sym(n)
   integer, intent(in)  :: iunit
   logical              :: selective=.false. ! Selective dynamics
   logical              :: cartesian=.true.  ! Cartesian or direct

   real(wp) :: ddum,latvec(3)
   real(wp) xx(10),scalar
   real(wp) :: coord(3)
   character(len=:),allocatable :: line
   character(len=80) :: args(90),args2(90)

   integer i,j,nn,ntype,ntype2,atnum,i_dummy1,i_dummy2,ncheck

   integer :: iat, inum, idum, err

   lattice=0

   rewind(iunit)
   ncheck=0
   ntype=0
   ! first line contains the symbols of different atom types
   call getline(iunit,line,err)
   if (err.ne.0) then
      call env%error(4,"Could not read POSCAR")
      return
   endif
   if (debug) print'(">",a)',line
   call parse(line,' ',args,ntype)

   ! this line contains the global scaling factor,
   call getline(iunit,line,err)
   if (err.ne.0) then
      call env%error(4,"Could not read POSCAR")
      return
   endif
   if (debug) print'(">",a)',line
   read(line,*,iostat=err) ddum
   if (err.ne.0) then
      call env%error(4,"Could not read POSCAR")
      return
   endif
   ! the Ang->au conversion is included in the scaling factor
   scalar = ddum*aatoau

   ! reading the lattice constants
   do i=1,3
      call getline(iunit,line,err)
      if (err.ne.0) then
         call env%error(4,"Could not read lattice from POSCAR")
         return
      endif
      if (debug) print'("->",a)',line
      read(line,*,iostat=err) latvec
      if (err.ne.0) then
         call env%error(4,"Could not read lattice from POSCAR")
         return
      endif
      lattice(:,i) = latvec * scalar
   enddo
   ! Either here are the numbers of each element,
   ! or (>vasp.5.1) here are the element symbols
   call getline(iunit,line,err)
   if (err.ne.0) then
      call env%error(4,"Could not read POSCAR")
      return
   endif
   if (debug) print'(">",a)',line
   ! try to read first element
   read(line,*,iostat=err) idum
   ! CONTCAR files have additional Element line here since vasp.5.1
   if (err.ne.0) then
      call parse(line,' ',args,ntype)
      call getline(iunit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) then
         call env%error(4,"Could not read POSCAR")
         return
      endif
   endif
   call parse(line,' ',args2,nn)
   if (nn.ne.ntype) then
      call env%error(4, 'Error reading number of atomtypes')
      return
   endif
   ncheck=0
   do i=1,nn
      read(args2(i),*,iostat=err) inum
      call elem(args(i),iat)
      if (iat < 1 .or. inum < 1) then
         call env%error(4, 'Error: unknown element.')
         return
      endif
      do j=1,inum
         ncheck=ncheck+1
         sym(ncheck) = args(i)(1:2)
         at(ncheck)=iat
      enddo
   enddo
   if (n.ne.ncheck) then
      call env%error(4,'Error reading Number of Atoms')
      return
   endif

   call getline(iunit,line,err)
   if (err.ne.0) then
      call env%error(4,"Could not read POSCAR")
      return
   endif
   if (debug) print'(">",a)',line
   line=adjustl(line)
   if (line(:1).eq.'s' .or. line(:1).eq.'S') then
      selective=.true.
      call getline(iunit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) then
         call env%error(4,"Could not read POSCAR")
         return
      endif
      line=adjustl(line)
   endif

   cartesian=(line(:1).eq.'c' .or. line(:1).eq.'C' .or. &
      &       line(:1).eq.'k' .or. line(:1).eq.'K')
   do i=1,n
      call getline(iunit,line,err)
      if (err.ne.0) then
         call env%error(4,"Could not read geometry from POSCAR")
         return
      endif
      if (debug) print'("-->",a)',line
      read(line,*,iostat=err) coord
      if (err.ne.0) then
         call env%error(4,"Could not read geometry from POSCAR")
         return
      endif

      if (cartesian) then
         xyz(:,i)=coord*scalar
      else
         xyz(:,i)=matmul(lattice,coord)
      endif

   enddo

end subroutine get_coord

!> read the coordinates from POSCAR
subroutine get_atomnumber(iunit,n)
   use mctc_econv
   use mctc_strings
   use mctc_systools

   implicit none

   integer, intent(out) :: n
   integer, intent(in)  :: iunit
   logical              :: selective=.false. ! Selective dynamics
   logical              :: cartesian=.true.  ! Cartesian or direct

   real(wp) :: ddum,latvec(3)
   real(wp) xx(10),scalar
   real(wp) :: coord(3)
   character(len=:),allocatable :: line
   character(len=80) :: args(90),args2(90)

   integer i,j,nn,ntype,ntype2,atnum

   integer :: iat, inum, idum, err

   rewind(iunit)
   ntype=0
   ! first line contains the symbols of different atom types
   call getline(iunit,line,err)
   if (err.ne.0) then
      call env%error(4,"Could not read POSCAR")
      return
   endif
   if (debug) print'(">",a)',line
   call parse(line,' ',args,ntype)

   ! this line contains the global scaling factor,
   call getline(iunit,line,err)
   if (err.ne.0) then
      call env%error(4,"Could not read POSCAR")
      return
   endif
   if (debug) print'(">",a)',line
   read(line,*,iostat=err) ddum
   if (err.ne.0) then
      call env%error(4,"Could not read POSCAR")
      return
   endif
   ! the Ang->au conversion is included in the scaling factor
   scalar = ddum*aatoau

   ! reading the lattice constants
   do i=1,3
      call getline(iunit,line,err)
      if (err.ne.0) then
         call env%error(4,"Could not read lattice from POSCAR")
         return
      endif
      if (debug) print'("->",a)',line
   enddo
   ! Either here are the numbers of each element,
   ! or (>vasp.5.1) here are the element symbols
   call getline(iunit,line,err)
   if (err.ne.0) then
      call env%error(4,"Could not read POSCAR")
      return
   endif
   if (debug) print'(">",a)',line
   ! try to read first element
   read(line,*,iostat=err) idum
   ! CONTCAR files have additional Element line here since vasp.5.1
   if (err.ne.0) then
      call parse(line,' ',args,ntype)
      call getline(iunit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) then
         call env%error(4,"Could not read POSCAR")
         return
      endif
   endif
   call parse(line,' ',args2,nn)
   if (nn.ne.ntype) then
      call env%error(4, 'Error reading number of atomtypes')
      return
   endif
   n=0
   do i=1,nn
      read(args2(i),*,iostat=err) inum
      call elem(args(i),iat)
      if (iat < 1 .or. inum < 1) then
         call env%error(4, 'Error: unknown element.')
         return
      endif
      do j=1,inum
         n=n+1
      enddo
   enddo

end subroutine get_atomnumber

end subroutine read_poscar

! ------------------------------------------------------------------------
! Turbomole's coord file for riper
! ------------------------------------------------------------------------
!> read geometry as coord from iunit
subroutine read_coord(iunit,mol,env)
   use iso_fortran_env, wp => real64
   use class_molecule
   use mctc_readin, getline => strip_line
   use mctc_environment
   use pbc_tools
   implicit none
   type(mctc_logger),intent(inout) :: env
   integer,intent(in) :: iunit !< file handle
   type(molecule),intent(inout) :: mol
   character(len=:),allocatable :: line
   integer :: err
   integer :: idum
   logical,parameter :: debug = .false.
   integer :: natoms
   integer :: i
   integer :: npbc
   logical :: exist

   character,parameter :: flag = '$'
   character(len=*),parameter :: flag_end = flag//'end'

   exist  = .false.
   npbc   = -1
   natoms = -1

!  we need to read this file twice, for a lot of reasons
   rewind(iunit)

!  read first line before the readloop starts, I have to do this
!  to avoid using backspace on iunit (dammit Turbomole format)
   if (debug) print'("first scan")'
   call getline(iunit,line,err)
   scanfile: do
      !  check if there is a $ in the *first* column
      if (index(line,flag).eq.1) then
         if (index(line,'coord').eq.2) then
            if (debug) print'("*>",a)',line
            call get_atomnumber(natoms,iunit,line,err)
         elseif (index(line,'periodic').eq.2) then
            if (debug) print'("*>",a)',line
            call get_periodic(npbc,iunit,line,err)
         else
!           get a new line
            if (debug) print'(" >",a)',line
            call getline(iunit,line,err)
         endif
      else ! not a keyword -> ignore
         if (debug) print'(" >",a)',line
         call getline(iunit,line,err)
      endif
   !  check for end of file, which I will tolerate as alternative to $end
      if (err.eq.iostat_end) exit scanfile
      if (index(line,flag_end).ne.0) exit scanfile
   enddo scanfile

   if (natoms.eq.-1) then
      call env%error(4,"coord keyword was not found")
      return
   endif
   if (npbc.eq.-1) npbc = 0

   if (natoms.lt.1) then ! this will catch an empty coord data group
      call env%error(4,'Found no atoms, cannot work without atoms!')
      return
   endif

   call mol%allocate(natoms)
   mol%npbc = npbc
   do i = 1, npbc
      mol%pbc(i) = .true.
   enddo
   if (mol%npbc == 0) then
      ! this might be wrong
      mol%lattice = reshape([0.0_wp,0.0_wp,0.0_wp,&
         &                   0.0_wp,0.0_wp,0.0_wp,&
         &                   0.0_wp,0.0_wp,0.0_wp],[3,3])
      mol%rec_lat = reshape([1.0_wp,0.0_wp,0.0_wp,&
         &                   0.0_wp,1.0_wp,0.0_wp,&
         &                   0.0_wp,0.0_wp,1.0_wp],[3,3])
   endif

   rewind(iunit)

!  read first line before the readloop starts, I have to do this
!  to avoid using backspace on iunit (dammit Turbomole format)
   if (debug) print'("one''n''half read")'
   if (npbc > 0) then
   call getline(iunit,line,err)
   readlat: do
      !  check if there is a $ in the *first* column
      if (index(line,flag).eq.1) then
         if (index(line,'lattice').eq.2) then
            if (exist) &
               call env%warning(4,"Multiple definitions of lattice present!")
            exist = .true.
            if (debug) print'("*>",a)',line
            call get_lattice(iunit,line,err,npbc,mol%lattice)
            call dlat_to_cell(mol%lattice,mol%cellpar)
            call dlat_to_rlat(mol%lattice,mol%rec_lat)
            mol%volume = dlat_to_dvol(mol%lattice)
         elseif (index(line,'cell').eq.2) then
            if (exist) &
               call env%warning(4,"Multiple definitions of lattice present!")
            exist = .true.
            if (debug) print'("*>",a)',line
            call get_cell(iunit,line,err,npbc,mol%cellpar)
            call cell_to_dlat(mol%cellpar,mol%lattice)
            call cell_to_rlat(mol%cellpar,mol%rec_lat)
            mol%volume = cell_to_dvol(mol%cellpar)
         else
!           get a new line
            if (debug) print'(" >",a)',line
            call getline(iunit,line,err)
         endif
      else ! not a keyword -> ignore
         if (debug) print'(" >",a)',line
         call getline(iunit,line,err)
      endif
   !  check for end of file, which I will tolerate as alternative to $end
      if (err.eq.iostat_end) exit readlat
      if (index(line,flag_end).ne.0) exit readlat
   enddo readlat

   rewind(iunit)
   endif

!  read first line before the readloop starts, I have to do this
!  to avoid using backspace on iunit (dammit Turbomole format)
   if (debug) print'("second read")'
   call getline(iunit,line,err)
   readfile: do
      !  check if there is a $ in the *first* column
      if (index(line,flag).eq.1) then
         if (index(line,'coord').eq.2) then
            if (debug) print'("*>",a)',line
            call get_coord(iunit,line,err,mol)
         else
!           get a new line
            if (debug) print'(" >",a)',line
            call getline(iunit,line,err)
         endif
      else ! not a keyword -> ignore
         if (debug) print'(" >",a)',line
         call getline(iunit,line,err)
      endif
   !  check for end of file, which I will tolerate as alternative to $end
      if (err.eq.iostat_end) exit readfile
      if (index(line,flag_end).ne.0) exit readfile
   enddo readfile

   if (.not.exist .and. npbc > 0) then
      call env%error(4,"There is no definition of the lattice present!")
      return
   endif

contains

!> count the number of lines in the coord block -> number of atoms
subroutine get_atomnumber(ncount,iunit,line,err)
   implicit none
   integer,intent(in) :: iunit          !< file handle
   character(len=:),allocatable :: line !< string buffer
   integer,intent(out) :: ncount        !< number of readin lines
   integer,intent(out) :: err           !< status of iunit
   ncount = 0
   do
      call getline(iunit,line,err)
      if (err.eq.iostat_end) return
      if (index(line,flag).ne.0) return
      if (debug) print'(" ->",a)',line

      if (line.eq.'') cycle ! skip empty lines
      ncount = ncount + 1   ! but count non-empty lines first
   enddo
end subroutine get_atomnumber

!> get the dimensionality of the system
subroutine get_periodic(npbc,iunit,line,err)
   implicit none
   integer,intent(in) :: iunit          !< file handle
   character(len=:),allocatable :: line !< string buffer
   integer,intent(out) :: npbc          !< number of periodic dimensions
   integer,intent(out) :: err           !< status of iunit
   integer :: idum
   if (get_value(line(10:),idum)) then
      if (idum.eq.0 .or. idum.eq.3) then
         npbc = idum
      else
         call env%error(4,"This kind of periodicity is not implemented")
         return
      endif
   else
      call env%error(4,"Could not read periodicity from '"//line//"'")
      return
   endif
   call getline(iunit,line,err)
end subroutine get_periodic

!> read the lattice vectors from the data group
subroutine get_lattice(iunit,line,err,npbc,lat)
   use mctc_econv
   implicit none
   integer,intent(in) :: iunit          !< file handle
   character(len=:),allocatable :: line !< string buffer
   integer,intent(out) :: err           !< status of iunit
   integer,intent(in)  :: npbc          !< number of periodic dimensions
   real(wp),intent(out) :: lat(3,3)     !< lattice vectors
   integer :: ncount
   real(wp) :: lvec(3)
   real(wp) :: conv
   if (npbc == 0) then
      call env%warning(4,"lattice data group is present but not periodic?")
      return
   endif
   if (len(line).le.8) then
      conv = 1.0_wp ! default is bohr
   elseif (index(line(9:),'bohr').gt.0) then
      conv = 1.0_wp
   elseif (index(line(9:),'angs').gt.0) then
      conv = aatoau
   else ! fall back to default
      call env%warning(4,"Could not make sense out of unit '"//line(7:)//"'")
      conv = 1.0_wp
   endif
   lat = 0.0_wp
   ncount = 0
   do
      call getline(iunit,line,err)
      if (err.eq.iostat_end) return
      if (index(line,flag).ne.0) return
      if (debug) print'(" ->",a)',line

      if (line.eq.'') cycle ! skip empty lines
      ncount = ncount + 1   ! but count non-empty lines first
      if (ncount.gt.npbc) then
         call env%error(4,"Input error in lattice data group")
         return
      endif
      read(line,*,iostat=err) lvec(1:npbc)
      if (err.ne.0) then
         call env%error(4,"Input error in lattice data group")
         return
      endif
      lat(1:npbc,ncount) = lvec(1:npbc)*conv
   enddo
end subroutine get_lattice

!> read the cell parameters and transform to lattice
subroutine get_cell(iunit,line,err,npbc,cellpar)
   use mctc_constants
   use mctc_econv
   implicit none
   integer,intent(in)   :: iunit        !< file handle
   character(len=:),allocatable :: line !< string buffer
   integer,intent(out)  :: err          !< status of iunit
   integer,intent(in)   :: npbc         !< number of periodic dimensions
   real(wp),intent(out) :: cellpar(6)   !< cell parameter
   real(wp) :: conv,vol
   if (npbc == 0) then
      call env%warning(4,"cell data group is present but not periodic?")
      return
   endif

   associate( alen => cellpar(1), blen => cellpar(2), clen => cellpar(3), &
              alp  => cellpar(4), bet  => cellpar(5), gam  => cellpar(6) )

   alen = 1.0_wp; blen = 1.0_wp; clen = 1.0_wp
   alp = 90.0_wp; bet = 90.0_wp; gam = 90.0_wp

   if (len(line).le.5) then
      conv = 1.0_wp ! default is bohr
   elseif (index(line(6:),'bohr').gt.0) then
      conv = 1.0_wp
   elseif (index(line(6:),'angs').gt.0) then
      conv = aatoau
   else ! fall back to default
      call env%warning(4,"Could not make sense out of unit '"//line(7:)//"'")
      conv = 1.0_wp
   endif
   call getline(iunit,line,err)
   if (err.eq.iostat_end) return
   if (index(line,flag).ne.0) return
   if (debug) print'(" ->",a)',line
   if (npbc == 3) &
      read(line,*,iostat=err) alen,blen,clen,alp,bet,gam
   if (npbc == 2) &
      read(line,*,iostat=err) alen,blen,             gam
   if (npbc == 1) &
      read(line,*,iostat=err) alen
   if (err.ne.0) then
      call env%error(4,"Could not read cell from '"//line//"'")
      return
   endif

   if(alen < 0.0d0 .or. blen < 0.0d0 .or. clen < 0.0d0 .or. &
      alp < 0.0d0 .or. bet < 0.0d0 .or. gam < 0.0d0) then
      call env%error(4,"Negative cell parameters make no sense!")
      return
   endif

   alen = alen*conv
   blen = blen*conv
   clen = clen*conv
   alp  = alp*pi/180.0_wp
   bet  = bet*pi/180.0_wp
   gam  = gam*pi/180.0_wp

   end associate

end subroutine get_cell

!> read the coordinates from coord data group
subroutine get_coord(iunit,line,err,mol)
   use mctc_econv
   use mctc_strings
   use class_molecule
   implicit none
   integer,intent(in) :: iunit            !< file handle
   character(len=:),allocatable :: line   !< string buffer
   integer,intent(out) :: err             !< status of iunit
   type(molecule),intent(inout) :: mol !< molecular structure information
   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv
   integer  :: ncount,iat,icoord
   real(wp) :: conv,ddum,xyz(3)
   logical  :: frac
   if (len(line).le.6) then
      frac = .false.
      conv = 1.0_wp ! default is bohr
   elseif (index(line(7:),'bohr').gt.0) then
      frac = .false.
      conv = 1.0_wp
   elseif (index(line(7:),'angs').gt.0) then
      frac = .false.
      conv = aatoau
   elseif (index(line(7:),'frac').gt.0) then
      frac = .true.
      conv = 1.0_wp
   else ! fall back to default
      call env%warning(4,"Could not make sense out of unit '"//line(7:)//"'")
      frac = .false.
      conv = 1.0_wp
   endif

   ncount = 0
   do
      call getline(iunit,line,err)
      if (err.eq.iostat_end) return
      if (index(line,flag).ne.0) return
      if (debug) print'(" ->",a)',line

      if (line.eq.'') cycle ! skip empty lines
      ncount = ncount + 1   ! but count non-empty lines first
      if (ncount.gt.mol%n) then
         call env%error(4,"Internal error while reading coord data group")
         return
      endif
      call parse(line,' ',argv,narg)
      if (narg.lt.4) then
         call env%error(4,"not enough arguments in line in coord data group")
         return
      endif
      call elem(argv(4),iat)
      if (iat.eq.0) then
         call env%error(4,"invalid element input in line in coord data group")
         return
      endif
      mol%sym(ncount) = trim(argv(4))
      mol%at(ncount) = iat
      do icoord = 1, 3
         if (get_value(trim(argv(icoord)),ddum)) then
            xyz(icoord) = conv*ddum
         else
            call env%error(4,"invalid coordinate input in line in coord data group")
            return
         endif
      enddo
      ! in case of fractional coordinates we perform a gemv with the lattice
      if (frac) then
         mol%xyz(:,ncount) = matmul(mol%lattice,xyz)
      else
         mol%xyz(:,ncount) = xyz
      endif
   enddo
end subroutine get_coord

end subroutine read_coord

! ------------------------------------------------------------------------
! cyrstal's unit 34
! ------------------------------------------------------------------------
subroutine read_fort34(iunit,mol,env)
   use iso_fortran_env, wp => real64
   use mctc_environment
   use class_molecule
   implicit none
   type(mctc_logger),intent(inout) :: env
   integer,intent(in) :: iunit !< file handle
   type(molecule),intent(inout) :: mol
   call env%error(4,"Currently not implemented, we are very sorry -- your dev team")
end subroutine read_fort34

subroutine elem(key1, nat)
   implicit double precision (a-h,o-z), integer(i-n)
   character(len=*) :: key1
   character(len=2) :: e
   character(len=2),parameter :: elemnt(118) = (/ &
      & 'h ','he', &
      & 'li','be','b ','c ','n ','o ','f ','ne', &
      & 'na','mg','al','si','p ','s ','cl','ar', &
      & 'k ','ca', &
      & 'sc','ti','v ','cr','mn','fe','co','ni','cu','zn', &
      &           'ga','ge','as','se','br','kr', &
      & 'rb','sr', &
      & 'y ','zr','nb','mo','tc','ru','rh','pd','ag','cd', &
      &           'in','sn','sb','te','i ','xe', &
      & 'cs','ba','la', &
      & 'ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb', &
      & 'lu','hf','ta','w ','re','os','ir','pt','au','hg', &
      &           'tl','pb','bi','po','at','rn', &
      & 'fr','ra','ac', &
      & 'th','pa','u ','np','pu','am','cm','bk','cf','es','fm','md','no', &
      & 'lr','rf','db','sg','bh','hs','mt','ds','rg','cn', &
      &           'nh','fl','mc','lv','ts','og' /)

   nat=0
   e='  '
   do i=1,len(key1)
      if(key1(i:i).ne.' ') L=i
   enddo
   k=1
   do j=1,l           
      if (k.gt.2)exit
      n=ichar(key1(J:J))
      if(len_trim(e).ge.1 .and. n.eq.ichar(' '))exit!break if space after elem. symbol
      if(len_trim(e).ge.1 .and. n.eq.9)exit !break if tab after elem. symbol
      if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
         e(k:k)=char(n+ichar('a')-ichar('A'))
         k=k+1
      endif
      if(n.ge.ichar('a') .and. n.le.ichar('z') )then
         e(k:k)=key1(j:j)
         k=k+1
      endif
   enddo

   do i=1,118
      if(e.eq.elemnt(i))then
         nat=i
         return
      endif
   enddo
end subroutine elem

end module geometry_reader
