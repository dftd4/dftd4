!! ------------------------------------------------------------------------
!  Purpose:
!  read in a geometry from a provided file, can read XMOL and Turbomole
!
!  Input:
!     fname  - file name with geometry
!
!  Output:
!     mol    - molecular structure
!! ------------------------------------------------------------------------
subroutine get_geometry(mol,fname)
   use iso_fortran_env, only : wp => real64

!! ------------------------------------------------------------------------
!  general purpose library
!! ------------------------------------------------------------------------
   use mctc_systools
   use mctc_readin
   use mctc_econv

!! ------------------------------------------------------------------------
!  Class definition
!! ------------------------------------------------------------------------
   use class_molecule

   implicit none

   type(molecule) :: mol
   character(len=*),intent(in)  :: fname

   character(len=:),allocatable :: line
   real(wp)          :: floats(3)
   character(len=80) :: strings(3)
   integer  :: nfloats,nstrings
   real(wp) :: scalfactor
   integer  :: nat,iat
   integer  :: ich ! file handle
   integer  :: err
   logical  :: exist

   inquire(file=fname,exist=exist)
   if (.not.exist) call raise('E',"'"//fname//"' was not found for readin")

   open(newunit=ich,file=fname)

   ! skip blank lines
   do
      call getline(ich,line,iostat=err)
      if (err.ne.0) call raise('E',"'"//fname//"' does not contain anything")
      if (line.ne."") exit
   enddo

   call readline(line,floats,strings,nfloats,nstrings)

   if (nfloats.eq.1 .and. floats(1).gt.0) then
      call getline(ich,line,iostat=err)
      if (err.ne.0) call raise('E',"'"//fname//"' does not contain anything")
      nat = int(floats(1))
   else
      nat = 0
      count_lines: do
         call getline(ich,line,iostat=err)
         if(index(line,'$').ne.0) exit count_lines
         if(err.eq.iostat_end) exit count_lines
         call readline(line,floats,strings,nfloats,nstrings)
         if(nfloats.ne.3) cycle count_lines
         iat = elem(strings(1))
         if(iat.le.0) cycle count_lines
         nat=nat+1
      enddo count_lines
   endif

   call mol%allocate(nat,.false.)

   ! reread file
   rewind(ich)
   ! skip blank lines, again
   do
      call getline(ich,line,iostat=err)
      if (line.ne."") exit
   enddo

   call readline(line,floats,strings,nfloats,nstrings)

   if (nfloats.eq.1 .and. floats(1).gt.0) then
      scalfactor = aatoau
      call getline(ich,line,iostat=err)
   else if (index(line,'$coord').gt.0) then
      scalfactor = 1.0_wp
   else
      call raise('E','Coordinate format not recognized')
   endif
   nat = 0
   read_lines: do
      call getline(ich,line,iostat=err)
      if(index(line,'$').ne.0) exit read_lines
      if(err.eq.iostat_end) exit read_lines
      call readline(line,floats,strings,nfloats,nstrings)
      if(nfloats.ne.3) cycle read_lines
      iat = elem(strings(1))
      if(iat.le.0) cycle read_lines
      nat=nat+1
      if (nat.gt.mol%nat) &
         call raise('E','Number of atoms does match provided geometry')
      mol%xyz(:,nat)=floats(1:3)*scalfactor
      mol%sym(nat)=trim(strings(1))
      mol%at(nat)=iat
   enddo read_lines

contains

!! ------------------------------------------------------------------------
!  Purpose:
!  reads a line cuts the at blanks and tabstops
!  returns all floats and strings in order of occurence
!
!  Input:
!     line    - string to be parsed
!
!  Output:
!     floats   - floating point numbers found in line
!     strings  - strings found in line
!     nfloats  - number of floats found
!     nstrings - number of strings found
!! ------------------------------------------------------------------------
subroutine readline(line,floats,strings,nfloats,nstrings)
   use iso_fortran_env, wp => real64
   implicit none
   real(wp),intent(out)          :: floats(3)
   character(len=*),intent(in)   :: line
   character(len=80),intent(out) :: strings(3)
   integer,intent(out)           :: nstrings
   integer,intent(out)           :: nfloats

   real(wp) :: num
   character(len=80) :: stmp,str
   character(len=1)  :: digit
   integer  :: i,ty,cs,cf

   stmp=''
   cs=1
   cf=1
   strings=''
   do i=1,len_trim(line)
      digit=line(i:i)
      ! should exclude tabstops and blanks, 9 is ascii code for tab
      if(digit.ne.' '.and.digit.ne.char(9)) then
         stmp=trim(stmp)//trim(digit)
      elseif(stmp.ne.'')then
         ! get type of string, 0=number, 1=character
         call checktype(stmp,num,str,ty)
         if(ty.eq.0) then
            floats(cf)=num
            cf=cf+1
         elseif(ty.eq.1) then
            strings(cs)=str
            cs=cs+1
         else
            call raise('S',"(readline) problem in checktype, must abort")
            exit
         endif
         stmp=''
      endif
      if(i.eq.len_trim(line)) then ! special case: end of line
         call checktype(stmp,num,str,ty)
         if(ty.eq.0) then
            floats(cf)=num
            cf=cf+1
         elseif(ty.eq.1) then
            strings(cs)=str
            cs=cs+1
         else
            call raise('S',"(readline) problem in checktype, must abort")
            exit
         endif
         stmp=''
      endif
   enddo
   nstrings=cs-1
   nfloats=cf-1
end subroutine readline

!! ------------------------------------------------------------------------
!  Purpose:
!  this checks the type of the string and returns it cast to real or as string.
!
!  Input:
!     field - string to be checked
!
!  Output:
!     str   - trimmed string, if string was found
!     num   - floating point number if float was found
!     ty    - 0 if float, 1 if string, 99 if failed
!! ------------------------------------------------------------------------
subroutine checktype(field,num,str,ty)
   use iso_fortran_env, only : wp => real64
   implicit none
   character(len=*),intent(in)  :: field
   character(len=*),intent(out) :: str
   real(wp),intent(out)         :: num
   integer :: i,e,ty
   logical :: is_num

   ty=99
   str=''
   is_num=.false.
   ! cast string on real and get error code; 0 means success.
   read(field,'(F10.5)',IOSTAT=e) num
   if(e.eq.0)is_num=.true.
   if(is_num)then
      if(index(field,'.').ne.0) then ! check for integer/real
         read(field,'(F30.16)')num
         ty=0
      else
         ! if integer, add .0 to string; otherwise cast to real does not work
         str=trim(field)//'.0'
         read(str,'(F30.16)')num
         str=''
         ty=0
      endif
   else
      str=field
      ty=1
   endif
end subroutine checktype

!! ------------------------------------------------------------------------
!  Purpose:
!  translates element symbol into ordinal number
!
!  Input:
!     string - arbitary length string with element symbol
!
!  Outout:
!     nat    - ordinal number
!! ------------------------------------------------------------------------
pure elemental function elem(string) result(nat)
   implicit none
   character(len=*),intent(in)  :: string
   integer :: nat
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
   integer :: i,j,k,l,n

   nat=0
   e='  '
   do i=1,len(string)
      if(string(i:i).ne.' ') L=i
   enddo
   k=1
   do j=1,l
      if (k.gt.2)exit
      n=ichar(string(J:J))
      ! break if space after elem. symbol
      if(len_trim(e).ge.1 .and. n.eq.ichar(' '))exit
      ! break if tab after elem. symbol
      if(len_trim(e).ge.1 .and. n.eq.9)exit
      if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
         e(k:k)=char(n+ichar('a')-ichar('A'))
         k=k+1
      endif
      if(n.ge.ichar('a') .and. n.le.ichar('z') )then
         e(k:k)=string(j:j)
         k=k+1
      endif
   enddo

   do i=1,118
      if(e.eq.elemnt(i))then
         nat=i
         return
      endif
   enddo
end function elem

end subroutine get_geometry

!> implements readers for different periodic geometry types
!
!  currently Turbomole coord format and VASP POSCAR are supported
module geometry_reader
   use iso_fortran_env, wp => real64
   use mctc_systools
   implicit none

   public :: read_geometry
   public :: read_coord
   public :: read_poscar
   public :: read_xmol
   private

   integer,private,parameter :: p_str_length = 48
   integer,private,parameter :: p_arg_length = 24

contains

!> interface to read in a geometry from a file
subroutine read_geometry(fname,mol)
   use iso_fortran_env, wp => real64
   use class_molecule
   implicit none
   character(len=*),intent(in)    :: fname !< geometry file name
   type(molecule),  intent(inout) :: mol   !< molecular structure information
   integer :: ifile
   integer :: filetype
   integer :: err
   logical :: exist

   character(len=:),allocatable :: line
   logical :: has_tmdg, named_poscar, named_coord
   integer :: is

   inquire(file=fname,exist=exist)
   if (.not.exist) call raise('E',"Could not find '"//fname//"'")
   open(newunit=ifile,file=fname,iostat=err,action='read',status='old')
   if (err.ne.0)   call raise('E',"Could not open '"//fname//"'")

   call getline(ifile,line,err)
   if (err.ne.0)   call raise('E',"'"//fname//"' seems empty to me.")

   is = index(fname,'/') + 1
   has_tmdg = index(line,'$') > 0
   named_poscar = index(fname(is:),'POSCAR') > 0
   named_coord = index(fname(is:),'coord') > 0

   if (has_tmdg .or. named_coord) then
      call read_coord(ifile,mol)
   else ! VASP's POSCAR input was default
      call read_poscar(ifile,mol)
   endif

   close(ifile)
end subroutine read_geometry

subroutine read_xmol(iunit,mol)
   use iso_fortran_env, wp => real64
   use mctc_econv
   use class_molecule
   use pbc_tools
   implicit none
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
   if (err.ne.0) call raise('E',"Could not read number of atoms, check format!")

   if (n.lt.1) &
      call raise('E','Found no atoms, cannot work without atoms!')

   call mol%allocate(n)
   mol%npbc = 0 ! Xmol is always molecular (there are extensions to this...)
   mol%pbc = .false.

   ! drop next record
   read(iunit,'(a)')

   n = 0
   do
      call getline(iunit,line,err)
      if (is_iostat_end(err)) exit
      if (err.ne.0) call raise('E',"Could not read geometry from Xmol file")
      if (debug) print'(">",a)',line
      read(line,*,iostat=err) chdum, xyz(1), xyz(2), xyz(3)
      if (debug) print'("->",a)',chdum
      if (debug) print'("->",3g0)',xyz

      iat = elem(line)
      if (debug) print'("->",g0)',iat
      if (iat > 0) then
         n = n+1
         if (n > mol%nat) call raise('E',"Atom number missmatch in Xmol file")
         mol%at(n) = iat
         mol%sym(n) = line(1:2)
         mol%xyz(:,n) = xyz*conv
      endif
   enddo

   if (n.ne.mol%nat) call raise('E',"Atom number missmatch in Xmol file")

end subroutine read_xmol

!> read geometry as POSCAR from iunit
subroutine read_poscar(iunit,mol)
   use iso_fortran_env, wp => real64
   use class_molecule
   use pbc_tools
   implicit none
   logical, parameter :: debug = .false.
   integer,intent(in) :: iunit !< file handle
   type(molecule),intent(inout) :: mol
   integer :: n

   call get_atomnumber(iunit,n)

   if (n.lt.1) &
      call raise('E','Found no atoms, cannot work without atoms!')

   call mol%allocate(n)
   mol%npbc = 3 ! VASP is always 3D
   mol%pbc = .true.

   call get_coord(iunit,mol%lattice,mol%nat,mol%xyz,mol%at,mol%sym)
   call dlat_to_cell(mol%lattice,mol%cellpar)
   call dlat_to_rlat(mol%lattice,mol%rec_lat)
   mol%volume = dlat_to_dvol(mol%lattice)

   call xyz_to_abc(mol%nat,mol%lattice,mol%xyz,mol%abc,mol%pbc)

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
   if (err.ne.0) call raise('E',"Could not read POSCAR")
   if (debug) print'(">",a)',line
   call parse(line,' ',args,ntype)

   ! this line contains the global scaling factor,
   call getline(iunit,line,err)
   if (err.ne.0) call raise('E',"Could not read POSCAR")
   if (debug) print'(">",a)',line
   read(line,*,iostat=err) ddum
   if (err.ne.0) call raise('E',"Could not read POSCAR")
   ! the Ang->au conversion is included in the scaling factor
   scalar = ddum*aatoau

   ! reading the lattice constants
   do i=1,3
      call getline(iunit,line,err)
      if (err.ne.0) call raise('E',"Could not read lattice from POSCAR")
      if (debug) print'("->",a)',line
      read(line,*,iostat=err) latvec
      if (err.ne.0) call raise('E',"Could not read lattice from POSCAR")
      lattice(:,i) = latvec * scalar
   enddo
   ! Either here are the numbers of each element,
   ! or (>vasp.5.1) here are the element symbols
   call getline(iunit,line,err)
   if (err.ne.0) call raise('E',"Could not read POSCAR")
   if (debug) print'(">",a)',line
   ! try to read first element
   read(line,*,iostat=err) idum
   ! CONTCAR files have additional Element line here since vasp.5.1
   if (err.ne.0) then
      call parse(line,' ',args,ntype)
      call getline(iunit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) call raise('E',"Could not read POSCAR")
   endif
   call parse(line,' ',args2,nn)
   if (nn.ne.ntype) call raise('E', 'Error reading number of atomtypes')
   ncheck=0
   do i=1,nn
      read(args2(i),*,iostat=err) inum
      iat = elem(args(i))
      if (iat < 1 .or. inum < 1) call raise('E', 'Error: unknown element.')
      do j=1,inum
         ncheck=ncheck+1
         sym(ncheck) = args(i)(1:2)
         at(ncheck)=iat
      enddo
   enddo
   if (n.ne.ncheck) call raise('E','Error reading Number of Atoms')

   call getline(iunit,line,err)
   if (err.ne.0) call raise('E',"Could not read POSCAR")
   if (debug) print'(">",a)',line
   line=adjustl(line)
   if (line(:1).eq.'s' .or. line(:1).eq.'S') then
      selective=.true.
      call getline(iunit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) call raise('E',"Could not read POSCAR")
      line=adjustl(line)
   endif

   cartesian=(line(:1).eq.'c' .or. line(:1).eq.'C' .or. &
      &       line(:1).eq.'k' .or. line(:1).eq.'K')
   do i=1,n
      call getline(iunit,line,err)
      if (err.ne.0) call raise('E',"Could not read geometry from POSCAR")
      if (debug) print'("-->",a)',line
      read(line,*,iostat=err) coord
      if (err.ne.0) call raise('E',"Could not read geometry from POSCAR")

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
   if (err.ne.0) call raise('E',"Could not read POSCAR")
   if (debug) print'(">",a)',line
   call parse(line,' ',args,ntype)

   ! this line contains the global scaling factor,
   call getline(iunit,line,err)
   if (err.ne.0) call raise('E',"Could not read POSCAR")
   if (debug) print'(">",a)',line
   read(line,*,iostat=err) ddum
   if (err.ne.0) call raise('E',"Could not read POSCAR")
   ! the Ang->au conversion is included in the scaling factor
   scalar = ddum*aatoau

   ! reading the lattice constants
   do i=1,3
      call getline(iunit,line,err)
      if (err.ne.0) call raise('E',"Could not read lattice from POSCAR")
      if (debug) print'("->",a)',line
   enddo
   ! Either here are the numbers of each element,
   ! or (>vasp.5.1) here are the element symbols
   call getline(iunit,line,err)
   if (err.ne.0) call raise('E',"Could not read POSCAR")
   if (debug) print'(">",a)',line
   ! try to read first element
   read(line,*,iostat=err) idum
   ! CONTCAR files have additional Element line here since vasp.5.1
   if (err.ne.0) then
      call parse(line,' ',args,ntype)
      call getline(iunit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) call raise('E',"Could not read POSCAR")
   endif
   call parse(line,' ',args2,nn)
   if (nn.ne.ntype) call raise('E', 'Error reading number of atomtypes')
   n=0
   do i=1,nn
      read(args2(i),*,iostat=err) inum
      iat = elem(args(i))
      if (iat < 1 .or. inum < 1) call raise('E', 'Error: unknown element.')
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
subroutine read_coord(iunit,mol)
   use iso_fortran_env, wp => real64
   use class_molecule
   use mctc_readin, getline => strip_line
   use pbc_tools
   implicit none
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

!  we need to read this file ~twice~ trice, for a lot of reasons
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

   if (natoms.eq.-1) &
      call raise('E',"coord keyword was not found")
   if (npbc.eq.-1) npbc = 0

   if (natoms.lt.1) & ! this will catch an empty coord data group
      call raise('E','Found no atoms, cannot work without atoms!')

   call mol%allocate(natoms,npbc > 0)
   mol%npbc = npbc
   do i = 1, npbc
      mol%pbc(i) = .true.
   enddo

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
               call raise('S',"Multiple definitions of lattice present!")
            exist = .true.
            if (debug) print'("*>",a)',line
            call get_lattice(iunit,line,err,npbc,mol%lattice)
            call dlat_to_cell(mol%lattice,mol%cellpar)
            call dlat_to_rlat(mol%lattice,mol%rec_lat)
            mol%volume = dlat_to_dvol(mol%lattice)
         elseif (index(line,'cell').eq.2) then
            if (exist) &
               call raise('S',"Multiple definitions of lattice present!")
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

   if (.not.exist .and. npbc > 0) &
      call raise('E',"There is no definition of the lattice present!")

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
         call raise('E',"This kind of periodicity is not implemented")
      endif
   else
      call raise('E',"Could not read periodicity from '"//line//"'")
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
      call raise('S',"lattice data group is present but not periodic?")
      return
   endif
   if (len(line).le.8) then
      conv = 1.0_wp ! default is bohr
   elseif (index(line(9:),'bohr').gt.0) then
      conv = 1.0_wp
   elseif (index(line(9:),'angs').gt.0) then
      conv = aatoau
   else ! fall back to default
      call raise('S',"Could not make sense out of unit '"//line(7:)//"'")
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
      if (ncount.gt.npbc) &
         call raise('E',"Input error in lattice data group")
      read(line,*,iostat=err) lvec(1:npbc)
      if (err.ne.0) &
         call raise('E',"Input error in lattice data group")
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
      call raise('S',"cell data group is present but not periodic?")
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
      call raise('S',"Could not make sense out of unit '"//line(7:)//"'")
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
   if (err.ne.0) &
      call raise('E',"Could not read cell from '"//line//"'")

   if(alen < 0.0d0 .or. blen < 0.0d0 .or. clen < 0.0d0 .or. &
      alp < 0.0d0 .or. bet < 0.0d0 .or. gam < 0.0d0) &
      call raise('E',"Negative cell parameters make no sense!")

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
   integer,intent(in) :: iunit          !< file handle
   character(len=:),allocatable :: line !< string buffer
   integer,intent(out) :: err           !< status of iunit
   type(molecule),intent(inout) :: mol  !< molecular structure information
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
      call raise('S',"Could not make sense out of unit '"//line(7:)//"'")
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
      if (ncount.gt.mol%nat) &
         call raise('E',"Internal error while reading coord data group")
      call parse(line,' ',argv,narg)
      if (narg.lt.4) &
         call raise('E',"not enough arguments in line in coord data group")
      iat = elem(argv(4))
      if (iat.eq.0) &
         call raise('E',"invalid element input in line in coord data group")
      mol%sym(ncount) = trim(argv(4))
      mol%at(ncount) = iat
      do icoord = 1, 3
         if (get_value(trim(argv(icoord)),ddum)) then
            xyz(icoord) = conv*ddum
         else
            call raise('E',"invalid coordinate input in line in coord data group")
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

!subroutine read_xmol(iunit,mol)
!   use iso_fortran_env, wp => real64
!   use class_molecule
!   use mctc_readin, getline => strip_line
!   use mctc_strings
!   use mctc_econv
!   implicit none
!   integer,intent(in) :: iunit !< file handle
!   type(molecule),intent(inout) :: mol
!   character(len=:),allocatable :: line
!   integer :: err
!   integer :: idum
!   logical,parameter :: debug = .false.
!   integer :: natoms
!   integer :: iat,icoord
!   integer :: ncount
!   logical :: exist
!   integer  :: narg
!   real(wp) :: ddum,xyz(3)
!   character(len=p_str_length),dimension(p_arg_length) :: argv
!
!   rewind(iunit)
!
!   call getline(iunit,line,err)
!   if (debug) print'(" ->",a)',line
!   read(line,*,iostat=err) natoms
!   if (err.ne.0) then
!      call raise('E',"Could not read number of atoms")
!   endif
!   read(iunit,'(a)')
!
!   call mol%allocate(natoms,.false.)
!
!   ncount = 0
!   do
!      call getline(iunit,line,err)
!      if (err.eq.iostat_end) exit
!      if (debug) print'(" ->",a)',line
!
!      if (line.eq.'') cycle ! skip empty lines
!      ncount = ncount + 1   ! but count non-empty lines first
!      if (ncount.gt.mol%nat) &
!         call raise('E',"Internal error while reading coord data group")
!      call parse(line,' ',argv,narg)
!      if (narg.lt.4) &
!         call raise('E',"not enough arguments in line in coord data group")
!      iat = elem(argv(1))
!      if (iat.eq.0) &
!         call raise('E',"invalid element input in line in coord data group")
!      mol%sym(ncount) = trim(argv(4))
!      mol%at(ncount) = iat
!      do icoord = 1, 3
!         if (get_value(trim(argv(icoord+1)),ddum)) then
!            xyz(icoord) = aatoau*ddum
!         else
!            call raise('E',"invalid coordinate input in line in coord data group")
!         endif
!      enddo
!
!   enddo
!
!end subroutine read_xmol

! ------------------------------------------------------------------------
!  Purpose:
!  translates element symbol into ordinal number
!
!  Input:
!     string - arbitary length string with element symbol
!
!  Outout:
!     nat    - ordinal number
! ------------------------------------------------------------------------
pure elemental function elem(string) result(nat)
   implicit none
   character(len=*),intent(in)  :: string
   integer :: nat
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
   integer :: i,j,k,l,n

   nat=0
   e='  '
   do i=1,len(string)
      if(string(i:i).ne.' ') L=i
   enddo
   k=1
   do j=1,l
      if (k.gt.2)exit
      n=ichar(string(J:J))
      ! break if space after elem. symbol
      if(len_trim(e).ge.1 .and. n.eq.ichar(' '))exit
      ! break if tab after elem. symbol
      if(len_trim(e).ge.1 .and. n.eq.9)exit
      if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
         e(k:k)=char(n+ichar('a')-ichar('A'))
         k=k+1
      endif
      if(n.ge.ichar('a') .and. n.le.ichar('z') )then
         e(k:k)=string(j:j)
         k=k+1
      endif
   enddo

   do i=1,118
      if(e.eq.elemnt(i))then
         nat=i
         return
      endif
   enddo
end function elem

end module geometry_reader
