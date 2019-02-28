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
