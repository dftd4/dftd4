!> @brief wrapper for IO functions to work with allocatable characters
module mctc_systools

   character,parameter :: space = ' '
   character,parameter :: colon = ':'
   character,parameter :: slash = '/'

contains

!> @brief reads a line from unit into an allocatable character
subroutine getline(unit,line,iostat)
   use iso_fortran_env, only : iostat_eor
   integer,intent(in) :: unit
   character(len=:),allocatable,intent(out) :: line
   integer,intent(out),optional :: iostat

   integer,parameter  :: buffersize=144
   character(len=buffersize) :: buffer
   integer :: size
   integer :: err

   line = ''
   do
      read(unit,'(a)',advance='no',iostat=err,size=size)  &
      &    buffer
      if (err.gt.0) then
         if (present(iostat)) iostat=err
         return ! an error occurred
      endif
      line = line // buffer(:size)
      if (err.lt.0) then
         if (err.eq.iostat_eor) err = 0
         if (present(iostat)) iostat=err
         return
      endif
   enddo

end subroutine getline

!> @brief searches for an file in a given path variable
subroutine rdpath(path,arg,fname,ex)
   implicit none
   character(len=*),intent(in)  :: arg  !< file to be found in path
   character(len=*),intent(in)  :: path !< path variable
   character(len=:),allocatable,intent(out) :: fname !< contains path on success
   logical,intent(out),optional :: ex   !< exit status

!* temporary variables
   character(len=:),allocatable :: scratch1
   character(len=:),allocatable :: scratch2
   character(len=:),allocatable :: fpath
   logical :: exist
   integer :: i

   scratch1 = path
   do
      i = index(scratch1,colon)
      if (i.eq.0) then
         scratch2 = scratch1
      else
         scratch2 = scratch1(:i-1)
         scratch1 = scratch1(i+1:)
      endif
      fpath = scratch2//slash//arg
      inquire(file=fpath,exist=exist)
      if (exist) exit
!     print*,fpath,space,scratch1,exist
      if (i.eq.0) exit
   enddo

!  print*,fpath,exist

   if (exist) fname = fpath
   if (present(ex)) ex = exist
   
end subroutine rdpath

!> @brief reads a command line argument in an allocatable character
subroutine rdarg(i,arg,iostat)
   integer,intent(in) :: i !< number of argument
   character(len=:),allocatable,intent(out) :: arg !< contains argument on exit
   integer,intent(out),optional :: iostat !< exit status
   integer :: l,err
   if (allocated(arg)) deallocate(arg)
   call get_command_argument(i,length=l,status=err)
   if (err.ne.0) then
      if (present(iostat)) then
         iostat = err
         return
      else
         call raise('E','Command argument corrupted')
      endif
   endif
   allocate( character(len=l) :: arg, stat=err )
   if (err.ne.0) then
      if (present(iostat)) then
         iostat = err
         return
      else
         call raise('E','could not be allocated')
      endif
   endif
   call get_command_argument(i,arg,status=err)
   if (err.ne.0) then
      if (present(iostat)) then
         iostat = err
         return
      else
         call raise('E','Command argument corrupted')
      endif
   endif
   if (present(iostat)) iostat=0
end subroutine rdarg

!> @brief reads a system cariable in an allocatable character
subroutine rdvar(name,var,iostat)
   character(len=*),intent(in) :: name
   character(len=:),allocatable,intent(out) :: var
   integer,intent(out),optional :: iostat
   integer :: l,err
   if (allocated(var)) deallocate(var)
   call get_environment_variable(name,length=l,status=err)
   if (err.ne.0) then
      if (present(iostat)) then
         iostat = err
         return
      else
         call raise('E','System variable unassigned')
      endif
   endif
   allocate( character(len=l) :: var, stat=err )
   if (err.ne.0) then
      if (present(iostat)) then
         iostat = err
         return
      else
         call raise('E','could not be allocated')
      endif
   endif
   call get_environment_variable(name,var,status=err)
   if (err.ne.0) then
      if (present(iostat)) then
         iostat = err
         return
      else
         call raise('E','System variable corrupted')
      endif
   endif
   if (present(iostat)) iostat=0
end subroutine rdvar

end module mctc_systools
