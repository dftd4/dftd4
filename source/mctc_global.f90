module mctc_global
   character(len=:),allocatable :: name
   character(len=:),allocatable,target :: msgbuffer
   integer :: msgid
   integer :: maxmsg = 100

   type :: errormsg
      character(len=:),allocatable :: msg
      integer :: len
   endtype errormsg

   type(errormsg),allocatable :: errorbuffer(:)

contains

subroutine init_errorbuffer
   implicit none
   if (allocated(errorbuffer)) deallocate ( errorbuffer )
   allocate ( errorbuffer(maxmsg) )
   msgid = 0
   msgbuffer = ' '
end subroutine init_errorbuffer

end module mctc_global

