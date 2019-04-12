!> error handler for MCTC library
subroutine raise(mode,message)
   use iso_fortran_env
   use mctc_global
   character,       intent(in) :: mode    !< kind of operation
   character(len=*),intent(in) :: message !< string containing error description
   select case(mode)
   case('S','s') ! save to message buffer
      if (allocated(errorbuffer)) then
         call save_warning
      else
         call warning
      endif
   case('F','f') ! flush message buffer
      if (msgid.gt.0) then
         write(output_unit,'(72(''#''))')
         write(output_unit,'(''# WARNING!'',1x,a)') message
         do i = 1, msgid
            write(output_unit,'(''#  -'',1x,a)')  &
            &  errorbuffer(i) % msg
         enddo
         write(output_unit,'(72(''#''))')
         ! after we are done with printing, clear the buffer
         call init_errorbuffer
      endif
   case('W','w') ! print warning directly
      call warning
   case('E','e')
      call error
   end select
contains
subroutine save_warning
   msgid = msgid + 1
   errorbuffer(msgid) % msg = message
   errorbuffer(msgid) % len = len(message)
end subroutine save_warning
subroutine warning
   write(output_unit,'(''#WARNING!'',1x,a)') message
end subroutine warning
subroutine error
   write(output_unit,'(''#ERROR!'',1x,a)')   message
   call terminate(1)
end subroutine error
end subroutine raise

subroutine terminate(signal)
   use iso_fortran_env, only : error_unit
   use mctc_global, only : name
   integer,intent(in) :: signal
   integer,parameter  :: p_exit_success = 0
   integer,parameter  :: p_exit_failure = 1
   integer,parameter  :: p_exit_external = -1
   if (.not.allocated(name)) name = 'program'
   select case(signal)
   case(p_exit_success)
      write(error_unit,'(  ''normal termination of'',1x,a)') name
      stop
   case(p_exit_external)
      write(error_unit,'(''external termination of'',1x,a)') name
      error stop
   case default
      write(error_unit,'(''abnormal termination of'',1x,a)') name
      error stop
   end select
end subroutine terminate

subroutine wsigint
   use iso_fortran_env, only : error_unit
   write(error_unit,'(''recieved SIGINT, terminating...'')')
   call terminate(-1)
end subroutine wsigint

subroutine wsigterm
   use iso_fortran_env, only : error_unit
   write(error_unit,'(''recieved SIGTERM, terminating...'')')
   call terminate(-1)
end subroutine wsigterm
