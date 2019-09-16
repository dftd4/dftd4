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

module mctc_environment
   use iso_fortran_env
   implicit none
   public :: mctc_logger
   public :: mctc_argparser
   private

   type string
      character(len=:), allocatable :: raw
   end type string

   type mctc_logger
      integer :: iunit = output_unit
      integer :: length = 0
      logical :: sane = .true.
      type(string), allocatable :: log(:)
   contains
      procedure :: reallocate => reallocate_log
      procedure :: error => raise_error
      procedure :: warning => raise_warning
      procedure :: checkpoint => check_state
      procedure :: write => write_log
      procedure :: reset => reset_log
   end type mctc_logger

   type mctc_argparser
      integer :: n = 0
      type(string), allocatable :: arg(:)
   contains
      procedure :: new => new_argparser
      generic :: get => get_argument_as_character, &
         &              get_argument_as_integer, &
         &              get_argument_as_real, &
         &              get_argument_as_logical
      procedure, private :: get_argument_as_character
      procedure, private :: get_argument_as_integer
      procedure, private :: get_argument_as_real
      procedure, private :: get_argument_as_logical
   end type mctc_argparser

contains

subroutine new_argparser(self)
   class(mctc_argparser), intent(inout) :: self
   character(len=:), allocatable :: arg
   integer :: iarg,len
   self%n = command_argument_count()
   allocate( self%arg(0:self%n) )
   do iarg = 0, self%n
      call get_command_argument(iarg,length=len)
      allocate(character(len=len) :: arg)
      call get_command_argument(iarg,arg)
      call move_alloc(arg,self%arg(iarg)%raw)
   enddo
end subroutine new_argparser

logical function get_argument_as_character(self,iarg,arg) result(status)
   class(mctc_argparser), intent(inout) :: self
   character(len=:), allocatable, intent(out) :: arg
   integer, intent(in) :: iarg
   status = iarg >= 0 .and. iarg <= self%n
   if (status) arg = self%arg(iarg)%raw
end function get_argument_as_character

logical function get_argument_as_integer(self,iarg,ival) result(status)
   class(mctc_argparser), intent(inout) :: self
   character(len=:), allocatable :: arg
   integer, intent(out) :: ival
   integer, intent(in) :: iarg
   integer :: idum
   integer :: err
   status = self%get(iarg,arg)
   if (status) then
      read(arg,*,iostat=err) idum
      status = err /= 0
   endif
   if (status) ival = idum
end function get_argument_as_integer

logical function get_argument_as_logical(self,iarg,lval) result(status)
   class(mctc_argparser), intent(inout) :: self
   character(len=:), allocatable :: arg
   logical, intent(out) :: lval
   integer, intent(in) :: iarg
   logical :: ldum
   integer :: err
   status = self%get(iarg,arg)
   if (status) then
      read(arg,*,iostat=err) ldum
      status = err /= 0
   endif
   if (status) lval = ldum
end function get_argument_as_logical

logical function get_argument_as_real(self,iarg,dval) result(status)
   use iso_fortran_env, only: wp => real64
   class(mctc_argparser), intent(inout) :: self
   character(len=:), allocatable :: arg
   real(wp), intent(out) :: dval
   integer, intent(in) :: iarg
   real(wp) :: ddum
   integer  :: err
   status = self%get(iarg,arg)
   if (status) then
      read(arg,*,iostat=err) ddum
      status = err /= 0
   endif
   if (status) dval = ddum
end function get_argument_as_real

subroutine reallocate_log(self)
   class(mctc_logger), intent(inout) :: self
   type(string), allocatable :: tmp(:)
   integer :: current, new_size
   if (allocated(self%log)) then
      if (self%length > size(self%log,1)) then
         current = size(self%log,1)
         new_size = current + current/2 + 1
         allocate(tmp(new_size))
         tmp(:current) = self%log
         deallocate(self%log)
         call move_alloc(tmp,self%log)
      endif
   else
      allocate(self%log(10))
   endif
end subroutine reallocate_log

subroutine deallocate_log(self)
   type(mctc_logger), intent(inout) :: self
   if (allocated(self%log)) deallocate(self%log)
   self%sane = .true.
end subroutine

subroutine raise_error(self,code,message)
   class(mctc_logger), intent(inout) :: self
   integer, intent(in) :: code
   character(len=*), intent(in) :: message
   integer :: pos
   self%sane = .false.
   call self%reallocate
   pos = self%length+1
   self%log(pos) = string(message)
   self%length = pos
end subroutine raise_error

subroutine raise_warning(self,code,message)
   class(mctc_logger), intent(inout) :: self
   integer, intent(in) :: code
   character(len=*), intent(in) :: message
   integer :: pos
   call self%reallocate
   pos = self%length+1
   self%log(pos) = string(message)
   self%length = pos
end subroutine raise_warning

subroutine check_state(self)
   class(mctc_logger), intent(inout) :: self
   if (.not.self%sane) then
      call self%write("Fatal error encountered must stop! Error backtrace:")
      call terminate(1)
   endif
end subroutine check_state

subroutine write_log(self,message)
   class(mctc_logger), intent(inout) :: self
   character(len=*), intent(in) :: message
   integer :: i
   if (self%length > 0) then
      write(self%iunit,'(72("#"))')
      if (self%sane) then
         write(self%iunit,'("# WARNING!",1x,a)') message
      else
         write(self%iunit,'("# ERROR!",1x,a)') message
      endif
      do i = 1, self%length
         write(self%iunit,'("#  -",1x,a)')  &
            &  self%log(i)%raw
      enddo
      write(self%iunit,'(72("#"))')
   endif
end subroutine write_log

subroutine reset_log(self)
   class(mctc_logger), intent(inout) :: self
   self%length = 0
end subroutine reset_log

end module mctc_environment
