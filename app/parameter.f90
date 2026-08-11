! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Loading of the command-line parameter database.
module dftd4_app_parameter
   use dftd4, only : load_parameters
   use mctc_env, only : error_type, fatal_error
   use mctc_env_system, only : is_windows
   implicit none
   private

   public :: load_parameter_database

contains


!> Load the parameter database for the command-line application.
subroutine load_parameter_database(error)
   !> Error handling.
   type(error_type), allocatable, intent(out) :: error

   character(len=4096) :: value, paths(7)
   character :: separator
   integer :: length, stat, pos, npath, ipath

   npath = 0
   call get_environment_variable("DFTD4_PARAMETER_FILE", value, length=length, status=stat)
   if (stat == 0 .and. length > 0) call add_path(paths, npath, value(:length))
   call add_path(paths, npath, "assets/parameters.toml")
   call add_path(paths, npath, "parameters.toml")

   call get_command_argument(0, value, length=length, status=stat)
   if (stat == 0 .and. length > 0) then
      pos = scan(value(:length), "/\", back=.true.)
      if (pos > 0) then
         separator = "/"
         if (is_windows()) separator = "\"
         call add_path(paths, npath, value(:pos)//".."//separator//"share"//separator//"dftd4"//separator//"parameters.toml")
         call add_path(paths, npath, value(:pos)//".."//separator//"assets"//separator//"parameters.toml")
         call add_path(paths, npath, value(:pos)//".."//separator//".."//separator//"assets"//separator//"parameters.toml")
         call add_path(paths, npath, value(:pos)//".."//separator//".."//separator//".."//separator// &
            & "assets"//separator//"parameters.toml")
      end if
   end if

   do ipath = 1, npath
      call load_parameters(trim(paths(ipath)), error)
      if (.not.allocated(error)) return
   end do

   if (allocated(error)) deallocate(error)
   call fatal_error(error, "Could not locate the DFT-D4 parameter database")
end subroutine load_parameter_database


!> Append a non-empty candidate path to the search list.
subroutine add_path(paths, npath, path)
   character(len=*), intent(inout) :: paths(:)
   integer, intent(inout) :: npath
   character(len=*), intent(in) :: path

   if (len_trim(path) == 0 .or. npath == size(paths)) return
   npath = npath + 1
   paths(npath) = trim(path)
end subroutine add_path


end module dftd4_app_parameter
