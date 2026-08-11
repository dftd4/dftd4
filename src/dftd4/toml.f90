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

!> TOML-backed database for DFT-D4 damping parameters.
module dftd4_toml
   use dftd4_damping_rational, only : rational_damping_param
   use dftd4_damping, only : damping_param
   use mctc_env, only : error_type, fatal_error
   use tomlf, only : toml_array, toml_error, toml_key, toml_parse, toml_table, &
      & get_value, len
   implicit none
   private

   public :: param_database


   !> Individual damping parameter record.
   type :: param_record
      !> Functional name identifying this record.
      character(len=:), allocatable :: key
      !> Parameter scheme identifying this record.
      character(len=:), allocatable :: id
      !> Actual rational damping parameters.
      type(rational_damping_param) :: param
      !> Name of the damping function.
      character(len=:), allocatable :: damping
      !> Reference to the publication defining the parameters.
      character(len=:), allocatable :: doi
   end type param_record

   !> Damping parameter database.
   type :: param_database
      !> Default parameter schemes.
      type(param_record), allocatable :: defaults(:)
      !> Parameter records indexed by functional and scheme.
      type(param_record), allocatable :: records(:)
      !> Mask selecting the default schemes.
      logical, allocatable :: mask(:)
   contains
      !> Read damping parameter data.
      generic :: load => load_from_file, load_from_unit, load_from_toml
      procedure, private :: load_from_file
      procedure, private :: load_from_unit
      procedure, private :: load_from_toml
      !> Retrieve a damping parameter from the database.
      procedure :: get
   end type param_database


contains


!> Read damping parameter data from a file.
subroutine load_from_file(self, file, error)
   !> Damping parameter database.
   class(param_database), intent(inout) :: self
   !> File name.
   character(len=*), intent(in) :: file
   !> Error handling.
   type(error_type), allocatable, intent(out) :: error

   integer :: unit
   logical :: exist

   inquire(file=file, exist=exist)
   if (.not.exist) then
      call fatal_error(error, "Could not find parameter file '"//trim(file)//"'")
      return
   end if

   open(file=file, newunit=unit, status="old", action="read")
   call self%load(unit, error)
   close(unit)
end subroutine load_from_file


!> Read damping parameter data from a formatted unit.
subroutine load_from_unit(self, unit, error)
   !> Damping parameter database.
   class(param_database), intent(inout) :: self
   !> Input unit.
   integer, intent(in) :: unit
   !> Error handling.
   type(error_type), allocatable, intent(out) :: error

   type(toml_error), allocatable :: parse_error
   type(toml_table), allocatable :: table

   call toml_parse(table, unit, parse_error)
   if (allocated(parse_error)) then
      allocate(error)
      call move_alloc(parse_error%message, error%message)
      return
   end if

   call self%load(table, error)
end subroutine load_from_unit


!> Read damping parameter data from a TOML table.
subroutine load_from_toml(self, table, error)
   !> Damping parameter database.
   class(param_database), intent(inout) :: self
   !> Parsed TOML data.
   type(toml_table), intent(inout) :: table
   !> Error handling.
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   if (allocated(self%defaults)) deallocate(self%defaults)
   if (allocated(self%records)) deallocate(self%records)
   if (allocated(self%mask)) deallocate(self%mask)

   call get_value(table, "default", child)
   if (.not.associated(child)) then
      call fatal_error(error, "Missing 'default' table in parameter file")
      return
   end if
   call load_default(self, child, error)
   if (allocated(error)) return

   call get_value(table, "parameter", child)
   if (.not.associated(child)) then
      call fatal_error(error, "Missing 'parameter' table in parameter file")
      return
   end if
   call load_parameter(self, child, error)
end subroutine load_from_toml


!> Read the default parameter schemes from the TOML table.
subroutine load_default(self, table, error)
   !> Damping parameter database.
   type(param_database), intent(inout) :: self
   !> TOML default table.
   type(toml_table), intent(inout) :: table
   !> Error handling.
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child, child2
   type(toml_array), pointer :: children
   type(toml_key), allocatable :: keys(:)
   type(param_record) :: stub
   character(len=:), allocatable :: val
   integer :: ik

   ! The TOML defaults contain values common to all D4 records.  Initialize
   ! the otherwise required rational fields so that omitted optional entries
   ! still have deterministic values while the default records are built.
   stub%param = rational_damping_param(s6=1.0, s8=0.0, s9=1.0, a1=0.0, a2=0.0, alp=16.0)

   call get_value(table, "d4", children)
   call get_value(table, "parameter", child)
   if (.not.associated(children) .or. .not.associated(child)) then
      call fatal_error(error, "Missing D4 defaults in parameter file")
      return
   end if
   call get_value(child, "d4", child2)
   if (.not.associated(child2)) then
      call fatal_error(error, "Missing D4 default parameters in parameter file")
      return
   end if

   call child2%get_keys(keys)
   call resize(self%defaults, size(keys))
   do ik = 1, size(keys)
      call get_value(child2, keys(ik)%key, child)
      if (.not.associated(child)) then
         call fatal_error(error, "Invalid D4 default parameter '"//keys(ik)%key//"'")
         return
      end if
      call load_record(self%defaults(ik), child, stub, error)
      self%defaults(ik)%key = ""
      if (allocated(error)) return
   end do

   allocate(self%mask(size(keys)), source=.false.)
   do ik = 1, len(children)
      call get_value(children, ik, val)
      associate(id => get_record(self%defaults, "", val))
         if (id > 0) self%mask(id) = .true.
      end associate
   end do
end subroutine load_default


!> Read functional parameter records from the TOML table.
subroutine load_parameter(self, table, error)
   !> Damping parameter database.
   type(param_database), intent(inout) :: self
   !> TOML parameter table.
   type(toml_table), intent(inout) :: table
   !> Error handling.
   type(error_type), allocatable, intent(out) :: error

   type(param_record) :: stub
   type(toml_key), allocatable :: keys(:), list(:)
   type(toml_table), pointer :: child, child2
   integer :: nr, ik, iv, id

   stub%param = rational_damping_param(s6=1.0, s8=0.0, s9=1.0, a1=0.0, a2=0.0, alp=16.0)
   nr = 0
   call table%get_keys(keys)
   call resize(self%records, size(keys))

   records: do ik = 1, size(keys)
      call get_value(table, keys(ik)%key, child)
      if (.not.associated(child)) then
         call fatal_error(error, "Missing D4 parameters for '"//keys(ik)%key//"'")
         return
      end if
      call get_value(child, "d4", child2)
      if (.not.associated(child2)) then
         call fatal_error(error, "Missing D4 parameters for '"//keys(ik)%key//"'")
         return
      end if

      call child2%get_keys(list)
      if (nr + size(list) > size(self%records)) call resize(self%records)
      do iv = 1, size(list)
         id = get_record(self%defaults, "", list(iv)%key)
         call get_value(child2, list(iv)%key, child)
         if (id > 0) then
            call load_record(self%records(iv + nr), child, self%defaults(id), error)
         else
            call load_record(self%records(iv + nr), child, stub, error)
         end if
         self%records(iv + nr)%key = keys(ik)%key
         if (allocated(error)) exit records
      end do
      nr = nr + size(list)
   end do records

   if (allocated(error)) return
   call resize(self%records, nr)
end subroutine load_parameter


!> Deserialize one rational damping record from a TOML table.
subroutine load_record(record, table, default, error)
   !> Parameter record.
   type(param_record), intent(inout) :: record
   !> TOML record.
   type(toml_table), intent(inout) :: table
   !> Default values.
   type(param_record), intent(in) :: default
   !> Error handling.
   type(error_type), allocatable, intent(out) :: error

   call table%get_key(record%id)
   call get_value(table, "damping", record%damping, default%damping)
   call get_value(table, "doi", record%doi, default%doi)
   call get_value(table, "s6", record%param%s6, default%param%s6)
   call get_value(table, "s8", record%param%s8, default%param%s8)
   call get_value(table, "s9", record%param%s9, default%param%s9)
   call get_value(table, "a1", record%param%a1, default%param%a1)
   call get_value(table, "a2", record%param%a2, default%param%a2)
   call get_value(table, "alp", record%param%alp, default%param%alp)
end subroutine load_record


!> Retrieve a rational damping parameter from the database.
subroutine get(self, param, method, scheme)
   !> Damping parameter database.
   class(param_database), intent(inout) :: self
   !> Damping parameters.
   class(damping_param), allocatable, intent(out) :: param
   !> Functional name.
   character(len=*), intent(in) :: method
   !> Optional parameter scheme, e.g. `bj-eeq-atm`.
   character(len=*), intent(in), optional :: scheme

   integer :: ir, id

   if (.not.allocated(self%records)) return

   if (present(scheme)) then
      ir = get_record(self%records, method, scheme)
   else
      ir = 0
      do id = 1, size(self%defaults)
         if (self%mask(id)) then
            ir = get_record(self%records, method, self%defaults(id)%id)
            if (ir > 0) exit
         end if
      end do
   end if
   if (ir == 0) return

   associate(record => self%records(ir))
      select case(record%damping)
      case("bj", "rational")
         block
            type(rational_damping_param), allocatable :: tmp
            allocate(tmp)
            tmp = record%param
            call move_alloc(tmp, param)
         end block
      end select
   end associate
end subroutine get


!> Find a record by functional and scheme identifiers.
pure function get_record(records, key, id) result(pos)
   !> Records to search.
   type(param_record), intent(in) :: records(:)
   !> Functional name.
   character(len=*), intent(in) :: key
   !> Scheme identifier.
   character(len=*), intent(in) :: id
   !> Record position, or zero if not found.
   integer :: pos

   integer :: ii

   pos = 0
   do ii = 1, size(records)
      if (allocated(records(ii)%key) .and. allocated(records(ii)%id)) then
         if (records(ii)%key == key .and. records(ii)%id == id) then
            pos = ii
            exit
         end if
      end if
   end do
end function get_record


!> Reallocate a list of parameter records.
pure subroutine resize(var, n)
   !> Array to resize.
   type(param_record), allocatable, intent(inout) :: var(:)
   !> Final size; grow by a factor when omitted.
   integer, intent(in), optional :: n

   type(param_record), allocatable :: tmp(:)
   integer :: this_size, new_size
   integer, parameter :: initial_size = 16

   if (allocated(var)) then
      this_size = size(var, 1)
      call move_alloc(var, tmp)
   else
      this_size = initial_size
   end if

   if (present(n)) then
      new_size = n
   else
      new_size = this_size + this_size / 2 + 1
   end if

   allocate(var(new_size))
   if (allocated(tmp)) then
      this_size = min(size(tmp, 1), size(var, 1))
      var(:this_size) = tmp(:this_size)
      deallocate(tmp)
   end if
end subroutine resize


end module dftd4_toml
