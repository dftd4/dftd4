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

module class_results
   use iso_fortran_env, wp => real64
   implicit none

   public :: dftd_results
   private

   type :: dftd_results
      real(wp), allocatable :: energy
      real(wp), allocatable :: gradient(:,:)
      real(wp), allocatable :: stress(:,:)
      real(wp), allocatable :: lattice_gradient(:,:)
      real(wp), allocatable :: charges(:)
      real(wp), allocatable :: dipole_moment(:)
      real(wp), allocatable :: polarizibilities(:)
      real(wp), allocatable :: c6_coefficients(:,:)
      real(wp), allocatable :: hessian(:,:)
   contains
      procedure :: allocate => allocate_results
      procedure :: deallocate => deallocate_results
      procedure :: write_json
      procedure :: write_toml
   end type dftd_results

contains

subroutine allocate_results(self,n)
   class(dftd_results), intent(inout) :: self
   integer, intent(in) :: n
   call self%deallocate
   self%energy = 0.0_wp
   allocate(self%gradient(3,n),         source = 0.0_wp)
   allocate(self%stress(3,3),           source = 0.0_wp)
   allocate(self%lattice_gradient(3,3), source = 0.0_wp)
   allocate(self%charges(n),            source = 0.0_wp)
   allocate(self%dipole_moment(3),      source = 0.0_wp)
   allocate(self%polarizibilities(n),   source = 0.0_wp)
   allocate(self%c6_coefficients(n,n),  source = 0.0_wp)
   allocate(self%hessian(3*n,3*n),      source = 0.0_wp)
end subroutine allocate_results

subroutine deallocate_results(self)
   class(dftd_results), intent(inout) :: self
   if (allocated(self%energy))           deallocate(self%energy)
   if (allocated(self%stress))           deallocate(self%stress)
   if (allocated(self%lattice_gradient)) deallocate(self%lattice_gradient)
   if (allocated(self%gradient))         deallocate(self%gradient)
   if (allocated(self%charges))          deallocate(self%charges)
   if (allocated(self%dipole_moment))    deallocate(self%dipole_moment)
   if (allocated(self%polarizibilities)) deallocate(self%polarizibilities)
   if (allocated(self%c6_coefficients))  deallocate(self%c6_coefficients)
   if (allocated(self%hessian))          deallocate(self%hessian)
end subroutine deallocate_results

subroutine write_json(self, unit, indentation)
   include 'dftd4_version.fh'
   class(dftd_results), intent(in) :: self
   integer, intent(in) :: unit
   character(len=*), intent(in), optional :: indentation
   character(len=*), parameter :: jsonkey = "('""',a,'"":',1x)"
   real(wp), allocatable :: array(:)

   write(unit, '("{")', advance='no')
   if (present(indentation)) write(unit, '(/,a)', advance='no') repeat(indentation, 1)
   write(unit, jsonkey, advance='no') 'version'
   write(unit, '(1x,a)', advance='no') '"'//version//'"'
   if (allocated(self%energy)) then
      write(unit, '(",")', advance='no')
      if (present(indentation)) write(unit, '(/,a)', advance='no') repeat(indentation, 1)
      write(unit, jsonkey, advance='no') 'energy'
      write(unit, '(1x,es25.16)', advance='no') self%energy
   end if
   if (allocated(self%stress)) then
      write(unit, '(",")', advance='no')
      if (present(indentation)) write(unit, '(/,a)', advance='no') repeat(indentation, 1)
      write(unit, jsonkey, advance='no') 'stress tensor'
      array = reshape(self%stress, [product(shape(self%stress))])
      call write_json_array(unit, array, indentation)
   end if
   if (allocated(self%lattice_gradient)) then
      write(unit, '(",")', advance='no')
      if (present(indentation)) write(unit, '(/,a)', advance='no') repeat(indentation, 1)
      write(unit, jsonkey, advance='no') 'lattice gradient'
      array = reshape(self%lattice_gradient, [product(shape(self%lattice_gradient))])
      call write_json_array(unit, array, indentation)
   end if
   if (allocated(self%gradient)) then
      write(unit, '(",")', advance='no')
      if (present(indentation)) write(unit, '(/,a)', advance='no') repeat(indentation, 1)
      write(unit, jsonkey, advance='no') 'gradient'
      array = reshape(self%gradient, [product(shape(self%gradient))])
      call write_json_array(unit, array, indentation)
   end if
   if (allocated(self%charges)) then
      write(unit, '(",")', advance='no')
      if (present(indentation)) write(unit, '(/,a)', advance='no') repeat(indentation, 1)
      write(unit, jsonkey, advance='no') 'partial charges'
      call write_json_array(unit, self%charges, indentation)
   end if
   if (allocated(self%dipole_moment)) then
      write(unit, '(",")', advance='no')
      if (present(indentation)) write(unit, '(/,a)', advance='no') repeat(indentation, 1)
      write(unit, jsonkey, advance='no') 'dipole moment'
      call write_json_array(unit, self%dipole_moment, indentation)
   end if
   if (allocated(self%polarizibilities)) then
      write(unit, '(",")', advance='no')
      if (present(indentation)) write(unit, '(/,a)', advance='no') repeat(indentation, 1)
      write(unit, jsonkey, advance='no') 'polarizibilities'
      call write_json_array(unit, self%polarizibilities, indentation)
   end if
   if (allocated(self%c6_coefficients)) then
      write(unit, '(",")', advance='no')
      if (present(indentation)) write(unit, '(/,a)', advance='no') repeat(indentation, 1)
      write(unit, jsonkey, advance='no') 'c6 coefficients'
      array = reshape(self%c6_coefficients, [product(shape(self%c6_coefficients))])
      call write_json_array(unit, array, indentation)
   end if
   if (allocated(self%hessian)) then
      write(unit, '(",")', advance='no')
      if (present(indentation)) write(unit, '(/,a)', advance='no') repeat(indentation, 1)
      write(unit, jsonkey, advance='no') 'hessian'
      array = reshape(self%hessian, [product(shape(self%hessian))])
      call write_json_array(unit, array, indentation)
   end if
   if (present(indentation)) write(unit, '(/)', advance='no')
   write(unit, '("}")')

end subroutine write_json

subroutine write_json_array(unit, array, indent)
   integer, intent(in) :: unit
   real(wp), intent(in) :: array(:)
   character(len=*), intent(in), optional :: indent
   integer :: i
   write(unit, '("[")', advance='no')
   do i = 1, size(array)
      if (present(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 2)
      write(unit, '(1x,es25.16)', advance='no') array(i)
      if (i /= size(array)) write(unit, '(",")', advance='no')
   end do
   if (present(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
   write(unit, '("]")', advance='no')
end subroutine write_json_array

subroutine write_toml(self, unit, indentation)
   include 'dftd4_version.fh'
   class(dftd_results), intent(in) :: self
   integer, intent(in) :: unit
   character(len=*), intent(in), optional :: indentation
   character(len=*), parameter :: tomlkey = "('""',a,'"" = ')"
   real(wp), allocatable :: array(:)

   write(unit, tomlkey, advance='no') 'version'
   write(unit, '(a)') '"'//version//'"'
   if (allocated(self%energy)) then
      write(unit, tomlkey, advance='no') 'energy'
      write(unit, '(es25.16)') self%energy
   end if
   if (allocated(self%stress)) then
      write(unit, tomlkey, advance='no') 'stress tensor'
      array = reshape(self%stress, [product(shape(self%stress))])
      call write_json_array(unit, array, indentation)
      write(unit, '(a)')
   end if
   if (allocated(self%lattice_gradient)) then
      write(unit, tomlkey, advance='no') 'lattice gradient'
      array = reshape(self%lattice_gradient, [product(shape(self%lattice_gradient))])
      call write_json_array(unit, array, indentation)
      write(unit, '(a)')
   end if
   if (allocated(self%gradient)) then
      write(unit, tomlkey, advance='no') 'gradient'
      array = reshape(self%gradient, [product(shape(self%gradient))])
      call write_json_array(unit, array, indentation)
      write(unit, '(a)')
   end if
   if (allocated(self%charges)) then
      write(unit, tomlkey, advance='no') 'partial charges'
      call write_json_array(unit, self%charges, indentation)
      write(unit, '(a)')
   end if
   if (allocated(self%dipole_moment)) then
      write(unit, tomlkey, advance='no') 'dipole moment'
      call write_json_array(unit, self%dipole_moment, indentation)
      write(unit, '(a)')
   end if
   if (allocated(self%polarizibilities)) then
      write(unit, tomlkey, advance='no') 'polarizibilities'
      call write_json_array(unit, self%polarizibilities, indentation)
      write(unit, '(a)')
   end if
   if (allocated(self%c6_coefficients)) then
      write(unit, tomlkey, advance='no') 'c6 coefficients'
      array = reshape(self%c6_coefficients, [product(shape(self%c6_coefficients))])
      call write_json_array(unit, array, indentation)
      write(unit, '(a)')
   end if
   if (allocated(self%hessian)) then
      if (present(indentation)) write(unit, '(/,a)', advance='no') repeat(indentation, 1)
      write(unit, tomlkey, advance='no') 'hessian'
      array = reshape(self%hessian, [product(shape(self%hessian))])
      call write_json_array(unit, array, indentation)
      write(unit, '(a)')
   end if

end subroutine write_toml

end module class_results
