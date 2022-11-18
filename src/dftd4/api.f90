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
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Definition of the public C-API of dftd4
!>
!>```c
!>{!./include/dftd4.h!}
!>```
module dftd4_api
   use iso_c_binding
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_structure, only : structure_type, new
   use dftd4_cutoff, only : realspace_cutoff
   use dftd4_damping, only : damping_param
   use dftd4_damping_rational, only : rational_damping_param
   use dftd4_disp, only : get_dispersion, get_pairwise_dispersion, get_properties
   use dftd4_model, only : d4_model, new_d4_model
   use dftd4_numdiff, only: get_dispersion_hessian
   use dftd4_param, only : get_rational_damping
   use dftd4_utils, only : wrap_to_central_cell
   use dftd4_version, only : get_dftd4_version
   implicit none
   private

   public :: get_version_api

   public :: vp_error
   public :: new_error_api, check_error_api, get_error_api, delete_error_api

   public :: vp_structure
   public :: new_structure_api, delete_structure_api, update_structure_api

   public :: vp_model
   public :: new_d4_model_api, custom_d4_model_api, delete_model_api

   public :: vp_param
   public :: new_rational_damping_api , load_rational_damping_api
   public :: delete_param_api

   public :: get_dispersion_api, get_pairwise_dispersion_api, get_properties_api

   !> Namespace for C routines
   character(len=*), parameter :: namespace = "dftd4_"

   !> Void pointer to error handle
   type :: vp_error
      !> Actual payload
      type(error_type), allocatable :: ptr
   end type vp_error

   !> Void pointer to molecular structure data
   type :: vp_structure
      !> Actual payload
      type(structure_type) :: ptr
   end type vp_structure

   !> Void pointer to dispersion model
   type :: vp_model
      !> Actual payload
      type(d4_model) :: ptr
   end type vp_model

   !> Void pointer to damping parameters
   type :: vp_param
      !> Actual payload
      class(damping_param), allocatable :: ptr
   end type vp_param


   logical, parameter :: debug = .false.


contains


!> Obtain library version as major * 10000 + minor + 100 + patch
function get_version_api() result(version) &
      & bind(C, name=namespace//"get_version")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_version_api
   integer(c_int) :: version
   integer :: major, minor, patch

   call get_dftd4_version(major, minor, patch)
   version = 10000_c_int * major + 100_c_int * minor + patch

end function get_version_api


!> Create new error handle object
function new_error_api() &
      & result(verror) &
      & bind(C, name=namespace//"new_error")
   !DEC$ ATTRIBUTES DLLEXPORT :: new_error_api
   type(vp_error), pointer :: error
   type(c_ptr) :: verror

   if (debug) print'("[Info]",1x, a)', "new_error"

   allocate(error)
   verror = c_loc(error)

end function new_error_api


!> Delete error handle object
subroutine delete_error_api(verror) &
      & bind(C, name=namespace//"delete_error")
   !DEC$ ATTRIBUTES DLLEXPORT :: delete_error_api
   type(c_ptr), intent(inout) :: verror
   type(vp_error), pointer :: error

   if (debug) print'("[Info]",1x, a)', "delete_error"

   if (c_associated(verror)) then
      call c_f_pointer(verror, error)

      deallocate(error)
      verror = c_null_ptr
   end if

end subroutine delete_error_api


!> Check error handle status
function check_error_api(verror) result(status) &
      & bind(C, name=namespace//"check_error")
   !DEC$ ATTRIBUTES DLLEXPORT :: check_error_api
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   integer(c_int) :: status

   if (debug) print'("[Info]",1x, a)', "check_error"

   if (c_associated(verror)) then
      call c_f_pointer(verror, error)

      if (allocated(error%ptr)) then
         status = 1
      else
         status = 0
      end if
   else
      status = 2
   end if

end function check_error_api


!> Get error message from error handle
subroutine get_error_api(verror, charptr, buffersize) &
      & bind(C, name=namespace//"get_error")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_error_api
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   character(kind=c_char), intent(inout) :: charptr(*)
   integer(c_int), intent(in), optional :: buffersize
   integer :: max_length

   if (debug) print'("[Info]",1x, a)', "get_error"

   if (c_associated(verror)) then
      call c_f_pointer(verror, error)

      if (present(buffersize)) then
         max_length = buffersize
      else
         max_length = huge(max_length) - 2
      end if

      if (allocated(error%ptr)) then
         call f_c_character(error%ptr%message, charptr, max_length)
      end if
   end if

end subroutine get_error_api


!> Create new molecular structure data (quantities in Bohr)
function new_structure_api(verror, natoms, numbers, positions, c_charge, &
      & c_lattice, c_periodic) result(vmol) &
      & bind(C, name=namespace//"new_structure")
   !DEC$ ATTRIBUTES DLLEXPORT :: new_structure_api
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   integer(c_int), value, intent(in) :: natoms
   integer(c_int), intent(in) :: numbers(natoms)
   real(c_double), intent(in) :: positions(3, natoms)
   real(c_double), intent(in), optional :: c_charge
   real(wp), allocatable :: charge
   real(c_double), intent(in), optional :: c_lattice(3, 3)
   real(wp), allocatable :: lattice(:, :)
   logical(c_bool), intent(in), optional :: c_periodic(3)
   logical, allocatable :: periodic(:)
   type(vp_structure), pointer :: mol
   type(c_ptr) :: vmol

   if (debug) print'("[Info]",1x, a)', "new_structure"

   vmol = c_null_ptr

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (present(c_lattice)) then
      allocate(lattice(3, 3))
      lattice(:, :) = c_lattice
   end if
   if (present(c_periodic)) then
      allocate(periodic(3))
      periodic(:) = c_periodic
   end if
   if (present(c_charge)) then
      charge = c_charge
   end if

   allocate(mol)
   call new(mol%ptr, numbers, positions, lattice=lattice, periodic=periodic, &
      charge=charge)
   vmol = c_loc(mol)

   call wrap_to_central_cell(mol%ptr%xyz, mol%ptr%lattice, mol%ptr%periodic)

   call verify_structure(error%ptr, mol%ptr)

end function new_structure_api


!> Delete molecular structure data
subroutine delete_structure_api(vmol) &
      & bind(C, name=namespace//"delete_structure")
   !DEC$ ATTRIBUTES DLLEXPORT :: delete_structure_api
   type(c_ptr), intent(inout) :: vmol
   type(vp_structure), pointer :: mol

   if (debug) print'("[Info]",1x, a)', "delete_structure"

   if (c_associated(vmol)) then
      call c_f_pointer(vmol, mol)

      deallocate(mol)
      vmol = c_null_ptr
   end if

end subroutine delete_structure_api


!> Update coordinates and lattice parameters (quantities in Bohr)
subroutine update_structure_api(verror, vmol, positions, lattice) &
      & bind(C, name=namespace//"update_structure")
   !DEC$ ATTRIBUTES DLLEXPORT :: update_structure_api
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   real(c_double), intent(in) :: positions(3, *)
   real(c_double), intent(in), optional :: lattice(3, 3)

   if (debug) print'("[Info]",1x, a)', "update_structure"

   if (.not.c_associated(verror)) then
      return
   end if
   call c_f_pointer(verror, error)

   if (.not.c_associated(vmol)) then
      call fatal_error(error%ptr, "Molecular structure data is missing")
      return
   end if
   call c_f_pointer(vmol, mol)

   if (mol%ptr%nat <= 0 .or. mol%ptr%nid <= 0 .or. .not.allocated(mol%ptr%num) &
      & .or. .not.allocated(mol%ptr%id) .or. .not.allocated(mol%ptr%xyz)) then
      call fatal_error(error%ptr, "Invalid molecular structure data provided")
      return
   end if

   mol%ptr%xyz(:, :) = positions(:3, :mol%ptr%nat)
   if (present(lattice)) then
      mol%ptr%lattice(:, :) = lattice(:3, :3)
   end if

   call wrap_to_central_cell(mol%ptr%xyz, mol%ptr%lattice, mol%ptr%periodic)

   call verify_structure(error%ptr, mol%ptr)

end subroutine update_structure_api


!> Create new D4 dispersion model
function new_d4_model_api(verror, vmol) &
      & result(vdisp) &
      & bind(C, name=namespace//"new_d4_model")
   !DEC$ ATTRIBUTES DLLEXPORT :: new_d4_model_api
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr) :: vdisp
   type(vp_model), pointer :: disp

   if (debug) print'("[Info]",1x, a)', "new_d4_model"

   vdisp = c_null_ptr

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vmol)) then
      call fatal_error(error%ptr, "Molecular structure data is missing")
      return
   end if
   call c_f_pointer(vmol, mol)

   allocate(disp)
   call new_d4_model(disp%ptr, mol%ptr)
   vdisp = c_loc(disp)

end function new_d4_model_api


!> Create new D4 dispersion model
function custom_d4_model_api(verror, vmol, ga, gc, wf) &
      & result(vdisp) &
      & bind(C, name=namespace//"custom_d4_model")
   !DEC$ ATTRIBUTES DLLEXPORT :: custom_d4_model_api
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr) :: vdisp
   type(vp_model), pointer :: disp
   real(c_double), value, intent(in) :: ga
   real(c_double), value, intent(in) :: gc
   real(c_double), value, intent(in) :: wf

   if (debug) print'("[Info]",1x, a)', "custom_d4_model"

   vdisp = c_null_ptr

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vmol)) then
      call fatal_error(error%ptr, "Molecular structure data is missing")
      return
   end if
   call c_f_pointer(vmol, mol)

   allocate(disp)
   call new_d4_model(disp%ptr, mol%ptr, ga=ga, gc=gc, wf=wf)
   vdisp = c_loc(disp)

end function custom_d4_model_api


!> Delete dispersion model
subroutine delete_model_api(vdisp) &
      & bind(C, name=namespace//"delete_model")
   !DEC$ ATTRIBUTES DLLEXPORT :: delete_model_api
   type(c_ptr), intent(inout) :: vdisp
   type(vp_model), pointer :: disp

   if (debug) print'("[Info]",1x, a)', "delete_model"

   if (c_associated(vdisp)) then
      call c_f_pointer(vdisp, disp)

      deallocate(disp)
      vdisp = c_null_ptr
   end if

end subroutine delete_model_api


!> Create new rational damping parameters
function new_rational_damping_api(verror, s6, s8, s9, a1, a2, alp) &
      & result(vparam) &
      & bind(C, name=namespace//"new_rational_damping")
   !DEC$ ATTRIBUTES DLLEXPORT :: new_rational_damping_api
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   real(c_double), value, intent(in) :: s6
   real(c_double), value, intent(in) :: s8
   real(c_double), value, intent(in) :: s9
   real(c_double), value, intent(in) :: a1
   real(c_double), value, intent(in) :: a2
   real(c_double), value, intent(in) :: alp
   type(c_ptr) :: vparam
   type(rational_damping_param), allocatable :: tmp
   type(vp_param), pointer :: param

   if (debug) print'("[Info]",1x, a)', "new_rational_damping"

   vparam = c_null_ptr

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   allocate(tmp)
   tmp = rational_damping_param(s6=s6, s8=s8, s9=s9, a1=a1, a2=a2, alp=alp)

   allocate(param)
   call move_alloc(tmp, param%ptr)
   vparam = c_loc(param)

end function new_rational_damping_api


!> Load rational damping parameters from internal storage
function load_rational_damping_api(verror, charptr, atm) &
      & result(vparam) &
      & bind(C, name=namespace//"load_rational_damping")
   !DEC$ ATTRIBUTES DLLEXPORT :: load_rational_damping_api
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   character(kind=c_char), intent(in) :: charptr(*)
   logical(c_bool), value, intent(in) :: atm
   character(len=:, kind=c_char), allocatable :: method
   type(c_ptr) :: vparam
   type(vp_param), pointer :: param
   real(wp), allocatable :: s9

   if (debug) print'("[Info]",1x, a)', "load_rational_damping"

   vparam = c_null_ptr

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   call c_f_character(charptr, method)

   if (atm) s9 = 1.0_wp
   allocate(param)
   call get_rational_damping(method, param%ptr, s9)
   if (.not.allocated(param%ptr)) then
      call fatal_error(error%ptr, "Functional '"//method//"' not known")
      return
   end if

   vparam = c_loc(param)

end function load_rational_damping_api


!> Delete damping parameters
subroutine delete_param_api(vparam) &
      & bind(C, name=namespace//"delete_param")
   !DEC$ ATTRIBUTES DLLEXPORT :: delete_param_api
   type(c_ptr), intent(inout) :: vparam
   type(vp_param), pointer :: param

   if (debug) print'("[Info]",1x, a)', "delete_param"

   if (c_associated(vparam)) then
      call c_f_pointer(vparam, param)

      deallocate(param)
      vparam = c_null_ptr
   end if

end subroutine delete_param_api


!> Calculate dispersion
subroutine get_dispersion_api(verror, vmol, vdisp, vparam, &
      & energy, c_gradient, c_sigma) &
      & bind(C, name=namespace//"get_dispersion")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_dispersion_api
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vdisp
   type(vp_model), pointer :: disp
   type(c_ptr), value :: vparam
   type(vp_param), pointer :: param
   real(c_double), intent(out) :: energy
   real(c_double), intent(out), optional :: c_gradient(3, *)
   real(wp), allocatable :: gradient(:, :)
   real(c_double), intent(out), optional :: c_sigma(3, 3)
   real(wp), allocatable :: sigma(:, :)
   logical :: has_grad, has_sigma


   if (debug) print'("[Info]",1x, a)', "get_dispersion"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vmol)) then
      call fatal_error(error%ptr, "Molecular structure data is missing")
      return
   end if
   call c_f_pointer(vmol, mol)

   if (.not.c_associated(vdisp)) then
      call fatal_error(error%ptr, "Dispersion model is missing")
      return
   end if
   call c_f_pointer(vdisp, disp)

   if (.not.c_associated(vparam)) then
      call fatal_error(error%ptr, "Damping parameters are missing")
      return
   end if
   call c_f_pointer(vparam, param)

   if (.not.allocated(param%ptr)) then
      call fatal_error(error%ptr, "Damping parameters are not initialized")
      return
   end if

   has_grad = present(c_gradient) 
   if (has_grad) then
      gradient = c_gradient(:3, :mol%ptr%nat)
   endif

   has_sigma = present(c_sigma) 
   if (has_sigma) then
      sigma = c_sigma(:3, :3)
   ! Still needs to be passed into dispersion subroutines,
   ! just won't be returned through the API. 
   ! Would need to refactor disperision
   ! subroutines to make sigma truly optional. 
   else if (has_grad) then
      allocate(sigma(3,3)) 
   endif

   ! Evaluate energy, gradient (optional), and 
   ! sigma (optional) analytically
   call get_dispersion(mol%ptr, disp%ptr, param%ptr, realspace_cutoff(), &
      & energy, gradient, sigma)

   if (has_grad) then
      c_gradient(:3, :mol%ptr%nat) = gradient
   endif

   if (has_sigma) then
      c_sigma(:3, :3) = sigma
   endif

end subroutine get_dispersion_api

!> Calculate hessian numerically
subroutine get_numerical_hessian_api(verror, vmol, vdisp, & 
                                   & vparam, c_hessian) &
      & bind(C, name=namespace//"get_numerical_hessian")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_numerical_hessian_api
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vdisp
   type(vp_model), pointer :: disp
   type(c_ptr), value :: vparam
   type(vp_param), pointer :: param
   real(c_double), intent(out) :: c_hessian(*)
   real(wp), allocatable :: hessian(:, :, :, :)
   integer :: nat_sq


   if (debug) print'("[Info]",1x, a)', "get_numerical_hessian"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vmol)) then
      call fatal_error(error%ptr, "Molecular structure data is missing")
      return
   end if
   call c_f_pointer(vmol, mol)
   nat_sq = mol%ptr%nat*mol%ptr%nat

   if (.not.c_associated(vdisp)) then
      call fatal_error(error%ptr, "Dispersion model is missing")
      return
   end if
   call c_f_pointer(vdisp, disp)

   if (.not.c_associated(vparam)) then
      call fatal_error(error%ptr, "Damping parameters are missing")
      return
   end if
   call c_f_pointer(vparam, param)

   if (.not.allocated(param%ptr)) then
      call fatal_error(error%ptr, "Damping parameters are not initialized")
      return
   end if

   ! Evaluate hessian numerically 
   hessian = reshape(c_hessian(:9*nat_sq), &
                    &(/3, mol%ptr%nat, 3, mol%ptr%nat/))
   call get_dispersion_hessian(mol%ptr, disp%ptr, param%ptr, &
    & realspace_cutoff(), hessian)
   c_hessian(:9*nat_sq) = reshape(hessian, (/9*nat_sq/))

end subroutine get_numerical_hessian_api

!> Calculate pairwise representation of dispersion energy
subroutine get_pairwise_dispersion_api(verror, vmol, vdisp, vparam, &
      & c_pair_energy2, c_pair_energy3) &
      & bind(C, name=namespace//"get_pairwise_dispersion")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_pairwise_dispersion_api
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vdisp
   type(vp_model), pointer :: disp
   type(c_ptr), value :: vparam
   type(vp_param), pointer :: param
   type(c_ptr), value, intent(in) :: c_pair_energy2
   real(wp), pointer :: pair_energy2(:, :)
   type(c_ptr), value, intent(in) :: c_pair_energy3
   real(wp), pointer :: pair_energy3(:, :)

   if (debug) print'("[Info]",1x, a)', "get_pairwise_dispersion"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vmol)) then
      call fatal_error(error%ptr, "Molecular structure data is missing")
      return
   end if
   call c_f_pointer(vmol, mol)

   if (.not.c_associated(vdisp)) then
      call fatal_error(error%ptr, "Dispersion model is missing")
      return
   end if
   call c_f_pointer(vdisp, disp)

   if (.not.c_associated(vparam)) then
      call fatal_error(error%ptr, "Damping parameters are missing")
      return
   end if
   call c_f_pointer(vparam, param)

   if (.not.allocated(param%ptr)) then
      call fatal_error(error%ptr, "Damping parameters are not initialized")
      return
   end if

   call c_f_pointer(c_pair_energy2, pair_energy2, [mol%ptr%nat, mol%ptr%nat])
   call c_f_pointer(c_pair_energy3, pair_energy3, [mol%ptr%nat, mol%ptr%nat])

   call get_pairwise_dispersion(mol%ptr, disp%ptr, param%ptr, realspace_cutoff(), &
      & pair_energy2, pair_energy3)

end subroutine get_pairwise_dispersion_api


!> Calculate dispersion
subroutine get_properties_api(verror, vmol, vdisp, &
      & c_cn, c_charges, c_c6, c_alpha) &
      & bind(C, name=namespace//"get_properties")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_properties_api
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vdisp
   type(vp_model), pointer :: disp
   real(c_double), intent(out), optional :: c_cn(*)
   real(wp), allocatable :: cn(:)
   real(c_double), intent(out), optional :: c_charges(*)
   real(wp), allocatable :: charges(:)
   real(c_double), intent(out), optional :: c_c6(*)
   real(wp), allocatable :: c6(:, :)
   real(c_double), intent(out), optional :: c_alpha(*)
   real(wp), allocatable :: alpha(:)

   if (debug) print'("[Info]",1x, a)', "get_properties"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vmol)) then
      call fatal_error(error%ptr, "Molecular structure data is missing")
      return
   end if
   call c_f_pointer(vmol, mol)

   if (.not.c_associated(vdisp)) then
      call fatal_error(error%ptr, "Dispersion model is missing")
      return
   end if
   call c_f_pointer(vdisp, disp)

   allocate(cn(mol%ptr%nat), charges(mol%ptr%nat), alpha(mol%ptr%nat), &
      & c6(mol%ptr%nat, mol%ptr%nat))
   call get_properties(mol%ptr, disp%ptr, realspace_cutoff(), cn, charges, c6, alpha)

   if (present(c_cn)) then
      c_cn(:size(cn)) = cn
   end if

   if (present(c_charges)) then
      c_charges(:size(charges)) = charges
   end if

   if (present(c_c6)) then
      c_c6(:size(c6)) = reshape(c6, [size(c6)])
   end if

   if (present(c_alpha)) then
      c_alpha(:size(alpha)) = alpha
   end if

end subroutine get_properties_api


subroutine f_c_character(rhs, lhs, len)
   character(kind=c_char), intent(out) :: lhs(*)
   character(len=*), intent(in) :: rhs
   integer, intent(in) :: len
   integer :: length
   length = min(len-1, len_trim(rhs))

   lhs(1:length) = transfer(rhs(1:length), lhs(1:length))
   lhs(length+1:length+1) = c_null_char

end subroutine f_c_character


subroutine c_f_character(rhs, lhs)
   character(kind=c_char), intent(in) :: rhs(*)
   character(len=:, kind=c_char), allocatable, intent(out) :: lhs

   integer :: ii

   do ii = 1, huge(ii) - 1
      if (rhs(ii) == c_null_char) then
         exit
      end if
   end do
   allocate(character(len=ii-1) :: lhs)
   lhs = transfer(rhs(1:ii-1), lhs)

end subroutine c_f_character


!> Cold fusion check
subroutine verify_structure(error, mol)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(in) :: mol
   integer :: iat, jat, stat
   stat = 0
   do iat = 1, mol%nat
      do jat = 1, iat - 1
         if (norm2(mol%xyz(:, jat) - mol%xyz(:, iat)) < 1.0e-9_wp) stat = stat + 1
      end do
   end do
   if (stat > 0) then
      call fatal_error(error, "Too close interatomic distances found")
   end if
end subroutine verify_structure


end module dftd4_api
