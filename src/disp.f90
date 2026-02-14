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

!> High-level wrapper to obtain the dispersion energy for a DFT-D4 calculation
module dftd4_disp
   use, intrinsic :: iso_fortran_env, only : error_unit
   use dftd4_blas, only : d4_gemv
   use dftd4_cutoff, only : realspace_cutoff, get_lattice_points
   use dftd4_damping, only : damping_param
   use dftd4_data, only : get_covalent_rad
   use dftd4_model, only : dispersion_model
   use dftd4_ncoord, only : get_coordination_number, add_coordination_number_derivs
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa
   use multicharge, only : get_charges
   implicit none
   private

   public :: get_dispersion, get_properties, get_pairwise_dispersion


contains


!> Wrapper to handle the evaluation of dispersion energy and derivatives
subroutine get_dispersion(mol, disp, param, cutoff, energy, gradient, sigma)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_dispersion

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Dispersion energy
   real(wp), intent(out) :: energy

   !> Dispersion gradient
   real(wp), intent(out), contiguous, optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(out), contiguous, optional :: sigma(:, :)
 
   logical :: grad
   integer :: mref
   real(wp), allocatable :: cn(:)
   real(wp), allocatable :: q(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: gwvec(:, :, :), gwdcn(:, :, :), gwdq(:, :, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: dEdcn(:), dEdq(:), energies(:)
   real(wp), allocatable :: lattr(:, :)
   type(error_type), allocatable :: error

   mref = maxval(disp%ref)
   grad = present(gradient).or.present(sigma)

   if (.not. allocated(disp%mchrg)) then
      write(error_unit, '("[Error]:", 1x, a)') "Not supported for non-self-consistent D4 version"
      error stop
   end if

   allocate(cn(mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, cn)

   allocate(q(mol%nat))
   if (grad) allocate(dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat))
   call get_charges(disp%mchrg, mol, error, q, dqdr, dqdL)
   if(allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   allocate(gwvec(mref, mol%nat, disp%ncoup))
   if (grad) allocate(gwdcn(mref, mol%nat, disp%ncoup), gwdq(mref, mol%nat, disp%ncoup))
   call disp%weight_references(mol, cn, q, gwvec, gwdcn, gwdq)

   allocate(c6(mol%nat, mol%nat))
   if (grad) allocate(dc6dcn(mol%nat, mol%nat), dc6dq(mol%nat, mol%nat))
   call disp%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   allocate(energies(mol%nat))
   energies(:) = 0.0_wp
   if (grad) then
      allocate(dEdcn(mol%nat), dEdq(mol%nat))
      dEdcn(:) = 0.0_wp
      dEdq(:) = 0.0_wp
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
   end if

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
   call param%get_dispersion2(mol, lattr, cutoff%disp2, disp%r4r2, &
      & c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, sigma)
   if (grad) then
      call d4_gemv(dqdr, dEdq, gradient, beta=1.0_wp)
      call d4_gemv(dqdL, dEdq, sigma, beta=1.0_wp)
   end if

   q(:) = 0.0_wp
   call disp%weight_references(mol, cn, q, gwvec, gwdcn, gwdq)
   call disp%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call param%get_dispersion3(mol, lattr, cutoff%disp3, disp%r4r2, &
      & c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, sigma)
   if (grad) then
      call add_coordination_number_derivs(mol, lattr, cutoff%cn, &
         & disp%rcov, disp%en, dEdcn, gradient, sigma)
   end if

   energy = sum(energies)

end subroutine get_dispersion


!> Wrapper to handle the evaluation of properties related to this dispersion model
subroutine get_properties(mol, disp, cutoff, cn, q, c6, alpha)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_properties

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: disp

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Coordination number
   real(wp), intent(out) :: cn(:)

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: q(:)

   !> C6 coefficients
   real(wp), intent(out) :: c6(:, :)

   !> Static polarizabilities
   real(wp), intent(out) :: alpha(:)

   integer :: mref
   real(wp), allocatable :: gwvec(:, :, :), lattr(:, :)
   type(error_type), allocatable :: error

   if (.not. allocated(disp%mchrg)) then
      write(error_unit, '("[Error]:", 1x, a)') "Not supported for non-self-consistent D4 version"
      error stop
   end if

   mref = maxval(disp%ref)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, cn)

   call get_charges(disp%mchrg, mol, error, q)
   if(allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   allocate(gwvec(mref, mol%nat, disp%ncoup))
   call disp%weight_references(mol, cn, q, gwvec)

   call disp%get_atomic_c6(mol, gwvec, c6=c6)
   call disp%get_polarizabilities(mol, gwvec, alpha=alpha)

end subroutine get_properties


!> Wrapper to handle the evaluation of pairwise representation of the dispersion energy
subroutine get_pairwise_dispersion(mol, disp, param, cutoff, energy2, energy3)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_pairwise_dispersion

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Pairwise representation of additive dispersion energy
   real(wp), intent(out) :: energy2(:, :)

   !> Pairwise representation of non-additive dispersion energy
   real(wp), intent(out) :: energy3(:, :)

   integer :: mref
   real(wp), allocatable :: cn(:), q(:), gwvec(:, :, :), c6(:, :), lattr(:, :)
   type(error_type), allocatable :: error

   if (.not. allocated(disp%mchrg)) then
      write(error_unit, '("[Error]:", 1x, a)') "Not supported for non-self-consistent D4 version"
      error stop
   end if

   mref = maxval(disp%ref)

   allocate(cn(mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, cn)

   allocate(q(mol%nat))
   call get_charges(disp%mchrg, mol, error, q)
   if(allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   allocate(gwvec(mref, mol%nat, disp%ncoup))
   call disp%weight_references(mol, cn, q, gwvec)

   allocate(c6(mol%nat, mol%nat))
   call disp%get_atomic_c6(mol, gwvec, c6=c6)

   energy2(:, :) = 0.0_wp
   energy3(:, :) = 0.0_wp
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
   call param%get_pairwise_dispersion2(mol, lattr, cutoff%disp2, disp%r4r2, &
      & c6, energy2)

   q(:) = 0.0_wp
   call disp%weight_references(mol, cn, q, gwvec)
   call disp%get_atomic_c6(mol, gwvec, c6=c6)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call param%get_pairwise_dispersion3(mol, lattr, cutoff%disp3, disp%r4r2, &
      & c6, energy3)

end subroutine get_pairwise_dispersion


end module dftd4_disp
