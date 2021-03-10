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
   use dftd4_blas, only : d4_gemv
   use dftd4_charge, only : get_charges
   use dftd4_cutoff, only : realspace_cutoff, get_lattice_points
   use dftd4_damping, only : damping_param
   use dftd4_data, only : get_covalent_rad
   use dftd4_model, only : d4_model
   use dftd4_ncoord, only : get_coordination_number
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa
   implicit none
   private

   public :: get_dispersion, get_properties, get_pairwise_dispersion


contains


!> Wrapper to handle the evaluation of dispersion energy and derivatives
subroutine get_dispersion(mol, disp, param, cutoff, energy, gradient, sigma)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d4_model), intent(in) :: disp

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
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: q(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: gwvec(:, :), gwdcn(:, :), gwdq(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: dEdcn(:), dEdq(:), energies(:)
   real(wp), allocatable :: lattr(:, :)

   mref = maxval(disp%ref)
   grad = present(gradient).and.present(sigma)

   allocate(cn(mol%nat))
   if (grad) allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, &
      & cn, dcndr, dcndL)

   allocate(q(mol%nat))
   if (grad) allocate(dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat))
   call get_charges(mol, q, dqdr, dqdL)

   allocate(gwvec(mref, mol%nat))
   if (grad) allocate(gwdcn(mref, mol%nat), gwdq(mref, mol%nat))
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
      call d4_gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
      call d4_gemv(dcndL, dEdcn, sigma, beta=1.0_wp)
   end if

   energy = sum(energies)

end subroutine get_dispersion


!> Wrapper to handle the evaluation of properties related to this dispersion model
subroutine get_properties(mol, disp, cutoff, cn, q, c6, alpha)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d4_model), intent(in) :: disp

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Coordination number
   real(wp), intent(out) :: cn(:)

   !> Atomic partial charges
   real(wp), intent(out) :: q(:)

   !> C6 coefficients
   real(wp), intent(out) :: c6(:, :)

   !> Static polarizibilities
   real(wp), intent(out) :: alpha(:)

   integer :: mref
   real(wp), allocatable :: gwvec(:, :), lattr(:, :)

   mref = maxval(disp%ref)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, cn)

   call get_charges(mol, q)

   allocate(gwvec(mref, mol%nat))
   call disp%weight_references(mol, cn, q, gwvec)

   call disp%get_atomic_c6(mol, gwvec, c6=c6)
   call disp%get_polarizibilities(mol, gwvec, alpha=alpha)

end subroutine get_properties


!> Wrapper to handle the evaluation of pairwise representation of the dispersion energy
subroutine get_pairwise_dispersion(mol, disp, param, cutoff, energy2, energy3)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d4_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Pairwise representation of additive dispersion energy
   real(wp), intent(out) :: energy2(:, :)

   !> Pairwise representation of non-additive dispersion energy
   real(wp), intent(out) :: energy3(:, :)

   integer :: mref
   real(wp), allocatable :: cn(:), q(:), gwvec(:, :), c6(:, :), lattr(:, :)

   mref = maxval(disp%ref)

   allocate(cn(mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, cn)

   allocate(q(mol%nat))
   call get_charges(mol, q)

   allocate(gwvec(mref, mol%nat))
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
