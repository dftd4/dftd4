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
   use dftd4_damping_rational, only : rational_damping_param
   use dftd4_data, only : get_covalent_rad
   use dftd4_model, only : dispersion_model, d4_model, d4s_model
   use dftd4_ncoord, only : get_coordination_number, add_coordination_number_derivs, &
      & add_coordination_number_hessian
   use dftd4_partition, only : work_partition
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa
   use multicharge, only : get_charges
   implicit none
   private

   public :: get_dispersion, get_dispersion2, get_dispersion3_hessian, get_properties, get_pairwise_dispersion


contains


!> Wrapper to handle the evaluation of dispersion energy and derivatives
subroutine get_dispersion(mol, disp, param, cutoff, energy, gradient, sigma, partition)
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

   !> Optional externally assigned work partition
   type(work_partition), intent(in), optional :: partition

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
   call param%get_dispersion2(mol, lattr, cutoff%disp2, cutoff%width2, &
      & disp%r4r2, c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, &
      & sigma, partition)
   if (grad) then
      call d4_gemv(dqdr, dEdq, gradient, beta=1.0_wp)
      call d4_gemv(dqdL, dEdq, sigma, beta=1.0_wp)
   end if

   q(:) = 0.0_wp
   call disp%weight_references(mol, cn, q, gwvec, gwdcn, gwdq)
   call disp%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call param%get_dispersion3(mol, lattr, cutoff%disp3, cutoff%width3, &
      & disp%r4r2, c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, &
      & sigma, partition)
   if (grad) then
      call add_coordination_number_derivs(mol, lattr, cutoff%cn, &
         & disp%rcov, disp%en, dEdcn, gradient, sigma)
   end if

   energy = sum(energies)

end subroutine get_dispersion


!> Wrapper to handle the evaluation of the two-body dispersion contribution
subroutine get_dispersion2(mol, disp, param, cutoff, energy, gradient, sigma, partition)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Two-body dispersion energy
   real(wp), intent(out) :: energy

   !> Two-body dispersion gradient
   real(wp), intent(out), contiguous, optional :: gradient(:, :)

   !> Two-body dispersion virial
   real(wp), intent(out), contiguous, optional :: sigma(:, :)

   !> Optional externally assigned work partition
   type(work_partition), intent(in), optional :: partition

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
      write(error_unit, '("[Error]:", 1x, a)') &
         & "Not supported for non-self-consistent D4 version"
      error stop
   end if

   allocate(cn(mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, cn)

   allocate(q(mol%nat))
   if (grad) allocate(dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat))
   call get_charges(disp%mchrg, mol, error, q, dqdr, dqdL)
   if (allocated(error)) then
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
   call param%get_dispersion2(mol, lattr, cutoff%disp2, cutoff%width2, &
      & disp%r4r2, c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, &
      & sigma, partition)
   if (grad) then
      call d4_gemv(dqdr, dEdq, gradient, beta=1.0_wp)
      call d4_gemv(dqdL, dEdq, sigma, beta=1.0_wp)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
      call add_coordination_number_derivs(mol, lattr, cutoff%cn, &
         & disp%rcov, disp%en, dEdcn, gradient, sigma)
   end if

   energy = sum(energies)

end subroutine get_dispersion2


!> Analytical Hessian of the ATM three-body contribution
subroutine get_dispersion3_hessian(error, mol, disp, param, cutoff, hessian, partition)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_dispersion3_hessian

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> ATM dispersion Hessian
   real(wp), intent(out) :: hessian(:, :)

   !> Optional externally assigned work partition
   type(work_partition), intent(in), optional :: partition

   integer :: mref, nat, ndim, ii, jj, icn, jcn
   real(wp), allocatable, target :: dcndr(:, :, :)
   real(wp), allocatable :: cn(:), dcndL(:, :, :)
   real(wp), allocatable :: q(:)
   real(wp), allocatable :: gwvec(:, :, :), gwdcn(:, :, :), gwd2cn(:, :, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: d2c6dcn2(:, :), d2c6dcnij(:, :)
   real(wp), allocatable :: dEdcn(:), dEdcndr(:, :), dEdcndcn(:, :), cnwork(:, :)
   real(wp), allocatable :: drt(:, :), dEdcndrt(:, :)
   real(wp), pointer :: dr(:, :)
   real(wp) :: tmp
   real(wp), allocatable :: lattr(:, :)

   nat = mol%nat
   ndim = 3*nat
   mref = maxval(disp%ref)

   ! mctc-lib currently evaluates Cartesian CN derivatives only when both
   ! dcndr and dcndL are requested.  Keep dcndL even though the ATM Hessian
   ! itself only consumes dcndr.
   allocate(cn(nat), dcndr(3, nat, nat), dcndL(3, 3, nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, cn, dcndr, dcndL)

   ! The D4 ATM contribution is evaluated with charge-independent C6
   ! coefficients, matching the q=0 branch in get_dispersion.
   allocate(q(nat), source=0.0_wp)
   allocate(gwvec(mref, nat, disp%ncoup), gwdcn(mref, nat, disp%ncoup), &
      & gwd2cn(mref, nat, disp%ncoup))
   allocate(c6(nat, nat), dc6dcn(nat, nat), d2c6dcn2(nat, nat), d2c6dcnij(nat, nat))

   select type(disp)
   type is (d4_model)
      call disp%weight_references_hessian(mol, cn, q, gwvec, gwdcn, gwd2cn)
      call disp%get_atomic_c6_hessian(mol, gwvec, gwdcn, gwd2cn, c6, &
         & dc6dcn, d2c6dcn2, d2c6dcnij)
   type is (d4s_model)
      call disp%weight_references_hessian(mol, cn, q, gwvec, gwdcn, gwd2cn)
      call disp%get_atomic_c6_hessian(mol, gwvec, gwdcn, gwd2cn, c6, &
         & dc6dcn, d2c6dcn2, d2c6dcnij)
   class default
      call fatal_error(error, "Analytical ATM Hessian is only implemented for D4 and D4S models")
      return
   end select

   hessian(:, :) = 0.0_wp
   allocate(dEdcn(nat), source=0.0_wp)
   allocate(dEdcndr(ndim, nat), source=0.0_wp)
   allocate(dEdcndcn(nat, nat), source=0.0_wp)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   select type(param)
   class is (rational_damping_param)
      call param%get_dispersion3_hessian(mol, lattr, cutoff%disp3, cutoff%width3, &
         & disp%r4r2, c6, dc6dcn, d2c6dcn2, d2c6dcnij, hessian, dEdcn, &
         & dEdcndr, dEdcndcn, partition)
   class default
      call fatal_error(error, "Analytical ATM Hessian is only implemented for rational damping")
      return
   end select

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call add_coordination_number_hessian(mol, lattr, cutoff%cn, disp%rcov, disp%en, &
      & dEdcn, hessian)

   dr(1:ndim, 1:nat) => dcndr
   allocate(cnwork(nat, ndim), drt(nat, ndim), dEdcndrt(nat, ndim))

   ! Assemble the CN chain rule explicitly under OpenMP.  Intrinsic MATMUL is
   ! not guaranteed to use the program's OpenMP team (and is often serial with
   ! a single-threaded BLAS), which otherwise leaves a sizeable serial tail in
   ! analytical ATM Hessian calculations.
   !$omp parallel default(none) &
   !$omp shared(nat, ndim, dEdcndr, dEdcndcn, dr, drt, dEdcndrt, cnwork, hessian) &
   !$omp private(ii, jj, icn, jcn, tmp)

   ! Store the Cartesian-by-CN matrices with CN as the contiguous dimension.
   ! Both contractions below then stream through contiguous columns instead of
   ! repeatedly traversing the second dimension of Fortran arrays.
   !$omp do collapse(2) schedule(static)
   do ii = 1, ndim
      do icn = 1, nat
         drt(icn, ii) = dr(ii, icn)
         dEdcndrt(icn, ii) = dEdcndr(ii, icn)
      end do
   end do
   !$omp end do

   ! cnwork = dEdcndcn * transpose(dr).  Keeping jcn as the accumulation
   ! order preserves the scalar result while the icn loop is vectorizable.
   !$omp do schedule(static)
   do jj = 1, ndim
      cnwork(:, jj) = 0.0_wp
      do jcn = 1, nat
         tmp = drt(jcn, jj)
         !$omp simd
         do icn = 1, nat
            cnwork(icn, jj) = cnwork(icn, jj) + dEdcndcn(icn, jcn)*tmp
         end do
      end do
   end do
   !$omp end do

   ! H += dEdcndr*dr^T + dr*dEdcndr^T + dr*dEdcndcn*dr^T
   !$omp do collapse(2) schedule(static)
   do jj = 1, ndim
      do ii = 1, ndim
         tmp = 0.0_wp
         !$omp simd reduction(+:tmp)
         do icn = 1, nat
            tmp = tmp + dEdcndrt(icn, ii)*drt(icn, jj) &
               & + drt(icn, ii)*dEdcndrt(icn, jj) &
               & + drt(icn, ii)*cnwork(icn, jj)
         end do
         hessian(ii, jj) = hessian(ii, jj) + tmp
      end do
   end do
   !$omp end do
   !$omp end parallel

end subroutine get_dispersion3_hessian


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
   call param%get_pairwise_dispersion2(mol, lattr, cutoff%disp2, cutoff%width2, &
      & disp%r4r2, c6, energy2)

   q(:) = 0.0_wp
   call disp%weight_references(mol, cn, q, gwvec)
   call disp%get_atomic_c6(mol, gwvec, c6=c6)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call param%get_pairwise_dispersion3(mol, lattr, cutoff%disp3, cutoff%width3, &
      & disp%r4r2, c6, energy3)

end subroutine get_pairwise_dispersion


end module dftd4_disp
