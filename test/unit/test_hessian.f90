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

module test_hessian
   use dftd4, only : d4_model, d4s_model, dispersion_model, d4_qmod, &
      & new_d4_model, new_d4s_model, rational_damping_param, realspace_cutoff, &
      & work_partition, new_work_partition
   use dftd4_disp, only : get_dispersion2_hessian
   use dftd4_numdiff, only : get_dispersion_hessian, get_dispersion_hessian_numdiff
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   implicit none
   private
   public :: collect_hessian

contains

subroutine collect_hessian(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
      & new_unittest("C6 second derivatives", test_c6), &
      & new_unittest("molecular Hessian", test_molecular), &
      & new_unittest("periodic Hessian", test_periodic) &
      & ]
end subroutine collect_hessian

subroutine test_c6(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(d4_model) :: d4
   type(d4s_model) :: d4s

   call new(mol, [8, 1, 6], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
      & 1.8_wp, 0.2_wp, 0.1_wp, -0.4_wp, 2.5_wp, 0.3_wp], [3, 3]))
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call check_c6(error, mol, d4)
   if (allocated(error)) return
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call check_c6(error, mol, d4s)
end subroutine test_c6

! Differentiate the existing first-derivative API, including diagonal pairs.
subroutine check_c6(error, mol, model)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(in) :: mol
   class(dispersion_model), intent(in) :: model
   real(wp), parameter :: step = 1.0e-5_wp, tol = 2.0e-7_wp
   real(wp) :: cn(3), q(3), c6(3, 3), dc(3, 3), dcc(3, 3), dccij(3, 3)
   real(wp) :: dq(3, 3), dqq(3, 3), dqqij(3, 3), dcq(3, 3), dcqij(3, 3)
   real(wp) :: cp(3, 3), qp(3, 3), cm(3, 3), qm(3, 3), ec, eq, diff, scale
   real(wp), allocatable :: gw(:, :, :), gc(:, :, :), gcc(:, :, :)
   real(wp), allocatable :: gq(:, :, :), gqq(:, :, :), gcq(:, :, :)
   integer :: mref, k, variable, branch, i, j

   mref = maxval(model%ref)
   allocate(gw(mref, 3, model%ncoup), gc(mref, 3, model%ncoup), gcc(mref, 3, model%ncoup), &
      & gq(mref, 3, model%ncoup), gqq(mref, 3, model%ncoup), gcq(mref, 3, model%ncoup))
   do branch = 1, 2
      cn = [1.2_wp, 0.4_wp, 2.3_wp]
      q = [-0.4_wp, 0.1_wp, 0.3_wp]
      if (branch == 2) then
         ! Exercise the CN normalization fallback and constant charge-scaling branch.
         cn(2) = 1.0e4_wp
         q(1) = -model%zeff(mol%id(1)) - 1.0_wp
      end if
      select type(model)
      type is (d4_model)
         call model%weight_references_hessian(mol, cn, q, gw, gc, gcc, gq, gqq, gcq)
      type is (d4s_model)
         call model%weight_references_hessian(mol, cn, q, gw, gc, gcc, gq, gqq, gcq)
      class default
         return
      end select
      call model%get_atomic_c6_hessian(mol, gw, gc, gcc, c6, dc, dcc, dccij, &
         & gq, gqq, gcq, dq, dqq, dqqij, dcq, dcqij)
      do variable = 1, 2
         do k = 1, 3
            if (variable == 1) cn(k) = cn(k) + step
            if (variable == 2) q(k) = q(k) + step
            call model%weight_references(mol, cn, q, gw, gc, gq)
            call model%get_atomic_c6(mol, gw, gc, gq, c6, cp, qp)
            if (variable == 1) cn(k) = cn(k) - 2*step
            if (variable == 2) q(k) = q(k) - 2*step
            call model%weight_references(mol, cn, q, gw, gc, gq)
            call model%get_atomic_c6(mol, gw, gc, gq, c6, cm, qm)
            if (variable == 1) cn(k) = cn(k) + step
            if (variable == 2) q(k) = q(k) + step
            diff = 0.0_wp
            scale = 1.0_wp
            do j = 1, 3
               do i = 1, 3
                  ec = 0.0_wp
                  eq = 0.0_wp
                  if (variable == 1) then
                     if (i == k) ec = ec + dcc(i, j)
                     if (j == k) ec = ec + dccij(i, j)
                     if (i == k) eq = eq + dcq(i, j)
                     if (j == k) eq = eq + dcqij(j, i)
                  else
                     if (i == k) ec = ec + dcq(i, j)
                     if (j == k) ec = ec + dcqij(i, j)
                     if (i == k) eq = eq + dqq(i, j)
                     if (j == k) eq = eq + dqqij(i, j)
                  end if
                  diff = max(diff, abs(ec-(cp(i, j)-cm(i, j))/(2*step)), &
                     & abs(eq-(qp(i, j)-qm(i, j))/(2*step)))
                  scale = max(scale, abs(ec), abs(eq))
               end do
            end do
            call check(error, diff/scale, 0.0_wp, thr=tol)
            if (allocated(error)) return
         end do
      end do
   end do
end subroutine check_c6

subroutine test_molecular(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "01")
   call check_models(error, mol)
end subroutine test_molecular

subroutine test_periodic(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   real(wp) :: lattice(3, 3)

   lattice = reshape([8.0_wp, 0.0_wp, 0.0_wp, 1.3_wp, 9.0_wp, 0.0_wp, &
      & -0.6_wp, 0.9_wp, 8.5_wp], [3, 3])
   call new(mol, [30, 8, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
      & 3.0_wp, 0.1_wp, 0.2_wp, 0.5_wp, 2.0_wp, 0.3_wp], [3, 3]), &
      & charge=1.0_wp, lattice=lattice, periodic=[.true., .true., .true.])
   call check_models(error, mol)
   if (allocated(error)) return
   call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), &
      & charge=-1.0_wp, lattice=lattice, periodic=[.true., .true., .true.])
   call check_models(error, mol)
end subroutine test_periodic

subroutine check_models(error, mol)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(in) :: mol
   type(d4_model) :: d4
   type(d4s_model) :: d4s
   integer :: kind, qmod

   do kind = 1, 2
      qmod = d4_qmod%eeq
      if (kind == 2) qmod = d4_qmod%eeqbc
      call new_d4_model(error, d4, mol, qmod=qmod)
      if (allocated(error)) return
      call check_hessian(error, mol, d4)
      if (allocated(error)) return
      call new_d4s_model(error, d4s, mol, qmod=qmod)
      if (allocated(error)) return
      call check_hessian(error, mol, d4s)
      if (allocated(error)) return
   end do
end subroutine check_models

subroutine check_hessian(error, mol, model)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(in) :: mol
   class(dispersion_model), intent(in) :: model
   type(rational_damping_param) :: param
   type(realspace_cutoff) :: cutoff
   type(work_partition) :: partition
   real(wp), allocatable :: hess(:, :), part_hess(:, :), sum_hess(:, :)
   real(wp), allocatable :: analytic(:, :, :, :), numeric(:, :, :, :)
   integer :: n, part, ic

   n = mol%nat
   param = rational_damping_param(s6=1.0_wp, s8=0.95948085_wp, s9=0.0_wp, &
      & a1=0.38574991_wp, a2=4.80688534_wp)
   cutoff = realspace_cutoff(cn=20.0_wp, disp2=18.0_wp, width2=4.0_wp, disp3=12.0_wp, width3=2.0_wp)
   allocate(hess(3*n, 3*n), part_hess(3*n, 3*n), sum_hess(3*n, 3*n))
   allocate(analytic(3, n, 3, n), numeric(3, n, 3, n))
   call get_dispersion2_hessian(error, mol, model, param, cutoff, hess)
   if (allocated(error)) return
   call get_dispersion_hessian_numdiff(mol, model, param, cutoff, numeric, error=error)
   if (allocated(error)) return
   call check(error, maxval(abs(hess-reshape(numeric, [3*n, 3*n]))), 0.0_wp, thr=2.0e-8_wp)
   if (allocated(error)) return
   call check(error, maxval(abs(hess-transpose(hess))), 0.0_wp, thr=1.0e-12_wp)
   if (allocated(error)) return
   do ic = 1, 3
      call check(error, maxval(abs(sum(hess(:, ic::3), dim=2))), 0.0_wp, thr=1.0e-11_wp)
      if (allocated(error)) return
   end do

   ! The small periodic cells also exercise empty partitions and self images.
   if (n <= 3) then
      sum_hess = 0.0_wp
      do part = 0, 7
         call new_work_partition(error, partition, part, 8)
         if (allocated(error)) return
         call get_dispersion2_hessian(error, mol, model, param, cutoff, part_hess, partition)
         if (allocated(error)) return
         sum_hess = sum_hess + part_hess
      end do
      call check(error, maxval(abs(sum_hess-hess)), 0.0_wp, thr=1.0e-12_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(part_hess)), 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
   end if

   ! Verify that the public dispatcher includes the existing analytic ATM term.
   param%s9 = 1.0_wp
   call get_dispersion_hessian(mol, model, param, cutoff, analytic, error=error)
   if (allocated(error)) return
   call get_dispersion_hessian_numdiff(mol, model, param, cutoff, numeric, error=error)
   if (allocated(error)) return
   call check(error, maxval(abs(analytic-numeric)), 0.0_wp, thr=2.0e-8_wp)
end subroutine check_hessian

end module test_hessian
