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

!> Implementation of the Axilrod-Teller-Muto triple dipole dispersion
!> contribution with a modified zero (Chai--Head-Gordon) damping together
!> with the critical radii from the rational (Becke--Johnson) damping.
module dftd4_damping_atm
   use dftd4_cutoff, only : smooth_cutoff, smooth_cutoff_r2
   use dftd4_partition, only : work_partition, owns_pair
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none
   private

   public :: get_atm_dispersion, get_atm_dispersion_hessian

   real(wp), parameter :: third = 1.0_wp / 3.0_wp

contains


!> Analytical Hessian of the ATM three-body dispersion expression
subroutine get_atm_dispersion_hessian(mol, trans, cutoff, width, s9, a1, a2, alp, r4r2, &
      & c6, dc6dcn, d2c6dcn2, d2c6dcnij, hessian, dEdcn, dEdcndr, dEdcndcn, partition)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Width of smooth cutoff
   real(wp), intent(in) :: width

   !> Scaling for dispersion coefficients
   real(wp), intent(in) :: s9

   !> Scaling parameter for critical radius
   real(wp), intent(in) :: a1

   !> Offset parameter for critical radius
   real(wp), intent(in) :: a2

   !> Exponent of zero damping function
   real(wp), intent(in) :: alp

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivatives of the C6 w.r.t. the coordination number
   real(wp), intent(in) :: dc6dcn(:, :), d2c6dcn2(:, :), d2c6dcnij(:, :)

   !> Second derivative of the energy w.r.t. the Cartesian coordinates
   real(wp), intent(inout) :: hessian(:, :)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout) :: dEdcn(:)

   !> Mixed derivative w.r.t. coordination number and Cartesian coordinates
   real(wp), intent(inout) :: dEdcndr(:, :)

   !> Second derivative w.r.t. the coordination numbers
   real(wp), intent(inout) :: dEdcndcn(:, :)

   !> Work partition of the atom pairs, absent selects the complete work
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, kat, izp, jzp, kzp, jtr, ktr
   integer :: ipair, npair, ntrans, ip, iq, il, im, ia, ib, ic, jc, ie, if_, nent
   integer :: at(3), pat(2, 3), ent_atom(6), ent_pair(6)
   real(wp) :: vec(3, 3), u(3), cutoff2, triple, r0ij, r0ik, r0jk, r0, alp3, aexp, cval
   real(wp) :: swp(3), dswp(3), d2swp(3), sval, sp(3), spq(3, 3)
   real(wp) :: ww, wp_(3), wpq(3, 3), lval(3), nval, np(3), npq(3, 3)
   real(wp) :: gv, gd, gdd, hv, hd, hdd, ang, angp(3), angpq(3, 3)
   real(wp) :: tval, tp(3), tpq(3, 3), fd, fdp(3), fdpq(3, 3)
   real(wp) :: av, ap(3), apq(3, 3), kv, kp(3), kpq(3, 3)
   real(wp) :: c6t(3), qv, qm(3), qmn(3, 3), pref, ec(3), ecc(3, 3), euc(3, 3)
   real(wp) :: ent_coef(6), grad(3, 3, 3), tmp
   real(wp) :: gk(3, 3, 3), egr(3, 3, 3), dg(3, 3), pq

   ! Thread-private arrays for reduction
   real(wp), allocatable :: hessian_local(:, :), dEdcn_local(:)
   real(wp), allocatable :: dEdcndr_local(:, :), dEdcndcn_local(:, :)
   integer, parameter :: sigp(3, 3) = reshape(&
      & [-1, 1, 0, -1, 0, 1, 0, -1, 1], [3, 3])
   real(wp), parameter :: cl(3, 3) = reshape(&
      & [1.0_wp, 1.0_wp, -1.0_wp, -1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, -1.0_wp, 1.0_wp], &
      & [3, 3])

   if (abs(s9) < epsilon(1.0_wp)) return
   cutoff2 = cutoff*cutoff
   alp3 = alp / 3.0_wp
   aexp = 0.5_wp * alp3
   npair = mol%nat*(mol%nat + 1)/2
   ntrans = size(trans, 2)

   !$omp parallel default(none) &
   !$omp shared(mol, trans, cutoff, width, s9, a1, a2, r4r2, c6, dc6dcn, &
   !$omp& d2c6dcn2, d2c6dcnij, hessian, dEdcn, dEdcndr, dEdcndcn, partition, &
   !$omp& cutoff2, alp3, aexp, npair, ntrans) &
   !$omp private(ipair, iat, jat, kat, izp, jzp, kzp, jtr, ktr, ip, iq, il, im, ia, ib, &
   !$omp& ic, jc, ie, if_, nent, at, pat, ent_atom, ent_pair, vec, u, triple, &
   !$omp& r0ij, r0ik, r0jk, r0, cval, swp, dswp, d2swp, sval, sp, spq, ww, wp_, &
   !$omp& wpq, lval, nval, np, npq, gv, gd, gdd, hv, hd, hdd, ang, angp, angpq, &
   !$omp& tval, tp, tpq, fd, fdp, fdpq, av, ap, apq, kv, kp, kpq, c6t, qv, qm, &
   !$omp& qmn, pref, ec, ecc, euc, ent_coef, grad, tmp, gk, egr, dg, pq, &
   !$omp& hessian_local, dEdcn_local, dEdcndr_local, dEdcndcn_local)
   allocate(hessian_local(size(hessian, 1), size(hessian, 2)), source=0.0_wp)
   allocate(dEdcn_local(size(dEdcn, 1)), source=0.0_wp)
   allocate(dEdcndr_local(size(dEdcndr, 1), size(dEdcndr, 2)), source=0.0_wp)
   allocate(dEdcndcn_local(size(dEdcndcn, 1), size(dEdcndcn, 2)), source=0.0_wp)
   ! Schedule pair/translation combinations rather than only atom pairs.  This
   ! retains O(N^2) tasks for molecules and exposes the lattice-image work to
   ! OpenMP for small periodic unit cells.
   !$omp do collapse(2) schedule(guided, 1)
   do ipair = 1, npair
      do jtr = 1, ntrans
            iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp) + 1.0_wp) - 1.0_wp))
            if (iat*(iat + 1)/2 < ipair) iat = iat + 1
            jat = ipair - iat*(iat - 1)/2
            if (.not.owns_pair(partition, iat, jat)) cycle
            izp = mol%id(iat)
            jzp = mol%id(jat)
            vec(:, 1) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            u(1) = sum(vec(:, 1)**2)
            if (u(1) > cutoff2 .or. u(1) < epsilon(1.0_wp)) cycle
            c6t(1) = c6(jat, iat)
            if (abs(c6t(1)) < epsilon(1.0_wp)) cycle
            r0ij = a1 * sqrt(3.0_wp*r4r2(jzp)*r4r2(izp)) + a2
            call smooth_cutoff_r2(u(1), cutoff, width, swp(1), dswp(1), d2swp(1))

            do kat = 1, jat
               kzp = mol%id(kat)
               triple = triple_scale(iat, jat, kat)
               c6t(2) = c6(kat, iat)
               c6t(3) = c6(kat, jat)
               if (any(abs(c6t(2:3)) < epsilon(1.0_wp))) cycle

               r0ik = a1 * sqrt(3.0_wp*r4r2(kzp)*r4r2(izp)) + a2
               r0jk = a1 * sqrt(3.0_wp*r4r2(kzp)*r4r2(jzp)) + a2
               r0 = r0ij*r0ik*r0jk
               cval = 6.0_wp * r0**alp3

               ! These C6/CN quantities do not depend on the lattice image of
               ! the third atom.  Prepare them once per atom triple.
               qv = sqrt(abs(c6t(1)*c6t(2)*c6t(3)))
               do ip = 1, 3
                  qm(ip) = 0.5_wp*qv/c6t(ip)
               end do
               do ip = 1, 3
                  do iq = 1, 3
                     if (ip == iq) then
                        qmn(ip, iq) = -0.25_wp*qv/(c6t(ip)*c6t(ip))
                     else
                        qmn(ip, iq) = 0.25_wp*qv/(c6t(ip)*c6t(iq))
                     end if
                  end do
               end do
               pref = s9*triple
               pq = pref*qv

               at(1) = iat; at(2) = jat; at(3) = kat
               pat(1, 1) = iat; pat(2, 1) = jat
               pat(1, 2) = iat; pat(2, 2) = kat
               pat(1, 3) = jat; pat(2, 3) = kat

               nent = 0
               do ip = 1, 3
                  if (pat(1, ip) /= pat(2, ip)) then
                     nent = nent + 1
                     ent_atom(nent) = pat(1, ip)
                     ent_pair(nent) = ip
                     ent_coef(nent) = dc6dcn(pat(1, ip), pat(2, ip))
                     nent = nent + 1
                     ent_atom(nent) = pat(2, ip)
                     ent_pair(nent) = ip
                     ent_coef(nent) = dc6dcn(pat(2, ip), pat(1, ip))
                  else
                     nent = nent + 1
                     ent_atom(nent) = pat(1, ip)
                     ent_pair(nent) = ip
                     ent_coef(nent) = 2.0_wp*dc6dcn(pat(1, ip), pat(1, ip))
                  end if
               end do

               do ktr = 1, ntrans
                  vec(:, 2) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, iat)
                  u(2) = sum(vec(:, 2)**2)
                  if (u(2) > cutoff2 .or. u(2) < epsilon(1.0_wp)) cycle
                  vec(:, 3) = vec(:, 2) - vec(:, 1)
                  u(3) = sum(vec(:, 3)**2)
                  if (u(3) > cutoff2 .or. u(3) < epsilon(1.0_wp)) cycle

                  ! switching function and its derivatives w.r.t. the squared distances
                  do ip = 2, 3
                     call smooth_cutoff_r2(u(ip), cutoff, width, swp(ip), dswp(ip), d2swp(ip))
                  end do
                  sval = swp(1)*swp(2)*swp(3)
                  sp(1) = dswp(1)*swp(2)*swp(3)
                  sp(2) = swp(1)*dswp(2)*swp(3)
                  sp(3) = swp(1)*swp(2)*dswp(3)
                  spq(1, 1) = d2swp(1)*swp(2)*swp(3)
                  spq(2, 2) = swp(1)*d2swp(2)*swp(3)
                  spq(3, 3) = swp(1)*swp(2)*d2swp(3)
                  spq(1, 2) = dswp(1)*dswp(2)*swp(3)
                  spq(2, 1) = spq(1, 2)
                  spq(1, 3) = dswp(1)*swp(2)*dswp(3)
                  spq(3, 1) = spq(1, 3)
                  spq(2, 3) = swp(1)*dswp(2)*dswp(3)
                  spq(3, 2) = spq(2, 3)

                  ! product of the squared distances
                  ww = u(1)*u(2)*u(3)
                  wp_(1) = u(2)*u(3)
                  wp_(2) = u(1)*u(3)
                  wp_(3) = u(1)*u(2)
                  wpq(:, :) = 0.0_wp
                  wpq(1, 2) = u(3); wpq(2, 1) = u(3)
                  wpq(1, 3) = u(2); wpq(3, 1) = u(2)
                  wpq(2, 3) = u(1); wpq(3, 2) = u(1)

                  ! triple product entering the angular term
                  lval(1) = u(1) + u(3) - u(2)
                  lval(2) = u(1) - u(3) + u(2)
                  lval(3) = -u(1) + u(3) + u(2)
                  nval = lval(1)*lval(2)*lval(3)
                  do ip = 1, 3
                     np(ip) = cl(1, ip)*lval(2)*lval(3) + cl(2, ip)*lval(1)*lval(3) &
                        & + cl(3, ip)*lval(1)*lval(2)
                  end do
                  do ip = 1, 3
                     do iq = 1, 3
                        tmp = 0.0_wp
                        do il = 1, 3
                           do im = 1, 3
                              if (il == im) cycle
                              ! remaining index of the product
                              tmp = tmp + cl(il, ip)*cl(im, iq)*lval(6 - il - im)
                           end do
                        end do
                        npq(ip, iq) = tmp
                     end do
                  end do

                  gv = ww**(-2.5_wp)
                  gd = -2.5_wp * ww**(-3.5_wp)
                  gdd = 8.75_wp * ww**(-4.5_wp)
                  hv = ww**(-1.5_wp)
                  hd = -1.5_wp * ww**(-2.5_wp)
                  hdd = 3.75_wp * ww**(-3.5_wp)

                  ang = 0.375_wp*nval*gv + hv
                  do ip = 1, 3
                     angp(ip) = 0.375_wp*(np(ip)*gv + nval*gd*wp_(ip)) + hd*wp_(ip)
                  end do
                  do ip = 1, 3
                     do iq = 1, 3
                        angpq(ip, iq) = 0.375_wp*(npq(ip, iq)*gv &
                           & + np(ip)*gd*wp_(iq) + np(iq)*gd*wp_(ip) &
                           & + nval*(gdd*wp_(ip)*wp_(iq) + gd*wpq(ip, iq))) &
                           & + hdd*wp_(ip)*wp_(iq) + hd*wpq(ip, iq)
                     end do
                  end do

                  ! zero damping function
                  tval = cval * ww**(-aexp)
                  do ip = 1, 3
                     tp(ip) = -aexp*tval*wp_(ip)/ww
                  end do
                  do ip = 1, 3
                     do iq = 1, 3
                        tpq(ip, iq) = aexp*(aexp + 1.0_wp)*tval*wp_(ip)*wp_(iq)/(ww*ww) &
                           & - aexp*tval*wpq(ip, iq)/ww
                     end do
                  end do
                  fd = 1.0_wp/(1.0_wp + tval)
                  do ip = 1, 3
                     fdp(ip) = -tp(ip)*fd*fd
                  end do
                  do ip = 1, 3
                     do iq = 1, 3
                        fdpq(ip, iq) = -tpq(ip, iq)*fd*fd + 2.0_wp*tp(ip)*tp(iq)*fd**3
                     end do
                  end do

                  av = ang*fd
                  do ip = 1, 3
                     ap(ip) = angp(ip)*fd + ang*fdp(ip)
                  end do
                  do ip = 1, 3
                     do iq = 1, 3
                        apq(ip, iq) = angpq(ip, iq)*fd + angp(ip)*fdp(iq) &
                           & + angp(iq)*fdp(ip) + ang*fdpq(ip, iq)
                     end do
                  end do

                  kv = sval*av
                  do ip = 1, 3
                     kp(ip) = sp(ip)*av + sval*ap(ip)
                  end do
                  do ip = 1, 3
                     do iq = 1, 3
                        kpq(ip, iq) = spq(ip, iq)*av + sp(ip)*ap(iq) &
                           & + sp(iq)*ap(ip) + sval*apq(ip, iq)
                     end do
                  end do

                  do ia = 1, 3
                     do ip = 1, 3
                        grad(:, ia, ip) = 2.0_wp*sigp(ia, ip)*vec(:, ip)
                     end do
                  end do

                  ! Cartesian second derivatives at fixed coordination number
                  ! contract with the pair gradients once instead of per component pair
                  do ia = 1, 3
                     do ic = 1, 3
                        do iq = 1, 3
                           tmp = 0.0_wp
                           do ip = 1, 3
                              tmp = tmp + kpq(ip, iq)*grad(ic, ia, ip)
                           end do
                           gk(ic, ia, iq) = pq*tmp
                        end do
                     end do
                  end do
                  do ia = 1, 3
                     do ib = 1, 3
                        tmp = 0.0_wp
                        do ip = 1, 3
                           tmp = tmp + kp(ip)*sigp(ia, ip)*sigp(ib, ip)
                        end do
                        dg(ia, ib) = 2.0_wp*pq*tmp
                     end do
                  end do

                  do ia = 1, 3
                     do ib = 1, 3
                        do ic = 1, 3
                           do jc = 1, 3
                              tmp = 0.0_wp
                              do iq = 1, 3
                                 tmp = tmp + gk(ic, ia, iq)*grad(jc, ib, iq)
                              end do
                              if (ic == jc) tmp = tmp + dg(ia, ib)
                              hessian_local(3*(at(ia)-1)+ic, 3*(at(ib)-1)+jc) = &
                                 & hessian_local(3*(at(ia)-1)+ic, 3*(at(ib)-1)+jc) + tmp
                           end do
                        end do
                     end do
                  end do

                  ! derivatives with respect to the C6 coefficients
                  do ip = 1, 3
                     ec(ip) = pref*qm(ip)*kv
                     do iq = 1, 3
                        ecc(ip, iq) = pref*qmn(ip, iq)*kv
                        euc(iq, ip) = pref*qm(ip)*kp(iq)
                     end do
                  end do

                  ! contract the mixed CN/Cartesian derivatives once per pair
                  do ia = 1, 3
                     do ic = 1, 3
                        do im = 1, 3
                           tmp = 0.0_wp
                           do ip = 1, 3
                              tmp = tmp + euc(ip, im)*grad(ic, ia, ip)
                           end do
                           egr(ic, ia, im) = tmp
                        end do
                     end do
                  end do

                  do ie = 1, nent
                     dEdcn_local(ent_atom(ie)) = dEdcn_local(ent_atom(ie)) &
                        & + ec(ent_pair(ie))*ent_coef(ie)
                     do if_ = 1, nent
                        dEdcndcn_local(ent_atom(ie), ent_atom(if_)) = &
                           & dEdcndcn_local(ent_atom(ie), ent_atom(if_)) &
                           & + ecc(ent_pair(ie), ent_pair(if_))*ent_coef(ie)*ent_coef(if_)
                     end do
                     do ia = 1, 3
                        do ic = 1, 3
                           dEdcndr_local(3*(at(ia)-1)+ic, ent_atom(ie)) = &
                              & dEdcndr_local(3*(at(ia)-1)+ic, ent_atom(ie)) &
                              & + egr(ic, ia, ent_pair(ie))*ent_coef(ie)
                        end do
                     end do
                  end do

                  ! second derivative of the C6 coefficients w.r.t. the coordination numbers
                  do ip = 1, 3
                     ia = pat(1, ip)
                     ib = pat(2, ip)
                     if (ia /= ib) then
                        dEdcndcn_local(ia, ia) = dEdcndcn_local(ia, ia) + ec(ip)*d2c6dcn2(ia, ib)
                        dEdcndcn_local(ib, ib) = dEdcndcn_local(ib, ib) + ec(ip)*d2c6dcn2(ib, ia)
                        dEdcndcn_local(ia, ib) = dEdcndcn_local(ia, ib) + ec(ip)*d2c6dcnij(ia, ib)
                        dEdcndcn_local(ib, ia) = dEdcndcn_local(ib, ia) + ec(ip)*d2c6dcnij(ia, ib)
                     else
                        dEdcndcn_local(ia, ia) = dEdcndcn_local(ia, ia) + ec(ip) &
                           & * (2.0_wp*d2c6dcn2(ia, ia) + 2.0_wp*d2c6dcnij(ia, ia))
                     end if
                  end do
               end do
            end do
      end do
   end do
   !$omp end do nowait
   ! Merge a completed thread while other threads can still process their tail
   ! pairs.  The nowait above overlaps this reduction with useful ATM work.
   !$omp critical (get_atm_dispersion_hessian_)
   hessian(:, :) = hessian(:, :) + hessian_local(:, :)
   dEdcn(:) = dEdcn(:) + dEdcn_local(:)
   dEdcndr(:, :) = dEdcndr(:, :) + dEdcndr_local(:, :)
   dEdcndcn(:, :) = dEdcndcn(:, :) + dEdcndcn_local(:, :)
   !$omp end critical (get_atm_dispersion_hessian_)
   deallocate(hessian_local, dEdcn_local, dEdcndr_local, dEdcndcn_local)
   !$omp end parallel

end subroutine get_atm_dispersion_hessian


!> Evaluation of the dispersion energy expression
subroutine get_atm_dispersion(mol, trans, cutoff, width, s9, a1, a2, alp, r4r2, &
      & c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient, sigma, partition)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Width of smooth cutoff
   real(wp), intent(in) :: width

   !> Scaling for dispersion coefficients
   real(wp), intent(in) :: s9

   !> Scaling parameter for critical radius
   real(wp), intent(in) :: a1

   !> Offset parameter for critical radius
   real(wp), intent(in) :: a2

   !> Exponent of zero damping function
   real(wp), intent(in) :: alp

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in), optional :: dc6dcn(:, :)

   !> Derivative of the C6 w.r.t. the partial charges
   real(wp), intent(in), optional :: dc6dq(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)

   !> Derivative of the energy w.r.t. the partial charges
   real(wp), intent(inout), optional :: dEdq(:)

   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   !> Optional externally assigned work partition
   type(work_partition), intent(in), optional :: partition

   logical :: grad

   if (abs(s9) < epsilon(1.0_wp)) return
   grad = present(dc6dcn) .and. present(dEdcn) .and. present(dc6dq) &
      & .and. present(dEdq) .and. present(gradient) .and. present(sigma)

   if (grad) then
      call get_atm_dispersion_derivs(mol, trans, cutoff, width, s9, a1, a2, &
         & alp, r4r2, c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient, sigma, partition)
   else
      call get_atm_dispersion_energy(mol, trans, cutoff, width, s9, a1, a2, &
         & alp, r4r2, c6, energy, partition)
   end if

end subroutine get_atm_dispersion


!> Evaluation of the dispersion energy expression
subroutine get_atm_dispersion_energy(mol, trans, cutoff, width, s9, a1, a2, alp, r4r2, &
      & c6, energy, partition)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Width of smooth cutoff
   real(wp), intent(in) :: width

   !> Scaling for dispersion coefficients
   real(wp), intent(in) :: s9

   !> Scaling parameter for critical radius
   real(wp), intent(in) :: a1

   !> Offset parameter for critical radius
   real(wp), intent(in) :: a2

   !> Exponent of zero damping function
   real(wp), intent(in) :: alp

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Work partition of the atom pairs, absent selects the complete work
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, kat, izp, jzp, kzp, jtr, ktr
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, rij, rjk, rik
   real(wp) :: c6ij, c6jk, c6ik, triple
   real(wp) :: r0ij, r0jk, r0ik, r0, r1, r2, r3, r5, rr, fdmp, ang
   real(wp) :: cutoff2, c9, dE, alp3, swij, swjk, swik, dswdr, sw

   ! Thread-private arrays for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:)

   cutoff2 = cutoff*cutoff
   alp3 = alp / 3.0_wp

   !$omp parallel default(none) &
   !$omp shared(mol, trans, c6, s9, a1, a2, alp3, r4r2, cutoff2, cutoff, width, partition) &
   !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
   !$omp& r2ij, r2jk, r2ik, rij, rjk, rik, c6ij, c6jk, c6ik, triple, &
   !$omp& r0ij, r0jk, r0ik, r0, r1, r2, r3, r5, rr, fdmp, ang, c9, dE, &
   !$omp& swij, swjk, swik, dswdr, sw) &
   !$omp shared(energy) &
   !$omp private(energy_local)
   allocate(energy_local(size(energy, 1)), source=0.0_wp)
   !$omp do schedule(dynamic)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         if (.not.owns_pair(partition, iat, jat)) cycle
         jzp = mol%id(jat)
         c6ij = c6(jat, iat)
         r0ij = a1 * sqrt(3*r4r2(jzp)*r4r2(izp)) + a2
         do jtr = 1, size(trans, 2)
            vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
            if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
            rij = sqrt(r2ij)
            call smooth_cutoff(rij, cutoff, width, swij, dswdr)
            do kat = 1, jat
               kzp = mol%id(kat)
               c6ik = c6(kat, iat)
               c6jk = c6(kat, jat)
               c9 = -s9 * sqrt(abs(c6ij*c6ik*c6jk))
               r0ik = a1 * sqrt(3*r4r2(kzp)*r4r2(izp)) + a2
               r0jk = a1 * sqrt(3*r4r2(kzp)*r4r2(jzp)) + a2
               r0 = r0ij * r0ik * r0jk
               triple = triple_scale(iat, jat, kat)
               do ktr = 1, size(trans, 2)
                  vik(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, iat)
                  r2ik = vik(1)*vik(1) + vik(2)*vik(2) + vik(3)*vik(3)
                  if (r2ik > cutoff2 .or. r2ik < epsilon(1.0_wp)) cycle
                  rik = sqrt(r2ik)
                  call smooth_cutoff(rik, cutoff, width, swik, dswdr)

                  ! vjk(:) = mol%xyz(:, kat) + trans(:, ktr)
                  !          - mol%xyz(:, jat) - trans(:, jtr)
                  vjk(:) = vik(:) - vij(:)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                  rjk = sqrt(r2jk)
                  call smooth_cutoff(rjk, cutoff, width, swjk, dswdr)
                  sw = swij * swik * swjk
                  if (sw <= 0.0_wp) cycle

                  r2 = r2ij*r2ik*r2jk
                  r1 = sqrt(r2)
                  r3 = r2 * r1
                  r5 = r3 * r2

                  fdmp = 1.0_wp / (1.0_wp + 6.0_wp * (r0 / r1)**alp3)
                  ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik) &
                     & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3

                  rr = ang*fdmp

                  dE = rr * c9 * triple * third * sw
                  energy_local(iat) = energy_local(iat) - dE
                  energy_local(jat) = energy_local(jat) - dE
                  energy_local(kat) = energy_local(kat) - dE
               end do
            end do
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_atm_dispersion_energy_)
   energy(:) = energy(:) + energy_local(:)
   !$omp end critical (get_atm_dispersion_energy_)
   deallocate(energy_local)
   !$omp end parallel

end subroutine get_atm_dispersion_energy


!> Evaluation of the dispersion energy expression
subroutine get_atm_dispersion_derivs(mol, trans, cutoff, width, s9, a1, a2, alp, r4r2, &
      & c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient, sigma, partition)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Width of smooth cutoff
   real(wp), intent(in) :: width

   !> Scaling for dispersion coefficients
   real(wp), intent(in) :: s9

   !> Scaling parameter for critical radius
   real(wp), intent(in) :: a1

   !> Offset parameter for critical radius
   real(wp), intent(in) :: a2

   !> Exponent of zero damping function
   real(wp), intent(in) :: alp

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in) :: dc6dcn(:, :)

   !> Derivative of the C6 w.r.t. the partial charges
   real(wp), intent(in) :: dc6dq(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout) :: dEdcn(:)

   !> Derivative of the energy w.r.t. the partial charges
   real(wp), intent(inout) :: dEdq(:)

   !> Dispersion gradient
   real(wp), intent(inout) :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout) :: sigma(:, :)

   !> Work partition of the atom pairs, absent selects the complete work
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, kat, izp, jzp, kzp, jtr, ktr
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, rij, rjk, rik
   real(wp) :: c6ij, c6jk, c6ik, triple
   real(wp) :: r0ij, r0jk, r0ik, r0, r1, r2, r3, r5, rr, fdmp, dfdmp, ang, dang
   real(wp) :: cutoff2, alp3, c9, dE, dE0, dE_third
   real(wp) :: dGij(3), dGjk(3), dGik(3), dS(3, 3)
   real(wp) :: swij, swjk, swik, dswijdr, dswjkdr, dswikdr, sw

   ! Thread-private arrays for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:)
   real(wp), allocatable :: dEdcn_local(:)
   real(wp), allocatable :: dEdq_local(:)
   real(wp), allocatable :: gradient_local(:, :)
   real(wp), allocatable :: sigma_local(:, :)

   cutoff2 = cutoff*cutoff
   alp3 = alp / 3.0_wp

   !$omp parallel default(none) &
   !$omp shared(mol, trans, c6, s9, a1, a2, alp, alp3, r4r2, cutoff2, &
   !$omp& cutoff, width, dc6dcn, dc6dq, partition) &
   !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
   !$omp& r2ij, r2jk, r2ik, rij, rjk, rik, c6ij, c6jk, c6ik, triple, &
   !$omp& r0ij, r0jk, r0ik, r0, r1, r2, r3, r5, rr, fdmp, dfdmp, ang, dang, &
   !$omp& c9, dE, dE0, dE_third, dGij, dGjk, dGik, dS, swij, swjk, swik, &
   !$omp& dswijdr, dswjkdr, dswikdr, sw) &
   !$omp shared(energy, gradient, sigma, dEdcn, dEdq) &
   !$omp private(energy_local, gradient_local, sigma_local, dEdcn_local, &
   !$omp& dEdq_local)
   allocate(energy_local(size(energy, 1)), source=0.0_wp)
   allocate(dEdcn_local(size(dEdcn, 1)), source=0.0_wp)
   allocate(dEdq_local(size(dEdq, 1)), source=0.0_wp)
   allocate(gradient_local(size(gradient, 1), size(gradient, 2)), source=0.0_wp)
   allocate(sigma_local(size(sigma, 1), size(sigma, 2)), source=0.0_wp)
   !$omp do schedule(dynamic)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         if (.not.owns_pair(partition, iat, jat)) cycle
         jzp = mol%id(jat)
         c6ij = c6(jat, iat)
         r0ij = a1 * sqrt(3*r4r2(jzp)*r4r2(izp)) + a2
         do jtr = 1, size(trans, 2)
            vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
            if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
            rij = sqrt(r2ij)
            call smooth_cutoff(rij, cutoff, width, swij, dswijdr)
            do kat = 1, jat
               kzp = mol%id(kat)
               c6ik = c6(kat, iat)
               c6jk = c6(kat, jat)
               c9 = -s9 * sqrt(abs(c6ij*c6ik*c6jk))
               r0ik = a1 * sqrt(3*r4r2(kzp)*r4r2(izp)) + a2
               r0jk = a1 * sqrt(3*r4r2(kzp)*r4r2(jzp)) + a2
               r0 = r0ij * r0ik * r0jk
               triple = triple_scale(iat, jat, kat)
               do ktr = 1, size(trans, 2)
                  vik(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, iat)
                  r2ik = vik(1)*vik(1) + vik(2)*vik(2) + vik(3)*vik(3)
                  if (r2ik > cutoff2 .or. r2ik < epsilon(1.0_wp)) cycle
                  rik = sqrt(r2ik)
                  call smooth_cutoff(rik, cutoff, width, swik, dswikdr)

                  ! vjk(:) = mol%xyz(:, kat) + trans(:, ktr)
                  !          - mol%xyz(:, jat) - trans(:, jtr)
                  vjk(:) = vik(:) - vij(:)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                  rjk = sqrt(r2jk)
                  call smooth_cutoff(rjk, cutoff, width, swjk, dswjkdr)
                  sw = swij * swik * swjk
                  if (sw <= 0.0_wp) cycle

                  r2 = r2ij*r2ik*r2jk
                  r1 = sqrt(r2)
                  r3 = r2 * r1
                  r5 = r3 * r2

                  fdmp = 1.0_wp / (1.0_wp + 6.0_wp * (r0 / r1)**alp3)
                  ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)&
                     & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3

                  rr = ang*fdmp

                  dfdmp = -2.0_wp * alp * (r0 / r1)**alp3 * fdmp**2

                  ! d/drij
                  dang = -0.375_wp * (r2ij**3 + r2ij**2 * (r2jk + r2ik)&
                     & + r2ij * (3.0_wp * r2jk**2 + 2.0_wp * r2jk*r2ik&
                     & + 3.0_wp * r2ik**2)&
                     & - 5.0_wp * (r2jk - r2ik)**2 * (r2jk + r2ik)) / r5
                  dE0 = rr * c9
                  dGij(:) = sw * c9 * (-dang*fdmp + ang*dfdmp) / r2ij * vij &
                     & - dE0 * dswijdr / rij * swik * swjk * vij

                  ! d/drik
                  dang = -0.375_wp * (r2ik**3 + r2ik**2 * (r2jk + r2ij)&
                     & + r2ik * (3.0_wp * r2jk**2 + 2.0_wp * r2jk * r2ij&
                     & + 3.0_wp * r2ij**2)&
                     & - 5.0_wp * (r2jk - r2ij)**2 * (r2jk + r2ij)) / r5
                  dGik(:) = sw * c9 * (-dang * fdmp + ang * dfdmp) / r2ik * vik &
                     & - dE0 * dswikdr / rik * swij * swjk * vik

                  ! d/drjk
                  dang = -0.375_wp * (r2jk**3 + r2jk**2*(r2ik + r2ij)&
                     & + r2jk * (3.0_wp * r2ik**2 + 2.0_wp * r2ik * r2ij&
                     & + 3.0_wp * r2ij**2)&
                     & - 5.0_wp * (r2ik - r2ij)**2 * (r2ik + r2ij)) / r5
                  dGjk(:) = sw * c9 * (-dang * fdmp + ang * dfdmp) / r2jk * vjk &
                     & - dE0 * dswjkdr / rjk * swij * swik * vjk

                  dE = dE0 * triple * sw
                  dE_third = dE * third
                  energy_local(iat) = energy_local(iat) - dE_third
                  energy_local(jat) = energy_local(jat) - dE_third
                  energy_local(kat) = energy_local(kat) - dE_third

                  gradient_local(:, iat) = gradient_local(:, iat) &
                     & - (dGij + dGik) * triple
                  gradient_local(:, jat) = gradient_local(:, jat) &
                     & + (dGij - dGjk) * triple
                  gradient_local(:, kat) = gradient_local(:, kat) &
                     & + (dGik + dGjk) * triple

                  dS(:, :) = spread(dGij, 1, 3) * spread(vij, 2, 3)&
                     & + spread(dGik, 1, 3) * spread(vik, 2, 3)&
                     & + spread(dGjk, 1, 3) * spread(vjk, 2, 3)

                  sigma_local(:, :) = sigma_local + dS * triple

                  dEdcn_local(iat) = dEdcn_local(iat) - dE * 0.5_wp &
                     & * (dc6dcn(iat, jat) / c6ij + dc6dcn(iat, kat) / c6ik)
                  dEdcn_local(jat) = dEdcn_local(jat) - dE * 0.5_wp &
                     & * (dc6dcn(jat, iat) / c6ij + dc6dcn(jat, kat) / c6jk)
                  dEdcn_local(kat) = dEdcn_local(kat) - dE * 0.5_wp &
                     & * (dc6dcn(kat, iat) / c6ik + dc6dcn(kat, jat) / c6jk)

                  dEdq_local(iat) = dEdq_local(iat) - dE * 0.5_wp &
                     & * (dc6dq(iat, jat) / c6ij + dc6dq(iat, kat) / c6ik)
                  dEdq_local(jat) = dEdq_local(jat) - dE * 0.5_wp &
                     & * (dc6dq(jat, iat) / c6ij + dc6dq(jat, kat) / c6jk)
                  dEdq_local(kat) = dEdq_local(kat) - dE * 0.5_wp &
                     & * (dc6dq(kat, iat) / c6ik + dc6dq(kat, jat) / c6jk)
               end do
            end do
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_atm_dispersion_derivs_)
   energy(:) = energy(:) + energy_local(:)
   dEdcn(:) = dEdcn(:) + dEdcn_local(:)
   dEdq(:) = dEdq(:) + dEdq_local(:)
   gradient(:, :) = gradient(:, :) + gradient_local(:, :)
   sigma(:, :) = sigma(:, :) + sigma_local(:, :)
   !$omp end critical (get_atm_dispersion_derivs_)
   deallocate(energy_local)
   deallocate(dEdcn_local)
   deallocate(dEdq_local)
   deallocate(gradient_local)
   deallocate(sigma_local)
   !$omp end parallel

end subroutine get_atm_dispersion_derivs


!> Logic exercise to distribute a triple energy to atomwise energies.
elemental function triple_scale(ii, jj, kk) result(triple)

   !> Atom indices
   integer, intent(in) :: ii, jj, kk

   !> Fraction of energy
   real(wp) :: triple

   if (ii == jj) then
      if (ii == kk) then
         ! ii'i" -> 1/6
         triple = 1.0_wp/6.0_wp
      else
         ! ii'j -> 1/2
         triple = 0.5_wp
      end if
   else
      if (ii /= kk .and. jj /= kk) then
         ! ijk -> 1 (full)
         triple = 1.0_wp
      else
         ! ijj' and iji' -> 1/2
         triple = 0.5_wp
      end if
   end if

end function triple_scale


end module dftd4_damping_atm
