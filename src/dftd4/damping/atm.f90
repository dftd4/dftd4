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
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none

   public :: get_atm_dispersion


contains


!> Evaluation of the dispersion energy expression
subroutine get_atm_dispersion(mol, trans, cutoff, s9, a1, a2, alp, r4r2, &
      & c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient, sigma)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

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

   logical :: grad

   if (abs(s9) < epsilon(1.0_wp)) return
   grad = present(dc6dcn) .and. present(dEdcn) .and. present(dc6dq) &
      & .and. present(dEdq) .and. present(gradient) .and. present(sigma)

   if (grad) then
      call get_atm_dispersion_derivs(mol, trans, cutoff, s9, a1, a2, alp, r4r2, &
         & c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient, sigma)
   else
      call get_atm_dispersion_energy(mol, trans, cutoff, s9, a1, a2, alp, r4r2, &
         & c6, energy)
   end if

end subroutine get_atm_dispersion


!> Evaluation of the dispersion energy expression
subroutine get_atm_dispersion_energy(mol, trans, cutoff, s9, a1, a2, alp, r4r2, &
      & c6, energy)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

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

   integer :: iat, jat, kat, izp, jzp, kzp, jtr, ktr
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple
   real(wp) :: r0ij, r0jk, r0ik, r0, r1, r2, r3, r5, rr, fdmp, ang
   real(wp) :: cutoff2, c9, dE

   cutoff2 = cutoff*cutoff

   !$omp parallel do schedule(runtime) default(none) reduction(+:energy) &
   !$omp shared(mol, trans, c6, s9, a1, a2, alp, r4r2, cutoff2) &
   !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
   !$omp& r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple, r0ij, r0jk, r0ik, r0, &
   !$omp& r1, r2, r3, r5, rr, fdmp, ang, c9, dE)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         c6ij = c6(jat, iat)
         r0ij = a1 * sqrt(3*r4r2(jzp)*r4r2(izp)) + a2
         do jtr = 1, size(trans, 2)
            vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
            if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
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
                  vjk(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, jat) &
                     & - trans(:, jtr)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                  r2 = r2ij*r2ik*r2jk
                  r1 = sqrt(r2)
                  r3 = r2 * r1
                  r5 = r3 * r2

                  fdmp = 1.0_wp / (1.0_wp + 6.0_wp * (r0 / r1)**(alp / 3.0_wp))
                  ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)&
                     & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3

                  rr = ang*fdmp

                  dE = rr * c9 * triple
                  energy(iat) = energy(iat) - dE/3
                  energy(jat) = energy(jat) - dE/3
                  energy(kat) = energy(kat) - dE/3
               end do
            end do
         end do
      end do
   end do

end subroutine get_atm_dispersion_energy


!> Evaluation of the dispersion energy expression
subroutine get_atm_dispersion_derivs(mol, trans, cutoff, s9, a1, a2, alp, r4r2, &
      & c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient, sigma)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

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

   integer :: iat, jat, kat, izp, jzp, kzp, jtr, ktr
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple
   real(wp) :: r0ij, r0jk, r0ik, r0, r1, r2, r3, r5, rr, fdmp, dfdmp, ang, dang
   real(wp) :: cutoff2, c9, dE, dGij(3), dGjk(3), dGik(3), dS(3, 3)

   cutoff2 = cutoff*cutoff

   !$omp parallel do schedule(runtime) default(none) &
   !$omp reduction(+:energy, gradient, sigma, dEdcn, dEdq) &
   !$omp shared(mol, trans, c6, s9, a1, a2, alp, r4r2, cutoff2, dc6dcn, dc6dq) &
   !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
   !$omp& r2ij, r2jk, r2ik, c6ij, c6jk, c6ik, triple, r0ij, r0jk, r0ik, r0, &
   !$omp& r1, r2, r3, r5, rr, fdmp, dfdmp, ang, dang, c9, dE, dGij, dGjk, &
   !$omp& dGik, dS)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         c6ij = c6(jat, iat)
         r0ij = a1 * sqrt(3*r4r2(jzp)*r4r2(izp)) + a2
         do jtr = 1, size(trans, 2)
            vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
            if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
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
                  vjk(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, jat) &
                     & - trans(:, jtr)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                  r2 = r2ij*r2ik*r2jk
                  r1 = sqrt(r2)
                  r3 = r2 * r1
                  r5 = r3 * r2

                  fdmp = 1.0_wp / (1.0_wp + 6.0_wp * (r0 / r1)**(alp / 3.0_wp))
                  ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)&
                     & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3

                  rr = ang*fdmp

                  dfdmp = -2.0_wp * alp * (r0 / r1)**(alp / 3.0_wp) * fdmp**2

                  ! d/drij
                  dang = -0.375_wp * (r2ij**3 + r2ij**2 * (r2jk + r2ik)&
                     & + r2ij * (3.0_wp * r2jk**2 + 2.0_wp * r2jk*r2ik&
                     & + 3.0_wp * r2ik**2)&
                     & - 5.0_wp * (r2jk - r2ik)**2 * (r2jk + r2ik)) / r5
                  dGij(:) = c9 * (-dang*fdmp + ang*dfdmp) / r2ij * vij

                  ! d/drik
                  dang = -0.375_wp * (r2ik**3 + r2ik**2 * (r2jk + r2ij)&
                     & + r2ik * (3.0_wp * r2jk**2 + 2.0_wp * r2jk * r2ij&
                     & + 3.0_wp * r2ij**2)&
                     & - 5.0_wp * (r2jk - r2ij)**2 * (r2jk + r2ij)) / r5
                  dGik(:) = c9 * (-dang * fdmp + ang * dfdmp) / r2ik * vik

                  ! d/drjk
                  dang = -0.375_wp * (r2jk**3 + r2jk**2*(r2ik + r2ij)&
                     & + r2jk * (3.0_wp * r2ik**2 + 2.0_wp * r2ik * r2ij&
                     & + 3.0_wp * r2ij**2)&
                     & - 5.0_wp * (r2ik - r2ij)**2 * (r2ik + r2ij)) / r5
                  dGjk(:) = c9 * (-dang * fdmp + ang * dfdmp) / r2jk * vjk

                  dE = rr * c9 * triple
                  energy(iat) = energy(iat) - dE/3
                  energy(jat) = energy(jat) - dE/3
                  energy(kat) = energy(kat) - dE/3

                  gradient(:, iat) = gradient(:, iat) - dGij - dGik
                  gradient(:, jat) = gradient(:, jat) + dGij - dGjk
                  gradient(:, kat) = gradient(:, kat) + dGik + dGjk

                  dS(:, :) = spread(dGij, 1, 3) * spread(vij, 2, 3)&
                     & + spread(dGik, 1, 3) * spread(vik, 2, 3)&
                     & + spread(dGjk, 1, 3) * spread(vjk, 2, 3)

                  sigma(:, :) = sigma + dS * triple

                  dEdcn(iat) = dEdcn(iat) - dE * 0.5_wp &
                     & * (dc6dcn(iat, jat) / c6ij + dc6dcn(iat, kat) / c6ik)
                  dEdcn(jat) = dEdcn(jat) - dE * 0.5_wp &
                     & * (dc6dcn(jat, iat) / c6ij + dc6dcn(jat, kat) / c6jk)
                  dEdcn(kat) = dEdcn(kat) - dE * 0.5_wp &
                     & * (dc6dcn(kat, iat) / c6ik + dc6dcn(kat, jat) / c6jk)

                  dEdq(iat) = dEdq(iat) - dE * 0.5_wp &
                     & * (dc6dq(iat, jat) / c6ij + dc6dq(iat, kat) / c6ik)
                  dEdq(jat) = dEdq(jat) - dE * 0.5_wp &
                     & * (dc6dq(jat, iat) / c6ij + dc6dq(jat, kat) / c6jk)
                  dEdq(kat) = dEdq(kat) - dE * 0.5_wp &
                     & * (dc6dq(kat, iat) / c6ik + dc6dq(kat, jat) / c6jk)
               end do
            end do
         end do
      end do
   end do

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
