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

!> Implementation of the rational (Becke--Johnson) damping function.
module dftd4_damping_rational
   use dftd4_damping, only : damping_param
   use dftd4_damping_atm, only : get_atm_dispersion
   use dftd4_data, only : get_r4r2_val
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none
   private

   public :: rational_damping_param


   !> Rational (Becke-Johnson) damping model
   type, extends(damping_param) :: rational_damping_param
      real(wp) :: s6 = 1.0_wp
      real(wp) :: s8
      real(wp) :: s9 = 1.0_wp
      real(wp) :: a1
      real(wp) :: a2
      real(wp) :: alp = 16.0_wp
   contains

      !> Evaluate pairwise dispersion energy expression
      procedure :: get_dispersion2

      !> Evaluate ATM three-body dispersion energy expression
      procedure :: get_dispersion3

   end type rational_damping_param


contains


!> Evaluation of the dispersion energy expression
subroutine get_dispersion2(self, mol, trans, cutoff, r4r2, c6, dc6dcn, dc6dq, &
      & energy, dEdcn, dEdq, gradient, sigma)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

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
   integer :: iat, jat, izp, jzp, jtr
   real(wp) :: vec(3), r2, cutoff2, r0ij, rrij, c6ij, t6, t8, d6, d8, edisp, gdisp
   real(wp) :: dE, dG(3), dS(3, 3)

   if (abs(self%s6) < epsilon(1.0_wp) .and. abs(self%s8) < epsilon(1.0_wp)) return
   grad = present(dc6dcn) .and. present(dEdcn) .and. present(dc6dq) &
      & .and. present(dEdq) .and. present(gradient) .and. present(sigma)
   cutoff2 = cutoff*cutoff

   if (grad) then
      !$omp parallel do schedule(runtime) default(none) &
      !$omp reduction(+:energy, gradient, sigma, dEdcn, dEdq) &
      !$omp shared(mol, self, c6, dc6dcn, dc6dq, trans, cutoff2, r4r2) &
      !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, c6ij, t6, t8, &
      !$omp& d6, d8, edisp, gdisp, dE, dG, dS)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            rrij = 3*r4r2(izp)*r4r2(jzp)
            r0ij = self%a1 * sqrt(rrij) + self%a2
            c6ij = c6(jat, iat)
            do jtr = 1, size(trans, 2)
               vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
               r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
               if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

               t6 = 1.0_wp/(r2**3 + r0ij**6)
               t8 = 1.0_wp/(r2**4 + r0ij**8)

               d6 = -6*r2**2*t6**2
               d8 = -8*r2**3*t8**2

               edisp = self%s6*t6 + self%s8*rrij*t8
               gdisp = self%s6*d6 + self%s8*rrij*d8

               dE = -c6ij*edisp * 0.5_wp
               dG(:) = -c6ij*gdisp*vec
               dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3) * 0.5_wp

               energy(iat) = energy(iat) + dE
               dEdcn(iat) = dEdcn(iat) - dc6dcn(iat, jat) * edisp
               dEdq(iat) = dEdq(iat) - dc6dq(iat, jat) * edisp
               sigma(:, :) = sigma + dS
               if (iat /= jat) then
                  energy(jat) = energy(jat) + dE
                  dEdcn(jat) = dEdcn(jat) - dc6dcn(jat, iat) * edisp
                  dEdq(jat) = dEdq(jat) - dc6dq(jat, iat) * edisp
                  gradient(:, iat) = gradient(:, iat) + dG
                  gradient(:, jat) = gradient(:, jat) - dG
                  sigma(:, :) = sigma + dS
               end if
            end do
         end do
      end do
   else
      !$omp parallel do schedule(runtime) default(none) reduction(+:energy) &
      !$omp shared(mol, self, c6, trans, cutoff2, r4r2) &
      !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, c6ij, &
      !$omp& t6, t8, edisp, dE)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            rrij = 3*r4r2(izp)*r4r2(jzp)
            r0ij = self%a1 * sqrt(rrij) + self%a2
            c6ij = c6(jat, iat)
            do jtr = 1, size(trans, 2)
               vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
               r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
               if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

               t6 = 1.0_wp/(r2**3 + r0ij**6)
               t8 = 1.0_wp/(r2**4 + r0ij**8)

               edisp = self%s6*t6 + self%s8*rrij*t8

               dE = -c6ij*edisp * 0.5_wp

               energy(iat) = energy(iat) + dE
               if (iat /= jat) then
                  energy(jat) = energy(jat) + dE
               end if
            end do
         end do
      end do
   end if

end subroutine get_dispersion2


!> Evaluation of the dispersion energy expression
subroutine get_dispersion3(self, mol, trans, cutoff, r4r2, c6, dc6dcn, dc6dq, &
      & energy, dEdcn, dEdq, gradient, sigma)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

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

   call get_atm_dispersion(mol, trans, cutoff, self%s9, self%a1, self%a2, &
      & self%alp, r4r2, c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, &
      & gradient, sigma)

end subroutine get_dispersion3


end module dftd4_damping_rational
