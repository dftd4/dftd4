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

!> Utility functions for the dispersion models 
module dftd4_model_utils
   use ieee_arithmetic, only : ieee_is_nan
   use mctc_env, only : wp
   implicit none
   
   public :: is_exceptional, weight_cn, zeta, dzeta, trapzd

contains


!> Check whether we are dealing with an exceptional value, NaN or Inf
elemental function is_exceptional(val)
   real(wp), intent(in) :: val
   logical :: is_exceptional
   is_exceptional = ieee_is_nan(val) .or. abs(val) > huge(val)
end function is_exceptional


elemental function weight_cn(wf,cn,cnref) result(cngw)
   real(wp),intent(in) :: wf, cn, cnref
   real(wp) :: cngw
   intrinsic :: exp
   cngw = exp ( -wf * ( cn - cnref )**2 )
end function weight_cn

!> charge scaling function
elemental function zeta(a, c, qref, qmod)
   real(wp), intent(in) :: a
   real(wp), intent(in) :: c
   real(wp), intent(in) :: qref
   real(wp), intent(in) :: qmod
   real(wp) :: zeta

   intrinsic :: exp

   if (qmod < 0.0_wp) then
      zeta = exp( a )
   else
      zeta = exp( a * ( 1.0_wp - exp( c * ( 1.0_wp - qref/qmod ) ) ) )
   endif

end function zeta

!> derivative of charge scaling function w.r.t. charge
elemental function dzeta(a, c, qref, qmod)
   real(wp), intent(in) :: a
   real(wp), intent(in) :: c
   real(wp), intent(in) :: qref
   real(wp), intent(in) :: qmod
   real(wp) :: dzeta

   intrinsic :: exp

   if (qmod < 0.0_wp) then
      dzeta = 0.0_wp
   else
      dzeta = - a * c * exp( c * ( 1.0_wp - qref/qmod ) ) &
         & * zeta(a,c,qref,qmod) * qref / ( qmod**2 )
   endif

end function dzeta

!> numerical Casimir--Polder integration
pure function trapzd(pol)
   real(wp), intent(in) :: pol(23)
   real(wp) :: trapzd

   real(wp), parameter :: freq(23) = [ &
      & 0.000001_wp, 0.050000_wp, 0.100000_wp, &
      & 0.200000_wp, 0.300000_wp, 0.400000_wp, &
      & 0.500000_wp, 0.600000_wp, 0.700000_wp, &
      & 0.800000_wp, 0.900000_wp, 1.000000_wp, &
      & 1.200000_wp, 1.400000_wp, 1.600000_wp, &
      & 1.800000_wp, 2.000000_wp, 2.500000_wp, &
      & 3.000000_wp, 4.000000_wp, 5.000000_wp, &
      & 7.500000_wp, 10.00000_wp]
   real(wp), parameter :: weights(23) = 0.5_wp * [ &
      &  ( freq (2) - freq (1) ),  &
      &  ( freq (2) - freq (1) ) + ( freq (3) - freq (2) ),  &
      &  ( freq (3) - freq (2) ) + ( freq (4) - freq (3) ),  &
      &  ( freq (4) - freq (3) ) + ( freq (5) - freq (4) ),  &
      &  ( freq (5) - freq (4) ) + ( freq (6) - freq (5) ),  &
      &  ( freq (6) - freq (5) ) + ( freq (7) - freq (6) ),  &
      &  ( freq (7) - freq (6) ) + ( freq (8) - freq (7) ),  &
      &  ( freq (8) - freq (7) ) + ( freq (9) - freq (8) ),  &
      &  ( freq (9) - freq (8) ) + ( freq(10) - freq (9) ),  &
      &  ( freq(10) - freq (9) ) + ( freq(11) - freq(10) ),  &
      &  ( freq(11) - freq(10) ) + ( freq(12) - freq(11) ),  &
      &  ( freq(12) - freq(11) ) + ( freq(13) - freq(12) ),  &
      &  ( freq(13) - freq(12) ) + ( freq(14) - freq(13) ),  &
      &  ( freq(14) - freq(13) ) + ( freq(15) - freq(14) ),  &
      &  ( freq(15) - freq(14) ) + ( freq(16) - freq(15) ),  &
      &  ( freq(16) - freq(15) ) + ( freq(17) - freq(16) ),  &
      &  ( freq(17) - freq(16) ) + ( freq(18) - freq(17) ),  &
      &  ( freq(18) - freq(17) ) + ( freq(19) - freq(18) ),  &
      &  ( freq(19) - freq(18) ) + ( freq(20) - freq(19) ),  &
      &  ( freq(20) - freq(19) ) + ( freq(21) - freq(20) ),  &
      &  ( freq(21) - freq(20) ) + ( freq(22) - freq(21) ),  &
      &  ( freq(22) - freq(21) ) + ( freq(23) - freq(22) ),  &
      &  ( freq(23) - freq(22) ) ]

   trapzd = sum(pol*weights)

end function trapzd


end module dftd4_model_utils
