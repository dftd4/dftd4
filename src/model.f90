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

!> Re-export of all dispersion models
module dftd4_model
   use dftd4_model_type, only : dispersion_model, d4_ref
   use dftd4_model_d4, only : d4_model, new_d4_model
   use dftd4_model_d4s, only : d4s_model, new_d4s_model
   implicit none
   private

   public :: dispersion_model, d4_ref
   public :: d4_model, new_d4_model
   public :: d4s_model, new_d4s_model

end module dftd4_model
