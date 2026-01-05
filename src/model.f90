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
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use dftd4_utils, only : lowercase
   use dftd4_model_type, only : dispersion_model, d4_qmod
   use dftd4_model_d4, only : d4_model, new_d4_model
   use dftd4_model_d4s, only : d4s_model, new_d4s_model
   implicit none
   private

   public :: dispersion_model, d4_qmod
   public :: d4_model, new_d4_model
   public :: d4s_model, new_d4s_model
   public :: new_dispersion_model


contains


!> Wrapper for creating a new dispersion model (D4 or D4S) from molecular 
!> structure input using a given model string. Defaults to D4 if no model
!> is specified.
subroutine new_dispersion_model(error, d4, mol, model, ga, gc, wf, qmod)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Dispersion model to be returned
   class(dispersion_model), allocatable, intent(out) :: d4

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dispersion model to be used
   character(len=*), intent(in), optional :: model

   !> Charge scaling height
   real(wp), intent(in), optional :: ga

   !> Charge scaling steepness
   real(wp), intent(in), optional :: gc

   !> Weighting factor for coordination number interpolation
   real(wp), intent(in), optional :: wf

   !> Charge model selection
   integer, intent(in), optional :: qmod

   character(len=:), allocatable :: mdl

   if (present(model)) then
      mdl = lowercase(trim(model)) 
   else
      mdl = "d4"
   end if

   if(mdl == "d4") then
      block 
         type(d4_model), allocatable :: tmp
         allocate(tmp)
         call new_d4_model(error, tmp, mol, ga=ga, gc=gc, wf=wf, qmod=qmod)
         call move_alloc(tmp, d4)
      end block 
   else if(mdl == "d4s") then
      block 
         type(d4s_model), allocatable :: tmp
         allocate(tmp)
         call new_d4s_model(error, tmp, mol, ga=ga, gc=gc, qmod=qmod)
         call move_alloc(tmp, d4)
      end block
   else
      call fatal_error(error, "Unknown model selected")
   end if

end subroutine new_dispersion_model


end module dftd4_model
