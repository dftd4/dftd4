! This file is part of dftd4.
!
! Copyright (C) 2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> provides definition of calculation options type
module class_set
   use iso_fortran_env, wp => real64
   use iso_c_binding
   use class_param, only : dftd_parameter
   implicit none

   public :: options
   public :: dftd_options
   public :: c_dftd_options
   public :: assignment(=)
   private

!> calculation setup
   type :: dftd_options
      sequence
      integer  :: lmbd = -1                 !< kind of non-additivity correction
      integer  :: refq = -1                 !< kind of charge model
      real(wp) :: wf = 0.0_wp               !< weighting factor
      real(wp) :: g_a = 0.0_wp              !< charge scale height
      real(wp) :: g_c = 0.0_wp              !< charge scale steepness
      logical  :: lmolpol = .false.         !< calculate molecular properties?
      logical  :: lenergy = .false.         !< calculate dispersion energy?
      logical  :: lgradient = .false.       !< calculate dispersion gradient?
      logical  :: lhessian = .false.        !< calculate dispersion hessian?
      logical  :: verbose = .false.         !< print more information
      logical  :: veryverbose = .false.     !< clutter the screen more
      logical  :: silent = .false.          !< clutter the screen less
   end type dftd_options

!> calculation setup
   type,bind(C) :: c_dftd_options
      integer(c_int)  :: lmbd = -1             !< kind of non-additivity correction
      integer(c_int)  :: refq = -1             !< kind of charge model
      real(c_double)  :: wf = 0.0_wp           !< weighting factor
      real(c_double)  :: g_a = 0.0_wp          !< charge scale height
      real(c_double)  :: g_c = 0.0_wp          !< charge scale steepness
      logical(c_bool) :: lmolpol = .false.     !< calculate molecular properties?
      logical(c_bool) :: lenergy = .false.     !< calculate dispersion energy?
      logical(c_bool) :: lgradient = .false.   !< calculate dispersion gradient?
      logical(c_bool) :: lhessian = .false.    !< calculate dispersion hessian?
      logical(c_bool) :: verbose = .false.     !< print more information
      logical(c_bool) :: veryverbose = .false. !< clutter the screen more
      logical(c_bool) :: silent = .false.      !< clutter the screen less
   end type c_dftd_options

   interface assignment(=)
      module procedure :: convert_dftd_options_c_to_f
      module procedure :: convert_dftd_options_f_to_c
   end interface


!> contains all calculation options for the DFT-D run
   type :: options
      character(len=:),allocatable :: fname !< geometry file
      character(len=:),allocatable :: func  !< functional name
      type(dftd_parameter) :: dparam        !< damping parameters
      logical  :: inparam = .false.         !< parameters given by command line
      integer  :: chrg = 0                  !< molecular charge
      logical  :: inchrg = .false.          !< charge given by command line
      integer  :: lmbd = 3                  !< kind of non-additivity correction
      integer  :: refq = -1                 !< kind of charge model
      real(wp) :: wf = 6.0_wp               !< weighting factor
      real(wp) :: g_a = 3.0_wp              !< charge scale height
      real(wp) :: g_c = 2.0_wp              !< charge scale steepness
      logical  :: lperiodic = .false.       !< use periodic boundary conditions
      logical  :: lorca = .false.           !< calculation for ORCA
      logical  :: ltmer = .false.           !< write TMER2++ output
      logical  :: lmolpol = .false.         !< calculate molecular properties?
      logical  :: lenergy = .false.         !< calculate dispersion energy?
      logical  :: lgradient = .false.       !< calculate dispersion gradient?
      logical  :: lhessian = .false.        !< calculate dispersion hessian?
      logical  :: verbose = .false.         !< print more information
      logical  :: veryverbose = .false.     !< clutter the screen more
      logical  :: silent = .false.          !< clutter the screen less
   contains
      procedure :: default_options
      procedure :: export => export_dftd_options
   end type

contains

subroutine default_options(self)
   implicit none
   class(options),intent(inout) :: self
   self%inparam     = .false.
   self%chrg        = 0
   self%inchrg      = .false.
   self%lmbd        = 3
   self%wf          = 6.0_wp
   self%g_a         = 3.0_wp
   self%g_c         = 2.0_wp
   self%lorca       = .false.
   self%ltmer       = .false.
   self%lmolpol     = .false.
   self%lenergy     = .false.
   self%lgradient   = .false.
   self%lhessian    = .false.
   self%verbose     = .false.
   self%veryverbose = .false.
   self%silent      = .false.
end subroutine default_options

pure function export_dftd_options(self) result(opt)
   implicit none
   class(options),intent(in) :: self
   type(dftd_options) :: opt

   opt%lmbd        = self%lmbd
   opt%refq        = self%refq
   opt%wf          = self%wf
   opt%g_a         = self%g_a
   opt%g_c         = self%g_c
   opt%lmolpol     = self%lmolpol
   opt%lenergy     = self%lenergy
   opt%lgradient   = self%lgradient
   opt%lhessian    = self%lhessian
   opt%verbose     = self%verbose
   opt%veryverbose = self%veryverbose
   opt%silent      = self%silent

end function export_dftd_options

pure elemental subroutine convert_dftd_options_c_to_f &
      (f_dopt,c_dopt)
   implicit none
   type(c_dftd_options),intent(in)  :: c_dopt
   type(dftd_options),  intent(out) :: f_dopt
   f_dopt%lmbd        = c_dopt%lmbd
   f_dopt%refq        = c_dopt%refq
   f_dopt%wf          = c_dopt%wf
   f_dopt%g_a         = c_dopt%g_a
   f_dopt%g_c         = c_dopt%g_c
   f_dopt%lmolpol     = c_dopt%lmolpol
   f_dopt%lenergy     = c_dopt%lenergy
   f_dopt%lgradient   = c_dopt%lgradient
   f_dopt%lhessian    = c_dopt%lhessian
   f_dopt%verbose     = c_dopt%verbose
   f_dopt%veryverbose = c_dopt%veryverbose
   f_dopt%silent      = c_dopt%silent
end subroutine convert_dftd_options_c_to_f

pure elemental subroutine convert_dftd_options_f_to_c &
      (c_dopt,f_dopt)
   implicit none
   type(c_dftd_options),intent(out) :: c_dopt
   type(dftd_options),  intent(in)  :: f_dopt
   c_dopt%lmbd        = f_dopt%lmbd
   c_dopt%refq        = f_dopt%refq
   c_dopt%wf          = f_dopt%wf
   c_dopt%g_a         = f_dopt%g_a
   c_dopt%g_c         = f_dopt%g_c
   c_dopt%lmolpol     = f_dopt%lmolpol
   c_dopt%lenergy     = f_dopt%lenergy
   c_dopt%lgradient   = f_dopt%lgradient
   c_dopt%lhessian    = f_dopt%lhessian
   c_dopt%verbose     = f_dopt%verbose
   c_dopt%veryverbose = f_dopt%veryverbose
   c_dopt%silent      = f_dopt%silent
end subroutine convert_dftd_options_f_to_c

end module class_set
