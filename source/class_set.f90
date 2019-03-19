!> @brief provides definition of calculation options type
module class_set
   use iso_fortran_env, wp => real64
   use class_param, only : dftd_parameter
   implicit none

   public :: options
   public :: dftd_options
   private

!> @brief calculation setup
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

!> @brief contains all calculation options for the DFT-D run
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

end module class_set
