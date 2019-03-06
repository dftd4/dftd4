!> @brief provides definition of calculation options type
module class_set
   use iso_fortran_env, wp => real64
   use class_param, only : dftd_parameter
   implicit none

   public :: options
   private

!> @brief contains all calculation options for the DFT-D run
   type :: options
      character(len=:),allocatable :: fname !< geometry file
      character(len=:),allocatable :: func  !< functional name
      type(dftd_parameter) :: dparam        !< damping parameters
      logical  :: inparam = .false.         !< parameters given by command line
      integer  :: chrg = 0                  !< molecular charge
      logical  :: inchrg = .false.          !< charge given by command line
      integer  :: lmbd = 3                  !< kind of non-additivity correction
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
   end type

end module class_set
