module class_set
   use iso_fortran_env, wp => real64
   use class_param, only : dftd_parameter
   implicit none

   public :: options
   private

   type :: options
      character(len=:),allocatable :: fname
      character(len=:),allocatable :: func
      type(dftd_parameter) :: dparam
      logical  :: inparam = .false.
      integer  :: chrg = 0
      logical  :: inchrg = .false.
      integer  :: lmbd = 3
      real(wp) :: wf = 6.0_wp
      real(wp) :: g_a = 3.0_wp
      real(wp) :: g_c = 2.0_wp
      logical  :: lorca = .false.
      logical  :: ltmer = .false.
      logical  :: lmolpol = .false.
      logical  :: lenergy = .false.
      logical  :: lgradient = .false.
      logical  :: lhessian = .false.
      logical  :: verbose = .false.
      logical  :: veryverbose = .false.
      logical  :: silent = .false.
   end type

end module class_set
