! ------------------------------------------------------------------------
!  Purpose:
!> @brief reads in command line arguments and parses them
! ------------------------------------------------------------------------
subroutine read_commandline_arguments(set)
   use iso_fortran_env, only : wp => real64
!$ use omp_lib

! ------------------------------------------------------------------------
!  general purpose library
! ------------------------------------------------------------------------
   use mctc_systools
   use mctc_readin

! ------------------------------------------------------------------------
!  Class definition
! ------------------------------------------------------------------------
   use class_set

! ------------------------------------------------------------------------
!  get parameters
! ------------------------------------------------------------------------
   use dftd4

   implicit none
   type(options) :: set !< options for calculation

   character(len=:),allocatable :: arg
   character(len=:),allocatable :: sec
   logical  :: getopts
   integer  :: skip
   integer  :: iarg
   integer  :: nargs
   logical  :: exist
   integer  :: nproc

   integer  :: idum
   logical  :: ldum
   real(wp) :: ddum

   nargs = command_argument_count()
   if (nargs.eq.0) then
      call help
      call terminate(1)
   endif

   getopts = .true.
   skip = 0
   do iarg = 1, nargs
      if (skip.gt.0) then
         skip = skip-1
         cycle
      endif
      call rdarg(iarg,arg)
      if (arg.eq.'--') then
         getopts=.false.
         cycle
      endif
      if (getopts) then
!        write(output_unit,'(i0,'':'',x,a)') iarg,arg ! debugging stuff

         if ((len(arg).gt.2).and. &
         &  (index(arg,'--').eq.0).and.(index(arg,'-').eq.1)) then
            call raise('S',"the use of '"//arg//"' is discouraged, "// &
            &              "please use '-"//arg//"' next time")
            arg = '-'//arg
         endif
         select case(arg)
! ------------------------------------------------------------------------
!  check the help, version and citation flag, exit if found
! ------------------------------------------------------------------------
         case('-h','--help')
            call help
            call terminate(0)
         case(     '--citation')
            call dftd4_citation
            call terminate(0)
         case(     '--license')
            call gpl_license
            call terminate(0)
         case(     '--version')
            call dftd4_header(.true.)
            call terminate(0)

! ------------------------------------------------------------------------
!  now for true options, like verbosity or printouts
! ------------------------------------------------------------------------
         case('-v','--verbose')
            set%verbose = .true.
         case('-V','--very-verbose')
            set%verbose = .true.
            set%veryverbose = .true.
         case('-s','--silent')
            set%silent = .true.
      !  OMP option
      !$ case('-P','--parallel')
      !$    skip = 1
      !$    call rdarg(iarg+1,sec)
      !$    if (get_value(sec,idum)) then
      !$       nproc = omp_get_num_threads()
      !$       if (idum.gt.nproc) &
      !$       & call raise('S','Process number higher than OMP_NUM_THREADS, '//&
      !$       &                'I hope you know what you are doing.')
      !$       call omp_set_num_threads(idum)
      !$    endif
         case('-c','--chrg')
            skip = 1
            call rdarg(iarg+1,sec)
            if (get_value(sec,idum)) then
               set%chrg = idum
            else
               call raise('E',"Could not read charge from '"//sec//"'")
            endif
            set%inchrg = .true.
         case('-f','--func')
            skip = 1
            call rdarg(iarg+1,sec)
            set%func = sec
         case(     '--molc6')
            set%lmolpol = .true.
         case(     '--param'); skip = 4
            call rdarg(iarg+1,arg); if (get_value(arg,ddum)) set%dparam%s6 = ddum
            call rdarg(iarg+2,arg); if (get_value(arg,ddum)) set%dparam%s8 = ddum
            call rdarg(iarg+3,arg); if (get_value(arg,ddum)) set%dparam%a1 = ddum
            call rdarg(iarg+4,arg); if (get_value(arg,ddum)) set%dparam%a2 = ddum
            ! sanity check
            if (set%dparam%s6.lt.0.0_wp.or.set%dparam%s6.gt.1.0_wp) then
               call raise('S','unphysical s6 chosen, use on your own risk')
            endif
            if (set%dparam%s8.lt.0.0_wp.or.set%dparam%s8.gt.3.0_wp) then
               call raise('S','unphysical s8 chosen, use on your own risk')
            endif
            if (set%dparam%a1.lt.0.0_wp.or.set%dparam%a1.gt.1.0_wp) then
               call raise('S','unphysical a1 chosen, use on your own risk')
            endif
            if (set%dparam%a2.lt.0.0_wp.or.set%dparam%a2.gt.7.0_wp) then
               call raise('S','unphysical a2 chosen, use on your own risk')
            endif
            set%inparam = .true.
         case(     '--s10'); skip = 1
            call rdarg(iarg+1,arg); if (get_value(arg,ddum)) set%dparam%s10 = ddum
         case(     '--zeta'); skip = 2
            call rdarg(iarg+1,arg); if (get_value(arg,ddum)) set%g_a = ddum
            call rdarg(iarg+2,arg); if (get_value(arg,ddum)) set%g_c = ddum
         case(     '--mbdscale'); skip = 1
            call rdarg(iarg+1,arg); if (get_value(arg,ddum)) set%dparam%s9 = ddum
         case(     '--wfactor'); skip = 1
            call rdarg(iarg+1,arg); if (get_value(arg,ddum)) set%wf = ddum
         case(     '--tmer'); set%ltmer = .true.
         case('-m','--mbd');            set%lmbd = p_mbd_rpalike
         case('-2','--nomb');           set%lmbd = p_mbd_none
         case('-3','--abc');            set%lmbd = p_mbd_approx_atm
         case('-g','--grad');          set%lgradient = .true.
         case(     '--hess');           set%lhessian = .true.
         case(     '--orca');            set%lorca = .true.

! ------------------------------------------------------------------------
!  no match => take as file name
! ------------------------------------------------------------------------
         case default
            inquire(file=arg,exist=exist)
            if (exist) then
               if (allocated(set%fname)) call raise('S',  &
               &  "There are multiple files provided. '"//set%fname//  &
               &  "' will be ignored in this run.")
               set%fname = arg
            else
               if (index(arg,'-').eq.1) then
                  call raise('S',"Unfortunately, '"//arg// &
                  &    "' is not supported in this program. Check with --help.")
               else
                  call raise('E',"You don't have a file named '"//arg//"' here")
               endif
            endif

         end select

      else ! getopts?
         inquire(file=arg,exist=exist)
         if (exist) then
            if (allocated(set%fname)) call raise('S',  &
            &  "There are multiple files provided. '"//set%fname//  &
            &  "' will be ignored in this run.")
            set%fname = arg
         else
            call raise('E',"You don't have a file named '"//arg//"' here.")
         endif
      endif ! getopts?
   enddo

   if (.not.allocated(set%fname)) then
      call raise('E','No geometry given, so there is nothing to do.')
   endif

end subroutine read_commandline_arguments
