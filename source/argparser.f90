! This file is part of dftd4.
!
! Copyright (C) 2017-2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> reads in command line arguments and parses them into options type
subroutine read_commandline_arguments(env,set)
   use iso_fortran_env, only : wp => real64, istdout => output_unit
!$ use omp_lib

! ------------------------------------------------------------------------
!  general purpose library
! ------------------------------------------------------------------------
   use mctc_environment

! ------------------------------------------------------------------------
!  Class definition
! ------------------------------------------------------------------------
   use class_set

! ------------------------------------------------------------------------
!  get parameters
! ------------------------------------------------------------------------
   use dftd4

   implicit none
   type(mctc_logger),intent(inout) :: env
   type(options),intent(inout) :: set !< options for calculation
   type(mctc_argparser) :: args

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

   call args%new
   if (args%n.eq.0) then
      call help(istdout)
      call env%error(2,"No command line arguments given")
      return
   endif

   getopts = .true.
   skip = 0
   do iarg = 1, args%n
      if (skip.gt.0) then
         skip = skip-1
         cycle
      endif
      ldum = args%get(iarg,arg)
      if (arg.eq.'--') then
         getopts=.false.
         cycle
      endif
      if (getopts) then
!        write(output_unit,'(i0,'':'',x,a)') iarg,arg ! debugging stuff

         if ((len(arg).gt.2).and. &
         &  (index(arg,'--').eq.0).and.(index(arg,'-').eq.1)) then
            call env%warning(2,"the use of '"//arg//"' is discouraged, "// &
            &              "please use '-"//arg//"' next time")
            arg = '-'//arg
         endif
         select case(arg)
! ------------------------------------------------------------------------
!  check the help, version and citation flag, exit if found
! ------------------------------------------------------------------------
         case('-h','--help')
            call help(istdout)
            call terminate(0)
         case(     '--citation')
            call dftd4_citation(istdout)
            call terminate(0)
         case(     '--license')
            call gpl_license(istdout)
            call terminate(0)
         case(     '--version')
            call dftd4_header(istdout,.true.)
            call terminate(0)

! ------------------------------------------------------------------------
!  now for true options, like verbosity or printouts
! ------------------------------------------------------------------------
         case('-v','--verbose')
            set%print_level = 2
         case('-V','--very-verbose')
            set%print_level = 3
         case('-s','--silent')
            set%print_level = 0
      !  OMP option
      !$ case('-P','--parallel')
      !$    skip = 1
      !$    if (args%get(iarg+1,idum)) then
      !$       nproc = omp_get_num_threads()
      !$       if (idum.gt.nproc) &
      !$       & call env%warning(2,'Process number higher than OMP_NUM_THREADS, '//&
      !$       &                'I hope you know what you are doing.')
      !$       call omp_set_num_threads(idum)
      !$    endif
         case('-c','--chrg')
            skip = 1
            if (args%get(iarg+1,idum)) then
               set%chrg = idum
            else
               call env%error(2,"Could not read charge from '"//sec//"'")
               return
            endif
            set%inchrg = .true.
         case('-f','--func')
            skip = 1
            if (args%get(iarg+1,sec)) then
               set%func = sec
            else
               call env%warning(2,"Could find functional name after '--func'")
            endif
         case(     '--molc6')
            set%lmolpol = .true.
         case(     '--param'); skip = 4
            if (args%get(iarg+1,ddum)) set%dparam%s6 = ddum
            if (args%get(iarg+2,ddum)) set%dparam%s8 = ddum
            if (args%get(iarg+3,ddum)) set%dparam%a1 = ddum
            if (args%get(iarg+4,ddum)) set%dparam%a2 = ddum
            ! sanity check
            if (set%dparam%s6.lt.0.0_wp.or.set%dparam%s6.gt.1.0_wp) then
               call env%warning(2,'unphysical s6 chosen, use on your own risk')
            endif
            if (set%dparam%s8.lt.0.0_wp.or.set%dparam%s8.gt.3.0_wp) then
               call env%warning(2,'unphysical s8 chosen, use on your own risk')
            endif
            if (set%dparam%a1.lt.0.0_wp.or.set%dparam%a1.gt.1.0_wp) then
               call env%warning(2,'unphysical a1 chosen, use on your own risk')
            endif
            if (set%dparam%a2.lt.0.0_wp.or.set%dparam%a2.gt.7.0_wp) then
               call env%warning(2,'unphysical a2 chosen, use on your own risk')
            endif
            set%inparam = .true.
         case(     '--s10'); skip = 1
            if (args%get(iarg+1,ddum)) then
               set%dparam%s10 = ddum
            else
               call env%warning(2,"could not read value of s10 from commandline")
            endif
         case(     '--zeta'); skip = 2
            if (args%get(iarg+1,ddum)) set%g_a = ddum
            if (args%get(iarg+2,ddum)) set%g_c = ddum
         case(     '--mbdscale'); skip = 1
            if (args%get(iarg+1,ddum)) set%dparam%s9 = ddum
         case(     '--wfactor'); skip = 1
            if (args%get(iarg+1,ddum)) set%wf = ddum
         case(     '--tmer'); set%ltmer = .true.
         case('-m','--mbd');   set%lmbd = p_mbd_rpalike
         case('-2','--nomb');  set%lmbd = p_mbd_none
         case('-3','--abc');   set%lmbd = p_mbd_approx_atm
         case('-g','--grad'); set%lgradient = .true.
         case(     '--hess');  set%lhessian = .true.
         case(     '--orca');   set%lorca = .true.
         case(     '--json'); set%json = .true.
         case(     '--toml'); set%toml = .true.

! ------------------------------------------------------------------------
!  no match => take as file name
! ------------------------------------------------------------------------
         case default
            inquire(file=arg,exist=exist)
            if (exist) then
               if (allocated(set%fname)) call env%warning(2,  &
               &  "There are multiple files provided. '"//set%fname//  &
               &  "' will be ignored in this run.")
               set%fname = arg
            else
               if (index(arg,'-').eq.1) then
                  call env%warning(2,"Unfortunately, '"//arg// &
                  &    "' is not supported in this program. Check with --help.")
               else
                  call env%error(2,"You don't have a file named '"//arg//"' here")
                  return
               endif
            endif

         end select

      else ! getopts?
         inquire(file=arg,exist=exist)
         if (exist) then
            if (allocated(set%fname)) call env%warning(2,  &
            &  "There are multiple files provided. '"//set%fname//  &
            &  "' will be ignored in this run.")
            set%fname = arg
         else
            call env%error(2,"You don't have a file named '"//arg//"' here.")
            return
         endif
      endif ! getopts?
   enddo

   if (.not.allocated(set%fname)) then
      call env%error(2,'No geometry given, so there is nothing to do.')
   endif

end subroutine read_commandline_arguments
