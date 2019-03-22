program dftd_tester
   use iso_fortran_env, istdout => output_unit, kdp => real64
!$ use omp_lib
! ------------------------------------------------------------------------
!  general purpose library
! ------------------------------------------------------------------------
   use mctc_global
   use mctc_timings
   use mctc_systools
   use mctc_econv
   use mctc_readin

! ------------------------------------------------------------------------
!  class definitions
! ------------------------------------------------------------------------
   use class_molecule
   use class_param
   use class_set

   implicit none

! ------------------------------------------------------------------------
!  local variables
! ------------------------------------------------------------------------
   integer  :: testid
   integer  :: idum,nargs
   character(len=:),allocatable :: arg,sec

!! ------------------------------------------------------------------------
!!  signal processing
!! ------------------------------------------------------------------------
!   external :: wSIGTERM
!   external :: wSIGINT
!!  here two important signal handlers are installed, it seems that
!!  FORTRAN by itself cannot handle signals in the way I expected it
!!  to do, but this will force it to die when I hit CTRL^C.
!   call signal(2,wSIGINT)
!   call signal(15,wSIGTERM)

! ------------------------------------------------------------------------
!  general setup
! ------------------------------------------------------------------------
!  initialize the timing system
   call start_timing_run
   call init_timing(10,verb=.true.) ! verbosity allows printing of cputime
   call start_timing(1)

!  initialize the messagebuffer for the error handler
   call init_errorbuffer

!  set this for mctc_global
   name = 'test'

   nargs = command_argument_count()
   if (nargs.lt.2) then
      call raise('E',"Please give the tester a test to run!")
   endif

   call rdarg(1,arg)
   call rdarg(2,sec)

   select case(arg)
   case('eeq_model')
      select case(sec)
      case('water'); call test_eeq_model_water
      case('ewald'); call test_eeq_model_ewald
      end select
   case('dftd4')
      select case(sec)
      case('properties'); call test_dftd4_properties
      case('energies');   call test_dftd4_energies
      case('api');        call test_dftd4_api
      case('periodic');   call test_dftd4_pbc
      case('pbc_disp');   call test_dftd4_pbc_energies
      end select
   case('geometry_reader')
      select case(sec)
      case('coord_3d');  call test_geometry_reader_file_coord_CaMgCO_3d
      case('coord_2d');  call test_geometry_reader_file_coord_C_2d
      case('coord_1d');  call test_geometry_reader_file_coord_C_1d
      case('coord_0d');  call test_geometry_reader_file_coord_general_0d
      case('xmol_0d');   call test_geometry_reader_file_xmol_water_0d
      case('poscar_3d'); call test_geometry_reader_file_poscar_sio2_3d
      end select
   case('pbc_tools')
      select case(sec)
      case('convert'); call test_pbc_tools_convert
      case('cutoff');  call test_pbc_tools_cutoff
      end select
   case('class_molecule')
      select case(sec)
      case('mic');  call test_class_molecule_mic_distances
      case('axis'); call test_class_molecule_axis_trafo
      end select
   end select

   ! falling through the tester is always an error
   call terminate(1)

end program dftd_tester
