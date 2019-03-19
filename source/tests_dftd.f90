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
   character(len=:),allocatable :: arg

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
   if (nargs.eq.0) then
      call raise('E',"Please give the tester a test to run!")
   endif

   call rdarg(1,arg)
   if (get_value(arg,idum)) then
      testid = idum
   else
      call raise('E',"Please give the tester a valid test number!")
   endif


   select case(testid)
   case default; stop 77 ! test not available
   case(1);      call test_1
   case(2);      call test_2
   case(3);      call test_3
   case(4);      call test_4
   end select

   ! falling through the tester is always an error
   call terminate(1)

contains

subroutine test_1
   use coordination_number
   use eeq_model
   implicit none
   type(molecule)       :: mol
   type(chrg_parameter) :: chrgeq

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [ 0.00000000000000_wp,  0.00000000000000_wp, -0.73578586109551_wp, &
      & 1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp, &
      &-1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp  &
      & ],shape(xyz))

   real(wp) :: es
   real(wp),allocatable :: cn(:),dcndr(:,:,:)
   real(wp),allocatable :: q(:),dqdr(:,:,:),ges(:,:)

   allocate( cn(nat),dcndr(3,nat,nat),q(nat),dqdr(3,nat,nat+1),ges(3,nat) )
   es  = 0.0_wp
   ges = 0.0_wp

   call mol%allocate(nat,.false.)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp
   
   call dncoord_erf(mol,cn,dcndr)
   ! test for correct CN and correct symmetry in CN
   call assert_close(cn(1),1.9895544818182_wp,thr)
   call assert_close(cn(2),cn(3),             thr)
   call assert_close(cn(1),cn(2)+cn(3),       thr)

   ! test CN derivative, correct summation of diagonals
   call assert_close(dcndr(3,1,1), 0.80973220519637E-01_wp,thr)
   call assert_close(dcndr(1,2,1), 0.52891177763104E-01_wp,thr)
   call assert_close(dcndr(3,1,3),-0.40486610259818E-01_wp,thr)

   call new_charge_model_2019(chrgeq,mol)
   call eeq_chrgeq(chrgeq,mol,cn,dcndr,q,dqdr,es,ges,.false.,.true.,.true.)

   ! test electrostatic energy
   call assert_close(es,-0.64308088326667E-01_wp,thr)

   ! test molecular gradient of ES, also check symmetry
   call assert_close(ges(3,1),-0.44053031985285E-01_wp,thr)
   call assert_close(ges(3,2),ges(3,3),                thr)
   call assert_close(ges(1,2), 0.18102071036002E-01_wp,thr)

   ! test for charge constraint
   call assert_close(sum(q),0.0_wp,            thr)
   call assert_close(q(1),-0.59582744708873_wp,thr)

   ! test derivative of partial charges
   call assert_close(dqdr(3,1,1),-0.41466763854730E+00_wp,thr)
   call assert_close(dqdr(1,2,1), 0.17196993180581E+00_wp,thr)
   call assert_close(dqdr(1,3,2), 0.18631992989121E-01_wp,thr)
   call assert_close(dqdr(2,1,3), 0.00000000000000E+00_wp,thr)

   call mol%deallocate

   ! done: everythings fine
   call terminate(0)
end subroutine test_1

subroutine test_2
   use dftd4
   implicit none
   type(molecule)       :: mol
   integer              :: ndim
   real(wp) :: molpol,molc6,molc8        ! molecular Polarizibility
   real(wp),allocatable :: gweights(:)   ! gaussian weights
   real(wp),allocatable :: refc6(:,:)    ! reference C6 coeffients
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: aw(:,:)

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [ 0.00000000000000_wp, 0.00000000000000_wp,-0.73578586109551_wp, &
      & 1.44183152868459_wp, 0.00000000000000_wp, 0.36789293054775_wp, &
      &-1.44183152868459_wp, 0.00000000000000_wp, 0.36789293054775_wp  &
      & ],shape(xyz))
   real(wp),parameter :: covcn(nat) = &
      [ 1.6105486019977_wp,  0.80527430099886_wp, 0.80527430099886_wp]
   real(wp),parameter :: q(nat) = &
      [-0.59582744708873_wp, 0.29791372354436_wp, 0.29791372354436_wp]
   real(wp),parameter :: g_a = 3.0_wp
   real(wp),parameter :: g_c = 2.0_wp
   real(wp),parameter :: wf  = 6.0_wp
   integer, parameter :: lmbd = p_mbd_approx_atm
   integer, parameter :: refqmode = p_refq_goedecker

   call mol%allocate(nat,.false.)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp

   call d4init(mol,g_a,g_c,refqmode,ndim)

   call assert_eq(ndim,8)

   allocate( gweights(ndim),refc6(ndim,ndim),&
             c6ab(mol%nat,mol%nat),aw(23,mol%nat) )

   call d4(mol,ndim,wf,g_a,g_c,covcn,gweights,refc6)
   call mdisp(mol,ndim,q,g_a,g_c,gweights,refc6,molc6,molc8,molpol,aw,c6ab)

   call assert_close(molpol,9.4271529107854_wp,thr)
   call assert_close(molc6, 44.521545727311_wp,thr)
   call assert_close(molc8, 798.69639423703_wp,thr)

   call assert_close(aw(1,1),6.7482856960122_wp,thr)
   call assert_close(aw(4,2),1.1637689328906_wp,thr)
   call assert_close(aw(7,2),aw(7,3),           thr)

   call assert_close(c6ab(1,2),c6ab(2,1),         thr)
   call assert_close(c6ab(1,1),24.900853125836_wp,thr)
   call assert_close(c6ab(1,3),4.1779697881925_wp,thr)
   call assert_close(c6ab(2,2),c6ab(2,3),         thr)

   call assert_close(sum(gweights),3.0_wp,               thr)
   call assert_close(gweights(2),0.18388886232894E-01_wp,thr)
   call assert_close(gweights(7),0.21400765336381E-01_wp,thr)

   call assert_close(refc6(5,1),10.282422024500_wp,thr)
   call assert_close(refc6(8,6),3.0374149102818_wp,thr)

   call mol%deallocate

   ! done: everythings fine
   call terminate(0)
end subroutine test_2

subroutine test_3
   use dftd4
   implicit none
   type(molecule)       :: mol
   integer  :: idum
   real(wp) :: energy

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [ 0.00000000000000_wp, 0.00000000000000_wp,-0.73578586109551_wp, &
      & 1.44183152868459_wp, 0.00000000000000_wp, 0.36789293054775_wp, &
      &-1.44183152868459_wp, 0.00000000000000_wp, 0.36789293054775_wp  &
      & ],shape(xyz))
   real(wp),parameter :: covcn(nat) = &
      [ 1.6105486019977_wp,  0.80527430099886_wp, 0.80527430099886_wp]
   real(wp),parameter :: q(nat) = &
      [-0.59582744708873_wp, 0.29791372354436_wp, 0.29791372354436_wp]
   real(wp),parameter :: g_a = 3.0_wp
   real(wp),parameter :: g_c = 2.0_wp
   real(wp),parameter :: wf  = 6.0_wp
   integer, parameter :: lmbd = p_mbd_approx_atm
   integer, parameter :: refqmode = p_refq_goedecker
   integer, parameter :: ndim = 8
   real(wp),parameter :: gweights(ndim) = &
      [ 0.15526686926080E-06_wp, 0.18388886232894E-01_wp, 0.89143504504233_wp, &
      & 0.90175913457907E-01_wp, 0.21400765336381E-01_wp, 0.97859923466362_wp, &
      & 0.21400765336381E-01_wp, 0.97859923466362_wp ]
   real(wp),parameter :: refc6(ndim,ndim) = reshape(&
      [ 0.0000000000000_wp,      0.0000000000000_wp,      0.0000000000000_wp,  &
      & 0.0000000000000_wp,      10.282422024500_wp,      6.7431228212696_wp,  &
      & 10.282422024500_wp,      6.7431228212696_wp,      0.0000000000000_wp,  &
      & 0.0000000000000_wp,      0.0000000000000_wp,      0.0000000000000_wp,  &
      & 12.052429296454_wp,      7.8894703511335_wp,      12.052429296454_wp,  &
      & 7.8894703511335_wp,      0.0000000000000_wp,      0.0000000000000_wp,  &
      & 0.0000000000000_wp,      0.0000000000000_wp,      13.246161891965_wp,  &
      & 8.6635841400632_wp,      13.246161891965_wp,      8.6635841400632_wp,  &
      & 0.0000000000000_wp,      0.0000000000000_wp,      0.0000000000000_wp,  &
      & 0.0000000000000_wp,      10.100325850238_wp,      6.6163452797181_wp,  &
      & 10.100325850238_wp,      6.6163452797181_wp,      10.282422024500_wp,  &
      & 12.052429296454_wp,      13.246161891965_wp,      10.100325850238_wp,  &
      & 0.0000000000000_wp,      0.0000000000000_wp,      7.6362416262742_wp,  &
      & 4.7593057612608_wp,      6.7431228212696_wp,      7.8894703511335_wp,  &
      & 8.6635841400632_wp,      6.6163452797181_wp,      0.0000000000000_wp,  &
      & 0.0000000000000_wp,      4.7593057612608_wp,      3.0374149102818_wp,  &
      & 10.282422024500_wp,      12.052429296454_wp,      13.246161891965_wp,  &
      & 10.100325850238_wp,      7.6362416262742_wp,      4.7593057612608_wp,  &
      & 0.0000000000000_wp,      0.0000000000000_wp,      6.7431228212696_wp,  &
      & 7.8894703511335_wp,      8.6635841400632_wp,      6.6163452797181_wp,  &
      & 4.7593057612608_wp,      3.0374149102818_wp,      0.0000000000000_wp,  &
      & 0.0000000000000_wp],     shape(refc6))
   type(dftd_parameter),parameter :: dparam_pwpb95 = dftd_parameter ( &
      &  s6=0.8200_wp, s8=-0.34639127_wp, a1=0.41080636_wp, a2=3.83878274_wp )
   type(dftd_parameter),parameter :: dparam_pbe    = dftd_parameter ( &
      &  s6=1.0000_wp, s8=0.95948085_wp, a1=0.38574991_wp, a2=4.80688534_wp )
   type(dftd_parameter),parameter :: dparam_random = dftd_parameter ( &
      &  s6=0.95_wp, s8=0.45_wp, s10=0.65_wp, s9=1.10_wp, a1=0.43_wp, a2=5.10_wp )

   call mol%allocate(nat,.false.)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp

   call d4init(mol,g_a,g_c,refqmode,idum)

   call assert_eq(idum,ndim)

   energy = +1.0_wp ! energy is intent(out)

   call edisp(mol,ndim,q,dparam_pwpb95,g_a,g_c,gweights,refc6,lmbd,energy)
   call assert_close(energy,-0.22526819184723E-03_wp,thr)

   call edisp(mol,ndim,q,dparam_pbe,g_a,g_c,gweights,refc6,lmbd,energy)
   call assert_close(energy,-0.19788865790096E-03_wp,thr)
   !-0.19558245089408E-03

   call edisp(mol,ndim,q,dparam_random,g_a,g_c,gweights,refc6,lmbd,energy)
   call assert_close(energy,-0.11213581758666E-03_wp,thr)

   call mol%deallocate

   ! done: everythings fine
   call terminate(0)
end subroutine test_3

subroutine test_4
   use dftd4
   implicit none
   type(molecule)       :: mol

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [ 0.00000000000000_wp,  0.00000000000000_wp, -0.73578586109551_wp, &
      & 1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp, &
      &-1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp  &
      & ],shape(xyz))
   type(dftd_parameter),parameter :: dparam_b2plyp = dftd_parameter ( &
      &  s6=0.6400_wp, s8=1.16888646_wp, a1=0.44154604_wp, a2=4.73114642_wp )
   type(dftd_parameter),parameter :: dparam_tpss   = dftd_parameter ( &
      &  s6=1.0000_wp, s8=1.76596355_wp, a1=0.42822303_wp, a2=4.54257102_wp )
   type(dftd_options),  parameter :: opt_1 = dftd_options ( &
      &  lmbd = p_mbd_approx_atm, refq = p_refq_goedecker, &
      &  wf = 6.0_wp, g_a = 3.0_wp, g_c = 2.0_wp, &
      &  lmolpol=.false., lenergy=.true., lgradient=.false., lhessian=.true., &
      &  verbose = .false., veryverbose = .false., silent = .false. )
   type(dftd_options),  parameter :: opt_2 = dftd_options ( &
      &  lmbd = p_mbd_approx_atm, refq = p_refq_goedecker, &
      &  wf = 6.0_wp, g_a = 3.0_wp, g_c = 2.0_wp, &
      &  lmolpol=.false., lenergy=.false., lgradient=.true., lhessian=.false., &
      &  verbose = .false., veryverbose = .false., silent = .true. )

   real(wp) :: energy
   real(wp),allocatable :: gradient(:,:),hessian(:,:)

   allocate( gradient(3,nat), hessian(3*nat,3*nat) )
   energy   = 0.0_wp
   gradient = 0.0_wp
   hessian  = 0.0_wp

   call mol%allocate(nat,.false.)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp
   
   call d4_calculation(istdout,opt_1,mol,dparam_tpss,energy,gradient,hessian)
   call assert_close(energy,-0.26682682254336E-03_wp,thr)

   call assert_close(hessian(3,1), 7.9334628241320_wp,thr)
   call assert_close(hessian(4,8),-3.2756224894310_wp,thr)
   call assert_close(hessian(5,3), 0.0000000000000_wp,thr)

   call d4_calculation(istdout,opt_2,mol,dparam_b2plyp,energy,gradient,hessian)
   call assert_close(energy,-0.13368190339570E-03_wp,thr)

   call assert_close(gradient(1,1), 0.00000000000000E+00_wp,thr)
   call assert_close(gradient(3,1), 0.39778648945254E-04_wp,thr)
   call assert_close(gradient(3,2),-0.19889324472627E-04_wp,thr)
   call assert_close(gradient(1,2),-gradient(1,3),          thr)

   call mol%deallocate

   ! done: everythings fine
   call terminate(0)
end subroutine test_4

subroutine assert_eq(val1,val2)
   use iso_fortran_env, istderr => error_unit
   integer,intent(in) :: val1,val2

   if (val1 /= val2) then
      write(istderr,'("assertion:",1x,g21.14," == ",g21.14,1x,"FAILED")') &
         val1,val2
      call terminate(1)
   endif
end subroutine assert_eq

subroutine assert_close(val1,val2,thr)
   use iso_fortran_env, istderr => error_unit
   real(wp),intent(in) :: val1,val2,thr
   real(wp) :: diff

   diff = val1 - val2
   if (abs(diff) > thr) then
      write(istderr,'("assertion:",1x,g21.14," == ",g21.14,1x,"FAILED")') &
         val1,val2
      call terminate(1)
   endif
end subroutine assert_close

end program dftd_tester
