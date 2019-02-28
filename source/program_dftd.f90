program dftd
   use iso_fortran_env, istdout => output_unit, kdp => real64
!$ use omp_lib
!! ------------------------------------------------------------------------
!  general purpose library
!! ------------------------------------------------------------------------
   use mctc_global
   use mctc_timings
   use mctc_econv

!! ------------------------------------------------------------------------
!  class definitions
!! ------------------------------------------------------------------------
   use class_molecule
   use class_set

!! ------------------------------------------------------------------------
!  interfaces
!! ------------------------------------------------------------------------
   use coordination_number
   use eeq_model
   use dftd4
   use dfuncpar

   implicit none

!! ------------------------------------------------------------------------
!  class declarations
!! ------------------------------------------------------------------------
   type(molecule)       :: mol
   type(options)        :: set
   type(dftd_parameter) :: dparam
   type(chrg_parameter) :: chrgeq

!! ------------------------------------------------------------------------
!  local variables
!! ------------------------------------------------------------------------
   integer  :: ndim                      ! matrix dimension
   integer  :: i,j,k,l,ii,jj
   integer  :: err
   real(wp) :: memory
   real(wp) :: energy                    ! dispersion energy
   real(wp) :: etmp,etwo,emany,er,el,es
   real(wp) :: molpol,molc6,molc8        ! molecular Polarizibility
   real(wp),allocatable :: gradient(:,:) ! nuclear gradient
   real(wp),allocatable :: hessian(:,:)  ! nuclear hessian
   real(wp),allocatable :: q(:)          ! partial charges
   real(wp),allocatable :: dqdr(:,:,:)   ! partial charges
   real(wp),allocatable :: covcn(:)      ! covalent coordination number
   real(wp),allocatable :: dcovcndr(:,:,:) ! covalent coordination number derivative
   real(wp),allocatable :: cn(:)         ! erf coordination number
   real(wp),allocatable :: dcndr(:,:,:)  ! erf coordination number derivative
   real(wp),allocatable :: gweights(:)   ! gaussian weights
   real(wp),allocatable :: refc6(:,:)    ! reference C6 coeffients
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: aw(:,:)
   real(wp),allocatable :: ges(:,:)
   real(wp),allocatable :: gr(:,:)
   real(wp),allocatable :: gl(:,:)
   real(wp),parameter   :: step = 1.0e-5_wp, step2 = 0.5_wp/step

!!! ------------------------------------------------------------------------
!!  signal processing
!!! ------------------------------------------------------------------------
!   external :: wSIGTERM
!   external :: wSIGINT
!!  here two important signal handlers are installed, it seems that
!!  FORTRAN by itself cannot handle signals in the way I expected it
!!  to do, but this will force it to die when I hit CTRL^C.
!   call signal(2,wSIGINT)
!   call signal(15,wSIGTERM)

!! ------------------------------------------------------------------------
!  general setup
!! ------------------------------------------------------------------------
!  initialize the timing system
   call start_timing_run
   call init_timing(10,verb=.true.) ! verbosity allows printing of cputime
   call start_timing(1)

!  initialize the messagebuffer for the error handler
   call init_errorbuffer

!  set this for mctc_global
   name = 'dftd4'

!! ------------------------------------------------------------------------
!  command line arguments
!! ------------------------------------------------------------------------
   call read_commandline_arguments(set)

!! ------------------------------------------------------------------------
!  get molecular geometry
!! ------------------------------------------------------------------------
   call get_geometry(mol,set%fname)
   if (set%inchrg) mol%chrg = set%chrg

!! ------------------------------------------------------------------------
!  Output:
!  Header, Citation and Licence
!! ------------------------------------------------------------------------
   call dftd4_header(set%verbose)
   if (.not.set%silent) then
   call dftd4_citation
   call prdate('S')
   write(istdout,'(a)')
   endif

!! ------------------------------------------------------------------------
!  Set defaults
!! ------------------------------------------------------------------------
   if (set%inparam) then
      dparam = set%dparam
      if (.not.(set%lgradient.or.set%lhessian)) set%lenergy = .true.
   else if (allocated(set%func)) then
      call d4par(set%func,dparam,set%lmbd)
      if (.not.(set%lgradient.or.set%lhessian)) set%lenergy = .true.
   else
      if (set%lenergy) &
         call raise('E','Dispersion energy requested but no parameters given')
      if (set%lgradient) &
         call raise('E','Dispersion gradient requested but no parameters given')
      if (set%lhessian) &
         call raise('E','Dispersion Hessian requested but no parameters given')
   endif
   if (.not.(set%lenergy.or.set%lgradient.or.set%lhessian)) set%lmolpol = .true.

!! ------------------------------------------------------------------------
!  Output: Initialization and Parameter setup
!! ------------------------------------------------------------------------
   call d4init(mol,set%g_a,set%g_c,p_refq_goedecker,ndim)
   memory = ((mol%nat)+(mol%nat)+(ndim)+(ndim*ndim) &
            +(3*mol%nat)+(mol%nat*mol%nat)+(23*mol%nat)&
            +(3*mol%nat*3*mol%nat)+(mol%nat)+(3*mol%nat) &
            +(3*mol%nat*mol%nat)+(3*mol%nat*(mol%nat+1)) &
            +(3*mol%nat*mol%nat)) * wp / (1024.0_wp**2)

   call generic_header(istdout,'Calculation Setup',49,10)
   write(istdout,'(3x,a," : ",a)')    'coordinate file     ', set%fname
   write(istdout,'(3x,a," : ",i6)')   'number of atoms     ', mol%nat
   write(istdout,'(3x,a," : ",i6)')   'charge              ', nint(mol%chrg)
   write(istdout,'(3x,a," : ",a)')    'non-additivity corr.', lmbd2string(set%lmbd)
   write(istdout,'(3x,a," : ",a)')    'charge model        ', refq2string(p_refq_goedecker)
   if (allocated(set%func)) &
   write(istdout,'(3x,a," : ",a)')    'functional          ', set%func
!$ write(istdout,'(3x,a," : ",i6)')   'omp threads         ',omp_get_num_threads()
   write(istdout,'(3x,a," : ",f6.1,1x,a)') &
      "memory needed (est.)",memory,"Mb"
   if (set%wf  /= 6.0_wp .or. set%verbose) &
   write(istdout,'(3x,a," : ",f6.1)')   'weighting factor    ', set%wf
   if (set%g_a /= 3.0_wp .or. set%verbose) &
   write(istdout,'(3x,a," : ",f6.1)')   'q-scale height      ', set%g_a
   if (set%g_c /= 2.0_wp .or. set%verbose) &
   write(istdout,'(3x,a," : ",f6.1)')   'q-scale steepness   ', set%g_c

   write(istdout,'(a)')

   if (.not.set%silent) &
   call prd4ref(mol)

   allocate( q(mol%nat),covcn(mol%nat),gweights(ndim),refc6(ndim,ndim),&
             gradient(3,mol%nat),c6ab(mol%nat,mol%nat),aw(23,mol%nat),&
             hessian(3*mol%nat,3*mol%nat),cn(mol%nat),ges(3,mol%nat), &
             dcndr(3,mol%nat,mol%nat),dqdr(3,mol%nat,mol%nat+1), &
             dcovcndr(3,mol%nat,mol%nat), stat = err )
   if (err /= 0) &
      call raise('E','Memory allocation failed')

   call dncoord_d4(mol,covcn,dcovcndr)
   call d4(mol,ndim,set%wf,set%g_a,set%g_c,covcn,gweights,refc6)

   call dncoord_erf(mol,cn,dcndr)

!! ------------------------------------------------------------------------
!  get partial charges
!! ------------------------------------------------------------------------
   if (set%verbose) &
   call eeq_header
   call new_charge_model(chrgeq,mol)
   call eeq_chrgeq(chrgeq,mol,cn,dcndr,q,dqdr,es,ges, &
                   .false.,.false.,.true.)
   if (set%verbose) &
   call print_chrgeq(istdout,chrgeq,mol,q,cn)

!! ------------------------------------------------------------------------
!  calculate properties
!! ------------------------------------------------------------------------
dispersion_properties: if (set%lmolpol) then
   call generic_header(istdout,'Molecular Properties',49,10)
   call mdisp(mol,ndim,q,set%g_a,set%g_c,gweights,refc6,molc6,molc8,molpol,aw,c6ab)
   call prmolc6(mol,molc6,molc8,molpol,covcn=covcn,q=q,c6ab=c6ab,alpha=aw(1,:))
endif dispersion_properties

if (.not.set%silent.and.(set%lenergy.or.set%lgradient.or.set%lhessian)) then
   call generic_header(istdout,'Damping Parameters',49,10)
   write(istdout,'(3x,a," : ",f10.4)')   's6                  ', dparam%s6
   write(istdout,'(3x,a," : ",f10.4)')   's8                  ', dparam%s8
   if (dparam%s10 /= 0.0_wp .or. set%verbose) &
   write(istdout,'(3x,a," : ",f10.4)')   's10                 ', dparam%s10
   if (dparam%s9  /= 1.0_wp .or. set%verbose) &
   write(istdout,'(3x,a," : ",f10.4)')   's9                  ', dparam%s9
   write(istdout,'(3x,a," : ",f10.4)')   'a1                  ', dparam%a1
   write(istdout,'(3x,a," : ",f10.4)')   'a2                  ', dparam%a2
   write(istdout,'(a)')
endif

!! ------------------------------------------------------------------------
!  calculate energy
!! ------------------------------------------------------------------------
dispersion_energy: if (set%lenergy) then
   call edisp(mol,ndim,q,dparam,set%g_a,set%g_c, &
              gweights,refc6,set%lmbd,energy,etwo=etwo,emany=emany)
endif dispersion_energy

!! ------------------------------------------------------------------------
!  calculate gradient
!! ------------------------------------------------------------------------
dispersion_gradient: if (set%lgradient) then
   call dispgrad(mol,ndim,q,dqdr,covcn,dcovcndr,dparam,set%wf,set%g_a,set%g_c, &
   &             refc6,set%lmbd,gradient,etmp)
   if (.not.set%lenergy) energy = etmp
!   allocate(gr(3,mol%nat),gl(3,mol%nat))
!   call generic_header(istdout,'numerical gradient',49,10)
!   do i = 1, mol%nat
!      do j = 1, 3
!         ii = 3*(i-1)+j
!         er = 0.0_wp
!         el = 0.0_wp
!         mol%xyz(j,i) = mol%xyz(j,i) + step
!         call dncoord_d4(mol,covcn,dcovcndr)
!         call dncoord_erf(mol,cn,dcndr)
!         call eeq_chrgeq(chrgeq,mol,cn,dcndr,q,dqdr,es,ges, &
!                               .false.,.false.,.true.)
!         call dispgrad(mol,ndim,q,dqdr,covcn,dcovcndr,dparam, &
!                       set%wf,set%g_a,set%g_c,refc6,set%lmbd,gr,er)
!         mol%xyz(j,i) = mol%xyz(j,i) - 2*step
!         call dncoord_d4(mol,covcn,dcovcndr)
!         call dncoord_erf(mol,cn,dcndr)
!         call eeq_chrgeq(chrgeq,mol,cn,dcndr,q,dqdr,es,ges, &
!                               .false.,.false.,.true.)
!         call dispgrad(mol,ndim,q,dqdr,covcn,dcovcndr,dparam, &
!                       set%wf,set%g_a,set%g_c,refc6,set%lmbd,gr,el)
!         mol%xyz(j,i) = mol%xyz(j,i) + step
!         gl (j,i) = (er-el)*step2
!      enddo
!   enddo

endif dispersion_gradient

!! ------------------------------------------------------------------------
!  calculate hessian
!! ------------------------------------------------------------------------
dispersion_hessian: if (set%lhessian) then
   allocate(gr(3,mol%nat),gl(3,mol%nat))
   do i = 1, mol%nat
      do j = 1, 3
         ii = 3*(i-1)+j
         gr = 0.0_wp
         gl = 0.0_wp
         mol%xyz(j,i) = mol%xyz(j,i) + step
         call dncoord_d4(mol,covcn,dcovcndr)
         call dncoord_erf(mol,cn,dcndr)
         call eeq_chrgeq(chrgeq,mol,cn,dcndr,q,dqdr,es,ges, &
                         .false.,.false.,.true.)
         call dispgrad(mol,ndim,q,dqdr,covcn,dcovcndr,dparam, &
                       set%wf,set%g_a,set%g_c,refc6,set%lmbd,gr,er)
         mol%xyz(j,i) = mol%xyz(j,i) - 2*step
         call dncoord_d4(mol,covcn,dcovcndr)
         call dncoord_erf(mol,cn,dcndr)
         call eeq_chrgeq(chrgeq,mol,cn,dcndr,q,dqdr,es,ges, &
                         .false.,.false.,.true.)
         call dispgrad(mol,ndim,q,dqdr,covcn,dcovcndr,dparam, &
                       set%wf,set%g_a,set%g_c,refc6,set%lmbd,gr,el)
         mol%xyz(j,i) = mol%xyz(j,i) + step
         do k = 1, mol%nat
            do l = 1, 3
               jj = 3*(k-1) + l
               hessian(jj,ii) = (gr(l,k)-gl(j,i))*step2
            enddo
         enddo
      enddo
   enddo
   deallocate(gr,gl)
endif dispersion_hessian

!! ------------------------------------------------------------------------
!  Output:
!! ------------------------------------------------------------------------
if (set%lenergy.or.set%lgradient.or.set%lhessian) &
   call generic_header(istdout,'Results',49,10)
   if (set%lenergy) then
      write(istdout,'('// &
      &      '1x,"Edisp  /kcal,au:",f11.4,1x,f12.8)') &
      &       energy*autokcal,energy
      if(set%verbose) &
      write(istdout,'('// &
      &      '1x,"E(2)   /kcal,au:",f11.4,1x,f12.8,'// &
      &      '/,1x,"Emany  /kcal,au:",f11.4,1x,f12.8)') &
      &       Etwo*autokcal,Etwo, &
      &       Emany*autokcal,Emany
      write(istdout,'(a)')
   endif
   if (set%ltmer.or.(.not.set%lorca.and.set%lenergy)) &
      call out_tmer('.EDISP',energy)

   if (set%lgradient) then
      if (set%verbose) then
         write(istdout,'(1x,a)') &
            "Dispersion gradient"
         call write_gradient(mol,istdout,gradient)
         write(istdout,'(a)')
      endif
      write(istdout,'(1x,"E(opt) /kcal,au:",f11.4,1x,f12.8)') &
      &       Etmp*autokcal,Etmp
      write(istdout,'(1x,"|G| =",1x,f34.10)') &
      &       norm2(gradient)
      write(istdout,'(a)')

      if (set%lorca) then
         call orca_gradient(mol,istdout,gradient)
      else
         call out_gradient(mol,'gradient',energy,gradient)
      endif
   endif

   if (set%lhessian) then
      call orca_hessian(mol,istdout,hessian)
   endif

   call raise('F','Some non-fatal runtime exceptions occurred, please check:')

!! ------------------------------------------------------------------------
!  Print timings
!! ------------------------------------------------------------------------
   if (.not.set%silent) then
   write(istdout,'(a)')
   call stop_timing_run
   call stop_timing(1)
   call prdate('E')
   call prtiming(1,'total')
   endif

   write(istdout,'(a)')
   call terminate(0)

end program dftd
