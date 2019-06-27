!> interface to the DFT-D4 module
subroutine d4_calculation(iunit,opt,mol,dparam,energy,gradient,hessian)
   use iso_fortran_env, wp => real64
!$ use omp_lib

! ------------------------------------------------------------------------
!  class definitions
! ------------------------------------------------------------------------
   use class_set
   use class_param
   use class_molecule

! ------------------------------------------------------------------------
!  interfaces
! ------------------------------------------------------------------------
   use coordination_number
   use eeq_model
   use dftd4
   use dfuncpar


   implicit none

   !> output unit, usually bound to STDOUT
   !  not used if opt%silent is set (maybe for fatal errors)
   integer, intent(in) :: iunit

! ------------------------------------------------------------------------
!  class declarations
! ------------------------------------------------------------------------
   !> molecule structure information
   type(molecule), intent(inout) :: mol
   !> calculation options for DFT-D4 method
   type(dftd_options), intent(in) :: opt
   !> damping parameter for the DFT-D4 method
   type(dftd_parameter),intent(in) :: dparam
   type(chrg_parameter)            :: chrgeq

! ------------------------------------------------------------------------
!  output variables
! ------------------------------------------------------------------------
   !> dispersion energy, always referenced
   real(wp),intent(out) :: energy
   !> nuclear gradient, only references if opt%lgradient is set
   real(wp),intent(out) :: gradient(3,mol%nat)
   !> nuclear hessian, only references if opt%lhessian is set
   real(wp),intent(out) :: hessian(3*mol%nat,3*mol%nat)

! ------------------------------------------------------------------------
!  local variables
! ------------------------------------------------------------------------
   integer  :: ndim                      ! matrix dimension
   integer  :: i,j,k,l,ii,jj
   integer  :: err
   real(wp) :: memory
   real(wp) :: etmp,etwo,emany,er,el,es
   real(wp) :: molpol,molc6,molc8        ! molecular Polarizibility
   real(wp) :: sigma(3,3)
   real(wp),allocatable :: q(:)          ! partial charges
   real(wp),allocatable :: dqdr(:,:,:)   ! partial charges
   real(wp),allocatable :: dqdL(:,:,:)   ! partial charges
   real(wp),allocatable :: covcn(:)      ! covalent coordination number
   real(wp),allocatable :: dcovcndr(:,:,:) ! covalent coordination number derivative
   real(wp),allocatable :: dcovcndL(:,:,:) ! covalent coordination number derivative
   real(wp),allocatable :: cn(:)         ! erf coordination number
   real(wp),allocatable :: dcndr(:,:,:)  ! erf coordination number derivative
   real(wp),allocatable :: dcndL(:,:,:)  ! erf coordination number derivative
   real(wp),allocatable :: gweights(:)   ! gaussian weights
   real(wp),allocatable :: refc6(:,:)    ! reference C6 coeffients
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: aw(:,:)
   real(wp),allocatable :: ges(:,:)
   real(wp),allocatable :: gr(:,:)
   real(wp),allocatable :: gl(:,:)
   real(wp),parameter   :: step = 1.0e-5_wp, step2 = 0.5_wp/step

! ------------------------------------------------------------------------
!  Output: Initialization and Parameter setup
! ------------------------------------------------------------------------
   call d4init(mol,opt%g_a,opt%g_c,p_refq_goedecker,ndim)
   memory = ((mol%nat)+(mol%nat)+(ndim)+(ndim*ndim) &
            +(mol%nat*mol%nat)+(23*mol%nat)+(mol%nat)+(3*mol%nat) &
            +(3*mol%nat*mol%nat)+(3*mol%nat*(mol%nat+1)) &
            +(3*mol%nat*mol%nat)) * wp / (1024.0_wp**2)

   call generic_header(iunit,'Calculation Setup',49,10)
   write(iunit,'(3x,a," : ",i6)')   'number of atoms     ', mol%nat
   write(iunit,'(3x,a," : ",i6)')   'charge              ', nint(mol%chrg)
   write(iunit,'(3x,a," : ",a)')    'non-additivity corr.', lmbd2string(opt%lmbd)
   write(iunit,'(3x,a," : ",a)')    'charge model        ', refq2string(p_refq_goedecker)
!$ write(iunit,'(3x,a," : ",i6)')   'omp threads         ',omp_get_num_threads()
   write(iunit,'(3x,a," : ",f6.1,1x,a)') &
      "memory needed (est.)",memory,"Mb"
   if (opt%wf  /= 6.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f6.1)')   'weighting factor    ', opt%wf
   if (opt%g_a /= 3.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f6.1)')   'q-scale height      ', opt%g_a
   if (opt%g_c /= 2.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f6.1)')   'q-scale steepness   ', opt%g_c

   write(iunit,'(a)')

   if (.not.opt%silent) &
   call prd4ref(mol)

   allocate( q(mol%nat),covcn(mol%nat),gweights(ndim),refc6(ndim,ndim),&
             c6ab(mol%nat,mol%nat),aw(23,mol%nat),cn(mol%nat), &
             dcndr(3,mol%nat,mol%nat),dqdr(3,mol%nat,mol%nat+1), &
             dcovcndr(3,mol%nat,mol%nat), stat = err )
   if (err /= 0) &
      call raise('E','Memory allocation failed')

   call dncoord_d4(mol,covcn,dcovcndr)
   call d4(mol,ndim,opt%wf,opt%g_a,opt%g_c,covcn,gweights,refc6)

   call dncoord_erf(mol,cn,dcndr)
   call dncoord_logcn(mol%nat,cn,dcndr,cn_max=8.0_wp)

! ------------------------------------------------------------------------
!  get partial charges
! ------------------------------------------------------------------------
   if (opt%verbose) &
   call eeq_header
   call new_charge_model(chrgeq,mol)
   call eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma, &
                   .false.,.false.,.true.)
   if (opt%verbose) &
   call print_chrgeq(iunit,chrgeq,mol,q,cn)

! ------------------------------------------------------------------------
!  calculate properties
! ------------------------------------------------------------------------
dispersion_properties: if (opt%lmolpol) then
   call generic_header(iunit,'Molecular Properties',49,10)
   call mdisp(mol,ndim,q,opt%g_a,opt%g_c,gweights,refc6,molc6,molc8,molpol,aw,c6ab)
   call prmolc6(mol,molc6,molc8,molpol,covcn=covcn,q=q,c6ab=c6ab,alpha=aw(1,:))
endif dispersion_properties

if (.not.opt%silent.and.(opt%lenergy.or.opt%lgradient.or.opt%lhessian)) then
   call generic_header(iunit,'Damping Parameters',49,10)
   write(iunit,'(3x,a," : ",f10.4)')   's6                  ', dparam%s6
   write(iunit,'(3x,a," : ",f10.4)')   's8                  ', dparam%s8
   if (dparam%s10 /= 0.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f10.4)')   's10                 ', dparam%s10
   if (dparam%s9  /= 1.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f10.4)')   's9                  ', dparam%s9
   write(iunit,'(3x,a," : ",f10.4)')   'a1                  ', dparam%a1
   write(iunit,'(3x,a," : ",f10.4)')   'a2                  ', dparam%a2
   write(iunit,'(a)')
endif

! ------------------------------------------------------------------------
!  calculate energy
! ------------------------------------------------------------------------
dispersion_energy: if (opt%lenergy) then
   call edisp(mol,ndim,q,dparam,opt%g_a,opt%g_c, &
              gweights,refc6,opt%lmbd,energy,etwo=etwo,emany=emany)
endif dispersion_energy

! ------------------------------------------------------------------------
!  calculate gradient
! ------------------------------------------------------------------------
dispersion_gradient: if (opt%lgradient) then
   call dispgrad(mol,ndim,q,dqdr,covcn,dcovcndr,dparam,opt%wf,opt%g_a,opt%g_c, &
   &             refc6,opt%lmbd,gradient,etmp)
   if (.not.opt%lenergy) energy = etmp
!   allocate(gr(3,mol%nat),gl(3,mol%nat))
!   call generic_header(iunit,'numerical gradient',49,10)
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
!                       opt%wf,opt%g_a,opt%g_c,refc6,opt%lmbd,gr,er)
!         mol%xyz(j,i) = mol%xyz(j,i) - 2*step
!         call dncoord_d4(mol,covcn,dcovcndr)
!         call dncoord_erf(mol,cn,dcndr)
!         call eeq_chrgeq(chrgeq,mol,cn,dcndr,q,dqdr,es,ges, &
!                               .false.,.false.,.true.)
!         call dispgrad(mol,ndim,q,dqdr,covcn,dcovcndr,dparam, &
!                       opt%wf,opt%g_a,opt%g_c,refc6,opt%lmbd,gr,el)
!         mol%xyz(j,i) = mol%xyz(j,i) + step
!         gl (j,i) = (er-el)*step2
!      enddo
!   enddo

endif dispersion_gradient

! ------------------------------------------------------------------------
!  calculate hessian
! ------------------------------------------------------------------------
dispersion_hessian: if (opt%lhessian) then
   allocate(gr(3,mol%nat),gl(3,mol%nat))
   do i = 1, mol%nat
      do j = 1, 3
         ii = 3*(i-1)+j
         gr = 0.0_wp
         gl = 0.0_wp
         mol%xyz(j,i) = mol%xyz(j,i) + step
         call dncoord_d4(mol,covcn,dcovcndr)
         call dncoord_erf(mol,cn,dcndr)
         call dncoord_logcn(mol%nat,cn,dcndr,cn_max=8.0_wp)
         call eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma, &
                         .false.,.false.,.true.)
         call dispgrad(mol,ndim,q,dqdr,covcn,dcovcndr,dparam, &
                       opt%wf,opt%g_a,opt%g_c,refc6,opt%lmbd,gr,er)
         mol%xyz(j,i) = mol%xyz(j,i) - 2*step
         call dncoord_d4(mol,covcn,dcovcndr)
         call dncoord_erf(mol,cn,dcndr)
         call dncoord_logcn(mol%nat,cn,dcndr,cn_max=8.0_wp)
         call eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma, &
                         .false.,.false.,.true.)
         call dispgrad(mol,ndim,q,dqdr,covcn,dcovcndr,dparam, &
                       opt%wf,opt%g_a,opt%g_c,refc6,opt%lmbd,gr,el)
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

end subroutine d4_calculation

!> interface to the DFT-D4 module for periodic boundary conditions
subroutine d4_pbc_calculation(iunit,opt,mol,dparam,energy,gradient,latgrad)
   use iso_fortran_env, wp => real64
!$ use omp_lib

! ------------------------------------------------------------------------
!  class definitions
! ------------------------------------------------------------------------
   use class_set
   use class_param
   use class_molecule

! ------------------------------------------------------------------------
!  interfaces
! ------------------------------------------------------------------------
   use coordination_number
   use eeq_model
   use dftd4
   use dfuncpar
   use pbc_tools

   implicit none

   !> output unit, usually bound to STDOUT
   !  not used if opt%silent is set (maybe for fatal errors)
   integer, intent(in) :: iunit

! ------------------------------------------------------------------------
!  class declarations
! ------------------------------------------------------------------------
   !> molecule structure information
   type(molecule), intent(in) :: mol
   !> calculation options for DFT-D4 method
   type(dftd_options), intent(in) :: opt
   !> damping parameter for the DFT-D4 method
   type(dftd_parameter),intent(in) :: dparam
   type(chrg_parameter)            :: chrgeq

! ------------------------------------------------------------------------
!  output variables
! ------------------------------------------------------------------------
   !> dispersion energy, always referenced
   real(wp),intent(out) :: energy
   !> nuclear gradient, only references if opt%lgradient is set
   real(wp),intent(out) :: gradient(3,mol%nat)
   !> nuclear hessian, only references if opt%lhessian is set
   real(wp),intent(out) :: latgrad(3,3)

! ------------------------------------------------------------------------
!  local variables
! ------------------------------------------------------------------------
   integer  :: ndim                      ! matrix dimension
   integer  :: i,j,k,l,ii,jj
   integer  :: err
   integer  :: rep_cn(3),rep_vdw(3),rep_mbd(3)
   real(wp) :: memory
   real(wp) :: etmp,etwo,emany,er,el,es
   real(wp) :: molpol,molc6,molc8        ! molecular Polarizibility
   real(wp) :: sigma(3,3),inv_lat(3,3)
   real(wp),allocatable :: q(:)          ! partial charges
   real(wp),allocatable :: dqdr(:,:,:)   ! partial charges
   real(wp),allocatable :: dqdL(:,:,:)   ! partial charges
   real(wp),allocatable :: covcn(:)      ! covalent coordination number
   real(wp),allocatable :: dcovcndr(:,:,:) ! covalent coordination number derivative
   real(wp),allocatable :: dcovcndL(:,:,:) ! covalent coordination number derivative
   real(wp),allocatable :: cn(:)         ! erf coordination number
   real(wp),allocatable :: dcndr(:,:,:)  ! erf coordination number derivative
   real(wp),allocatable :: dcndL(:,:,:)  ! erf coordination number derivative
   real(wp),allocatable :: gweights(:)   ! gaussian weights
   real(wp),allocatable :: refc6(:,:)    ! reference C6 coeffients
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: aw(:,:)
   real(wp),allocatable :: gtmp(:,:)
   real(wp),allocatable :: stmp(:,:)

   real(wp),parameter   :: rthr_cn  =  900.0_wp
   real(wp),parameter   :: rthr_mbd = 1600.0_wp
   real(wp),parameter   :: rthr_vdw = 4000.0_wp

! ------------------------------------------------------------------------
!  Output: Initialization and Parameter setup
! ------------------------------------------------------------------------
   call d4init(mol,opt%g_a,opt%g_c,p_refq_goedecker,ndim)
   memory = ((mol%nat)+(mol%nat)+(ndim)+(ndim*ndim) &
            +(mol%nat*mol%nat)+(23*mol%nat)+(mol%nat)+(3*mol%nat) &
            +(3*mol%nat*mol%nat)+(3*mol%nat*(mol%nat+1)) &
            +(3*3*mol%nat)+(3*3*(mol%nat+1))+(3*3*mol%nat) &
            +(3*mol%nat*mol%nat)) * wp / (1024.0_wp**2)

   call get_realspace_cutoff(mol%lattice,rthr_cn, rep_cn)
   call get_realspace_cutoff(mol%lattice,rthr_vdw,rep_vdw)
   call get_realspace_cutoff(mol%lattice,rthr_mbd,rep_mbd)

   call generic_header(iunit,'Calculation Setup',49,10)
   write(iunit,'(3x,a," : ",i6)')   'number of atoms     ', mol%nat
   write(iunit,'(3x,a," : ",i6)')   'charge              ', nint(mol%chrg)
   write(iunit,'(3x,a," : ",a)')    'non-additivity corr.', lmbd2string(opt%lmbd)
   write(iunit,'(3x,a," : ",a)')    'charge model        ', refq2string(p_refq_goedecker)
!$ write(iunit,'(3x,a," : ",i6)')   'omp threads         ',omp_get_num_threads()
   write(iunit,'(3x,a," : ",f6.1,1x,a)') &
      "memory needed (est.)",memory,"Mb"
   if (opt%wf  /= 6.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f6.1)')   'weighting factor    ', opt%wf
   if (opt%g_a /= 3.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f6.1)')   'q-scale height      ', opt%g_a
   if (opt%g_c /= 2.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f6.1)')   'q-scale steepness   ', opt%g_c
   write(iunit,'(3x,a," : ",i0,2(" × ",i0))') 'supercell for CN    ',2*rep_cn +1
   write(iunit,'(3x,a," : ",i0,2(" × ",i0))') 'supercell for disp. ',2*rep_vdw+1
   write(iunit,'(3x,a," : ",i0,2(" × ",i0))') 'supercell for ATM   ',2*rep_mbd+1

   write(iunit,'(a)')

   if (.not.opt%silent) &
   call prd4ref(mol)

   allocate( q(mol%nat),covcn(mol%nat),gweights(ndim),refc6(ndim,ndim),&
             c6ab(mol%nat,mol%nat),aw(23,mol%nat),cn(mol%nat), &
             dcndr(3,mol%nat,mol%nat),dcndL(3,3,mol%nat), &
             dqdr(3,mol%nat,mol%nat+1),dqdL(3,3,mol%nat+1), &
             dcovcndr(3,mol%nat,mol%nat), dcovcndL(3,3,mol%nat), &
             source = 0.0_wp, stat = err )
   if (err /= 0) &
      call raise('E','Memory allocation failed')
   sigma = 0.0_wp

   call pbc_dncoord_d4(mol,covcn,dcovcndr,dcovcndL,rep_cn,rthr_cn)
   call d4(mol,ndim,opt%wf,opt%g_a,opt%g_c,covcn,gweights,refc6)

   call pbc_derfcoord(mol,cn,dcndr,dcndL,rthr_cn)
   call dncoord_logcn(mol%nat,cn,dcndr,dcndL,cn_max=8.0_wp)

! ------------------------------------------------------------------------
!  get partial charges
! ------------------------------------------------------------------------
   if (opt%verbose) &
   call eeq_header
   call new_charge_model(chrgeq,mol)
   call eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL,es,gtmp,stmp, &
                   .false.,.false.,.true.)
   if (opt%verbose) &
   call print_chrgeq(iunit,chrgeq,mol,q,cn)

! ------------------------------------------------------------------------
!  calculate properties
! ------------------------------------------------------------------------
dispersion_properties: if (opt%lmolpol) then
   call generic_header(iunit,'Molecular Properties',49,10)
   call mdisp(mol,ndim,q,opt%g_a,opt%g_c,gweights,refc6,molc6,molc8,molpol,aw,c6ab)
   call prmolc6(mol,molc6,molc8,molpol,covcn=covcn,q=q,c6ab=c6ab,alpha=aw(1,:))
endif dispersion_properties

if (.not.opt%silent.and.(opt%lenergy.or.opt%lgradient.or.opt%lhessian)) then
   call generic_header(iunit,'Damping Parameters',49,10)
   write(iunit,'(3x,a," : ",f10.4)')   's6                  ', dparam%s6
   write(iunit,'(3x,a," : ",f10.4)')   's8                  ', dparam%s8
   if (dparam%s10 /= 0.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f10.4)')   's10                 ', dparam%s10
   if (dparam%s9  /= 1.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f10.4)')   's9                  ', dparam%s9
   write(iunit,'(3x,a," : ",f10.4)')   'a1                  ', dparam%a1
   write(iunit,'(3x,a," : ",f10.4)')   'a2                  ', dparam%a2
   write(iunit,'(a)')
endif

! ------------------------------------------------------------------------
!  calculate energy
! ------------------------------------------------------------------------
dispersion_energy: if (opt%lenergy) then
   call edisp_3d(mol,ndim,q,rep_vdw,rep_mbd,rthr_vdw,rthr_mbd,dparam, &
      &          opt%g_a,opt%g_c,gweights,refc6,opt%lmbd,energy, &
      &          etwo=etwo,embd=emany)
endif dispersion_energy

! ------------------------------------------------------------------------
!  calculate gradient
! ------------------------------------------------------------------------
dispersion_gradient: if (opt%lgradient) then
   call dispgrad_3d(mol,ndim,q,covcn,dcovcndr,dcovcndL,rep_vdw,rep_mbd, &
      &             rthr_vdw,rthr_mbd,dparam,opt%wf,opt%g_a,opt%g_c, &
      &             refc6,opt%lmbd,gradient,sigma,etmp,dqdr,dqdL)
   if (.not.opt%lenergy) energy = etmp
   inv_lat = mat_inv_3x3(mol%lattice)
   call sigma_to_latgrad(sigma,inv_lat,latgrad)
endif dispersion_gradient

end subroutine d4_pbc_calculation
