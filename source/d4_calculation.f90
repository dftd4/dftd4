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

!> interface to the DFT-D4 module
module dispersion_calculator
contains
subroutine d4_calculation(iunit,env,opt,mol,dparam,dresults)
   use iso_fortran_env, wp => real64
!$ use omp_lib

   use mctc_environment

! ------------------------------------------------------------------------
!  class definitions
! ------------------------------------------------------------------------
   use class_set
   use class_param
   use class_molecule
   use class_results

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
   type(mctc_logger), intent(inout) :: env
   !> molecule structure information
   type(molecule), intent(inout) :: mol
   !> calculation options for DFT-D4 method
   type(dftd_options), intent(in) :: opt
   !> damping parameter for the DFT-D4 method
   type(dftd_parameter),intent(in) :: dparam
   type(chrg_parameter)            :: chrgeq
   type(dispersion_model),allocatable :: dispm
! ------------------------------------------------------------------------
!  output variables
! ------------------------------------------------------------------------
   type(dftd_results), intent(out) :: dresults

! ------------------------------------------------------------------------
!  local variables
! ------------------------------------------------------------------------
   integer  :: ndim                      ! matrix dimension
   integer  :: i,j,k,l,ii,jj
   integer  :: err
   real(wp) :: memory
   real(wp) :: etmp,etwo,emany,er,el,es
   real(wp) :: molpol,molc6,molc8        ! molecular Polarizibility
   real(wp) :: sigma(3,3),stmp(3,3)
   real(wp) :: latgrad(3,3),inv_lat(3,3)
   !> dispersion energy, always referenced
   real(wp) :: energy
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
   real(wp),allocatable :: c8ab(:,:)
   real(wp),allocatable :: aw(:,:)
   real(wp),allocatable :: ges(:,:)
   real(wp),allocatable :: gr(:,:)
   real(wp),allocatable :: gl(:,:)
   real(wp),allocatable :: gradient(:,:)
   real(wp),allocatable :: hessian(:,:)
   real(wp),parameter   :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp),parameter   :: rthr_mbd = 1600.0_wp
   real(wp),parameter   :: rthr_vdw = 4000.0_wp
   logical,parameter :: voigt_mask(3,3) = reshape( &
      & [.true.,.true.,.true.,.false.,.true.,.true.,.false.,.false.,.true.],&
      & shape(voigt_mask))

   integer :: stat
   logical :: minpr, verb, debug

   allocate(dispm)

   minpr = opt%print_level > 0
   verb  = opt%print_level > 1
   debug = opt%print_level > 2

   if (verb) then
      call dftd4_header(iunit,debug)
      call dftd4_citation(iunit)
   endif

! ------------------------------------------------------------------------
!  Output: Initialization and Parameter setup
! ------------------------------------------------------------------------
   call dispm%new(mol%at,p_refq_goedecker,opt%g_a,opt%g_c)
   ndim = sum(dispm%atoms*dispm%nref)
   memory = ((mol%n)+(mol%n)+(ndim)+(ndim*ndim) &
      &     +(mol%n*mol%n)+(23*mol%n)+(mol%n)+(3*mol%n) &
      &     +(3*mol%n*mol%n)+(3*mol%n*(mol%n+1)) &
      &     +(3*mol%n*mol%n)) * wp / (1024.0_wp**2)

   if (minpr) then
   call generic_header(iunit,'Calculation Setup',49,10)
   write(iunit,'(3x,a," : ",i6)')   'number of atoms     ', mol%n
   write(iunit,'(3x,a," : ",i6)')   'charge              ', nint(mol%chrg)
   write(iunit,'(3x,a," : ",a)')    'non-additivity corr.', lmbd2string(opt%lmbd)
   write(iunit,'(3x,a," : ",a)')    'charge model        ', refq2string(p_refq_goedecker)
!$ write(iunit,'(3x,a," : ",i6)')   'omp threads         ',omp_get_num_threads()
   write(iunit,'(3x,a," : ",f6.1,1x,a)') &
      "memory needed (est.)",memory,"Mb"
   if (opt%wf  /= 6.0_wp .or. verb) &
   write(iunit,'(3x,a," : ",f6.1)')   'weighting factor    ', opt%wf
   if (opt%g_a /= 3.0_wp .or. verb) &
   write(iunit,'(3x,a," : ",f6.1)')   'q-scale height      ', opt%g_a
   if (opt%g_c /= 2.0_wp .or. verb) &
   write(iunit,'(3x,a," : ",f6.1)')   'q-scale steepness   ', opt%g_c

   write(iunit,'(a)')

   call prd4ref(iunit,mol,dispm)
   endif

   allocate( q(mol%n),covcn(mol%n),gweights(ndim),refc6(ndim,ndim),&
      &      c6ab(mol%n,mol%n),aw(23,mol%n),cn(mol%n),dqdL(3,3,mol%n+1), &
      &      dcndr(3,mol%n,mol%n),dcndL(3,3,mol%n),dqdr(3,mol%n,mol%n+1), &
      &      dcovcndr(3,mol%n,mol%n),dcovcndL(3,3,mol%n), &
      &      source = 0.0_wp, stat = err )
   if (err /= 0) then
      call env%error(3,'Memory allocation failed')
      return
   endif
   energy = 0.0_wp

   call pbc_dncoord_d4(mol,covcn,dcovcndr,dcovcndL)
   call d4(mol,dispm,ndim,opt%wf,covcn,gweights,refc6)

   call pbc_dncoord_erf(mol,cn,dcndr,dcndL)
   call dncoord_logcn(mol%n,cn,dcndr,cn_max=8.0_wp)

! ------------------------------------------------------------------------
!  get partial charges
! ------------------------------------------------------------------------
   if (verb) &
   call eeq_header(iunit)
   call new_charge_model(chrgeq,mol)
   stat = eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma, &
                   .false.,.false.,.true.)
   if (stat.ne.0) then
      call env%error(3,"EEQ model could not be solved")
      return
   endif
   if (verb) &
   call print_chrgeq(iunit,chrgeq,mol,q,cn)

   dresults%charges = q
   dresults%dipole_moment = matmul(mol%xyz,q)

! ------------------------------------------------------------------------
!  calculate properties
! ------------------------------------------------------------------------
dispersion_properties: if (opt%lmolpol) then
   if (minpr) &
   call generic_header(iunit,'Molecular Properties',49,10)
   call mdisp(mol,dispm,ndim,q,gweights,refc6,molc6,molc8,molpol,aw,c6ab)
   dresults%polarizibilities = aw(1,:)
   dresults%c6_coefficients = c6ab
   if (minpr) &
   call prmolc6(iunit,mol,molc6,molc8,molpol,covcn=covcn,q=q,c6ab=c6ab,alpha=aw(1,:))
endif dispersion_properties

if (minpr.and.(opt%lenergy.or.opt%lgradient.or.opt%lhessian)) then
   call generic_header(iunit,'Damping Parameters',49,10)
   write(iunit,'(3x,a," : ",f10.4)')   's6                  ', dparam%s6
   write(iunit,'(3x,a," : ",f10.4)')   's8                  ', dparam%s8
   if (dparam%s10 /= 0.0_wp .or. verb) &
   write(iunit,'(3x,a," : ",f10.4)')   's10                 ', dparam%s10
   if (dparam%s9  /= 1.0_wp .or. verb) &
   write(iunit,'(3x,a," : ",f10.4)')   's9                  ', dparam%s9
   write(iunit,'(3x,a," : ",f10.4)')   'a1                  ', dparam%a1
   write(iunit,'(3x,a," : ",f10.4)')   'a2                  ', dparam%a2
   write(iunit,'(a)')
endif

! ------------------------------------------------------------------------
!  calculate energy
! ------------------------------------------------------------------------
dispersion_energy: if (.not.opt%lgradient .and. opt%lenergy) then
   call edisp_3d(mol,dispm,ndim,q,rthr_vdw,rthr_mbd,dparam, &
                 gweights,refc6,opt%lmbd,energy,etwo=etwo,embd=emany)
   dresults%energy = energy
endif dispersion_energy

! ------------------------------------------------------------------------
!  calculate gradient
! ------------------------------------------------------------------------
dispersion_gradient: if (opt%lgradient) then
   allocate( gradient(3,mol%n), source = 0.0_wp )
   call dispgrad_3d(mol,dispm,ndim,q,covcn,dcovcndr,dcovcndL,rthr_vdw,rthr_mbd, &
      &             dparam,opt%wf,refc6,opt%lmbd,gradient,sigma,energy, &
      &             dqdr,dqdL)
   if (mol%npbc > 0) then
      inv_lat = mat_inv_3x3(mol%lattice)
      call sigma_to_latgrad(sigma,inv_lat,latgrad)
   endif

   dresults%energy = energy
   dresults%gradient = gradient
   if (mol%npbc > 0) then
      dresults%stress = sigma/mol%volume
      dresults%lattice_gradient = latgrad
   endif

endif dispersion_gradient

! ------------------------------------------------------------------------
!  calculate hessian
! ------------------------------------------------------------------------
dispersion_hessian: if (opt%lhessian) then
   allocate(hessian(3*mol%n,3*mol%n),gr(3,mol%n),gl(3,mol%n), source = 0.0_wp)
   do i = 1, mol%n
      do j = 1, 3
         ii = 3*(i-1)+j
         gr = 0.0_wp
         gl = 0.0_wp

         mol%xyz(j,i) = mol%xyz(j,i) + step
         call pbc_dncoord_erf(mol,cn,dcndr,dcndL)
         call dncoord_logcn(mol%n,cn,dcndr,dcndL,cn_max=8.0_wp)
         stat = eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma, &
            &            .false.,.false.,.true.)
         call pbc_dncoord_d4(mol,covcn,dcovcndr,dcovcndL)
         call dispgrad_3d(mol,dispm,ndim,q,covcn,dcovcndr,dcovcndL,rthr_vdw,rthr_mbd, &
            &             dparam,opt%wf,refc6,opt%lmbd, &
            &             gr,stmp,er,dqdr,dqdL)

         mol%xyz(j,i) = mol%xyz(j,i) - 2*step
         call pbc_dncoord_erf(mol,cn,dcndr,dcndL)
         call dncoord_logcn(mol%n,cn,dcndr,dcndL,cn_max=8.0_wp)
         stat = eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma, &
                         .false.,.false.,.true.)
         call pbc_dncoord_d4(mol,covcn,dcovcndr,dcovcndL)
         call dispgrad_3d(mol,dispm,ndim,q,covcn,dcovcndr,dcovcndL,rthr_vdw,rthr_mbd, &
            &             dparam,opt%wf,refc6,opt%lmbd, &
            &             gl,stmp,el,dqdr,dqdL)

         mol%xyz(j,i) = mol%xyz(j,i) + step
         do k = 1, mol%n
            do l = 1, 3
               jj = 3*(k-1) + l
               hessian(jj,ii) = (gr(l,k)-gl(l,k))*step2
            enddo
         enddo
      enddo
   enddo
   deallocate(gr,gl)
   dresults%hessian = hessian
endif dispersion_hessian

end subroutine d4_calculation

subroutine d3_calculation(iunit,env,opt,mol,dparam,dresults)
   use iso_fortran_env, wp => real64
!$ use omp_lib

   use mctc_environment

! ------------------------------------------------------------------------
!  class definitions
! ------------------------------------------------------------------------
   use class_set
   use class_param
   use class_molecule
   use class_results

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
   type(mctc_logger), intent(inout) :: env
   !> molecule structure information
   type(molecule), intent(inout) :: mol
   !> calculation options for DFT-D4 method
   type(dftd_options), intent(in) :: opt
   !> damping parameter for the DFT-D4 method
   type(dftd_parameter),intent(in) :: dparam
   type(chrg_parameter)            :: chrgeq
   type(dispersion_model),allocatable :: dispm
! ------------------------------------------------------------------------
!  output variables
! ------------------------------------------------------------------------
   type(dftd_results), intent(out) :: dresults

! ------------------------------------------------------------------------
!  local variables
! ------------------------------------------------------------------------
   integer  :: ndim                      ! matrix dimension
   integer  :: i,j,k,l,ii,jj
   integer  :: err
   real(wp) :: memory
   real(wp) :: etmp,etwo,emany,er,el,es
   real(wp) :: molpol,molc6,molc8        ! molecular Polarizibility
   real(wp) :: sigma(3,3),stmp(3,3)
   real(wp) :: latgrad(3,3),inv_lat(3,3)
   !> dispersion energy, always referenced
   real(wp) :: energy
   real(wp),allocatable :: q(:)          ! partial charges
   real(wp),allocatable :: dqdr(:,:,:)   ! partial charges
   real(wp),allocatable :: dqdL(:,:,:)   ! partial charges
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
   real(wp),allocatable :: gradient(:,:)
   real(wp),allocatable :: hessian(:,:)
   real(wp),parameter   :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp),parameter   :: rthr_mbd = 1600.0_wp
   real(wp),parameter   :: rthr_vdw = 4000.0_wp
   logical,parameter :: voigt_mask(3,3) = reshape( &
      & [.true.,.true.,.true.,.false.,.true.,.true.,.false.,.false.,.true.],&
      & shape(voigt_mask))

   integer :: stat
   logical :: minpr, verb, debug

   allocate(dispm)

   minpr = opt%print_level > 0
   verb  = opt%print_level > 1
   debug = opt%print_level > 2

   if (verb) then
      call dftd4_header(iunit,debug)
      call dftd4_citation(iunit)
   endif
   if (minpr) call d3_model_warning(iunit)

! ------------------------------------------------------------------------
!  Output: Initialization and Parameter setup
! ------------------------------------------------------------------------
   call dispm%new(mol%at)
   ndim = sum(dispm%atoms*dispm%nref)
   memory = ((mol%n)+(mol%n)+(ndim)+(ndim*ndim) &
      &     +(mol%n*mol%n)+(23*mol%n)+(mol%n)+(3*mol%n) &
      &     +(3*mol%n*mol%n)+(3*mol%n*(mol%n+1)) &
      &     +(3*mol%n*mol%n)) * wp / (1024.0_wp**2)

   if (minpr) then
   call generic_header(iunit,'Calculation Setup',49,10)
   write(iunit,'(3x,a," : ",i6)')   'number of atoms     ', mol%n
   write(iunit,'(3x,a," : ",a)')    'non-additivity corr.', lmbd2string(opt%lmbd)
!$ write(iunit,'(3x,a," : ",i6)')   'omp threads         ',omp_get_num_threads()
   write(iunit,'(3x,a," : ",f6.1,1x,a)') &
      "memory needed (est.)",memory,"Mb"
   if (opt%wf  /= 6.0_wp .or. verb) &
   write(iunit,'(3x,a," : ",f6.1)')   'weighting factor    ', opt%wf

   write(iunit,'(a)')

   call prd4ref(iunit,mol,dispm)
   endif

   allocate( q(mol%n),gweights(ndim),refc6(ndim,ndim),&
      &      c6ab(mol%n,mol%n),aw(23,mol%n),cn(mol%n),dqdL(3,3,mol%n+1), &
      &      dcndr(3,mol%n,mol%n),dcndL(3,3,mol%n),dqdr(3,mol%n,mol%n+1), &
      &      source = 0.0_wp, stat = err )
   if (err /= 0) then
      call env%error(3,'Memory allocation failed')
      return
   endif
   energy = 0.0_wp

   call pbc_dncoord_erf(mol,cn,dcndr,dcndL)
   call d4(mol,dispm,ndim,opt%wf,cn,gweights,refc6)

   q = 0.0_wp
   dqdr = 0.0_wp
   dqdL = 0.0_wp

! ------------------------------------------------------------------------
!  calculate properties
! ------------------------------------------------------------------------
dispersion_properties: if (opt%lmolpol) then
   if (minpr) &
   call generic_header(iunit,'Molecular Properties',49,10)
   call mdisp(mol,dispm,ndim,q,gweights,refc6,molc6,molc8,molpol,aw,c6ab)
   dresults%polarizibilities = aw(1,:)
   dresults%c6_coefficients = c6ab
   if (minpr) &
   call prmolc6(iunit,mol,molc6,molc8,molpol,cn=cn,q=q,c6ab=c6ab,alpha=aw(1,:))
endif dispersion_properties

if (minpr.and.(opt%lenergy.or.opt%lgradient.or.opt%lhessian)) then
   call generic_header(iunit,'Damping Parameters',49,10)
   write(iunit,'(3x,a," : ",f10.4)')   's6                  ', dparam%s6
   write(iunit,'(3x,a," : ",f10.4)')   's8                  ', dparam%s8
   if (dparam%s10 /= 0.0_wp .or. verb) &
   write(iunit,'(3x,a," : ",f10.4)')   's10                 ', dparam%s10
   if (dparam%s9  /= 1.0_wp .or. verb) &
   write(iunit,'(3x,a," : ",f10.4)')   's9                  ', dparam%s9
   write(iunit,'(3x,a," : ",f10.4)')   'a1                  ', dparam%a1
   write(iunit,'(3x,a," : ",f10.4)')   'a2                  ', dparam%a2
   write(iunit,'(a)')
endif

! ------------------------------------------------------------------------
!  calculate energy
! ------------------------------------------------------------------------
dispersion_energy: if (.not.opt%lgradient .and. opt%lenergy) then
   call edisp_3d(mol,dispm,ndim,q,rthr_vdw,rthr_mbd,dparam, &
                 gweights,refc6,opt%lmbd,energy,etwo=etwo,embd=emany)
   dresults%energy = energy
endif dispersion_energy

! ------------------------------------------------------------------------
!  calculate gradient
! ------------------------------------------------------------------------
dispersion_gradient: if (opt%lgradient) then
   allocate( gradient(3,mol%n), source = 0.0_wp )
   call dispgrad_3d(mol,dispm,ndim,q,cn,dcndr,dcndL,rthr_vdw,rthr_mbd, &
      &             dparam,opt%wf,refc6,opt%lmbd,gradient,sigma,energy, &
      &             dqdr,dqdL)
   if (mol%npbc > 0) then
      inv_lat = mat_inv_3x3(mol%lattice)
      call sigma_to_latgrad(sigma,inv_lat,latgrad)
   endif

   dresults%energy = energy
   dresults%gradient = gradient
   if (mol%npbc > 0) then
      dresults%stress = sigma/mol%volume
      dresults%lattice_gradient = latgrad
   endif

endif dispersion_gradient

! ------------------------------------------------------------------------
!  calculate hessian
! ------------------------------------------------------------------------
dispersion_hessian: if (opt%lhessian) then
   allocate(hessian(3*mol%n,3*mol%n),gr(3,mol%n),gl(3,mol%n), source = 0.0_wp)
   do i = 1, mol%n
      do j = 1, 3
         ii = 3*(i-1)+j
         gr = 0.0_wp
         gl = 0.0_wp

         mol%xyz(j,i) = mol%xyz(j,i) + step
         call pbc_dncoord_erf(mol,cn,dcndr,dcndL)
         call dispgrad_3d(mol,dispm,ndim,q,cn,dcndr,dcndL,rthr_vdw,rthr_mbd, &
            &             dparam,opt%wf,refc6,opt%lmbd, &
            &             gr,stmp,er,dqdr,dqdL)

         mol%xyz(j,i) = mol%xyz(j,i) - 2*step
         call pbc_dncoord_erf(mol,cn,dcndr,dcndL)
         call dispgrad_3d(mol,dispm,ndim,q,cn,dcndr,dcndL,rthr_vdw,rthr_mbd, &
            &             dparam,opt%wf,refc6,opt%lmbd, &
            &             gl,stmp,el,dqdr,dqdL)

         mol%xyz(j,i) = mol%xyz(j,i) + step
         do k = 1, mol%n
            do l = 1, 3
               jj = 3*(k-1) + l
               hessian(jj,ii) = (gr(l,k)-gl(l,k))*step2
            enddo
         enddo
      enddo
   enddo
   deallocate(gr,gl)
   dresults%hessian = hessian
endif dispersion_hessian

end subroutine d3_calculation
end module dispersion_calculator

