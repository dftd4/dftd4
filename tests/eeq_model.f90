!> @brief test calculation of EEQ charges
subroutine test_eeq_model_water
   use iso_fortran_env, wp => real64
   use assertion
   use class_molecule
   use class_param
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

   real(wp) :: es,sigma(3,3)
   real(wp),allocatable :: cn(:),dcndr(:,:,:),dcndL(:,:,:)
   real(wp),allocatable :: q(:),dqdr(:,:,:),dqdL(:,:,:),ges(:,:)
   integer :: stat

   allocate( cn(nat),dcndr(3,nat,nat),dcndL(3,3,nat),q(nat),dqdr(3,nat,nat+1),ges(3,nat) )
   es  = 0.0_wp
   ges = 0.0_wp

   call mol%allocate(nat,.false.)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp
   mol%pbc = .false.
   mol%npbc = 0
   mol%lattice = 0.0_wp
   
   call pbc_dncoord_erf(mol,cn,dcndr,dcndL)
   ! test for correct CN and correct symmetry in CN
   call assert_close(cn(1),1.9895544848535_wp,thr)
   call assert_close(cn(2),cn(3),             thr)
   call assert_close(cn(1),cn(2)+cn(3),       thr)

   ! test CN derivative, correct summation of diagonals
   call assert_close(dcndr(3,1,1),-0.80973198569003E-01_wp,thr)
   call assert_close(dcndr(1,2,1), 0.52891163425093E-01_wp,thr)
   call assert_close(dcndr(3,1,3),-0.40486599284501E-01_wp,thr)

   call new_charge_model_2019(chrgeq,mol)
   stat = eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma,.false.,.true.,.true.)
   call assert_eq(stat,0)

   ! test electrostatic energy
   call assert_close(es,-0.64308088326667E-01_wp,thr)

   ! test molecular gradient of ES, also check symmetry
   call assert_close(ges(3,1),-0.44053032330503E-01_wp,thr)
   call assert_close(ges(3,2),ges(3,3),                thr)
   call assert_close(ges(1,2), 0.18102071270235E-01_wp,thr)

   ! test for charge constraint
   call assert_close(sum(q),0.0_wp,            thr)
   call assert_close(q(1),-0.59582744684480_wp,thr)

   ! test derivative of partial charges
   call assert_close(dqdr(3,1,1),-0.41466764014389E+00_wp,thr)
   call assert_close(dqdr(1,2,1), 0.17196993288918E+00_wp,thr)
   call assert_close(dqdr(1,3,2), 0.18631993794394E-01_wp,thr)
   call assert_close(dqdr(2,1,3), 0.00000000000000E+00_wp,thr)

   call mol%deallocate

   ! done: everythings fine
   call terminate(afail)
end subroutine test_eeq_model_water

subroutine test_eeq_model_ewald
   use iso_fortran_env, wp => real64
   use assertion
   use class_molecule
   use class_param
   use eeq_model
   use coordination_number
   use pbc_tools
   implicit none
   real(wp),parameter :: thr = 1.0e-9_wp
   integer, parameter :: nat = 6
   integer, parameter :: at(nat) = [14,14,8,8,8,8]
   real(wp),parameter :: abc(3,nat) = reshape(&
      &[.095985472469032_wp, .049722204206931_wp, 0.10160624337938_wp, &
      & 0.54722204206931_wp, 0.52863628207623_wp, 0.38664208660311_wp, &
      & 0.29843937068984_wp, 0.39572194413818_wp, 0.20321248675876_wp, &
      & 0.23364982659922_wp, 0.85647058758674_wp, 0.31884968761485_wp, &
      & 0.72250232459952_wp, 0.65548544066844_wp, .056207709103487_wp, &
      & 0.70514214000043_wp, 0.28321754549582_wp, 0.36424822189074_wp],&
      & shape(abc))
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[ 8.7413053236641_wp,  0.0000000000000_wp,  0.0000000000000_wp,   &
      &  0.0000000000000_wp,  8.7413053236641_wp,  0.0000000000000_wp,   &
      &  0.0000000000000_wp,  0.0000000000000_wp,  8.7413053236641_wp],  &
      & shape(lattice))
   integer, parameter :: wsc_rep(3) = [1,1,1]

   type(molecule)       :: mol
   type(chrg_parameter) :: chrgeq
   real(wp)             :: energy
   real(wp)             :: sigma(3,3)
   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: dcndr(:,:,:)
   real(wp),allocatable :: dcndL(:,:,:)
   real(wp),allocatable :: q(:)
   real(wp),allocatable :: dqdr(:,:,:)
   real(wp),allocatable :: dqdL(:,:,:)
   real(wp),allocatable :: gradient(:,:)
   integer :: stat

   call mol%allocate(nat,.true.)
   mol%at   = at
   mol%abc  = abc
   mol%npbc = 3
   mol%pbc  = .true.
   mol%lattice = lattice
   mol%volume = dlat_to_dvol(lattice)
   call dlat_to_cell(lattice,mol%cellpar)
   call dlat_to_rlat(lattice,mol%rec_lat)
   call coord_trafo(nat,lattice,abc,mol%xyz)
   call mol%wrap_back

   allocate( cn(nat), dcndr(3,nat,nat), dcndL(3,3,nat), &
      &      q(nat), dqdr(3,nat,nat+1), dqdL(3,3,nat+1), &
      &      gradient(3,nat), source = 0.0_wp )
   cn    = 0.0_wp
   dcndr = 0.0_wp
   dcndL = 0.0_wp
   q     = 0.0_wp
   dqdr  = 0.0_wp
   dqdL  = 0.0_wp
   gradient = 0.0_wp
   energy = 0.0_wp
   sigma = 0.0_wp

   call generate_wsc(mol,mol%wsc,wsc_rep)

   call pbc_dncoord_erf(mol,cn,dcndr,dcndL,900.0_wp)

   call assert_close(cn(2),3.6864725130236_wp,thr)
   call assert_close(cn(5),1.0523558225297_wp,thr)
   call assert_close(cn(6),1.1699488478421_wp,thr)

   call assert_close(dcndr(2,3,2),-0.24281192795725E-02_wp,thr)
   call assert_close(dcndr(1,6,6),-0.42240789876965E+00_wp,thr)
   call assert_close(dcndr(1,1,2),-0.71913132896636E-01_wp,thr)

   call assert_close(dcndL(1,3,4),-0.36068395425405_wp,thr)
   call assert_close(dcndL(3,3,1),-0.92808088688242_wp,thr)
   call assert_close(dcndL(2,1,3),-0.44871260782914_wp,thr)

   call new_charge_model_2019(chrgeq,mol)

   stat = eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL,energy,gradient,sigma,&
      &            .false.,.true.,.true.)
   call assert_eq(stat,0)

   call assert_close(energy,-0.90576568382295E-01_wp,thr)

   call assert_close(gradient(2,1), 0.59224535975576E-02_wp,thr)
   call assert_close(gradient(1,3),-0.17291053212435E-02_wp,thr)
   call assert_close(gradient(3,5), 0.13109339409759E-02_wp,thr)
   call assert_close(gradient(1,6), 0.10491055487530E-01_wp,thr)

   call assert_close(sigma(1,1),-0.69556402379101E-01_wp,thr)
   call assert_close(sigma(2,3),-0.11063392011917E-01_wp,thr)
   call assert_close(sigma(3,1), 0.15351359672016E-01_wp,thr)

   call assert_close(sum(q),0.0_wp,            thr)
   call assert_close(q(1), 0.39808668429315_wp,thr)
   call assert_close(q(3),-0.16133275860565_wp,thr)
   call assert_close(q(4),-0.19939812633388_wp,thr)

   call assert_close(dqdr(1,4,2),-0.16619680140869E-01_wp,thr)
   call assert_close(dqdr(2,5,7),-0.71528076609028E-04_wp,thr)
   call assert_close(dqdr(1,2,2), 0.36177437391610E-01_wp,thr)
   call assert_close(dqdr(3,1,4),-0.39909176876716E-01_wp,thr)

   call assert_close(dqdL(2,3,2),-0.66131872796832E-02_wp,thr)
   call assert_close(dqdL(2,1,7), 0.66206886462644E-02_wp,thr)
   call assert_close(dqdL(3,2,1), 0.51752004657789E-01_wp,thr)
   call assert_close(dqdL(1,1,5),-0.29765138003366E-01_wp,thr)

   call mol%deallocate

   ! done
   call terminate(afail)

end subroutine test_eeq_model_ewald

subroutine test_eeq_numgrad
   use iso_fortran_env, wp => real64, istderr => error_unit
   use assertion
   use mctc_constants
   use class_molecule
   use class_param
   use coordination_number
   use eeq_model
   implicit none

   real(wp),parameter :: thr = 1.0e-9_wp
   integer, parameter :: nat = 12
   integer, parameter :: at(nat) = [9,1,6,6,6,6,6,6,9,1,1,9]
   real(wp),parameter :: xyz(3,nat) = reshape([&
      &-3.61707571665846_wp,-2.21326241566733_wp,-0.20914111039455_wp, &
      &-4.12085912923274_wp, 2.64361027892242_wp,-0.77360733295900_wp, &
      & 0.70571179547784_wp,-1.46804568562192_wp, 0.40260319330016_wp, &
      &-1.70284271199247_wp,-0.55881754958407_wp,-0.07245619686109_wp, &
      &-2.22721391221737_wp, 1.98624548249273_wp,-0.40689082126896_wp, &
      &-0.21016646823535_wp, 3.64828289401216_wp,-0.24136575031326_wp, &
      & 2.24393351740138_wp, 2.86731529508407_wp, 0.23334807696843_wp, &
      & 2.63467443489279_wp, 0.29590222227677_wp, 0.54871426109542_wp, &
      &-0.65510226134586_wp, 6.12711283170949_wp,-0.54503487549757_wp, &
      & 3.78023043231610_wp, 4.20058933929303_wp, 0.35763959771486_wp, &
      & 1.05737533649814_wp,-3.45504988937755_wp, 0.68084171455529_wp, &
      & 4.98699896397970_wp,-0.51782531551223_wp, 1.02289296328971_wp], shape(xyz))
   type(molecule)    :: mol
   type(chrg_parameter) :: chrgeq
   integer  :: i,j
   real(wp),parameter :: step = 1.0e-6_wp, step2 = 0.5_wp/step
   real(wp) :: er,el
   real(wp),allocatable :: numg(:,:)
   real(wp) :: es,gsolv,sigma(3,3)
   real(wp),allocatable :: cn(:),dcndr(:,:,:),dcndL(:,:,:)
   real(wp),allocatable :: q(:),grd(:,:)
   real(wp),allocatable :: dqdr(:,:,:)
   real(wp),allocatable :: dqdL(:,:,:)
   real(wp),allocatable :: numq(:,:,:),ql(:),qr(:)
   integer :: stat

   allocate( cn(nat),q(nat), grd(3,nat), dcndr(3,nat,nat), &
      &      dqdr(3,nat,nat),dcndL(3,3,nat), source = 0.0_wp )
   es  = 0.0_wp

   call mol%allocate(nat)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp
   call generate_wsc(mol,mol%wsc)

   call pbc_dncoord_erf(mol,cn,dcndr,dcndL)
   call dncoord_logcn(mol%n,cn,dcndr,cn_max=8.0_wp)

   call new_charge_model_2019(chrgeq,mol)
   stat = eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL, &
      &            es,grd,sigma,.false.,.true.,.true.)
   call assert_eq(stat,0)

   write(*,*) matmul(xyz,q)

   allocate( numg(3,nat), source = 0.0_wp )
   allocate( numq(3,nat,nat), ql(nat), qr(nat), source = 0.0_wp )
   print'(/,"analytical gradient")'
   print *, grd
   do i = 1, nat
      do j = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(j,i) = mol%xyz(j,i) + step
         call pbc_ncoord_erf(mol,cn)
         call dncoord_logcn(mol%n,cn,cn_max=8.0_wp)
         stat = eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,qr,dqdr,dqdL, &
            &            er,grd,sigma,.false.,.false.,.false.)

         mol%xyz(j,i) = mol%xyz(j,i) - 2*step
         call pbc_ncoord_erf(mol,cn)
         call dncoord_logcn(mol%n,cn,cn_max=8.0_wp)
         stat = eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,ql,dqdr,dqdL, &
            &            el,grd,sigma,.false.,.false.,.false.)

         mol%xyz(j,i) = mol%xyz(j,i) + step
         numg(j,i) = step2 * (er-el)
         numq(j,i,:) = step2 * (qr-ql)
      enddo
   enddo

   print'(/,"numerical gradient")'
   print *, numg
   print'(/,"difference gradient")'
   print*,grd-numg
   print'(/,"difference response")'
   print*,dqdr(:,1,1)-numq(:,1,1)

   call assert_close(norm2(grd-numg),0.0_wp,1.0e-8_wp)
   call assert_close(norm2(dqdr-numq),0.0_wp,1.0e-8_wp)

   call terminate(afail)
end subroutine test_eeq_numgrad

subroutine test_peeq_numgrad
   use iso_fortran_env, wp => real64, istderr => error_unit
   use assertion
   use mctc_constants
   use class_molecule
   use class_param
   use coordination_number
   use eeq_model
   use pbc_tools
   implicit none

   real(wp),parameter :: thr = 1.0e-9_wp
   ! CaF2
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [9,9,20]
   real(wp),parameter :: abc(3,nat) = reshape(&
      &[0.24_wp, 0.26_wp, 0.27_wp, &
      & 0.75_wp, 0.75_wp, 0.75_wp, &
      & 0.00_wp, 0.00_wp, 0.00_wp], shape(abc))
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[5.9598811567890_wp,      2.1071361905157_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      6.3214085715472_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      0.0000000000000_wp,      7.2993338808807_wp],   &
      & shape(lattice))
   type(molecule)    :: mol
   type(chrg_parameter) :: chrgeq
   integer  :: i,j
   real(wp),parameter :: step = 1.0e-6_wp, step2 = 0.5_wp/step
   real(wp) :: er,el
   real(wp),allocatable :: numg(:,:)
   real(wp) :: es,gsolv,sigma(3,3)
   real(wp),allocatable :: cn(:),dcndr(:,:,:),dcndL(:,:,:)
   real(wp),allocatable :: q(:),grd(:,:)
   real(wp),allocatable :: dqdr(:,:,:)
   real(wp),allocatable :: dqdL(:,:,:)
   real(wp),allocatable :: numq(:,:,:),ql(:),qr(:)
   integer :: stat

   call mol%allocate(nat)
   mol%at   = at
   mol%npbc = 3
   mol%pbc  = .true.
   mol%lattice = lattice
   call coord_trafo(nat,lattice,abc,mol%xyz)
   call mol%update

   call generate_wsc(mol,mol%wsc)

   allocate( cn(nat),q(nat), grd(3,nat), dcndr(3,nat,nat), dcndL(3,3,nat), &
      &      dqdr(3,nat,nat), dqdL(3,3,nat), source = 0.0_wp )
   es  = 0.0_wp
   sigma = 0.0_wp

   call pbc_dncoord_erf(mol,cn,dcndr,dcndL)
   call dncoord_logcn(mol%n,cn,dcndr,cn_max=8.0_wp)

   call new_charge_model_2019(chrgeq,mol)
   stat = eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL, &
      &            es,grd,sigma,.false.,.true.,.true.)
   call assert_eq(stat,0)
   print*,es

   allocate( numg(3,nat), source = 0.0_wp )
   allocate( numq(3,nat,nat), ql(nat), qr(nat), source = 0.0_wp )
   print'(/,"analytical gradient")'
   print *, grd
   do i = 1, nat
      do j = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(j,i) = mol%xyz(j,i) + step
         call pbc_ncoord_erf(mol,cn)
         call dncoord_logcn(mol%n,cn,cn_max=8.0_wp)
         stat = eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,qr,dqdr,dqdL, &
            &            er,grd,sigma,.false.,.false.,.false.)

         mol%xyz(j,i) = mol%xyz(j,i) - 2*step
         call pbc_ncoord_erf(mol,cn)
         call dncoord_logcn(mol%n,cn,cn_max=8.0_wp)
         stat = eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,ql,dqdr,dqdL, &
            &            el,grd,sigma,.false.,.false.,.false.)

         mol%xyz(j,i) = mol%xyz(j,i) + step
         numg(j,i) = step2 * (er-el)
         numq(j,i,:) = step2 * (qr-ql)
      enddo
   enddo

   print'(/,"numerical gradient")'
   print *, numg
   print'(/,"difference gradient")'
   print*,grd-numg
   print'(/,"difference response")'
   print*,dqdr(:,1,1)-numq(:,1,1)
   print*,dqdr(:,1,nat)-numq(:,1,nat)
   print*,dqdr(:,nat,nat/2)-numq(:,nat,nat/2)

   call assert_close(norm2(grd-numg),0.0_wp,1.0e-9_wp)
   call assert_close(norm2(dqdr-numq),0.0_wp,1.0e-8_wp)

   call terminate(afail)
end subroutine test_peeq_numgrad
