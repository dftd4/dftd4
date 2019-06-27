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

   allocate( cn(nat),dcndr(3,nat,nat),q(nat),dqdr(3,nat,nat+1),ges(3,nat) )
   es  = 0.0_wp
   ges = 0.0_wp

   call mol%allocate(nat,.false.)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp
   
   call dncoord_erf(mol,cn,dcndr)
   ! test for correct CN and correct symmetry in CN
   call assert_close(cn(1),1.9895544848535_wp,thr)
   call assert_close(cn(2),cn(3),             thr)
   call assert_close(cn(1),cn(2)+cn(3),       thr)

   ! test CN derivative, correct summation of diagonals
   call assert_close(dcndr(3,1,1), 0.80973198569003E-01_wp,thr)
   call assert_close(dcndr(1,2,1), 0.52891163425093E-01_wp,thr)
   call assert_close(dcndr(3,1,3),-0.40486599284501E-01_wp,thr)

   call new_charge_model_2019(chrgeq,mol)
   call eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma,.false.,.true.,.true.)

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

   call pbc_derfcoord(mol,cn,dcndr,dcndL,900.0_wp)

   call assert_close(cn(2),3.6864725130236_wp,thr)
   call assert_close(cn(5),1.0523558225297_wp,thr)
   call assert_close(cn(6),1.1699488478421_wp,thr)

   call assert_close(dcndr(2,3,2),0.24281192795725E-02_wp,thr)
   call assert_close(dcndr(1,6,6),0.42240789876965E+00_wp,thr)
   call assert_close(dcndr(1,1,2),0.71913132896636E-01_wp,thr)

   call new_charge_model_2019(chrgeq,mol)

   call eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL,energy,gradient,sigma,&
      &            .false.,.true.,.true.)

   call assert_close(energy,-0.90576568382295E-01_wp,thr)

   call assert_close(gradient(2,1), 0.91930618913995E-03_wp,thr)
   call assert_close(gradient(1,3),-0.74402915552075E-03_wp,thr)
   call assert_close(gradient(3,5),-0.62990375312383E-03_wp,thr)
   call assert_close(gradient(1,6), 0.48745476414615E-02_wp,thr)

   call assert_close(sum(q),0.0_wp,            thr)
   call assert_close(q(1), 0.39808668429315_wp,thr)
   call assert_close(q(3),-0.16133275860565_wp,thr)
   call assert_close(q(4),-0.19939812633388_wp,thr)

   call assert_close(dqdr(1,4,2), 0.15400805612228E-01_wp,thr)
   call assert_close(dqdr(2,5,7), 0.14771235997449E-02_wp,thr)
   call assert_close(dqdr(1,2,2), 0.21715273501275E-01_wp,thr)
   call assert_close(dqdr(3,1,4),-0.33298482841066E-01_wp,thr)

   call mol%deallocate

   ! done
   call terminate(afail)

end subroutine test_eeq_model_ewald

subroutine test_eeq_model_cell_gradient
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

   call pbc_derfcoord(mol,cn,dcndr,dcndL,900.0_wp)

   call assert_close(cn(3),1.9369733824363_wp,thr)

   call assert_close(dcndr(1,6,6),0.42240789876965E+00_wp,thr)

   call assert_close(dcndL(1,3,4),-0.36068395425405_wp,thr)
   call assert_close(dcndL(3,3,1),-0.92808088688242_wp,thr)
   call assert_close(dcndL(2,1,3),-0.44871260782914_wp,thr)

   call new_charge_model_2019(chrgeq,mol)

   call eeq_chrgeq(chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL,energy,gradient,sigma,&
      &            .false.,.true.,.true.)

   call assert_close(energy,-0.90576568382295E-01_wp,thr)

   call assert_close(gradient(2,1), 0.91930618913995E-03_wp,thr)

   call assert_close(sigma(1,1),-0.69556402225023E-01_wp,thr)
   call assert_close(sigma(2,3),-0.11063392011917E-01_wp,thr)
   call assert_close(sigma(3,1), 0.15351359672016E-01_wp,thr)

   call assert_close(sum(q),0.0_wp,            thr)
   call assert_close(q(1), 0.39808668429315_wp,thr)

   call assert_close(dqdr(3,1,4),-0.33298482841066E-01_wp,thr)

   call assert_close(dqdL(2,3,2),-0.66131872796832E-02_wp,thr)
   call assert_close(dqdL(2,1,7), 0.66206886462644E-02_wp,thr)
   call assert_close(dqdL(3,2,1), 0.51752004657789E-01_wp,thr)
   call assert_close(dqdL(1,1,5),-0.29765138003366E-01_wp,thr)

   call mol%deallocate

   ! done
   call terminate(afail)

end subroutine test_eeq_model_cell_gradient
