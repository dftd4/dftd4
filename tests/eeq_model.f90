!> @brief test calculation of EEQ charges
subroutine test_eeq_model_water
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
end subroutine test_eeq_model_water

subroutine test_eeq_model_ewald
   use iso_fortran_env, wp => real64, istdout => output_unit
   use assertion
   use class_molecule
   use class_param
   use eeq_model
   use coordination_number
   use pbc_tools
   implicit none
   real(wp),parameter :: thr = 1.0e-10_wp
   ! CaF2
!   integer, parameter :: nat = 3
!   integer, parameter :: at(nat) = [9,9,20]
!   real(wp),parameter :: abc(3,nat) = reshape(&
!      &[0.25_wp, 0.25_wp, 0.25_wp, &
!      & 0.75_wp, 0.75_wp, 0.75_wp, &
!      & 0.00_wp, 0.00_wp, 0.00_wp], shape(abc))
!   real(wp),parameter :: lattice(3,3) = reshape(&
!      &[5.9598811567890_wp,      2.1071361905157_wp,      3.6496669404404_wp,    &
!      & 0.0000000000000_wp,      6.3214085715472_wp,      3.6496669404404_wp,    &
!      & 0.0000000000000_wp,      0.0000000000000_wp,      7.2993338808807_wp],   &
!      & shape(lattice))
   integer, parameter :: nat = 6
   integer, parameter :: at(nat) = [14,14,8,8,8,8]
   real(wp),parameter :: abc(3,nat) = reshape(&
      &[.095985472469032_wp,     .049722204206931_wp,     0.10160624337938_wp, &
      & 0.54722204206931_wp,     0.52863628207623_wp,     0.38664208660311_wp, &
      & 0.29843937068984_wp,     0.39572194413818_wp,     0.20321248675876_wp, &
      & 0.23364982659922_wp,     0.85647058758674_wp,     0.31884968761485_wp, &
      & 0.72250232459952_wp,     0.65548544066844_wp,     .056207709103487_wp, &
      & 0.70514214000043_wp,     0.28321754549582_wp,     0.36424822189074_wp],&
      & shape(abc))
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[ 8.7413053236641_wp,      0.0000000000000_wp,      0.0000000000000_wp,   &
      &  0.0000000000000_wp,      8.7413053236641_wp,      0.0000000000000_wp,   &
      &  0.0000000000000_wp,      0.0000000000000_wp,      8.7413053236641_wp],  &
      & shape(lattice))
   integer, parameter :: wsc_rep(3) = [1,1,1]

   type(molecule)       :: mol
   type(chrg_parameter) :: chrgeq
   real(wp)             :: energy
   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: dcndr(:,:,:)
   real(wp),allocatable :: q(:)
   real(wp),allocatable :: dqdr(:,:,:)
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

   allocate( cn(nat), dcndr(3,nat,nat), q(nat), dqdr(3,nat,nat+1), &
      &      gradient(3,nat), source = 0.0_wp )
   cn    = 0.0_wp
   dcndr = 0.0_wp
   q     = 0.0_wp
   dqdr  = 0.0_wp
   gradient = 0.0_wp
   energy = 0.0_wp

   call generate_wsc(mol,mol%wsc,wsc_rep)

   call print_pbcsum(istdout,mol)

   call pbc_derfcoord(mol,cn,dcndr,900.0d0)

   call assert_close(cn(2),3.6864723770739_wp,thr)
   call assert_close(cn(5),1.0523557617324_wp,thr)
   call assert_close(cn(6),1.1699487724900_wp,thr)

   call assert_close(dcndr(2,3,2),0.24281208190040E-02_wp,thr)
   call assert_close(dcndr(1,6,6),0.42240778640017E+00_wp,thr)
   call assert_close(dcndr(1,1,2),0.71913063700324E-01_wp,thr)

   call new_charge_model_2019(chrgeq,mol)

   call eeq_chrgeq(chrgeq,mol,cn,dcndr,q,dqdr,energy,gradient,&
      &            .false.,.true.,.true.)

   call assert_close(energy,-0.90576573072013E-01_wp,thr)
   
   call assert_close(gradient(2,1), 0.59224555288324E-02_wp,thr)
   call assert_close(gradient(1,3),-0.17291063655794E-02_wp,thr)
   call assert_close(gradient(3,5), 0.13109348427124E-02_wp,thr)
   call assert_close(gradient(1,6), 0.10491053931694E-01_wp,thr)

   call assert_close(sum(q),0.0_wp,            thr)
   call assert_close(q(1), 0.39808669536607_wp,thr)
   call assert_close(q(3),-0.16133275992899_wp,thr)
   call assert_close(q(4),-0.19939813398061_wp,thr)

   call assert_close(dqdr(1,4,2),-0.16619684227607E-01_wp,thr)
   call assert_close(dqdr(2,5,7),-0.71529121856161E-04_wp,thr)
   call assert_close(dqdr(1,2,2), 0.36177439010511E-01_wp,thr)
   call assert_close(dqdr(3,1,4),-0.39909178885213E-01_wp,thr)

   call mol%deallocate

   ! done
   call terminate(0)

end subroutine test_eeq_model_ewald
