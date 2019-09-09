! This file is part of dftd4.
!
! Copyright (C) 2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> @brief contains all functions needed for calculating partial charges
module eeq_model
   use iso_fortran_env, wp => real64

! ------------------------------------------------------------------------
!  class definitions
! ------------------------------------------------------------------------
   use class_param, only : chrg_parameter

   implicit none

!> @brief π itself
   real(wp),private,parameter :: pi = 3.1415926535897932384626433832795029_wp
!> @brief √π
   real(wp),private,parameter :: sqrtpi = sqrt(pi)
!> @brief √(2/π)
   real(wp),private,parameter :: sqrt2pi = sqrt(2.0_wp/pi)


   !> select one charge model as default one
   interface new_charge_model
      module procedure :: new_charge_model_2019
   end interface new_charge_model

   interface eeq_chrgeq
      module procedure :: eeq_chrgeq_core
   end interface eeq_chrgeq

contains

!> load parameters for the extended electronegativity model
!
!  Fetch fitted parameters for electronegativity, CN scaling factor,
!  chemical hardness and atom radius for a given molecule type.
subroutine new_charge_model_2019(chrgeq,mol)
   use class_molecule
   implicit none
   type(chrg_parameter) :: chrgeq !< charge model for EEQ model
   type(molecule) :: mol          !< molecular structure data

   integer,parameter :: max_elem = 86
   real(wp) :: en(max_elem)
   real(wp) :: gamm(max_elem)
   real(wp) :: cnfak(max_elem)
   real(wp) :: alp(max_elem)

! ------------------------------------------------------------------------
!  PARAMETRISATION BY S. SPICHER (NEW)
! ------------------------------------------------------------------------
   parameter( en =(/&
    1.23695041_wp, 1.26590957_wp, 0.54341808_wp, 0.99666991_wp, 1.26691604_wp, &
    1.40028282_wp, 1.55819364_wp, 1.56866440_wp, 1.57540015_wp, 1.15056627_wp, &
    0.55936220_wp, 0.72373742_wp, 1.12910844_wp, 1.12306840_wp, 1.52672442_wp, &
    1.40768172_wp, 1.48154584_wp, 1.31062963_wp, 0.40374140_wp, 0.75442607_wp, &
    0.76482096_wp, 0.98457281_wp, 0.96702598_wp, 1.05266584_wp, 0.93274875_wp, &
    1.04025281_wp, 0.92738624_wp, 1.07419210_wp, 1.07900668_wp, 1.04712861_wp, &
    1.15018618_wp, 1.15388455_wp, 1.36313743_wp, 1.36485106_wp, 1.39801837_wp, &
    1.18695346_wp, 0.36273870_wp, 0.58797255_wp, 0.71961946_wp, 0.96158233_wp, &
    0.89585296_wp, 0.81360499_wp, 1.00794665_wp, 0.92613682_wp, 1.09152285_wp, &
    1.14907070_wp, 1.13508911_wp, 1.08853785_wp, 1.11005982_wp, 1.12452195_wp, &
    1.21642129_wp, 1.36507125_wp, 1.40340000_wp, 1.16653482_wp, 0.34125098_wp, &
    0.58884173_wp, 0.68441115_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
    0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
    0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
    0.56999999_wp, 0.87936784_wp, 1.02761808_wp, 0.93297476_wp, 1.10172128_wp, &
    0.97350071_wp, 1.16695666_wp, 1.23997927_wp, 1.18464453_wp, 1.14191734_wp, &
    1.12334192_wp, 1.01485321_wp, 1.12950808_wp, 1.30804834_wp, 1.33689961_wp, &
    1.27465977_wp /))
   parameter( gamm =(/&
   -0.35015861_wp, 1.04121227_wp, 0.09281243_wp, 0.09412380_wp, 0.26629137_wp, &
    0.19408787_wp, 0.05317918_wp, 0.03151644_wp, 0.32275132_wp, 1.30996037_wp, &
    0.24206510_wp, 0.04147733_wp, 0.11634126_wp, 0.13155266_wp, 0.15350650_wp, &
    0.15250997_wp, 0.17523529_wp, 0.28774450_wp, 0.42937314_wp, 0.01896455_wp, &
    0.07179178_wp,-0.01121381_wp,-0.03093370_wp, 0.02716319_wp,-0.01843812_wp, &
   -0.15270393_wp,-0.09192645_wp,-0.13418723_wp,-0.09861139_wp, 0.18338109_wp, &
    0.08299615_wp, 0.11370033_wp, 0.19005278_wp, 0.10980677_wp, 0.12327841_wp, &
    0.25345554_wp, 0.58615231_wp, 0.16093861_wp, 0.04548530_wp,-0.02478645_wp, &
    0.01909943_wp, 0.01402541_wp,-0.03595279_wp, 0.01137752_wp,-0.03697213_wp, &
    0.08009416_wp, 0.02274892_wp, 0.12801822_wp,-0.02078702_wp, 0.05284319_wp, &
    0.07581190_wp, 0.09663758_wp, 0.09547417_wp, 0.07803344_wp, 0.64913257_wp, &
    0.15348654_wp, 0.05054344_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
    0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
    0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
    0.11000000_wp,-0.02786741_wp, 0.01057858_wp,-0.03892226_wp,-0.04574364_wp, &
   -0.03874080_wp,-0.03782372_wp,-0.07046855_wp, 0.09546597_wp, 0.21953269_wp, &
    0.02522348_wp, 0.15263050_wp, 0.08042611_wp, 0.01878626_wp, 0.08715453_wp, &
    0.10500484_wp /))
   parameter( cnfak =(/&
    0.04916110_wp, 0.10937243_wp,-0.12349591_wp,-0.02665108_wp,-0.02631658_wp, &
    0.06005196_wp, 0.09279548_wp, 0.11689703_wp, 0.15704746_wp, 0.07987901_wp, &
   -0.10002962_wp,-0.07712863_wp,-0.02170561_wp,-0.04964052_wp, 0.14250599_wp, &
    0.07126660_wp, 0.13682750_wp, 0.14877121_wp,-0.10219289_wp,-0.08979338_wp, &
   -0.08273597_wp,-0.01754829_wp,-0.02765460_wp,-0.02558926_wp,-0.08010286_wp, &
   -0.04163215_wp,-0.09369631_wp,-0.03774117_wp,-0.05759708_wp, 0.02431998_wp, &
   -0.01056270_wp,-0.02692862_wp, 0.07657769_wp, 0.06561608_wp, 0.08006749_wp, &
    0.14139200_wp,-0.05351029_wp,-0.06701705_wp,-0.07377246_wp,-0.02927768_wp, &
   -0.03867291_wp,-0.06929825_wp,-0.04485293_wp,-0.04800824_wp,-0.01484022_wp, &
    0.07917502_wp, 0.06619243_wp, 0.02434095_wp,-0.01505548_wp,-0.03030768_wp, &
    0.01418235_wp, 0.08953411_wp, 0.08967527_wp, 0.07277771_wp,-0.02129476_wp, &
   -0.06188828_wp,-0.06568203_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
   -0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
   -0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
   -0.11000000_wp,-0.03585873_wp,-0.03132400_wp,-0.05902379_wp,-0.02827592_wp, &
   -0.07606260_wp,-0.02123839_wp, 0.03814822_wp, 0.02146834_wp, 0.01580538_wp, &
   -0.00894298_wp,-0.05864876_wp,-0.01817842_wp, 0.07721851_wp, 0.07936083_wp, &
    0.05849285_wp /))
   parameter( alp =(/&
    0.55159092_wp, 0.66205886_wp, 0.90529132_wp, 1.51710827_wp, 2.86070364_wp, &
    1.88862966_wp, 1.32250290_wp, 1.23166285_wp, 1.77503721_wp, 1.11955204_wp, &
    1.28263182_wp, 1.22344336_wp, 1.70936266_wp, 1.54075036_wp, 1.38200579_wp, &
    2.18849322_wp, 1.36779065_wp, 1.27039703_wp, 1.64466502_wp, 1.58859404_wp, &
    1.65357953_wp, 1.50021521_wp, 1.30104175_wp, 1.46301827_wp, 1.32928147_wp, &
    1.02766713_wp, 1.02291377_wp, 0.94343886_wp, 1.14881311_wp, 1.47080755_wp, &
    1.76901636_wp, 1.98724061_wp, 2.41244711_wp, 2.26739524_wp, 2.95378999_wp, &
    1.20807752_wp, 1.65941046_wp, 1.62733880_wp, 1.61344972_wp, 1.63220728_wp, &
    1.60899928_wp, 1.43501286_wp, 1.54559205_wp, 1.32663678_wp, 1.37644152_wp, &
    1.36051851_wp, 1.23395526_wp, 1.65734544_wp, 1.53895240_wp, 1.97542736_wp, &
    1.97636542_wp, 2.05432381_wp, 3.80138135_wp, 1.43893803_wp, 1.75505957_wp, &
    1.59815118_wp, 1.76401732_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
    1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
    1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
    1.63999999_wp, 1.47055223_wp, 1.81127084_wp, 1.40189963_wp, 1.54015481_wp, &
    1.33721475_wp, 1.57165422_wp, 1.04815857_wp, 1.78342098_wp, 2.79106396_wp, &
    1.78160840_wp, 2.47588882_wp, 2.37670734_wp, 1.76613217_wp, 2.66172302_wp, &
    2.82773085_wp /))

   integer :: i

   call chrgeq%allocate(mol%n)

   do i = 1, mol%n
      chrgeq%xi   (i) = en   (mol%at(i))
      chrgeq%gam  (i) = gamm (mol%at(i))
      chrgeq%kappa(i) = cnfak(mol%at(i))
      chrgeq%alpha(i) = alp  (mol%at(i))**2 ! square it
   enddo

end subroutine new_charge_model_2019

!> print a charge model and its results in a human readable form
subroutine print_chrgeq(iunit,chrgeq,mol,q,cn)
   use iso_fortran_env, wp => real64
   use class_molecule
   implicit none
   integer, intent(in) :: iunit
   type(molecule),intent(in) :: mol          !< structure information
   type(chrg_parameter),intent(in) :: chrgeq !< parametrisation
   real(wp),intent(in)    :: q(mol%n)      !< partial charges
   real(wp),intent(in)    :: cn(mol%n)     !< erf-CN
   real(wp) :: mmom(3)
   integer  :: i
   write(iunit,'(a)')
   write(iunit,'(7x,"   #   Z   ")',advance='no')
   write(iunit,'("         q")',advance='no')
   write(iunit,'("        CN")',advance='no')
   write(iunit,'("        EN")',advance='no')
   write(iunit,'("       Aii")',advance='no')
   write(iunit,'(a)')
   do i=1,mol%n
      write(iunit,'(i11,1x,i3,1x,a2)',advance='no') &
      &     i,mol%at(i),mol%sym(i)
      write(iunit,'(f10.3)',advance='no')q(i)
      write(iunit,'(f10.3)',advance='no')cn(i)
      write(iunit,'(f10.3)',advance='no')chrgeq%xi(i)-chrgeq%kappa(i)*sqrt(cn(i))
      write(iunit,'(f10.3)',advance='no')chrgeq%gam(i)+sqrt2pi/sqrt(chrgeq%alpha(i))
      write(iunit,'(a)')
   enddo
   mmom = 0.0_wp
   do i = 1, mol%n
      mmom = mmom + q(i)*mol%xyz(:,i)
   enddo
   write(iunit,'(a)')
   write(iunit,'(7x,a)')'dipole moment:'
   write(iunit,'(18x,a)')'x           y           z       tot (au)'
   write(iunit,'(7x,4f12.3)') mmom,norm2(mmom)
   write(iunit,'(a)')

end subroutine print_chrgeq

! ========================================================================
!  Purpose:
! ------------------------------------------------------------------------
!> calculate energies and properties from the EEQ model
!
!  calculates only energy and partial charges by default, but can
!  calculate derivatives of the energy and partial charges with
!  repect to the nuclear coordinates. Works in the space of the
!  Lagranian for coupled perturbed charges, allocatation should be
!  on n+1 not n.
!
!  @note For computational effiency this routine does NOT work with the
!        energy expression but with the Lagrangian (n+1 instead of n),
!        therefore, all intermediates live in the space of the Lagrangian.
!
!  Implemented by SAW in 2018, see focusing lab course report for details.
subroutine eeq_chrgeq_core &
      & (chrgeq,mol,cn,dcndr,dcndL,q,dqdr,dqdL,energy,gradient,sigma,&
      &  lverbose,lgrad,lcpq)
   use iso_fortran_env, wp => real64, istdout => output_unit
   use class_molecule
   use coordination_number
   implicit none

! ------------------------------------------------------------------------
!  Input
! ------------------------------------------------------------------------
   type(molecule),intent(in) :: mol                    !< structure information
   type(chrg_parameter),intent(in) :: chrgeq           !< parametrisation
   real(wp),intent(in)    :: cn(mol%n)               !< erf-CN
   real(wp),intent(in)    :: dcndr(3,mol%n,mol%n)  !< derivative of erf-CN
   real(wp),intent(in)    :: dcndL(3,3,mol%n)        !< lattice derivative of erf-CN
   logical, intent(in)    :: lverbose                  !< toggles printout
   logical, intent(in)    :: lgrad                     !< flag for gradient calculation
   logical, intent(in)    :: lcpq                      !< do partial charge derivative
! ------------------------------------------------------------------------
!  Output
! ------------------------------------------------------------------------
   real(wp),intent(out)   :: q(mol%n)                !< partial charges
   real(wp),intent(out)   :: dqdr(3,mol%n,mol%n+1) !< derivative of partial charges
   real(wp),intent(out)   :: dqdL(3,3,mol%n+1)       !< lattice derivative of partial charges
   real(wp),intent(inout) :: energy                    !< electrostatic energy
   real(wp),intent(inout) :: sigma(3,3)                !< sigma tensor of IES
   real(wp),intent(inout) :: gradient(3,mol%n)       !< molecular gradient of IES
!
! ------------------------------------------------------------------------
!  charge model
! ------------------------------------------------------------------------
   integer  :: m ! dimension of the Lagrangian
   real(wp),allocatable :: Amat(:,:)
   real(wp),allocatable :: Xvec(:)
   real(wp),allocatable :: Ainv(:,:)
   real(wp),allocatable :: dAmatdr(:,:,:)
   real(wp),allocatable :: dAmatdL(:,:,:)
   real(wp),allocatable :: dXvecdr(:,:,:)
   real(wp),allocatable :: dXvecdL(:,:,:)

! ------------------------------------------------------------------------
!  local variables
! ------------------------------------------------------------------------
   integer  :: i,j,k,l
   real(wp) :: r,rij(3),r2
   real(wp) :: gamij,gamij2
   real(wp) :: arg,tmp,dtmp
   real(wp) :: lambda
   real(wp) :: es

!! ------------------------------------------------------------------------
!  for periodic case
!! ------------------------------------------------------------------------
   integer  :: wscAt
   real(wp) :: gamii,riw(3),dAtmp(3),stmp(3,3),t(3)
   real(wp) :: cf ! convergence factor
   real(wp),parameter :: zero(3) = [0.0_wp,0.0_wp,0.0_wp]
   integer, parameter :: ewaldCutD(3) = [2,2,2]
   integer, parameter :: ewaldCutR(3) = [2,2,2]


! ------------------------------------------------------------------------
!  scratch variables
! ------------------------------------------------------------------------
   real(wp),allocatable :: xtmp(:)
   real(wp),allocatable :: atmp(:,:)
   real(wp),allocatable :: Xfac(:)
   real(wp),allocatable :: Afac(:,:)

! ------------------------------------------------------------------------
!  Lapack work variables
! ------------------------------------------------------------------------
   integer, allocatable :: ipiv(:)
   real(wp),allocatable :: temp(:)
   real(wp),allocatable :: work(:)
   integer  :: lwork
   integer  :: info
   real(wp) :: test(1)
   real(wp),external :: ddot

   ! quick return if possible
   if (mol%n == 1) then
      q = mol%chrg
      ! should always hold, even for extendend systems
      es = 0.0_wp
      return
   endif

! ------------------------------------------------------------------------
!  initizialization
! ------------------------------------------------------------------------
   m    = mol%n+1
   q    = 0.0_wp
   allocate( ipiv(m), source = 0 )
   allocate( Amat(m,m), Xvec(m), Xfac(m), source = 0.0_wp )

! ------------------------------------------------------------------------
!  set up the A matrix and X vector
! ------------------------------------------------------------------------
!  αi -> alpha(i), ENi -> xi(i), κi -> kappa(i), Jii -> gam(i)
!  γij = 1/√(αi+αj)
!  Xi  = -ENi + κi·√CNi
!  Aii = Jii + 2/√π·γii
!  Aij = erf(γij·Rij)/Rij = 2/√π·F0(γ²ij·R²ij)
! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Setup of the A matrix and the RHS X vector"

!  prepare some arrays
!$omp parallel default(none) &
!$omp shared(mol,cn,chrgeq) &
!$omp private(i,tmp) &
!$omp shared(Xvec,Xfac)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      tmp = chrgeq%kappa(i)/(sqrt(cn(i))+1e-14_wp)
      Xvec(i) = -chrgeq%xi(i) + tmp*cn(i)
      Xfac(i) = 0.5_wp*tmp
   enddo
!$omp enddo
!$omp endparallel

   if (mol%npbc > 0) then
   cf = sqrtpi/mol%volume**(1.0_wp/3.0_wp)
   ! build Ewald matrix
!$omp parallel default(none) &
!$omp private(i,j,wscAt,riw,gamii,gamij,t) &
!$omp shared(mol,chrgeq,cf) &
!$omp reduction(+:Amat)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      gamii = 1.0/sqrt(chrgeq%alpha(i)+chrgeq%alpha(i))
      Amat(i,i) = chrgeq%gam(i) + sqrt2pi/sqrt(chrgeq%alpha(i)) &
         ! reciprocal for 0th atom
         + eeq_ewald_3d_rec(zero,ewaldCutR,mol%rec_lat,mol%volume,cf) &
         ! direct for 0th atom
         + eeq_ewald_3d_dir(zero,ewaldCutD,mol%lattice,gamii,cf) &
         ! neutralizing background contribution
         - pi/mol%volume/cf**2
      do j = 1, i-1
         gamij = 1.0_wp/sqrt(chrgeq%alpha(i)+chrgeq%alpha(j))
         do wscAt = 1, mol%wsc%itbl(j,i)
            t = mol%wsc%lattr(:,wscAt,j,i)
            riw   = mol%xyz(:,i) - mol%xyz(:,j) - matmul(mol%lattice,t)
            Amat(i,j) = Amat(i,j) + mol%wsc%w(j,i) * ( &
               ! reciprocal lattice sums
               + eeq_ewald_3d_rec(riw,ewaldCutR,mol%rec_lat,mol%volume,cf) &
               ! direct lattice sums
               + eeq_ewald_3d_dir(riw,ewaldCutD,mol%lattice,gamij,cf) &
               ! neutralizing background contribution
               - pi/mol%volume/cf**2 )
         end do
         Amat(j,i) = Amat(i,j)
      end do
   end do
!$omp enddo
!$omp endparallel
   else
!$omp parallel default(none) &
!$omp shared(mol,chrgeq) &
!$omp private(i,j,r,gamij) &
!$omp shared(Xvec,Xfac,Amat)
!$omp do schedule(dynamic)
   ! prepare A matrix
   do i = 1, mol%n
      ! EN of atom i
      do j = 1, i-1
         r = norm2(mol%xyz(:,j) - mol%xyz(:,i))
         gamij = 1.0_wp/sqrt(chrgeq%alpha(i)+chrgeq%alpha(j))
         Amat(j,i) = erf(gamij*r)/r
         Amat(i,j) = Amat(j,i)
      enddo
      Amat(i,i) = chrgeq%gam(i) + sqrt2pi/sqrt(chrgeq%alpha(i))
   enddo
!$omp enddo
!$omp endparallel
   endif

   Amat(m,1:m) = 1.0_wp
   Amat(1:m,m) = 1.0_wp
   Amat(m,m  ) = 0.0_wp
   Xvec(m)     = mol%chrg

! ------------------------------------------------------------------------
!  invert the A matrix using a Bunch-Kaufman factorization
!  A⁻¹ = (L·D·L^T)⁻¹ = L^T·D⁻¹·L
! ------------------------------------------------------------------------
if (lcpq) then
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "A matrix inversion by Bunch-Kaufman factorization"
   allocate( Ainv(m,m), source = Amat )
   allocate( Xtmp(m),   source = 0.0_wp )

   ! assume work space query, set best value to test after first dsytrf call
   call dsytrf('L',m,Ainv,m,ipiv,test,-1,info)
   lwork=int(test(1))
   allocate( work(lwork), source = 0.0_wp )

   ! Bunch-Kaufman factorization A = L*D*L**T
   call dsytrf('L',m,Ainv,m,ipiv,work,lwork,info)
   if(info > 0)then
      call raise('E', '(eeq_inversion) DSYTRF failed')
   endif

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Ainv matrix
   ! Ainv matrix is overwritten with lower triangular part of A⁻¹
   call dsytri('L',m,Ainv,m,ipiv,work,info)
   if (info > 0) then
      call raise('E', '(eeq_inversion) DSYTRI failed')
   endif

   ! symmetrizes A⁻¹ matrix from lower triangular part of inverse matrix
   do i = 1, m
      do j = i+1, m
         Ainv(i,j)=Ainv(j,i)
      enddo
   enddo

   call dsymv('u',m,1.0_wp,Ainv,m,Xvec,1,0.0_wp,Xtmp,1)
else
! ------------------------------------------------------------------------
!  solve the linear equations to obtain partial charges
! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Solve the linear equations to obtain partial charges"
   ! generate temporary copy
   allocate( Atmp(m,m), source = Amat )
   allocate( Xtmp(m),   source = Xvec )

   ! assume work space query, set best value to test after first dsysv call
   call dsysv('u',m,1,Atmp,m,ipiv,Xtmp,m,test,-1,info)
   lwork = int(test(1))
   allocate( work(lwork), source = 0.0_wp )

   call dsysv('u',m,1,Atmp,m,ipiv,Xtmp,m,work,lwork,info)
   if(info > 0) call raise('E','(eeq_solve) DSYSV failed')

   !print'(3f20.14)',Xtmp

endif
   q = Xtmp(:mol%n)
   if(abs(sum(q)-mol%chrg) > 1.e-6_wp) &
      call raise('E','(eeq_solve) charge constrain error')
   lambda = Xtmp(m)
   if (lverbose) then
      write(istdout,'(72("-"))')
      write(istdout,'(1x,a,1x,f20.14)') &
         "lambda       :",lambda,&
         "total charge :",sum(q)
   endif

! ------------------------------------------------------------------------
!  calculate isotropic electrostatic (IES) energy
! ------------------------------------------------------------------------
!  E = ∑i (ENi - κi·√CNi)·qi + ∑i (Jii + 2/√π·γii)·q²i
!      + ½ ∑i ∑j,j≠i qi·qj·2/√π·F0(γ²ij·R²ij)
!    = q·(½A·q - X)
! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Isotropic electrostatic (IES) energy calculation"
   call dsymv('u',m,0.5_wp,Amat,m,Xtmp,1,-1.0_wp,Xvec,1)
   es = ddot(m,Xtmp,1,Xvec,1)
   energy = es + energy
   if (lverbose) then
      write(istdout,'(72("-"))')
      write(istdout,'(1x,a,1x,f20.14)') &
         "energy",es
   endif

! ------------------------------------------------------------------------
!  calculate molecular gradient of the IES energy
! ------------------------------------------------------------------------
!  dE/dRj -> g(:,j), ∂Xi/∂Rj -> -dcn(:,i,j), ½∂Aij/∂Rj -> dAmat(:,j,i)
!  dE/dR = (½∂A/∂R·q - ∂X/∂R)·q
!  ∂Aij/∂Rj = ∂Aij/∂Ri
! ------------------------------------------------------------------------
do_molecular_gradient: if (lgrad .or. lcpq) then
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "molecular gradient calculation"
   allocate( dAmatdr(3,mol%n,m), dXvecdr(3,mol%n,m), Afac(3,mol%n), &
      &      source = 0.0_wp )
   !allocate( dAmat(3,mol%n,m), Afac(3,mol%n), source = 0.0_wp )
   if (mol%npbc > 0) then
   allocate( dAmatdL(3,3,m), dXvecdL(3,3,m), source = 0.0_wp )
!$omp parallel default(none) &
!$omp shared(mol,chrgeq,dcndr,dcndL,Xfac,Xtmp,cf) &
!$omp private(i,j,wscAt,riw,gamij,dAtmp,stmp,t) &
!$omp shared(dXvecdr,dXvecdL) &
!$omp reduction(+:Afac,dAmatdr,dAmatdL)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      dXvecdr(:,:,i) = +dcndr(:,:,i)*Xfac(i)
      dXvecdL(:,:,i) = +dcndL(:,:,i)*Xfac(i)
      do j = 1, i-1
         ! over WSC partner
         gamij = 1.0_wp/sqrt(chrgeq%alpha(i) + chrgeq%alpha(j))
         do wscAt = 1, mol%wsc%itbl(j,i)
            t = mol%wsc%lattr(:,wscAt,j,i)
            riw = mol%xyz(:,i) - mol%xyz(:,j) - matmul(mol%lattice,t)
            call eeq_ewald_dx_3d_rec(riw,ewaldCutR,mol%rec_lat,mol%volume,cf, &
               &                     dAtmp,stmp)
            dAtmp = dAtmp*mol%wsc%w(j,i)
            stmp = stmp*mol%wsc%w(j,i)
            dAmatdr(:,i,j) = dAmatdr(:,i,j) + dAtmp(:)*Xtmp(i)
            dAmatdr(:,j,i) = dAmatdr(:,j,i) - dAtmp(:)*Xtmp(j)
            Afac(:,i)    =  Afac(:,i)   + dAtmp(:)*Xtmp(j)
            Afac(:,j)    =  Afac(:,j)   - dAtmp(:)*Xtmp(i)
            dAmatdL(:,:,i) = dAmatdL(:,:,i) + stmp(:,:)*Xtmp(j)
            dAmatdL(:,:,j) = dAmatdL(:,:,j) + stmp(:,:)*Xtmp(i)
            call eeq_ewald_dx_3d_dir(riw,ewaldCutD,mol%lattice,gamij,cf, &
               &                     dAtmp,stmp)
            dAtmp = dAtmp*mol%wsc%w(j,i)
            stmp = stmp*mol%wsc%w(j,i)
            dAmatdr(:,i,j) = dAmatdr(:,i,j) + dAtmp(:)*Xtmp(i)
            dAmatdr(:,j,i) = dAmatdr(:,j,i) - dAtmp(:)*Xtmp(j)
            Afac(:,i)    =  Afac(:,i)   + dAtmp(:)*Xtmp(j)
            Afac(:,j)    =  Afac(:,j)   - dAtmp(:)*Xtmp(i)
            dAmatdL(:,:,i) = dAmatdL(:,:,i) + stmp(:,:)*Xtmp(j)
            dAmatdL(:,:,j) = dAmatdL(:,:,j) + stmp(:,:)*Xtmp(i)
         enddo  ! k WSC partner
      enddo     ! j

      gamij = 1.0_wp/sqrt(chrgeq%alpha(i) + chrgeq%alpha(i))
      call eeq_ewald_dx_3d_rec(zero,ewaldCutR,mol%rec_lat,mol%volume,cf, &
         &                     dAtmp,stmp)
      dAmatdL(:,:,i) = dAmatdL(:,:,i) + stmp(:,:)*Xtmp(i)
      call eeq_ewald_dx_3d_dir(zero,ewaldCutD,mol%lattice,gamij,cf, &
         &                     dAtmp,stmp)
      dAmatdL(:,:,i) = dAmatdL(:,:,i) + stmp(:,:)*Xtmp(i)
      do j = 1,3
         dAmatdL(j,j,i) = dAmatdL(j,j,i) + cf/sqrtpi/3.0_wp*Xtmp(i)
      enddo
   enddo
!$omp enddo
!$omp endparallel
   else
!$omp parallel default(none) &
!$omp shared(mol,dcndr,Amat,Xfac,Xtmp,chrgeq) &
!$omp private(i,j,rij,r2,gamij,arg,dtmp) &
!$omp shared(dXvecdr,dAmatdr) &
!$omp reduction(+:Afac)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      dXvecdr(:,:,i) = +dcndr(:,:,i)*Xfac(i) ! merge dX and dA for speedup
      do j = 1, i-1
         rij = mol%xyz(:,i) - mol%xyz(:,j)
         r2 = sum(rij**2)
         gamij = 1.0_wp/sqrt(chrgeq%alpha(i) + chrgeq%alpha(j))
         arg = gamij**2*r2
         dtmp = 2.0_wp*gamij*exp(-arg)/(sqrtpi*r2)-Amat(j,i)/r2
         Afac(:,i) = +dtmp*rij*Xtmp(j) + Afac(:,i)
         Afac(:,j) = -dtmp*rij*Xtmp(i) + Afac(:,j)
         dAmatdr(:,i,j) = +dtmp*rij*Xtmp(i) !+ dcndr(:,i,j)*Xfac(i)
         dAmatdr(:,j,i) = -dtmp*rij*Xtmp(j) !+ dcndr(:,j,i)*Xfac(j)
      enddo
   enddo
!$omp enddo
!$omp endparallel
   endif
endif do_molecular_gradient

   if (lgrad) then
   call dgemv('n',3*mol%n,m,+1.0_wp,dAmatdr,3*mol%n,Xtmp,1,1.0_wp,gradient,1)
   call dgemv('n',3*mol%n,m,+1.0_wp,dXvecdr,3*mol%n,Xtmp,1,1.0_wp,gradient,1)
   if (mol%npbc > 0) then
      call dgemv('n',3*3,m,+0.5_wp,dAmatdL,3*3,Xtmp,1,1.0_wp,sigma,1)
      call dgemv('n',3*3,m,-1.0_wp,dXvecdL,3*3,Xtmp,1,1.0_wp,sigma,1)
   endif
   endif

! ------------------------------------------------------------------------
!  calculate gradient of the partial charge w.r.t. the nuclear coordinates
! ------------------------------------------------------------------------
do_partial_charge_derivative: if (lcpq) then
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "calculating the derivative of the partial charges"
   dqdr = 0.0_wp
   if (mol%npbc > 0) then
      dqdL = 0.0_wp
   endif
   do i = 1, mol%n
      dAmatdr(:,i,i) = Afac(:,i) + dAmatdr(:,i,i)
   enddo
   !call dsymm('r','l',3*mol%n,m,-1.0_wp,Ainv,m,dAmatdr,3*mol%n,1.0_wp,dqdr,3*mol%n)
   call dgemm('n','n',3*mol%n,m,m,-1.0_wp,dAmatdr,3*mol%n,Ainv,m,1.0_wp, &
      &       dqdr,3*mol%n)
   call dgemm('n','n',3*mol%n,m,m,-1.0_wp,dXvecdr,3*mol%n,Ainv,m,1.0_wp, &
      &       dqdr,3*mol%n)

   if (mol%npbc > 0) then
      call dgemm('n','n',3*3,m,m,-1.0_wp,dAmatdL,3*3,Ainv,m,1.0_wp,dqdL,3*3)
      call dgemm('n','n',3*3,m,m,+1.0_wp,dXvecdL,3*3,Ainv,m,1.0_wp,dqdL,3*3)
   endif
   !print'(/,"analytical gradient")'
   !print'(3f20.14)',dqdr(:,:,:mol%n)

endif do_partial_charge_derivative

! ------------------------------------------------------------------------
!  Clean up
! ------------------------------------------------------------------------
   if (allocated(Amat))    deallocate(Amat)
   if (allocated(dAmatdr)) deallocate(dAmatdr)
   if (allocated(dAmatdL)) deallocate(dAmatdL)
   if (allocated(Afac))    deallocate(Afac)
   if (allocated(Xvec))    deallocate(Xvec)
   if (allocated(dXvecdr)) deallocate(dXvecdr)
   if (allocated(dXvecdL)) deallocate(dXvecdL)
   if (allocated(Xfac))    deallocate(Xfac)
   if (allocated(Xtmp))    deallocate(Xtmp)
   if (allocated(Atmp))    deallocate(Atmp)
   if (allocated(temp))    deallocate(temp)
   if (allocated(work))    deallocate(work)
   if (allocated(ipiv))    deallocate(ipiv)

end subroutine eeq_chrgeq_core

!> direct space contribution to the Ewald matrix
pure function eeq_ewald_3d_dir(riw,rep,dlat,gamij,cf) result(Amat)
   use iso_fortran_env, wp => real64
   use mctc_constants
   implicit none
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   integer, intent(in) :: rep(3)    !< images to consider
   real(wp),intent(in) :: dlat(3,3) !< direct lattice
   real(wp),intent(in) :: gamij     !< interaction radius
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp) :: Amat                 !< element of interaction matrix
   real(wp),parameter :: eps = 1.0e-9_wp
   integer  :: dx,dy,dz
   real(wp) :: distiw,rij(3)
   real(wp) :: t(3)
   Amat = 0.0_wp
   do concurrent(dx = -rep(1):rep(1), &
         &       dy = -rep(2):rep(2), &
         &       dz = -rep(3):rep(3))
      t = [dx,dy,dz]
      rij = riw + matmul(dlat,t)
      distiw = norm2(rij)
      ! self-interaction correction
      if(distiw < eps) then
         Amat = Amat - cf/sqrtpi
      else
         Amat = Amat - erf(   cf*distiw)/distiw &
            &        + erf(gamij*distiw)/distiw
      end if
   end do
end function eeq_ewald_3d_dir

!> direct space contribution to the gradient and sigma tensor
pure subroutine eeq_ewald_dx_3d_dir(riw,rep,dlat,gamij,cf,dAmat,sigma)
   use iso_fortran_env, wp => real64
   use mctc_constants
   use pbc_tools
   implicit none
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   integer, intent(in) :: rep(3)    !< images to consider
   real(wp),intent(in) :: dlat(3,3) !< direct lattice
   real(wp),intent(in) :: gamij     !< interaction radius
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp),intent(out) :: dAmat(3) !< element of interaction matrix
   real(wp),intent(out) :: sigma(3,3)
   real(wp),parameter :: eps = 1.0e-9_wp
   integer  :: dx,dy,dz,i
   real(wp) :: distiw,rij(3),arga,argb
   real(wp) :: t(3),dtmp,stmp(3)
   dAmat = 0.0_wp
   sigma = 0.0_wp
   do concurrent(dx = -rep(1):rep(1), &
         &       dy = -rep(2):rep(2), &
         &       dz = -rep(3):rep(3))
      ! real contributions
      t = [dx,dy,dz]
      rij = riw + matmul(dlat,t)
      distiw = norm2(rij)
      if(distiw < eps) cycle
      arga = cf**2   *distiw**2
      stmp = + exp(-arga)/sqrtpi * cf * 2.0_wp / 3.0_wp
      do i = 1, 3
         sigma(i,i) = sigma(i,i) + stmp(i)! * rij(i)**2
      enddo
      argb = gamij**2*distiw**2
      dtmp = - 2*cf*exp(-arga)/(sqrtpi*distiw**2) &
         &   + erf(cf*distiw)/(distiw**3)           &
         &   + 2*gamij*exp(-argb)/(sqrtpi*distiw**2) &
         &   - erf(gamij*distiw)/(distiw**3)
      dAmat = dAmat + rij*dtmp
      sigma = sigma + dtmp*outer_prod_3x3(rij,rij)
   enddo

end subroutine eeq_ewald_dx_3d_dir

!> reciprocal lattice contribution to the Ewald matrix
pure function eeq_ewald_3d_rec(riw,rep,rlat,vol,cf) result(Amat)
   use iso_fortran_env, wp => real64
   use mctc_constants
   implicit none
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   integer, intent(in) :: rep(3)    !< images to consider
   real(wp),intent(in) :: rlat(3,3) !< reciprocal lattice
   real(wp),intent(in) :: vol       !< direct cell volume
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp) :: Amat                 !< element of interaction matrix
   integer  :: dx,dy,dz
   real(wp) :: rik2,rik(3)
   real(wp) :: t(3)
   Amat = 0.0_wp
   do concurrent(dx = -rep(1):rep(1), &
         &       dy = -rep(2):rep(2), &
         &       dz = -rep(3):rep(3))
      if (dx==0 .and. dy==0 .and. dz==0) cycle
      t = [dx,dy,dz]
      rik = matmul(rlat,t)
      rik2 = dot_product(rik,rik)
      Amat=Amat + cos(dot_product(rik,riw) ) * 4.0_wp*pi/vol &
         * exp(-rik2/(4.0_wp*cf**2))/rik2
   end do
end function eeq_ewald_3d_rec

!> reciprocal lattice contribution to the gradient and sigma tensor
pure subroutine eeq_ewald_dx_3d_rec(riw,rep,rlat,vol,cf,dAmat,sigma)
   use iso_fortran_env, wp => real64
   use mctc_constants
   use pbc_tools
   implicit none
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   integer, intent(in) :: rep(3)    !< images to consider
   real(wp),intent(in) :: rlat(3,3) !< reciprocal lattice
   real(wp),intent(in) :: vol       !< direct cell volume
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp),intent(out) :: dAmat(3) !< element of interaction matrix
   real(wp),intent(out) :: sigma(3,3)
   integer  :: dx,dy,dz
   real(wp) :: rik2,rik(3)
   real(wp) :: t(3),dtmp,fpivol
   real(wp) :: expterm,arg
   integer  :: i
   dAmat = 0.0_wp
   sigma = 0.0_wp
   fpivol = 4.0_wp*pi/vol
   do concurrent(dx = -rep(1):rep(1), &
         &       dy = -rep(2):rep(2), &
         &       dz = -rep(3):rep(3))
      if (dx==0 .and. dy==0 .and. dz==0) cycle
      t = [dx,dy,dz]
      rik = matmul(rlat,t)
      rik2 = dot_product(rik,rik)
      expterm = exp(-rik2/(4.0_wp*cf**2))/rik2
      arg = dot_product(rik,riw)
      ! d/dx (sin**2 + cos**2) = -2*sin*cos - 2*cos*sin
      dtmp = -sin(arg) * expterm * fpivol
      dAmat = dAmat + rik*dtmp
      sigma = sigma + fpivol * expterm * cos(arg) * (&
         & reshape([-1.0_wp, 0.0_wp, 0.0_wp, &
         &           0.0_wp,-1.0_wp, 0.0_wp, &
         &           0.0_wp, 0.0_wp,-1.0_wp],shape(sigma)) &
         & * (1.0_wp + rik2/(4.0_wp*cf**2)*2.0_wp/3.0_wp) &
         & + (2.0_wp/rik2 + 0.5_wp/cf**2) * outer_prod_3x3(rik,rik))
   end do
end subroutine eeq_ewald_dx_3d_rec

end module eeq_model
