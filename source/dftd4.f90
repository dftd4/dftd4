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

!> provides DFT-D4 dispersion
module dftd4
   use iso_fortran_env, only : wp => real64
   use mctc_constants
!! ========================================================================
!  mix in the covalent coordination number from the ncoord module
!  also get the CN-Parameters to inline the CN-derivative in the gradient
   use coordination_number, only : kn,k1,k4,k5,k6
   use class_param, only : dftd_parameter
   use mctc_param, only : rcov => covalent_radius_d3, &
      &                   en => pauling_en, &
      &                   r4r2 => sqrt_z_r4_over_r2, &
      &                   gam => chemical_hardness, &
      &                   zeff => def2ecp_nuclear_charges
   implicit none

   real(wp) :: thopi,ootpi
   parameter ( thopi = 3._wp/pi )
   parameter ( ootpi = 0.5_wp/pi )

   integer,parameter :: p_refq_gfn2xtb          = 0
   integer,parameter :: p_refq_gasteiger        = 1
   integer,parameter :: p_refq_hirshfeld        = 2
   integer,parameter :: p_refq_periodic         = 3
   integer,parameter :: p_refq_gfn2xtb_gbsa_h2o = 4
   integer,parameter :: p_refq_goedecker        = 5

   integer,parameter :: p_mbd_none       = 0 !< just pair-wise dispersion
   integer,parameter :: p_mbd_rpalike    = 1 !< RPA-like (=MBD) non-additivity
   integer,parameter :: p_mbd_exact_atm  = 2 !< integrate C9 from polarizibilities
   integer,parameter :: p_mbd_approx_atm = 3 !< approximate C9 from C6

   integer, parameter :: max_elem = 118

   integer, dimension(max_elem)      :: refn ! for D4
   real(wp),dimension(7,max_elem)    :: refq
   real(wp),dimension(7,max_elem)    :: refh
   real(wp),dimension(7,max_elem)    :: dftq,pbcq,gffq,solq,clsq
   real(wp),dimension(7,max_elem)    :: dfth,pbch,gffh,solh,clsh
   real(wp),dimension(7,max_elem)    :: hcount
   real(wp),dimension(7,max_elem)    :: ascale
   real(wp),dimension(7,max_elem)    :: refcovcn
   real(wp),dimension(7,max_elem)    :: refcn
   integer, dimension(7,max_elem)    :: refsys
   real(wp),dimension(23,7,max_elem) :: alphaiw
   real(wp),dimension(17)       :: secq
   real(wp),dimension(17)       :: sscale
   real(wp),dimension(17)       :: seccn
   real(wp),dimension(17)       :: seccnd3
   real(wp),dimension(23,17)    :: secaiw

   include 'param_ref.inc'

   type :: dispersion_model
      real(wp) :: g_a = 0.0_wp
      real(wp) :: g_c = 0.0_wp
      integer, dimension(max_elem) :: atoms = 0
      integer, dimension(max_elem) :: nref = 0
      integer, dimension(7,max_elem) :: ncount = 0
      real(wp),dimension(7,max_elem) :: cn = 0.0_wp
      real(wp),dimension(7,max_elem) :: q = 0.0_wp
      real(wp),dimension(23,7,max_elem) :: alpha = 0.0_wp
      real(wp),dimension(7,7,max_elem,max_elem) :: c6 = 0.0_wp
   contains
      procedure :: allocate => allocate_dispersion_model
      procedure :: deallocate => deallocate_dispersion_model
      generic :: new => new_d4_model, new_d3_model
      procedure, private :: new_d4_model => initialize_d4_model
      procedure, private :: new_d3_model => initialize_d3_model
   end type dispersion_model

contains

subroutine allocate_dispersion_model(self)
   class(dispersion_model), intent(inout) :: self
   call self%deallocate
end subroutine allocate_dispersion_model

subroutine deallocate_dispersion_model(self)
   class(dispersion_model), intent(out) :: self
end subroutine deallocate_dispersion_model

subroutine initialize_d4_model(self,atoms,mode,g_a,g_c)
   class(dispersion_model), intent(inout) :: self
   integer, intent(in)  :: atoms(:)
   real(wp),intent(in)  :: g_a,g_c
   integer, intent(in)  :: mode

   integer  :: i,ia,is,icn,j,ii,jj
   integer  :: cncount(0:18)
   real(wp) :: alpha(23),iz,c6
   real(wp) :: tmp_hq(7,max_elem)

   intrinsic :: nint

   call self%allocate

   self%g_a = g_a
   self%g_c = g_c

   select case(mode)
   case(p_refq_hirshfeld,p_refq_periodic)
      self%q = dftq
      tmp_hq = dfth
   case(p_refq_gasteiger)
      self%q = gffq
      tmp_hq = gffh
   case(p_refq_goedecker)
      self%q = clsq
      tmp_hq = clsh
   case(p_refq_gfn2xtb_gbsa_h2o)
      self%q = solq
      tmp_hq = solh
   case default
      self%q = refq
      tmp_hq = refh
   end select

   self%atoms = 0
   self%nref = 0

   !> set up ncount und alpha, also obtain the dimension of the dispmat
   do i = 1, size(atoms,1)
      cncount = 0
      cncount(0) = 1
      ia = atoms(i)
      if (self%atoms(ia).eq.0) then
         self%nref(ia) = refn(ia)
         do j = 1, refn(ia)
            is = refsys(j,ia)
            iz = zeff(is)
            alpha = sscale(is)*secaiw(:,is) &
               &    * zeta(self%g_a,gam(is)*self%g_c,iz,tmp_hq(j,ia)+iz)
            icn = nint(refcn(j,ia))
            self%cn(j,ia) = refcovcn(j,ia)
            cncount(icn) = cncount(icn) + 1
            self%alpha(:,j,ia) = max(ascale(j,ia)*(alphaiw(:,j,ia) &
               &                     - hcount(j,ia)*alpha), 0.0_wp)
         enddo
         do j = 1, refn(ia)
            icn = cncount(nint(refcn(j,ia)))
            self%ncount(j,ia) = icn*(icn+1)/2
         enddo
      endif
      self%atoms(ia) = self%atoms(ia)+1
   enddo

   ! integrate C6 coefficients
   !$omp parallel default(none) private(i,j,ii,jj,alpha,c6) shared(self)
   !$omp do schedule(runtime)
   do i = 1, max_elem
      do j = 1, i
         if (self%atoms(i) > 0 .and. self%atoms(j) > 0) then
            do ii = 1, self%nref(i)
               do jj = 1, self%nref(j)
                  alpha = self%alpha(:,ii,i)*self%alpha(:,jj,j)
                  c6 = thopi * trapzd(alpha)
                  self%c6(jj,ii,j,i) = c6
                  self%c6(ii,jj,i,j) = c6
               enddo
            enddo
         endif
      enddo
   enddo
   !$omp end do
   !$omp end parallel
end subroutine initialize_d4_model

!> This routine is initializing a D3-like dispersion model (notice the absence
!  of DFT- in D3-like) using the polarizibilities from the D4 model
!  we published in JCP2019 and PCCP2019. So if you want to use this
!  D3-like model in the framework of the `dftd4`-program, you have to
!  cite those papers AND the D3(BJ) paper. Also you have to explain that
!  you are evalulating the D4 energy expression with your D3-like model
!  and cite the JCP2019 paper for this (citing it once is okay).
!
!  Note that, it is easier to use the D4 standard model directly and
!  just cite the JCP2019 and PCCP2019 paper.
!
!  Oh, and the D3-like model is still not as good the D4 standard model,
!  I already pushed it to the limit.                                 --- SAW
subroutine initialize_d3_model(self,atoms)
   class(dispersion_model), intent(inout) :: self
   integer, intent(in)  :: atoms(:)

   integer  :: i,ia,is,j,ii,jj
   real(wp) :: alpha(23),c6

   intrinsic :: nint

   call self%allocate

   self%g_a = 0.0_wp
   self%g_c = 0.0_wp
   self%q = 0.0_wp
   self%ncount = 1

   self%atoms = 0
   self%nref = 0

   !> set up ncount und alpha, also obtain the dimension of the dispmat
   do i = 1, size(atoms,1)
      ia = atoms(i)
      if (self%atoms(ia).eq.0) then
         self%nref(ia) = refn(ia)
         do j = 1, refn(ia)
            is = refsys(j,ia)
            alpha = sscale(is)*secaiw(:,is)
            self%cn(j,ia) = refcn(j,ia)
            self%alpha(:,j,ia) = max(ascale(j,ia)*(alphaiw(:,j,ia) &
               &                     - hcount(j,ia)*alpha), 0.0_wp)
         enddo
      endif
      self%atoms(ia) = self%atoms(ia)+1
   enddo

   ! integrate C6 coefficients
   !$omp parallel default(none) private(i,j,ii,jj,alpha,c6) shared(self)
   !$omp do schedule(runtime)
   do i = 1, max_elem
      do j = 1, i
         if (self%atoms(i) > 0 .and. self%atoms(j) > 0) then
            do ii = 1, self%nref(i)
               do jj = 1, self%nref(j)
                  alpha = self%alpha(:,ii,i)*self%alpha(:,jj,j)
                  c6 = thopi * trapzd(alpha)
                  self%c6(jj,ii,j,i) = c6
                  self%c6(ii,jj,i,j) = c6
               enddo
            enddo
         endif
      enddo
   enddo
   !$omp end do
   !$omp end parallel
end subroutine initialize_d3_model

!> prints molecular properties in a human readable format
!
!  molecular polarizibilities and molecular C6/C8 coefficients are
!  printed together with the used partial charge, coordination number
!  atomic C6 and static polarizibilities.
subroutine prmolc6(id,mol,molc6,molc8,molpol,  &
           &       cn,covcn,q,qlmom,c6ab,alpha,rvdw,hvol)
   use mctc_econv
   use class_molecule
   implicit none
   integer, intent(in)  :: id
   type(molecule),intent(in) :: mol !< molecular structure information
   real(wp),intent(in)  :: molc6    !< molecular C6 coefficient in au
   real(wp),intent(in)  :: molc8    !< molecular C8 coefficient in au
   real(wp),intent(in)  :: molpol   !< molecular static dipole polarizibility
   real(wp),intent(in),optional :: cn(mol%n)
   real(wp),intent(in),optional :: covcn(mol%n)
   real(wp),intent(in),optional :: q(mol%n)
   real(wp),intent(in),optional :: qlmom(3,mol%n)
   real(wp),intent(in),optional :: c6ab(mol%n,mol%n)
   real(wp),intent(in),optional :: alpha(mol%n)
   real(wp),intent(in),optional :: rvdw(mol%n)
   real(wp),intent(in),optional :: hvol(mol%n)
   integer :: i
   if(present(cn).or.present(covcn).or.present(q).or.present(c6ab) &
   &   .or.present(alpha).or.present(rvdw).or.present(hvol)) then
   write(id,'(a)')
   write(id,'(7x,''   #   Z   '')',advance='no')
   if(present(cn))   write(id,'(''        CN'')',advance='no')
   if(present(covcn))write(id,'(''     covCN'')',advance='no')
   if(present(q))    write(id,'(''         q'')',advance='no')
   if(present(qlmom))write(id,   '(''   n(s)'')',advance='no')
   if(present(qlmom))write(id,   '(''   n(p)'')',advance='no')
   if(present(qlmom))write(id,   '(''   n(d)'')',advance='no')
   if(present(c6ab)) write(id,'(''      C6AA'')',advance='no')
   if(present(alpha))write(id,'(''      α(0)'')',advance='no')
   if(present(rvdw)) write(id,'(''    RvdW/Å'')',advance='no')
   if(present(hvol)) write(id,'(''    relVol'')',advance='no')
   write(id,'(a)')
   do i=1,mol%n
      write(id,'(i11,1x,i3,1x,a2)',advance='no') &
      &     i,mol%at(i),mol%sym(i)
      if(present(cn))   write(id,'(f10.3)',advance='no')cn(i)
      if(present(covcn))write(id,'(f10.3)',advance='no')covcn(i)
      if(present(q))    write(id,'(f10.3)',advance='no')q(i)
      if(present(qlmom))write(id, '(f7.3)',advance='no')qlmom(1,i)
      if(present(qlmom))write(id, '(f7.3)',advance='no')qlmom(2,i)
      if(present(qlmom))write(id, '(f7.3)',advance='no')qlmom(3,i)
      if(present(c6ab)) write(id,'(f10.3)',advance='no')c6ab(i,i)
      if(present(alpha))write(id,'(f10.3)',advance='no')alpha(i)
      if(present(rvdw)) write(id,'(f10.3)',advance='no')rvdw(i)*autoaa
      if(present(hvol)) write(id,'(f10.3)',advance='no')hvol(i)
      write(id,'(a)')
   enddo
   endif
   write(id,'(/,12x,''Mol. C6AA /au·bohr⁶  :'',f18.6,'// &
   &         '/,12x,''Mol. C8AA /au·bohr⁸  :'',f18.6,'// &
   &         '/,12x,''Mol. α(0) /au        :'',f18.6,/)') &
   &          molc6,molc8,molpol
end subroutine prmolc6

!> calculate molecular dispersion properties
!
!  calculates molecular C6/C8 coefficients and molecular static
!  polarizibility, optionally return dipole polarizibilities
!  partitioned to atoms, C6AA coefficients for each atom,
!  radii derived from the polarizibilities and relative
!  volumes relative to the atom
subroutine mdisp(mol,dispm,ndim,q, &
           &     gw,refc6,molc6,molc8,molpol,aout,cout,rout,vout)
   use class_molecule
   implicit none
   type(dispersion_model),intent(in) :: dispm
   type(molecule),intent(in) :: mol   !< molecular structure information
   integer, intent(in)  :: ndim       !< dimension of reference systems
   real(wp),intent(in)  :: q(mol%n) !< partial charges
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: refc6(ndim,ndim)
   real(wp),intent(out) :: molc6  !< molecular C6 coefficient in au
   real(wp),intent(out) :: molc8  !< molecular C8 coefficient in au
   real(wp),intent(out) :: molpol !< molecular static dipole polarizibility
   real(wp),intent(out),optional :: aout(23,mol%n)
   real(wp),intent(out),optional :: cout(mol%n,mol%n)
   real(wp),intent(out),optional :: rout(mol%n)
   real(wp),intent(out),optional :: vout(mol%n)

   integer  :: i,ii,ia,j,ja,k
   integer, allocatable :: itbl(:,:)
   real(wp) :: oth,iz
   real(wp),allocatable :: zetvec(:)
   real(wp),allocatable :: rvdw(:)
   real(wp),allocatable :: phv(:)
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: c8ab(:,:)
   real(wp),allocatable :: aw(:,:)
   parameter (oth=1._wp/3._wp)

   allocate( zetvec(ndim),rvdw(mol%n),phv(mol%n),c6ab(mol%n,mol%n), &
   &         c8ab(mol%n,mol%n), aw(23,mol%n),  source = 0.0_wp )
   allocate( itbl(7,mol%n), source = 0 )

   molc6  = 0._wp
   molc8  = 0._wp
   molpol = 0._wp

   k = 0
   do i = 1, mol%n
      ia = mol%at(i)
      iz = zeff(ia)
      do ii = 1, dispm%nref(ia)
         k = k+1
         itbl(ii,i) = k
         zetvec(k) = gw(k) * zeta(dispm%g_a,gam(ia)*dispm%g_c,dispm%q(ii,ia)+iz,q(i)+iz)
         aw(:,i) = aw(:,i) + zetvec(k) * dispm%alpha(:,ii,ia)
      enddo
!     van-der-Waals radius, alpha = 4/3 pi r**3 <=> r = (3/(4pi) alpha)**(1/3)
      rvdw(i) = (0.25_wp*thopi*aw(1,i))**oth
!     pseudo-hirshfeld volume
      phv(i) = aw(1,i)/dispm%alpha(1,1,ia)
      c6ab(i,i) = thopi * trapzd(aw(:,i)**2)
      c8ab(i,i) = 3*r4r2(ia)**2*c6ab(i,i)
      molpol = molpol + aw(1,i)
      molc6  = molc6  + c6ab(i,i)
      molc8 = molc8 + 3*r4r2(ia)**2*c6ab(i,i)
      do j = 1, i-1
         ja = mol%at(j)
         c6ab(j,i) = thopi * trapzd(aw(:,i)*aw(:,j))
         c6ab(i,j) = c6ab(j,i)
         c8ab(j,i) = 3*r4r2(ia)*r4r2(ja)*c6ab(j,i)
         c8ab(i,j) = c8ab(j,i)
         molc6 = molc6 + 2*c6ab(j,i)
         molc8 = molc8 + 6*r4r2(ia)*r4r2(ja)*c6ab(j,i)
      enddo
   enddo

   if (present(aout)) aout = aw
   if (present(vout)) vout = phv
   if (present(rout)) rout = rvdw
   if (present(cout)) cout = c6ab

end subroutine mdisp

subroutine prd4ref(iunit,mol,dispm)
   use class_molecule
   implicit none
   type(dispersion_model),intent(in) :: dispm
   integer, intent(in) :: iunit
   type(molecule),intent(in) :: mol

   integer :: i,ii,ia
   logical :: printed(118)

   printed = .false.

   write(iunit,'(a)')
   do i = 1, mol%n
      ia = mol%at(i)
      if (printed(ia)) cycle
      write(iunit,'(a3,1x,a10,1x,a14,1x,a15)') mol%sym(i), &
         &   'q(ref)','cov. CN(ref)','α(AIM,ref)'
      do ii = 1, dispm%nref(ia)
         write(iunit,'(4x,f10.7,5x,f10.7,5x,f10.7)') &
            &      dispm%q(ii,ia),dispm%cn(ii,ia),dispm%alpha(1,ii,ia)
      enddo
      printed(ia) = .true.
   enddo
   write(iunit,'(a)')
end subroutine prd4ref

!> charge scaling function
pure elemental function zeta(a,c,qref,qmod)
   implicit none
   real(wp),intent(in) :: qmod,qref
   real(wp),intent(in) :: a,c
   real(wp)            :: zeta

   intrinsic :: exp

   if (qmod.lt.0._wp) then
      zeta = exp( a )
   else
      zeta = exp( a * ( 1._wp - exp( c * ( 1._wp - qref/qmod ) ) ) )
   endif

end function zeta

!> derivative of charge scaling function w.r.t. charge
pure elemental function dzeta(a,c,qref,qmod)
   implicit none
   real(wp),intent(in) :: qmod,qref
   real(wp),intent(in) :: a,c
   real(wp)            :: dzeta

   intrinsic :: exp

   if (qmod.lt.0._wp) then
      dzeta = 0._wp
   else
      dzeta = - a * c * exp( c * ( 1._wp - qref/qmod ) ) &
      &           * zeta(a,c,qref,qmod) * qref / ( qmod**2 )
   endif

end function dzeta

!> numerical Casimir--Polder integration
pure function trapzd(pol)
   implicit none
   real(wp),intent(in) :: pol(23)
   real(wp)            :: trapzd

   real(wp),parameter  :: freq(23) = (/ &
&   0.000001_wp,0.050000_wp,0.100000_wp, &
&   0.200000_wp,0.300000_wp,0.400000_wp, &
&   0.500000_wp,0.600000_wp,0.700000_wp, &
&   0.800000_wp,0.900000_wp,1.000000_wp, &
&   1.200000_wp,1.400000_wp,1.600000_wp, &
&   1.800000_wp,2.000000_wp,2.500000_wp, &
&   3.000000_wp,4.000000_wp,5.000000_wp, &
&   7.500000_wp,10.00000_wp /)
!  just precalculate all weights and get the job done
   real(wp),parameter :: weights(23) = 0.5_wp * (/ &
&  ( freq (2) - freq (1) ),  &
&  ( freq (2) - freq (1) ) + ( freq (3) - freq (2) ),  &
&  ( freq (3) - freq (2) ) + ( freq (4) - freq (3) ),  &
&  ( freq (4) - freq (3) ) + ( freq (5) - freq (4) ),  &
&  ( freq (5) - freq (4) ) + ( freq (6) - freq (5) ),  &
&  ( freq (6) - freq (5) ) + ( freq (7) - freq (6) ),  &
&  ( freq (7) - freq (6) ) + ( freq (8) - freq (7) ),  &
&  ( freq (8) - freq (7) ) + ( freq (9) - freq (8) ),  &
&  ( freq (9) - freq (8) ) + ( freq(10) - freq (9) ),  &
&  ( freq(10) - freq (9) ) + ( freq(11) - freq(10) ),  &
&  ( freq(11) - freq(10) ) + ( freq(12) - freq(11) ),  &
&  ( freq(12) - freq(11) ) + ( freq(13) - freq(12) ),  &
&  ( freq(13) - freq(12) ) + ( freq(14) - freq(13) ),  &
&  ( freq(14) - freq(13) ) + ( freq(15) - freq(14) ),  &
&  ( freq(15) - freq(14) ) + ( freq(16) - freq(15) ),  &
&  ( freq(16) - freq(15) ) + ( freq(17) - freq(16) ),  &
&  ( freq(17) - freq(16) ) + ( freq(18) - freq(17) ),  &
&  ( freq(18) - freq(17) ) + ( freq(19) - freq(18) ),  &
&  ( freq(19) - freq(18) ) + ( freq(20) - freq(19) ),  &
&  ( freq(20) - freq(19) ) + ( freq(21) - freq(20) ),  &
&  ( freq(21) - freq(20) ) + ( freq(22) - freq(21) ),  &
&  ( freq(22) - freq(21) ) + ( freq(23) - freq(22) ),  &
&  ( freq(23) - freq(22) ) /)

!!  do average between trap(1)-trap(22) .and. trap(2)-trap(23)
!   tmp1 = 0.5_wp * ( &
!&  ( freq (2) - freq (1) ) * ( pol (2) + pol (1) )+ &
!&  ( freq (4) - freq (3) ) * ( pol (4) + pol (3) )+ &
!&  ( freq (6) - freq (5) ) * ( pol (6) + pol (5) )+ &
!&  ( freq (8) - freq (7) ) * ( pol (8) + pol (7) )+ &
!&  ( freq(10) - freq (9) ) * ( pol(10) + pol (9) )+ &
!&  ( freq(12) - freq(11) ) * ( pol(12) + pol(11) )+ &
!&  ( freq(14) - freq(13) ) * ( pol(14) + pol(13) )+ &
!&  ( freq(16) - freq(15) ) * ( pol(16) + pol(15) )+ &
!&  ( freq(18) - freq(17) ) * ( pol(18) + pol(17) )+ &
!&  ( freq(20) - freq(19) ) * ( pol(20) + pol(19) )+ &
!&  ( freq(22) - freq(21) ) * ( pol(22) + pol(21) ))
!   tmp2 = 0.5_wp * ( &
!&  ( freq (3) - freq (2) ) * ( pol (3) + pol (2) )+ &
!&  ( freq (5) - freq (4) ) * ( pol (5) + pol (4) )+ &
!&  ( freq (7) - freq (6) ) * ( pol (7) + pol (6) )+ &
!&  ( freq (9) - freq (8) ) * ( pol (9) + pol (8) )+ &
!&  ( freq(11) - freq(10) ) * ( pol(11) + pol(10) )+ &
!&  ( freq(13) - freq(12) ) * ( pol(13) + pol(12) )+ &
!&  ( freq(15) - freq(14) ) * ( pol(15) + pol(14) )+ &
!&  ( freq(17) - freq(16) ) * ( pol(17) + pol(16) )+ &
!&  ( freq(19) - freq(18) ) * ( pol(19) + pol(18) )+ &
!&  ( freq(21) - freq(20) ) * ( pol(21) + pol(20) )+ &
!&  ( freq(23) - freq(22) ) * ( pol(23) + pol(22) ))

   trapzd = sum(pol*weights)

end function trapzd

!> @brief coordination number gaussian weight
pure elemental function cngw(wf,cn,cnref)
   implicit none
   real(wp),intent(in) :: wf,cn,cnref
   real(wp)            :: cngw ! CN-gaussian-weight

   intrinsic :: exp

   cngw = exp ( -wf * ( cn - cnref )**2 )

end function cngw

!> derivative of gaussian weight w.r.t. coordination number
pure elemental function dcngw(wf,cn,cnref)
!  use dftd4, only : cngw
   implicit none
   real(wp),intent(in) :: wf,cn,cnref
   real(wp) :: dcngw

   dcngw = 2*wf*(cnref-cn)*cngw(wf,cn,cnref)

end function dcngw

!> BJ damping function ala DFT-D3(BJ)
!
!  f(n,rab) = sn*rab**n/(rab**n + R0**n)  w/ R0 = a1*sqrt(C6/C8)+a2
!  see: https://doi.org/10.1002/jcc.21759
pure elemental function fdmpr_bj(n,r,c) result(fdmp)
   implicit none
   integer, intent(in)  :: n  !< order
   real(wp),intent(in)  :: r  !< distance
   real(wp),intent(in)  :: c  !< critical radius
   real(wp) :: fdmp
   fdmp = 1.0_wp / ( r**n + c**n )
end function fdmpr_bj
!> derivative of BJ damping function ala DFT-D3(BJ)
pure elemental function fdmprdr_bj(n,r,c) result(dfdmp)
   implicit none
   integer, intent(in)  :: n  !< order
   real(wp),intent(in)  :: r  !< distance
   real(wp),intent(in)  :: c  !< critical radius
   real(wp) :: dfdmp
   dfdmp = -n*r**(n-1) * fdmpr_bj(n,r,c)**2
end function fdmprdr_bj

!> original DFT-D3(0) damping
!
!  f(n,rab) = sn/(1+6*(4/3*R0/rab)**alp)  w/ R0 of unknown origin
pure elemental function fdmpr_zero(n,r,c,alp) result(fdmp)
   implicit none
   integer, intent(in)  :: n   !< order
   real(wp),intent(in)  :: r   !< distance
   real(wp),intent(in)  :: c   !< critical radius
   integer, intent(in)  :: alp !< exponent
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: fdmp
   fdmp = 1.0_wp / (r**n*(1 + six * (c/r)**(n+alp)))
end function fdmpr_zero
!> derivative of original DFT-D3(0) damping
pure elemental function fdmprdr_zero(n,r,c,alp) result(dfdmp)
   implicit none
   integer, intent(in)  :: n   !< order
   real(wp),intent(in)  :: r   !< distance
   real(wp),intent(in)  :: c   !< critical radius
   integer, intent(in)  :: alp !< exponent
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: dfdmp
   dfdmp = -( n*r**(n-1)*(1+six*(c/r)**(alp)) &
             - alp*r**n/c*six*(c/r)**(alp-1) ) &
           * fdmpr_zero(n,r,c,alp)**2
!  fdmp = 1.0_wp / (r**n*(1 + 6.0_wp * (c/r)**(n+alp)))
end function fdmprdr_zero

!> fermi damping function from TS and MBD methods
!
!  f(n,rab) = sn/(1+exp[-alp*(rab/R0-1)]) w/ R0 as experimenal vdW-Radii
pure elemental function fdmpr_fermi(n,r,c,alp) result(fdmp)
   implicit none
   integer, intent(in)  :: n   !< order
   real(wp),intent(in)  :: r   !< distance
   real(wp),intent(in)  :: c   !< critical radius
   integer, intent(in)  :: alp !< steepness
   real(wp) :: fdmp
   fdmp = 1.0_wp / (r**n*(1.0_wp+exp(-alp*(r/c - 1.0))))
end function fdmpr_fermi
!> derivative of fermi damping function from TS and MBD methods
pure elemental function fdmprdr_fermi(n,r,c,alp) result(dfdmp)
   implicit none
   integer, intent(in)  :: n   !< order
   real(wp),intent(in)  :: r   !< distance
   real(wp),intent(in)  :: c   !< critical radius
   integer, intent(in)  :: alp !< steepness
   real(wp) :: dfdmp
   dfdmp = -(-alp/c*r**n*exp(-alp*(r/c - 1.0)) &
             + n*r**(n-1)*(1.0_wp+exp(-alp*(r/c - 1.0)))) &
             * fdmpr_fermi(n,r,c,alp)**2
end function fdmprdr_fermi

!> optimized power zero damping (M. Head-Gordon)
!
!  f(n,rab) = sn*rab**(n+alp)/(rab**(n+alp) + R0**(n+alp))
!  see: https://dx.doi.org/10.1021/acs.jpclett.7b00176
pure elemental function fdmpr_op(n,r,c,alp) result(fdmp)
   implicit none
   integer, intent(in)  :: n   !< order
   real(wp),intent(in)  :: r   !< distance
   real(wp),intent(in)  :: c   !< critical radius
   integer, intent(in)  :: alp !< optimized power
   real(wp) :: fdmp
   fdmp = r**alp / (r**(n+alp)*c**(n+alp))
end function fdmpr_op
!> derivative optimized power zero damping (M. Head-Gordon)
pure elemental function fdmprdr_op(n,r,c,alp) result(dfdmp)
   implicit none
   integer, intent(in)  :: n   !< order
   real(wp),intent(in)  :: r   !< distance
   real(wp),intent(in)  :: c   !< critical radius
   integer, intent(in)  :: alp !< optimized power
   real(wp) :: dfdmp
   dfdmp = (alp*r*(alp-1) - (n+alp)*r**alp*r**(n+alp-1)) &
           * fdmpr_op(n,r,c,alp)**2
!  fdmp = r**alp / (r**(n+alp)*c**(n+alp))
end function fdmprdr_op

!> Sherrill's M-zero damping function
!
!  f(n,rab) = sn/(1+6*(4/3*R0/rab+a2*R0)**(-alp))
!  see: https://dx.doi.org/10.1021/acs.jpclett.6b00780
pure elemental function fdmpr_zerom(n,r,c,rsn,alp) result(fdmp)
   implicit none
   integer, intent(in)  :: n   !< order
   real(wp),intent(in)  :: r   !< distance
   real(wp),intent(in)  :: c   !< critical radius
   real(wp),intent(in)  :: rsn !< offset for critical radius
   integer, intent(in)  :: alp !< exponent
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: fdmp
   fdmp = 1.0_wp / (r**n*(1 + six * (r/c+rsn*c)**(-alp)))
end function fdmpr_zerom
!> derivative of Sherrill's M-zero damping function
pure elemental function fdmprdr_zerom(n,r,c,rsn,alp) result(dfdmp)
   implicit none
   integer, intent(in)  :: n   !< order
   real(wp),intent(in)  :: r   !< distance
   real(wp),intent(in)  :: c   !< critical radius
   real(wp),intent(in)  :: rsn !< offset for critical radius
   integer, intent(in)  :: alp !< exponent
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: dfdmp
   dfdmp = -( n*r**(n-1)*(1+six*(r/c+rsn*c)**(-alp)) &
              - alp*r**n/c*six*(r/c+rsn*c)**(-alp-1) ) &
           * fdmpr_zerom(n,r,c,rsn,alp)**2
end function fdmprdr_zerom

!> basic D4 calculation
!
!  obtain gaussian weights for the references systems based on
!  input CN and integrate reference C6 coefficients
subroutine d4(mol,dispm,ndim,wf,covcn,gw,refc6)
   use class_molecule
   implicit none
   type(dispersion_model),intent(in) :: dispm
   type(molecule),intent(in) :: mol !< molecular structure information
   integer, intent(in)  :: ndim     !< calculation dimension
   real(wp),intent(in)  :: wf
   real(wp),intent(in)  :: covcn(mol%n)  !< covalent coordination number
   real(wp),intent(out) :: gw(ndim)
   real(wp),intent(out) :: refc6(ndim,ndim)

   integer  :: i,ia,ii,iii,j,jj,ja,k,l
   integer,allocatable :: itbl(:,:)
   real(wp) :: twf,norm,aiw(23)

   intrinsic :: maxval

   allocate( itbl(7,mol%n), source = 0 )

   gw = 0._wp
   refc6 = 0._wp

   k = 0
   do i = 1, mol%n
      do ii = 1, dispm%nref(mol%at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, mol%n
      ia = mol%at(i)
      norm = 0.0_wp
      do ii = 1, dispm%nref(ia)
         do iii = 1, dispm%ncount(ii,ia)
            twf = iii*wf
            norm = norm + cngw(twf,covcn(i),dispm%cn(ii,ia))
         enddo
      enddo
      norm = 1._wp / norm
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         do iii = 1, dispm%ncount(ii,ia)
            twf = iii*wf
            gw(k) = gw(k) + cngw(twf,covcn(i),dispm%cn(ii,ia)) * norm
         enddo
!    --- okay, if we run out of numerical precision, gw(k) will be NaN.
!        In case it is NaN, it will not match itself! So we can rescue
!        this exception. This can only happen for very high CNs.
         if (gw(k).ne.gw(k)) then
            if (maxval(dispm%cn(:dispm%nref(ia),ia)).eq.dispm%cn(ii,ia)) then
               gw(k) = 1.0_wp
            else
               gw(k) = 0.0_wp
            endif
         endif
         ! diagonal terms for image-cell interaction under PBC
         do jj = 1, ii
            l = itbl(jj,i)
            aiw = dispm%alpha(:,ii,ia)*dispm%alpha(:,jj,ia)
            refc6(l,k) = dispm%c6(jj,ii,ia,ia)
            refc6(k,l) = refc6(l,k)
         enddo
         ! off-diagonal terms
         do j = 1, i-1
            ja = mol%at(j)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               aiw = dispm%alpha(:,ii,ia)*dispm%alpha(:,jj,ja)
               refc6(l,k) = dispm%c6(jj,ii,ja,ia)
               refc6(k,l) = refc6(l,k)
            enddo
         enddo
      enddo
   enddo

end subroutine d4

!> build a dispersion matrix (mainly for SCC calculations)
pure subroutine build_dispmat(mol,dispm,ndim,par,refc6,dispmat)
   use class_molecule
   implicit none
   type(dispersion_model),intent(in) :: dispm
   type(molecule),intent(in) :: mol !< molecular structure information
   integer, intent(in)  :: ndim
   type(dftd_parameter),intent(in)  :: par
   real(wp),intent(in)  :: refc6(ndim,ndim)
   real(wp),intent(out) :: dispmat(ndim,ndim)

   integer  :: i,ii,ia,j,jj,ja,k,l
   integer, allocatable :: itbl(:,:)
   real(wp) :: refc8,refc10,r,oor6,oor8,oor10,cutoff
   real(wp), parameter :: rthr = 72.0_wp ! slightly larger than in gradient

   allocate( itbl(7,mol%n), source = 0 )

   dispmat = 0.0_wp

   k = 0
   do i = 1, mol%n
      do ii = 1, dispm%nref(mol%at(i))
         k = k + 1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, mol%n
      ia = mol%at(i)
      do j = 1, i-1
         ja = mol%at(j)
         cutoff = par%a1*sqrt(3._wp*r4r2(mol%at(i))*r4r2(mol%at(j)))+par%a2
         r = norm2(mol%xyz(:,j)-mol%xyz(:,i))
         if (r.gt.rthr) cycle
!        oor6  = 1.0_wp/(r**6  + cutoff**6 )
!        oor8  = 1.0_wp/(r**8  + cutoff**8 )
!        oor10 = 1.0_wp/(r**10 + cutoff**10)
         oor6  = fdmpr_bj( 6,r,cutoff)
         oor8  = fdmpr_bj( 8,r,cutoff)
         oor10 = fdmpr_bj(10,r,cutoff)
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               refc8 = 3.0_wp * r4r2(ia)*r4r2(ja) * refc6(k,l)
               refc10 = 49.0_wp/40.0_wp * refc8**2/refc6(k,l)
               dispmat(k,l) = &
               &  - par%s6 * ( refc6(k,l) * oor6 ) &
               &  - par%s8 * ( refc8      * oor8 ) &
               &  - par%s8 * ( refc10     * oor8 )
               dispmat(l,k) = dispmat(k,l)
            enddo
         enddo
      enddo
   enddo

end subroutine build_dispmat

!> build a weighted dispersion matrix (mainly for SCC calculations)
subroutine build_wdispmat(mol,dispm,ndim,par,refc6,gw,wdispmat)
   use class_molecule
   implicit none
   type(dispersion_model),intent(in) :: dispm
   type(molecule),intent(in) :: mol !< molecular structure information
   integer, intent(in)  :: ndim
   type(dftd_parameter),intent(in)  :: par
   real(wp),intent(in)  :: refc6(ndim,ndim)
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(out) :: wdispmat(ndim,ndim)

   integer :: i,ii,ia,j,jj,ja,k,l
   integer, allocatable :: itbl(:,:)
   real(wp) :: refc8,refc10,cutoff,oor6,oor8,oor10,r,gwgw,r4r2ij
   real(wp), parameter :: rthr = 72.0_wp ! slightly larger than in gradient
   real(wp), parameter :: gwcut = 1.0e-7_wp

   allocate( itbl(7,mol%n), source = 0 )

   wdispmat = 0.0_wp

   k = 0
   do i = 1, mol%n
      do ii = 1, dispm%nref(mol%at(i))
         k = k + 1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, mol%n
      ia = mol%at(i)
      do j = 1, i-1
         ja = mol%at(j)
         r4r2ij = 3.0_wp*r4r2(ia)*r4r2(ja)
         cutoff = par%a1*sqrt(r4r2ij)+par%a2
!        r2 = sum( (mol%xyz(:,j)-mol%xyz(:,i))**2 )
!        oor6  = 1.0_wp/(r2**3 + cutoff**6 )
!        oor8  = 1.0_wp/(r2**4 + cutoff**8 )
!        oor10 = 1.0_wp/(r2**5 + cutoff**10)
         r = norm2(mol%xyz(:,j)-mol%xyz(:,i))
         if (r.gt.rthr) cycle
         oor6  = fdmpr_bj( 6,r,cutoff)
         oor8  = fdmpr_bj( 8,r,cutoff)
         oor10 = fdmpr_bj(10,r,cutoff)
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               gwgw = gw(k)*gw(l)
               if (gwgw.lt.gwcut) cycle
               refc8  = r4r2ij * refc6(k,l)
               refc10 = 49.0_wp/40.0_wp * r4r2ij**2 * refc6(k,l)
               wdispmat(k,l) = gw(k)*gw(l) * ( &
               &  - par%s6  * ( refc6(k,l)  * oor6 ) &
               &  - par%s8  * ( refc8       * oor8 ) &
               &  - par%s10 * ( refc10      * oor10) )
               wdispmat(l,k) = wdispmat(k,l)
               if (abs(refc6(k,l)).lt.0.1_wp) then
               print*,ia,ja,refc6(k,l),i,ii,j,jj
               endif
            enddo
         enddo
      enddo
   enddo

end subroutine build_wdispmat

!> calculate contribution to the Fockian
subroutine disppot(mol,dispm,ndim,q,wdispmat,gw,hdisp)
   use class_molecule
   implicit none
   type(dispersion_model),intent(in) :: dispm
   type(molecule),intent(in) :: mol !< molecular structure information
   integer, intent(in)  :: ndim
   real(wp),intent(in)  :: q(mol%n)
   real(wp),intent(in)  :: wdispmat(ndim,ndim)
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(out) :: hdisp(mol%n)

   integer  :: i,ii,k,ia
   real(wp) :: iz
   real(wp),parameter   :: gw_cut = 1.0e-7_wp
   real(wp),allocatable :: zetavec(:)
   real(wp),allocatable :: zerovec(:)
   real(wp),allocatable :: dumvec(:)

   intrinsic :: sum,dble

   allocate( zetavec(ndim),zerovec(ndim),dumvec(ndim), source = 0._wp )

   zetavec = 0.0_wp
   zerovec = 0.0_wp
   dumvec  = 0.0_wp
   hdisp   = 0.0_wp

   k = 0
   do i = 1, mol%n
       ia = mol%at(i)
       iz = zeff(ia)
       do ii = 1, dispm%nref(ia)
          k = k + 1
          if (gw(k).lt.gw_cut) cycle
          zerovec(k) = dzeta(dispm%g_a,gam(ia)*dispm%g_c,dispm%q(ii,ia)+iz,q(i)+iz)
          zetavec(k) =  zeta(dispm%g_a,gam(ia)*dispm%g_c,dispm%q(ii,ia)+iz,q(i)+iz)
      enddo
   enddo
!  create vector -> dispmat(ndim,dnim) * zetavec(ndim) = dumvec(ndim)
   call dsymv('U',ndim,1._wp,wdispmat,ndim,zetavec,1,0._wp,dumvec,1)
!  call dgemv('N',ndim,ndim,1._wp,wdispmat,ndim,zetavec, &
!  &     1,0._wp,dumvec,1)
!  get atomic reference contributions
   k = 0
   do i = 1, mol%n
      ia = mol%at(i)
      hdisp(i) = sum(dumvec(k+1:k+dispm%nref(ia))*zerovec(k+1:k+dispm%nref(ia)))
      k = k + dispm%nref(ia)
   enddo

   deallocate(zetavec,zerovec,dumvec)

end subroutine disppot

!> calculate dispersion energy in SCC
function edisp_scc(mol,dispm,ndim,q,wdispmat,gw) result(ed)
   use class_molecule
   implicit none
   type(dispersion_model),intent(in) :: dispm
   type(molecule),intent(in) :: mol !< molecular structure information
   integer, intent(in)  :: ndim
   real(wp),intent(in)  :: q(mol%n)
   real(wp),intent(in)  :: wdispmat(ndim,ndim)
   real(wp),intent(in)  :: gw(ndim)
   real(wp) :: ed

   integer  :: i,ii,k,ia
   real(wp) :: iz
   real(wp),parameter   :: gw_cut = 1.0e-7_wp
   real(wp),allocatable :: zetavec(:)
   real(wp),allocatable :: dumvec(:)

   intrinsic :: sum,dble

   allocate( zetavec(ndim),dumvec(ndim), source = 0._wp )

   ed = 0.0_wp

   k = 0
   do i = 1, mol%n
       ia = mol%at(i)
       iz = zeff(ia)
       do ii = 1, dispm%nref(ia)
          k = k + 1
          if (gw(k).lt.gw_cut) cycle
          zetavec(k) =  zeta(dispm%g_a,gam(ia)*dispm%g_c,dispm%q(ii,ia)+iz,q(i)+iz)
      enddo
   enddo
!  create vector -> dispmat(ndim,dnim) * zetavec(ndim) = dumvec(ndim)
   call dsymv('U',ndim,0.5_wp,wdispmat,ndim,zetavec,1,0.0_wp,dumvec,1)
!  call dgemv('N',ndim,ndim,0.5_wp,wdispmat,ndim,zetavec, &
!  &           1,0.0_wp,dumvec,1)
   ed = dot_product(dumvec,zetavec)

   deallocate(zetavec,dumvec)

end function edisp_scc

!> calculate non additivity by RPA-like scheme
subroutine dispmb(mol,E,aw,oor6ab)
   use class_molecule
   implicit none
   type(molecule),intent(in) :: mol !< molecular structure information
   real(wp),intent(in)  :: aw(23,mol%n)
   real(wp),intent(in)  :: oor6ab(mol%n,mol%n)
   real(wp),intent(out) :: E

   integer  :: i,j,ii,jj,k
   integer  :: info
   real(wp) :: tau(3,3),spur(23),d_,d2,r(3),r2,alpha
   real(wp) :: two(23),atm(23),d3
   real(wp),allocatable :: T (:,:)
   real(wp),allocatable :: A (:,:)
   real(wp),allocatable :: AT(:,:)
   real(wp),allocatable :: F (:,:)
   real(wp),allocatable :: F_(:,:)
   real(wp),allocatable :: d (:)
   real(wp),allocatable :: w (:)

   intrinsic :: sum,sqrt,minval,log

   allocate( T(3*mol%n,3*mol%n),  A(3*mol%n,3*mol%n), AT(3*mol%n,3*mol%n), &
   &         F(3*mol%n,3*mol%n), F_(3*mol%n,3*mol%n),  d(3*mol%n), &
   &         w(12*mol%n), &
   &         source = 0.0_wp )

   spur = 0.0_wp

   do i = 1, 3*mol%n
      F(i,i) = 1.0_wp
   enddo

   do i = 1, mol%n
      do j  = 1, i-1
         r  = mol%xyz(:,j) - mol%xyz(:,i)
         r2 = sum(r**2)
         do ii = 1, 3
            tau(ii,ii) = (3*r(ii)*r(ii)-r2)/r2
            do jj = ii+1, 3
               tau(ii,jj) = (3*r(ii)*r(jj))/r2
               tau(jj,ii) = tau(ii,jj)
            enddo
         enddo
         tau = tau*sqrt(oor6ab(i,j))
         T(3*i-2:3*i,3*j-2:3*j) = tau
         T(3*j-2:3*j,3*i-2:3*i) = tau
      enddo
   enddo

   !call prmat(6,T,3*mol%n,3*mol%n,'T')

   do k = 1, 23
      A = 0.0_wp
      do i =  1, mol%n
         alpha = sqrt(aw(k,i))
         A(3*i-2,3*i-2) = alpha
         A(3*i-1,3*i-1) = alpha
         A(3*i  ,3*i  ) = alpha
      enddo

      AT  = 0.0_wp
      call dgemm('N','N',3*mol%n,3*mol%n,3*mol%n,1.0_wp,A,3*mol%n,T, &
  &             3*mol%n,0.0_wp,F_,3*mol%n)
      call dgemm('N','N',3*mol%n,3*mol%n,3*mol%n,1.0_wp,F_,3*mol%n,A, &
  &             3*mol%n,0.0_wp,AT,3*mol%n)

      F_ = F - AT

      d = 0.0_wp
      call dsyev('N','U',3*mol%n,F_,3*mol%n,d,w,12*mol%n,info)
      if (info.ne.0) then
         E = 0.0_wp
         return
      endif
      if (minval(d).le.0.0_wp) then
         E = 0.0_wp
         return
      endif

      call dgemm('N','N',3*mol%n,3*mol%n,3*mol%n,1.0_wp,AT,3*mol%n,AT, &
  &             3*mol%n,0.0_wp,F_,3*mol%n)
!     call dgemm('N','N',3*mol%n,3*mol%n,3*mol%n,1.0_wp,F_,3*mol%n,AT, &
! &             3*mol%n,0.0_wp,A,3*mol%n)

      d_ = 1.0_wp; d2 = 0.0_wp!; d3 = 0.0_wp
      do i = 1, 3*mol%n
         d_ = d_ * d(i)
         d2 = d2 - F_(i,i)
!        d3 = d3 - A(i,i)
      enddo
      spur(k) = log(d_) - d2*0.5
!     two(k) = d2/2.0_wp
!     atm(k) = d3/3.0_wp
   enddo

   E = trapzd(spur)*ooTPI
   !print*,'     full contribution', trapzd(spur)*ooTPI
   !print*,' manybody contribution', trapzd(spur-two)*ooTPI
   !print*,'  twobody contribution', trapzd(two)*ootpi
   !print*,'threebody contribution', trapzd(atm)*ootpi

   deallocate(T,A,AT,F,F_,d)
end subroutine dispmb

function refq2string(ref) result(string)
   implicit none
   integer,intent(in) :: ref
   character(len=:),allocatable :: string
   select case(ref)
   case default;                  string = 'unknown'
   case(p_refq_gfn2xtb);          string = 'GFN2'
   case(p_refq_gasteiger);        string = 'EEQ'
   case(p_refq_hirshfeld);        string = 'extern'
   case(p_refq_periodic);         string = 'EEQ'
   case(p_refq_gfn2xtb_gbsa_h2o); string = 'GFN2/GBSA'
   case(p_refq_goedecker);        string = 'EEQ'
   end select
end function refq2string

function lmbd2string(lmbd) result(string)
   implicit none
   integer,intent(in) :: lmbd
   character(len=:),allocatable :: string
   select case(lmbd)
   case default;           string = 'unknown'
   case(p_mbd_none);       string = 'none'
   case(p_mbd_rpalike);    string = 'RPA like'
   case(p_mbd_exact_atm);  string = 'ATM'
   case(p_mbd_approx_atm); string = 'ATM'
   end select
end function lmbd2string

!> compute D4 energy under periodic boundary conditions
subroutine edisp_3d(mol,dispm,ndim,q,r_thr,mbd_thr,par, &
      &             gw,refc6,mbd,ed,etwo,embd)
   use iso_fortran_env, only : wp => real64
   use class_molecule
   use pbc_tools, only: get_realspace_cutoff
   implicit none
   type(dispersion_model),intent(in) :: dispm
   type(molecule),intent(in) :: mol
   integer, intent(in)  :: ndim
   real(wp),intent(in)  :: q(mol%n)
   integer              :: rep(3),mbd_rep(3)
   real(wp),intent(in)  :: r_thr,mbd_thr
   integer              :: tx,ty,tz
   real(wp)             :: t(3)
   real(wp)             :: aiw(23)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: refc6(ndim,ndim)
   integer, intent(in)  :: mbd
   real(wp),intent(out) :: ed
   real(wp),intent(out),optional :: etwo
   real(wp),intent(out),optional :: embd

   integer  :: i,ii,iii,j,jj,k,l,ia,ja,ij
   integer, allocatable :: itbl(:,:)
   real(wp) :: iz
   real(wp) :: qmod,eabc
   real(wp) :: norm,dnorm
   real(wp) :: dexpw,expw
   real(wp) :: twf,tgw,r4r2ij
   real(wp) :: rij(3),r,r2,r4,r6,r8,R0
   real(wp) :: oor6,oor8,oor10,door6,door8,door10,cutoff
   real(wp) :: c8abns,disp,ddisp,x1,x2,x3
   real(wp) :: c6ii,c6ij,dic6ii,dic6ij,djc6ij,dizii,dizij,djzij
   real(wp) :: rcovij,expterm,den,dcndr

   real(wp) :: drdx(3),dtmp,gwk,dgwk
   real(wp),allocatable :: r2ab(:)
   real(wp),allocatable :: dc6dcn(:)
   real(wp),allocatable :: zetavec(:)
   real(wp),allocatable :: zerovec(:)
   real(wp) :: cn_thr,gw_thr

   parameter(cn_thr = 1600.0_wp)
   parameter(gw_thr=0.000001_wp)

   intrinsic :: present,sqrt,sum,maxval,exp,abs

   if (mol%npbc > 0) then
      call get_realspace_cutoff(mol%lattice,r_thr,rep)
      call get_realspace_cutoff(mol%lattice,mbd_thr,mbd_rep)
   endif
   where(.not.mol%pbc) rep = 0
   where(.not.mol%pbc) mbd_rep = 0

   !  print'(" * Allocating local memory")'
   allocate( zetavec(ndim),zerovec(ndim), source = 0.0_wp )
   allocate( itbl(7,mol%n), source = 0 )

   ed = 0.0_wp
   eabc = 0.0_wp

   k = 0
   do i = 1, mol%n
      do ii = 1, dispm%nref(mol%at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

!  print'(" * Entering first OMP section")'
!$OMP parallel default(none) &
!$omp private(i,ii,iii,ia,iz,k)  &
!$omp shared (dispm,mol,itbl,q) &
!$omp shared (gw,zetavec,zerovec,r_thr)
!$omp do
   do i = 1, mol%n
      ia = mol%at(i)
      iz = zeff(ia)
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         zetavec(k) = zeta(dispm%g_a,gam(ia)*dispm%g_c,dispm%q(ii,ia)+iz,q(i)+iz) * gw(k)
         ! NEW: q=0 for ATM
         zerovec(k) =  zeta(dispm%g_a,gam(ia)*dispm%g_c,dispm%q(ii,ia)+iz,iz) * gw(k)

      enddo
   enddo
!$omp end do
!$omp end parallel

!$OMP parallel default(none) &
!$omp private(i,j,ia,ja,ij,k,l,c6ii,c6ij,disp,  &
!$omp         rij,r2,r,r4r2ij,r0,oor6,oor8,oor10, &
!$omp         t,tx,ty,tz)  &
!$omp shared(mol,dispm,itbl,zetavec,refc6,par,rep,r_Thr) &
!$omp shared(r2ab) reduction(+:ed)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      ia = mol%at(i)
      ! temps
      c6ij   = 0.0_wp
      ! all refs
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         do jj = 1, dispm%nref(ia)
            l = itbl(jj,i)
            c6ij   = c6ij   +  zetavec(k) *  zetavec(l) * refc6(k,l)
         enddo
      enddo
      ! i in primitive cell with i in images
      r4r2ij = 3*r4r2(ia)*r4r2(ia)
      r0 = par%a1*sqrt(r4r2ij) + par%a2
      do tx = -rep(1), rep(1)
      do ty = -rep(2), rep(2)
      do tz = -rep(3), rep(3)
         ! cycle i with i interaction in same cell
         if (tx.eq.0.and.ty.eq.0.and.tz.eq.0) cycle
         rij = tx*mol%lattice(:,1) + ty*mol%lattice(:,2) + tz*mol%lattice(:,3)
         r2  = sum( rij**2 )
         if (r2.gt.r_thr) cycle
         r   = sqrt(r2)
         oor6 = 1._wp/(r2**3+r0**6)
         oor8 = 1._wp/(r2**4+r0**8)
         oor10 = 1._wp/(r2**5+r0**10)
         disp = par%s6*oor6 + par%s8*r4r2ij*oor8   &
            & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10
         ed = ed - c6ij*disp/2
      enddo ! tz
      enddo ! ty
      enddo ! tx
      ! over all j atoms
      do j = 1, i-1
         ja = mol%at(j)
         ! temps
         c6ij   = 0.0_wp
         ! all refs
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               c6ij   = c6ij   +  zetavec(k) *  zetavec(l) * refc6(k,l)
            enddo
         enddo
         r4r2ij = 3*r4r2(ia)*r4r2(ja)
         r0 = par%a1*sqrt(r4r2ij) + par%a2
         do tx = -rep(1), rep(1)
         do ty = -rep(2), rep(2)
         do tz = -rep(3), rep(3)
            t = tx*mol%lattice(:,1) + ty*mol%lattice(:,2) + tz*mol%lattice(:,3)
            rij = mol%xyz(:,i) - mol%xyz(:,j) + t
            r2 = sum(rij**2)
            if (r2.gt.r_thr) cycle
            r = sqrt(r2)
            oor6 = 1._wp/(r2**3+r0**6)
            oor8 = 1._wp/(r2**4+r0**8)
            oor10 = 1._wp/(r2**5+r0**10)
            disp = par%s6*oor6 + par%s8*r4r2ij*oor8 &
               & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10
            ed = ed - c6ij*disp
         enddo ! tz
         enddo ! ty
         enddo ! tx
      enddo
   enddo
!$omp enddo
!$omp end parallel

   eabc = 0.0_wp
   if (mbd > 0) then
      call abcappr_3d(mol,dispm,ndim,par,zerovec,refc6,itbl,mbd_rep,mbd_thr,eabc)
   endif

   !  print'(" * Dispersion all done, saving variables")'
   if (present(etwo)) etwo = ed
   if (present(embd)) embd = eabc
   ed = ed + eabc

end subroutine edisp_3d

!> calculates threebody dispersion gradient from C6 coefficients
subroutine abcappr_3d(mol,dispm,ndim,par,zvec,refc6,itbl,rep,r_thr,eabc)
   use class_molecule
   use pbc_tools
   implicit none
   type(dispersion_model),intent(in) :: dispm
   type(molecule),intent(in) :: mol !< molecular structure information
   integer, intent(in)  :: ndim
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: zvec(ndim)
   real(wp),intent(in)  :: refc6(ndim,ndim)
   integer, intent(in)  :: itbl(7,mol%n)
   integer, intent(in)  :: rep(3)
   real(wp),intent(in)  :: r_thr
   real(wp),intent(out) :: eabc

   integer  :: i,ii,ia,j,jj,ja,k,l
   integer  :: ij
   real(wp),allocatable :: c6ab(:)
   real(wp) :: r,r2ij
   real(wp) :: c6ij
   real(wp) :: gw_thr
   parameter(gw_thr=0.0001_wp)

   intrinsic :: present,sqrt

   allocate( c6ab(mol%n*(mol%n+1)/2), source = 0.0_wp )

   eabc = 0.0_wp

!$OMP parallel default(none) &
!$omp private(i,ia,j,ja,ij,r2ij,c6ij,k,l,r)  &
!$omp shared (mol,dispm,itbl,refc6,zvec) &
!$omp shared (c6ab)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      ia = mol%at(i)
      ij = i*(i-1)/2+i

      ! temps
      c6ij = 0.0_wp
      ! all refs
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         do jj = 1, dispm%nref(ia)
            l = itbl(jj,i)
            c6ij = c6ij + zvec(k)*zvec(l)*refc6(k,l)
         enddo
      enddo
      ! save
      c6ab(ij) = c6ij

      do j = 1, i-1
!        if(i.eq.j) cycle
         ja = mol%at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         !if(mol%dist(j,i)**2.gt.r_thr) cycle

         ! temps
         c6ij = 0.0_wp
         ! all refs
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(k)*zvec(l)*refc6(k,l)
            enddo
         enddo
         ! save
         c6ab(ij) = c6ij
      enddo
   enddo
!$omp enddo
!$omp end parallel

   call abcappr_3d_dftd3_like_style(mol%n,mol%at,mol%xyz,par,r_thr,rep, &
      &                             mol%lattice,c6ab,eabc)
   !call abcappr_3d_wsc(mol,par,[2,2,2],r_thr,c6ab,eabc)
   !call abcappr_3d_bvk(mol,par,rep,r_thr,c6ab,eabc)

end subroutine abcappr_3d

subroutine abcappr_3d_dftd3_like_style(nat,at,xyz,par,thr,rep,dlat,c6ab,eabc)
   implicit none
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in) :: thr
   integer, intent(in) :: rep(3)
   real(wp),intent(in) :: c6ab(nat*(nat+1)/2)
   real(wp),intent(in) :: dlat(3,3)

   real(wp),intent(inout) :: eabc

   real(wp),parameter :: six = 6.0_wp, oth = 1.0_wp/3.0_wp

   integer iat,jat,kat
   real(wp) x1
   real(wp) r2,r
   real(wp) fdmp,tmp1

   real(wp) rij(3),rik(3),rjk(3)
   real(wp), allocatable,dimension(:,:,:,:) ::  dc6dr  !d(E)/d(r_ij) derivative wrt. dist. iat-jat
   !dCN(jat)/d(r_ij)
   real(wp) :: r9ijk
   real(wp) vec(3)
   integer ij,ik,jk

   real(wp),dimension(3) ::ijvec,ikvec,jkvec,t,s,dumvec
   integer tx,ty,tz,sx,sy,sz
   real(wp) rij2,rik2,rjk2,c9,c6ij,c6ik,c6jk,rijk,rijk3
   real(wp) :: cij,cjk,cik,cijk
   real(wp) time1,time2,rijk2,dc9,dfdmp,dang,ang
   integer,dimension(3) :: repmin,repmax

   allocate(dc6dr(-rep(3):rep(3),-rep(2):rep(2), &
      &           -rep(1):rep(1),nat*(nat+1)/2))
   dc6dr = 0.0_wp


   !        write(*,*)'!!!!!!!!!!    THREEBODY  GRADIENT  !!!!!!!!!!'
   eabc=0.0_wp
   !        write(*,*)'thr:',sqrt(thr)


!$omp parallel default(none) &
!$omp shared(nat,xyz,c6ab,rep,par,at,dlat,thr) &
!$omp private(iat,jat,ij,ijvec,c6ij,kat,ik,jk,ikvec,jkvec,c6ik,c6jk,c9, &
!$omp&        cij,cik,cjk,cijk,tx,ty,tz,repmin,repmax,t,rij2,sx,sy,sz, &
!$omp&        s,rik2,dumvec,rjk2,rijk2,rijk,fdmp,rijk3,ang,r9ijk,r) &
!$omp reduction(+:eabc)
!$omp do schedule(dynamic)
   do iat=3,nat
      do jat=2,iat-1
         ij=iat*(iat-1)/2+jat
         ijvec=xyz(:,jat)-xyz(:,iat)

         c6ij=c6ab(ij)
         do kat=1,jat-1
            ik=iat*(iat-1)/2+kat
            jk=jat*(jat-1)/2+kat
            ikvec=xyz(:,kat)-xyz(:,iat)
            jkvec=xyz(:,kat)-xyz(:,jat)

            c6ik=c6ab(ik)
            c6jk=c6ab(jk)
            c9=-par%s9*sqrt(c6ij*c6ik*c6jk)
            cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
            cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
            cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
            cijk = cij*cjk*cik

            do tx=-rep(1), rep(1)
            do ty=-rep(2), rep(2)
            do tz=-rep(3), rep(3)
               repmin(1)=max(-rep(1),tx-rep(1))
               repmax(1)=min(+rep(1),tx+rep(1))
               repmin(2)=max(-rep(2),ty-rep(2))
               repmax(2)=min(+rep(2),ty+rep(2))
               repmin(3)=max(-rep(3),tz-rep(3))
               repmax(3)=min(+rep(3),tz+rep(3))
               t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
               rij2=SUM((ijvec+t)*(ijvec+t))
               if(rij2.gt.thr)cycle

               !rr0ij=sqrt(rij2)/r0ab(at(iat),at(jat))


               do sx=repmin(1), repmax(1)
               do sy=repmin(2), repmax(2)
               do sz=repmin(3), repmax(3)
                  s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
                  rik2=SUM((ikvec+s)*(ikvec+s))
                  if(rik2.gt.thr)cycle

                  dumvec=jkvec+s-t
                  rjk2=SUM(dumvec*dumvec)
                  if(rjk2.gt.thr)cycle
                  !rr0ik=sqrt(rik2)/r0ab(at(iat),at(kat))
                  !rr0jk=sqrt(rjk2)/r0ab(at(jat),at(kat))
                  rijk2=(rij2*rjk2*rik2)
                  ! first calculate the three components for the energy calculation fdmp
                  ! and ang
                  !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
                  !damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

                  rijk=sqrt(rijk2)
                  fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
                  rijk3=rijk*rijk2
                  ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                     *(-rij2+rjk2+rik2)/(rijk3*rijk2) &
                     +1.0_wp/(rijk3)

                  r9ijk=ang*fdmp
                  eabc=eabc-r9ijk*c9

               enddo !sz
               enddo !sy
               enddo !sx
            enddo !tz
            enddo !ty
            enddo !tx
         enddo !kat
      enddo !jat
   enddo !iat
!$omp enddo

   ! Now the interaction with jat=iat of the triples iat,iat,kat
!$omp do schedule(dynamic)
   DO iat=2,nat
      jat=iat
      ij=iat*(iat-1)/2+jat
      ijvec=0.0_wp

      c6ij=c6ab(ij)
      DO kat=1,iat-1
         jk=jat*(jat-1)/2+kat
         ik=jk

         c6ik=c6ab(ik)
         c6jk=c6ik
         ikvec=xyz(:,kat)-xyz(:,iat)
         jkvec=ikvec
         c9=-par%s9*sqrt(c6ij*c6ik*c6jk)
         cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
         cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
         cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
         cijk = cij*cjk*cik
         do tx=-rep(1), rep(1)
         do ty=-rep(2), rep(2)
         do tz=-rep(3), rep(3)
            repmin(1)=max(-rep(1),tx-rep(1))
            repmax(1)=min(+rep(1),tx+rep(1))
            repmin(2)=max(-rep(2),ty-rep(2))
            repmax(2)=min(+rep(2),ty+rep(2))
            repmin(3)=max(-rep(3),tz-rep(3))
            repmax(3)=min(+rep(3),tz+rep(3))
            IF (tx.eq.0 .and. ty.eq.0 .and. tz.eq.0) cycle
            t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
            dumvec=t
            rij2=SUM(dumvec*dumvec)
            if(rij2.gt.thr)cycle

            !rr0ij=sqrt(rij2)/r0ab(at(iat),at(jat))

            do sx=repmin(1), repmax(1)
            do sy=repmin(2), repmax(2)
            do sz=repmin(3), repmax(3)
               ! every result * 0.5

               s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
               dumvec=ikvec+s
               dumvec=dumvec*dumvec
               rik2=SUM(dumvec)
               if(rik2.gt.thr)cycle

               dumvec=jkvec+s-t
               dumvec=dumvec*dumvec
               rjk2=SUM(dumvec)
               if(rjk2.gt.thr)cycle
               !rr0ik=sqrt(rik2)/r0ab(at(iat),at(kat))
               !rr0jk=sqrt(rjk2)/r0ab(at(jat),at(kat))


               rijk2=(rij2*rjk2*rik2)
               !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
               !damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

               rijk=sqrt(rijk2)
               fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
               rijk3=rijk*rijk2
               ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                  *(-rij2+rjk2+rik2)/(rijk3*rijk2) &
                  +1.0_wp/(rijk3)


               r9ijk=ang*fdmp/2.0_wp   !factor 1/2 for doublecounting
               eabc=eabc-r9ijk*c9

            enddo !sz
            enddo !sy
            enddo !sx

         enddo !tz
         enddo !ty
         enddo !tx
      enddo !kat
   enddo !iat
!$omp enddo
   ! And now kat=jat, but cycling throug all imagecells without t=s. and jat>iat going though all cells    (iat,jat,jat)
   ! But this counts only 1/2

!$omp do schedule(dynamic)
   do iat=2,nat
      do jat=1,iat-1
         kat=jat
         ij=iat*(iat-1)/2+jat
         jk=jat*(jat-1)/2+kat
         ik=ij

         c6ij=c6ab(ij)
         c6ik=c6ij

         c6jk=c6ab(jk)
         ikvec=xyz(:,kat)-xyz(:,iat)
         ijvec=ikvec
         jkvec=0.0_wp
         cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
         cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
         cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
         cijk = cij*cjk*cik

         c9=-par%s9*sqrt(c6ij*c6ik*c6jk)
         do tx=-rep(1), rep(1)
         do ty=-rep(2), rep(2)
         do tz=-rep(3), rep(3)
            repmin(1)=max(-rep(1),tx-rep(1))
            repmax(1)=min(+rep(1),tx+rep(1))
            repmin(2)=max(-rep(2),ty-rep(2))
            repmax(2)=min(+rep(2),ty+rep(2))
            repmin(3)=max(-rep(3),tz-rep(3))
            repmax(3)=min(+rep(3),tz+rep(3))

            t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
            dumvec=ijvec+t
            dumvec=dumvec*dumvec
            rij2=SUM(dumvec)
            if(rij2.gt.thr)cycle

            !rr0ij=SQRT(rij2)/r0ab(at(iat),at(jat))

            do sx=repmin(1), repmax(1)
            do sy=repmin(2), repmax(2)
            do sz=repmin(3), repmax(3)
               ! every result * 0.5
               IF (tx.eq.sx .and. ty.eq.sy  &
                  .and. tz.eq.sz) cycle
               s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
               dumvec=ikvec+s
               dumvec=dumvec*dumvec
               rik2=SUM(dumvec)
               if(rik2.gt.thr)cycle
               !rr0ik=SQRT(rik2)/r0ab(at(iat),at(kat))

               dumvec=jkvec+s-t
               dumvec=dumvec*dumvec
               rjk2=SUM(dumvec)
               if(rjk2.gt.thr)cycle
               !rr0jk=SQRT(rjk2)/r0ab(at(jat),at(kat))

               !              if (rij*rjk*rik.gt.thr)cycle

               rijk2=(rij2*rjk2*rik2)
               !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
               !damp9=1./(1.+6._wp*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

               rijk=sqrt(rijk2)
               fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
               rijk3=rijk*rijk2
               ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                  *(-rij2+rjk2+rik2)/(rijk2*rijk3) &
                  +1.0_wp/(rijk3)
               r9ijk=ang*fdmp/2.0_wp   !factor 1/2 for doublecounting
               eabc=eabc-r9ijk*c9

            enddo !sz
            enddo !sy
            enddo !sx

         enddo !tz
         enddo !ty
         enddo !tx
      enddo !kat
   enddo !iat
!$omp enddo


   ! and finally the self interaction iat=jat=kat all

!$omp do schedule(dynamic)
   do iat=1,nat
      jat=iat
      kat=iat
      ijvec=0.0_wp
      ij=iat*(iat-1)/2+jat
      ik=iat*(iat-1)/2+kat
      jk=jat*(jat-1)/2+kat
      ikvec=ijvec
      jkvec=ikvec
      c6ij=c6ab(ij)
      c6ik=c6ij
      c6jk=c6ij
      c9=-par%s9*sqrt(c6ij*c6ij*c6ij)
      cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
      cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
      cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
      cijk = cij*cjk*cik

      do tx=-rep(1), rep(1)
      do ty=-rep(2), rep(2)
      do tz=-rep(3), rep(3)
         repmin(1)=max(-rep(1),tx-rep(1))
         repmax(1)=min(+rep(1),tx+rep(1))
         repmin(2)=max(-rep(2),ty-rep(2))
         repmax(2)=min(+rep(2),ty+rep(2))
         repmin(3)=max(-rep(3),tz-rep(3))
         repmax(3)=min(+rep(3),tz+rep(3))
         if ((tx.eq.0) .and.(ty.eq.0) .and.(tz.eq.0))cycle !IF iat and jat are the same then cycle
         t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
         dumvec=t
         dumvec=dumvec*dumvec
         rij2=SUM(dumvec)
         if(rij2.gt.thr)cycle
         !rr0ij=SQRT(rij2)/r0ab(at(iat),at(jat))

         do sx=repmin(1), repmax(1)
         do sy=repmin(2), repmax(2)
         do sz=repmin(3), repmax(3)
            if ((sx.eq.0) .and.( sy.eq.0) .and.( sz.eq.0))cycle !IF iat and kat are the same then cycle
            if ((sx.eq.tx) .and. (sy.eq.ty)  &
               .and. (sz.eq.tz)) cycle      !If kat and jat are the same then cycle

            ! every result * 1/6 becaues every triple is counted twice due to the two loops t and s going from -rep to rep -> *1/2
            !
            !plus 1/3 becaues every triple is three times in each unitcell
            s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
            dumvec=s
            dumvec=dumvec*dumvec
            rik2=SUM(dumvec)
            if(rik2.gt.thr)cycle
            !rr0ik=SQRT(rik2)/r0ab(at(iat),at(kat))

            dumvec=jkvec+s-t
            dumvec=dumvec*dumvec
            rjk2=SUM(dumvec)
            if(rjk2.gt.thr)cycle
            !rr0jk=SQRT(rjk2)/r0ab(at(jat),at(kat))

            rijk2=(rij2*rjk2*rik2)
            !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
            !damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

            rijk=sqrt(rijk2)
            fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
            rijk3=rijk*rijk2
            ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
               *(-rij2+rjk2+rik2)/(rijk2*rijk3) &
               +1.0_wp/(rijk3)
            r9ijk=ang*fdmp/6.0_wp
            eabc=eabc-c9*r9ijk

         enddo !sz
         enddo !sy
         enddo !sx
      enddo !tz
      enddo !ty
      enddo !tx


   enddo !iat
!$omp enddo
!$omp end parallel


end subroutine abcappr_3d_dftd3_like_style

!> compute D4 gradient under pbc
subroutine dispgrad_3d(mol,dispm,ndim,q,cn,dcndr,dcndL,r_thr,mbd_thr,par,wf, &
      &                refc6,mbd,g,sigma,eout,dqdr,dqdL,aout)
   use iso_fortran_env, only : wp => real64
   use class_molecule
   use pbc_tools
   implicit none
   type(dispersion_model),intent(in) :: dispm
   type(molecule),intent(in) :: mol
   integer, intent(in)  :: ndim
   real(wp),intent(in)  :: q(mol%n)
   real(wp),intent(in)  :: cn(mol%n)
   real(wp),intent(in)  :: dcndr(3,mol%n,mol%n)
   real(wp),intent(in)  :: dcndL(3,3,mol%n)
   integer  :: rep(3)
   integer  :: mbd_rep(3)
   real(wp),intent(in)  :: r_thr
   real(wp),intent(in)  :: mbd_thr
   integer              :: tx,ty,tz
   real(wp)             :: t(3)
   real(wp)             :: aiw(23)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: wf
   real(wp),intent(in)  :: refc6(ndim,ndim)
   integer, intent(in)  :: mbd
   real(wp),intent(inout)        :: g(3,mol%n)
   real(wp),intent(inout)        :: sigma(3,3)
   real(wp),intent(out),optional :: eout
   real(wp),intent(in), optional :: dqdr(3,mol%n,mol%n+1)
   real(wp),intent(in), optional :: dqdL(3,3,mol%n+1)
   real(wp),intent(out),optional :: aout(23,mol%n)


   integer  :: i,ii,iii,j,jj,k,l,ia,ja,ij
   integer, allocatable :: itbl(:,:)
   real(wp) :: iz
   real(wp) :: qmod,eabc,ed
   real(wp) :: norm,dnorm
   real(wp) :: dexpw,expw
   real(wp) :: twf,tgw,r4r2ij
   real(wp) :: rij(3),r,r2,r4,r6,r8,R0
   real(wp) :: oor6,oor8,oor10,door6,door8,door10,cutoff
   real(wp) :: c8abns,disp,ddisp,x1,x2,x3
   real(wp) :: c6ii,c6ij,dic6ii,dic6ij,djc6ij,dizii,dizij,djzij
   real(wp) :: rcovij,expterm,den

   real(wp) :: drdx(3),dtmp,gwk,dgwk
   real(wp),allocatable :: r2ab(:)
   real(wp),allocatable :: dc6dcn(:)
   real(wp),allocatable :: zvec(:)
   real(wp),allocatable :: dzvec(:)
   real(wp),allocatable :: gw(:)
   real(wp),allocatable :: dgw(:)
   real(wp),allocatable :: dc6dq(:)
   real(wp),allocatable :: dzdq(:)
   real(wp) :: cn_thr,gw_thr

   parameter(cn_thr = 1600.0_wp)
   parameter(gw_thr=0.000001_wp)

   intrinsic :: present,sqrt,sum,maxval,exp,abs

   if (mol%npbc > 0) then
      call get_realspace_cutoff(mol%lattice,r_thr,rep)
      call get_realspace_cutoff(mol%lattice,mbd_thr,mbd_rep)
   endif
   where(.not.mol%pbc) rep = 0
   where(.not.mol%pbc) mbd_rep = 0

   !  print'(" * Allocating local memory")'
   allocate( dc6dcn(mol%n),                  &
      &      r2ab(mol%n*(mol%n+1)/2),     &
      &      zvec(ndim),dzvec(ndim),  &
      &      gw(ndim),dgw(ndim),dc6dq(mol%n),dzdq(ndim),  &
      &         source = 0.0_wp )
   allocate( itbl(7,mol%n), source = 0 )

   ed = 0.0_wp
   eabc = 0.0_wp

   k = 0
   do i = 1, mol%n
      do ii = 1, dispm%nref(mol%at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

!  print'(" * Entering first OMP section")'
!$OMP parallel default(none) &
!$omp private(i,ii,iii,ia,iz,k,norm,dnorm,twf,tgw,dexpw,expw,gwk,dgwk)  &
!$omp shared (mol,dispm,itbl,wf,cn,q) &
!$omp shared (gw,dgw,zvec,dzvec,dzdq,r_thr)
!$omp do
   do i = 1, mol%n
      ia = mol%at(i)
      iz = zeff(ia)
      norm  = 0.0_wp
      dnorm = 0.0_wp
      do ii=1,dispm%nref(ia)
         do iii = 1, dispm%ncount(ii,ia)
            twf = iii*wf
            tgw = cngw(twf,cn(i),dispm%cn(ii,ia))
            norm  =  norm + tgw
            dnorm = dnorm + 2*twf*(dispm%cn(ii,ia)-cn(i))*tgw
         enddo
      enddo
      norm = 1._wp/norm
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         dexpw=0.0_wp
         expw=0.0_wp
         do iii = 1, dispm%ncount(ii,ia)
            twf = wf*iii
            tgw = cngw(twf,cn(i),dispm%cn(ii,ia))
            expw  =  expw + tgw
            dexpw = dexpw + 2*twf*(dispm%cn(ii,ia)-cn(i))*tgw
         enddo

         ! save
         gwk = expw*norm
         if (gwk.ne.gwk) then
            if (maxval(dispm%cn(:dispm%nref(ia),ia)).eq.dispm%cn(ii,ia)) then
               gwk = 1.0_wp
            else
               gwk = 0.0_wp
            endif
         endif
         zvec(k) = zeta(dispm%g_a,gam(ia)*dispm%g_c,dispm%q(ii,ia)+iz,q(i)+iz) * gwk
         ! NEW: q=0 for ATM
         gw(k) =  zeta(dispm%g_a,gam(ia)*dispm%g_c,dispm%q(ii,ia)+iz,iz) * gwk

         dgwk = norm*(dexpw-expw*dnorm*norm)
         if (dgwk.ne.dgwk) then
            dgwk = 0.0_wp
         endif
         dzvec(k) = zeta(dispm%g_a,gam(ia)*dispm%g_c,dispm%q(ii,ia)+iz,q(i)+iz) * dgwk
         dzdq(k) = dzeta(dispm%g_a,gam(ia)*dispm%g_c,dispm%q(ii,ia)+iz,q(i)+iz) * gwk
         ! NEW: q=0 for ATM
         dgw(k) = zeta(dispm%g_a,gam(ia)*dispm%g_c,dispm%q(ii,ia)+iz,iz) * dgwk
      enddo
   enddo
!$omp end do
!$omp end parallel

!$OMP parallel default(none) &
!$omp private(i,j,ia,ja,ij,k,l,c6ii,c6ij,dic6ii,dic6ij,djc6ij,disp,ddisp, &
!$omp         rij,r2,r,r4r2ij,r0,oor6,oor8,oor10,door6,door8,door10, &
!$omp         t,tx,ty,tz,dtmp,drdx,dizii,dizij,djzij)  &
!$omp shared(mol,dispm,itbl,zvec,dzvec,refc6,par,dzdq,rep,r_Thr) &
!$omp shared(r2ab) reduction(+:dc6dq,dc6dcn,ed,g,sigma)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      ia = mol%at(i)
      ! temps
      c6ij   = 0.0_wp
      dic6ij = 0.0_wp
      djc6ij = 0.0_wp
      dizij  = 0.0_wp
      djzij  = 0.0_wp
      ! all refs
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         do jj = 1, dispm%nref(ia)
            l = itbl(jj,i)
            c6ij   = c6ij   +  zvec(k) *  zvec(l) * refc6(k,l)
            dic6ij = dic6ij + dzvec(k) *  zvec(l) * refc6(k,l)
            djc6ij = djc6ij +  zvec(k) * dzvec(l) * refc6(k,l)
            dizij  = dizij  +  dzdq(k) *  zvec(l) * refc6(k,l)
            djzij  = djzij  +  zvec(k) *  dzdq(l) * refc6(k,l)
         enddo
      enddo
      ! i in primitive cell with i in images
      r4r2ij = 3*r4r2(ia)*r4r2(ia)
      r0 = par%a1*sqrt(r4r2ij) + par%a2
      do tx = -rep(1), rep(1)
      do ty = -rep(2), rep(2)
      do tz = -rep(3), rep(3)
         if (tx.eq.0.and.ty.eq.0.and.tz.eq.0) cycle
         ! cycle i with i interaction in same cell
         t = [tx,ty,tz]
         rij = matmul(mol%lattice,t)
         r2  = sum( rij**2 )
         if (r2.gt.r_thr) cycle
         r   = sqrt(r2)
         oor6 = 1._wp/(r2**3+r0**6)
         oor8 = 1._wp/(r2**4+r0**8)
         oor10 = 1._wp/(r2**5+r0**10)
         door6 = -6*r2**2*r*oor6**2
         door8 = -8*r2**3*r*oor8**2
         door10 = -10*r2**4*r*oor10**2
         disp = par%s6*oor6 + par%s8*r4r2ij*oor8   &
            & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10
         ddisp= par%s6*door6 + par%s8*r4r2ij*door8 &
            & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*door10
         ed = ed - c6ij*disp/2
         ! save this
         dtmp = c6ij*ddisp
         dc6dq(i)  = dc6dq(i)  + (dizij  + djzij )*disp/2
         dc6dcn(i) = dc6dcn(i) + (dic6ij + djc6ij)*disp/2
         drdx = rij/r
         sigma = sigma - dtmp * outer_prod_3x3(drdx,rij)/2
      enddo ! tz
      enddo ! ty
      enddo ! tx
      ! over all j atoms
      do j = 1, i-1
         ja = mol%at(j)
         ! temps
         c6ij   = 0.0_wp
         dic6ij = 0.0_wp
         djc6ij = 0.0_wp
         dizij  = 0.0_wp
         djzij  = 0.0_wp
         ! all refs
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               c6ij   = c6ij   +  zvec(k) *  zvec(l) * refc6(k,l)
               dic6ij = dic6ij + dzvec(k) *  zvec(l) * refc6(k,l)
               djc6ij = djc6ij +  zvec(k) * dzvec(l) * refc6(k,l)
               dizij  = dizij  +  dzdq(k) *  zvec(l) * refc6(k,l)
               djzij  = djzij  +  zvec(k) *  dzdq(l) * refc6(k,l)
            enddo
         enddo
         r4r2ij = 3*r4r2(ia)*r4r2(ja)
         r0 = par%a1*sqrt(r4r2ij) + par%a2
         do tx = -rep(1), rep(1)
         do ty = -rep(2), rep(2)
         do tz = -rep(3), rep(3)
            t = [tx,ty,tz]
            rij = mol%xyz(:,i) - mol%xyz(:,j) + matmul(mol%lattice,t)
            r2 = sum(rij**2)
            if (r2.gt.r_thr) cycle
            r = sqrt(r2)
            oor6 = 1._wp/(r2**3+r0**6)
            oor8 = 1._wp/(r2**4+r0**8)
            oor10 = 1._wp/(r2**5+r0**10)
            door6 = -6*r2**2*r*oor6**2
            door8 = -8*r2**3*r*oor8**2
            door10 = -10*r2**4*r*oor10**2
            disp = par%s6*oor6 + par%s8*r4r2ij*oor8 &
               & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10
            ddisp= par%s6*door6 + par%s8*r4r2ij*door8 &
               & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*door10
            ed = ed - c6ij*disp
            ! save this
            dc6dq(i)  = dc6dq(i)  + dizij  *disp
            dc6dq(j)  = dc6dq(j)  + djzij  *disp
            dc6dcn(i) = dc6dcn(i) + dic6ij *disp
            dc6dcn(j) = dc6dcn(j) + djc6ij *disp
            dtmp = c6ij*ddisp
            drdx = rij/r
            g(:,i) = g(:,i) - dtmp * drdx
            g(:,j) = g(:,j) + dtmp * drdx

            sigma = sigma - dtmp * outer_prod_3x3(drdx,rij)
         enddo ! tz
         enddo ! ty
         enddo ! tx
      enddo
   enddo
!$omp enddo
!$omp end parallel

   eabc = 0.0_wp
   if (mbd > 0) then
      call dabcappr_3d(mol,dispm,ndim,par,gw,dgw,refc6,itbl,mbd_rep,mbd_thr, &
         &             g,sigma,dc6dcn,eabc)
      !dc6dcn = 0.0_wp
   endif

   if(present(dqdr)) then
      ! handle dqdr  :: gradient is exact e-11
      call dgemv('n',3*mol%n,mol%n,-1.0_wp,dqdr,3*mol%n,dc6dq,1,1.0_wp,g,1)
   endif
   if (mol%npbc > 0) then
      if(present(dqdL)) then
         ! handle dqdL
         call dgemv('n',3*3,mol%n,-1.0_wp,dqdL,3*3,dc6dq,1,1.0_wp,sigma,1)
      endif
   endif

   ! always handle dcndr :: gradient is exact e-11
   call dgemv('n',3*mol%n,mol%n,-1.0_wp,dcndr,3*mol%n,dc6dcn,1,1.0_wp,g,1)
   if (mol%npbc > 0) then
      call dgemv('n',3*3,mol%n,-1.0_wp,dcndL,3*3,dc6dcn,1,1.0_wp,sigma,1)
   endif

   !  print'(" * Dispersion all done, saving variables")'
   if (present(eout)) eout = ed + eabc

   if (present(aout)) then
      aout = 0._wp
      do i = 1, mol%n
         ia = mol%at(i)
         do ii = 1, dispm%nref(ia)
            aout(:,i) = aout(:,i) + zvec(k) * dispm%alpha(:,ii,ia)
         enddo
      enddo
   endif
end subroutine dispgrad_3d

!> calculates threebody dispersion gradient from C6 coefficients
subroutine dabcappr_3d(mol,dispm,ndim,par,zvec,dzvec,refc6,itbl,rep,r_thr,g,sigma,dc6dcn,eabc)
   use class_molecule
   use pbc_tools
   implicit none
   type(dispersion_model),intent(in) :: dispm
   type(molecule),intent(in) :: mol !< molecular structure information
   integer, intent(in)  :: ndim
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: zvec(ndim)
   real(wp),intent(in)  :: dzvec(ndim)
   real(wp),intent(in)  :: refc6(ndim,ndim)
   integer, intent(in)  :: itbl(7,mol%n)
   integer, intent(in)  :: rep(3)
   real(wp),intent(in)  :: r_thr
   real(wp),intent(inout) :: g(3,mol%n)
   real(wp),intent(inout) :: sigma(3,3)
   real(wp),intent(inout) :: dc6dcn(mol%n)
   real(wp),intent(out)   :: eabc

   integer  :: i,ii,ia,j,jj,ja,k,kk,ka,l,m,n
   integer  :: ij,jk,ik
   real(wp),allocatable :: c6ab(:),dc6ab(:,:)
   real(wp) :: r2ij,r2jk,r2ik,r
   real(wp) :: cii,cij,cjk,cik,ciii,ciij,cijk
   real(wp) :: c9iii,c9iij,c9ijk,oor9ijk,rijk
   real(wp) :: rij(3),rjk(3),rik(3)
   real(wp) :: dijfdmp,dikfdmp,djkfdmp
   real(wp) :: dijatm,dikatm,djkatm
   real(wp) :: dijoor9ijk,djkoor9ijk,dikoor9ijk
   real(wp) :: c6ij,dic6ij,djc6ij
   real(wp) :: dic9iii,dic9iij,djc9iij,dic9ijk,djc9ijk,dkc9ijk
   real(wp),parameter :: zero(3) = [0.0_wp,0.0_wp,0.0_wp]
   real(wp) :: gw_thr
   parameter(gw_thr=0.0001_wp)

   intrinsic :: present,sqrt

   allocate( c6ab(mol%n*(mol%n+1)/2),dc6ab(mol%n,mol%n), source = 0.0_wp )

   eabc = 0.0_wp

!$OMP parallel default(none) &
!$omp private(i,ia,j,ja,ij,r2ij,c6ij,dic6ij,djc6ij,k,l,r)  &
!$omp shared (mol,dispm,itbl,refc6,zvec,dzvec) &
!$omp shared (c6ab,dc6ab)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      ia = mol%at(i)
      ij = i*(i-1)/2+i

      ! temps
      c6ij = 0.0_wp
      dic6ij = 0.0_wp
      djc6ij = 0.0_wp
      ! all refs
      do ii = 1, dispm%nref(ia)
         k = itbl(ii,i)
         do jj = 1, dispm%nref(ia)
            l = itbl(jj,i)
            c6ij = c6ij + zvec(k)*zvec(l)*refc6(k,l)
            dic6ij = dic6ij + dzvec(k)*zvec(l)*refc6(k,l)
            djc6ij = djc6ij + zvec(k)*dzvec(l)*refc6(k,l)
         enddo
      enddo
      ! save
      c6ab(ij) = c6ij
      dc6ab(i,i) = dic6ij! + djc6ij

      do j = 1, i-1
!        if(i.eq.j) cycle
         ja = mol%at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         !if(mol%dist(j,i)**2.gt.r_thr) cycle

         ! temps
         c6ij = 0.0_wp
         dic6ij = 0.0_wp
         djc6ij = 0.0_wp
         ! all refs
         do ii = 1, dispm%nref(ia)
            k = itbl(ii,i)
            do jj = 1, dispm%nref(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(k)*zvec(l)*refc6(k,l)
               dic6ij = dic6ij + dzvec(k)*zvec(l)*refc6(k,l)
               djc6ij = djc6ij + zvec(k)*dzvec(l)*refc6(k,l)
            enddo
         enddo
         ! save
         c6ab(ij) = c6ij
         dc6ab(i,j) = dic6ij
         dc6ab(j,i) = djc6ij
      enddo
   enddo
!$omp enddo
!$omp end parallel

   call dabcappr_3d_dftd3_like_style(mol%n,mol%at,mol%xyz,par,r_thr,rep,&
      &    mol%lattice,c6ab,dc6ab,eabc,dc6dcn,g,sigma)
   !call dabcappr_3d_wsc(mol,par,[2,2,2],r_thr,c6ab,dc6ab,g,dc6dcn,eabc)
   !call dabcappr_3d_bvk(mol,par,rep,r_thr,c6ab,dc6ab,g,dc6dcn,eabc)

end subroutine dabcappr_3d

subroutine dabcappr_3d_dftd3_like_style(nat,at,xyz,par,thr,rep,dlat,c6ab,dc6ab, &
      &    eabc,dc6dcn,g,sigma)
   use mctc_constants
   use pbc_tools
   implicit none
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in) :: thr
   integer, intent(in) :: rep(3)
   real(wp),intent(in) :: c6ab(nat*(nat+1)/2)
   real(wp),intent(in) :: dc6ab(nat,nat)    !dC6(iat,jat)/cCN(iat) in dc6ab(i,j) for ABC-grad
   real(wp),intent(in) :: dlat(3,3)

   real(wp),intent(inout) :: g(3,nat)
   real(wp),intent(inout) :: sigma(3,3)
   real(wp),intent(inout) :: eabc
   real(wp),intent(inout) :: dc6dcn(nat)

   real(wp),parameter :: six = 6.0_wp, oth = 1.0_wp/3.0_wp

   integer iat,jat,kat
   real(wp) x1
   real(wp) r2,r
   real(wp) fdmp,tmp1
   real(wp) :: rco,den,tmp,dtmp

   real(wp) rij(3),rik(3),rjk(3)
   !dCN(jat)/d(r_ij)
   real(wp) :: r9ijk
   real(wp) vec(3)
   real(wp) :: r3,g3(3,3)
   integer ij,ik,jk

   real(wp),dimension(3) ::ijvec,ikvec,jkvec,t,s,dumvec
   integer tx,ty,tz,sx,sy,sz
   real(wp) rij2,rik2,rjk2,c9,c6ij,c6ik,c6jk,rijk,rijk3
   real(wp) :: cij,cjk,cik,cijk
   real(wp) time1,time2,rijk2,dc9,dfdmp,dang,ang
   integer,dimension(3) :: repmin,repmax


   !        write(*,*)'!!!!!!!!!!    THREEBODY  GRADIENT  !!!!!!!!!!'
   eabc=0.0_wp
   !        write(*,*)'thr:',sqrt(thr)

!$omp parallel default(none) &
!$omp shared(nat,xyz,c6ab,rep,par,at,dlat,dc6ab,thr) &
!$omp private(iat,jat,ij,ijvec,c6ij,kat,ik,jk,ikvec,jkvec,c6ik,c6jk,c9, &
!$omp&        cij,cik,cjk,cijk,tx,ty,tz,repmin,repmax,t,rij2,sx,sy,sz, &
!$omp&        s,rik2,vec,rjk2,rijk2,rijk,fdmp,rijk3,ang,r9ijk,dfdmp, &
!$omp&        r,dang,tmp1,dc9,rij,rik,rjk,r3,g3) &
!$omp reduction(+:eabc,dc6dcn,sigma,g)
!$omp do schedule(dynamic)
   iAt_ijk: do iat=3,nat
      jAt_ijk: do jat=2,iat-1
         ij=iat*(iat-1)/2+jat
         ijvec=xyz(:,jat)-xyz(:,iat)

         c6ij=c6ab(ij)
         kAt_ijk: do kat=1,jat-1
            ik=iat*(iat-1)/2+kat
            jk=jat*(jat-1)/2+kat
            ikvec=xyz(:,kat)-xyz(:,iat)
            jkvec=xyz(:,kat)-xyz(:,jat)

            c6ik=c6ab(ik)
            c6jk=c6ab(jk)
            c9=-par%s9*sqrt(c6ij*c6ik*c6jk)
            cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
            cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
            cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
            cijk = cij*cjk*cik

            r3 = 0.0_wp
            g3 = 0.0_wp
            do tx=-rep(1), rep(1)
            do ty=-rep(2), rep(2)
            do tz=-rep(3), rep(3)
               repmin(1)=max(-rep(1),tx-rep(1))
               repmax(1)=min(+rep(1),tx+rep(1))
               repmin(2)=max(-rep(2),ty-rep(2))
               repmax(2)=min(+rep(2),ty+rep(2))
               repmin(3)=max(-rep(3),tz-rep(3))
               repmax(3)=min(+rep(3),tz+rep(3))
               t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
               rij=ijvec+t
               rij2=SUM(rij*rij)
               if(rij2.gt.thr)cycle


               do sx=repmin(1), repmax(1)
               do sy=repmin(2), repmax(2)
               do sz=repmin(3), repmax(3)
                  s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
                  rik = ikvec+s
                  rik2=SUM(rik*rik)
                  if(rik2.gt.thr)cycle

                  rjk = jkvec+s-t
                  rjk2=SUM(rjk*rjk)
                  if(rjk2.gt.thr)cycle
                  rijk2=(rij2*rjk2*rik2)
                  ! first calculate the three components for the energy calculation fdmp
                  ! and ang

                  rijk=sqrt(rijk2)
                  fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
                  rijk3=rijk*rijk2
                  ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                     *(-rij2+rjk2+rik2)/(rijk3*rijk2) &
                     +1.0_wp/(rijk3)

                  r9ijk=ang*fdmp
                  r3 = r3 + r9ijk
                  !
                  !start calculating the gradient components dfdmp, dang and dc9

                  !dfdmp is the same for all three distances
                  dfdmp = -(oth*six*par%alp*((cijk/rijk)**oth)**par%alp)*fdmp**2

                  !start calculating the derivatives of each part w.r.t. r_ij
                  r=sqrt(rij2)


                  dang=-0.375_wp*(rij2**3+rij2**2*(rjk2+rik2) &
                     +rij2*(3.0_wp*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2) &
                     -5.0*(rjk2-rik2)**2*(rjk2+rik2)) &
                     /(r*rijk3*rijk2)

                  tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)/r
                  g3(:,1) = g3(:,1) - tmp1*rij
                  g3(:,2) = g3(:,2) + tmp1*rij
                  sigma = sigma + outer_prod_3x3(rij,rij)*tmp1

                  !start calculating the derivatives of each part w.r.t. r_ik

                  r=sqrt(rik2)


                  dang=-0.375_wp*(rik2**3+rik2**2*(rjk2+rij2) &
                     +rik2*(3.0_wp*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2) &
                     -5.0*(rjk2-rij2)**2*(rjk2+rij2)) &
                     /(r*rijk3*rijk2)

                  tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)/r
                  !                 tmp1=-dc9
                  g3(:,1) = g3(:,1) - tmp1*rik
                  g3(:,3) = g3(:,3) + tmp1*rik
                  sigma = sigma + outer_prod_3x3(rik,rik)*tmp1

                  !
                  !start calculating the derivatives of each part w.r.t. r_jk

                  r=sqrt(rjk2)

                  dang=-0.375_wp*(rjk2**3+rjk2**2*(rik2+rij2) &
                     +rjk2*(3.0_wp*rik2**2+2.0*rik2*rij2+3.0*rij2**2) &
                     -5.0*(rik2-rij2)**2*(rik2+rij2)) &
                     /(r*rijk3*rijk2)

                  tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)/r
                  g3(:,2) = g3(:,2) - tmp1*rjk
                  g3(:,3) = g3(:,3) + tmp1*rjk
                  sigma = sigma + outer_prod_3x3(rjk,rjk)*tmp1

               enddo !sz
               enddo !sy
               enddo !sx
            enddo !tz
            enddo !ty
            enddo !tx
            eabc = eabc - r3*c9
            g(:,iat) = g(:,iat) + g3(:,1)
            g(:,jat) = g(:,jat) + g3(:,2)
            g(:,kat) = g(:,kat) + g3(:,3)
            !calculating the CN derivative dE_disp(ijk)/dCN(i)
            dc9=0.5_wp*c9*(dc6ab(iat,jat)/c6ij+dc6ab(iat,kat)/c6ik)
            dc6dcn(iat) = dc6dcn(iat) + r3*dc9
            dc9=0.5_wp*c9*(dc6ab(jat,iat)/c6ij+dc6ab(jat,kat)/c6jk)
            dc6dcn(jat) = dc6dcn(jat) + r3*dc9
            dc9=0.5_wp*c9*(dc6ab(kat,iat)/c6ik+dc6ab(kat,jat)/c6jk)
            dc6dcn(kat) = dc6dcn(kat) + r3*dc9
         enddo kAt_ijk
      enddo jAt_ijk
   enddo iAt_ijk
!$omp enddo

!$omp do schedule(dynamic)
   ! Now the interaction with jat=iat of the triples iat,iat,kat
   iAt_iik: do iat=2,nat
      jat=iat
      ij=iat*(iat-1)/2+jat
      ijvec=0.0_wp

      c6ij=c6ab(ij)
      kAt_iik: do kat=1,iat-1
         jk=jat*(jat-1)/2+kat
         ik=jk

         c6ik=c6ab(ik)
         c6jk=c6ik
         ikvec=xyz(:,kat)-xyz(:,iat)
         jkvec=ikvec
         c9=-par%s9*sqrt(c6ij*c6ik*c6jk)
         cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
         cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
         cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
         cijk = cij*cjk*cik

         g3 = 0.0_wp
         r3 = 0.0_wp
         do tx=-rep(1), rep(1)
         do ty=-rep(2), rep(2)
         do tz=-rep(3), rep(3)
            repmin(1)=max(-rep(1),tx-rep(1))
            repmax(1)=min(+rep(1),tx+rep(1))
            repmin(2)=max(-rep(2),ty-rep(2))
            repmax(2)=min(+rep(2),ty+rep(2))
            repmin(3)=max(-rep(3),tz-rep(3))
            repmax(3)=min(+rep(3),tz+rep(3))
            IF (tx.eq.0 .and. ty.eq.0 .and. tz.eq.0) cycle
            t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
            rij=t
            rij2=SUM(rij*rij)
            if(rij2.gt.thr)cycle

            do sx=repmin(1), repmax(1)
            do sy=repmin(2), repmax(2)
            do sz=repmin(3), repmax(3)
               ! every result * 0.5

               s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
               rik=ikvec+s
               rik2=SUM(rik*rik)
               if(rik2.gt.thr)cycle

               rjk=jkvec+s-t
               rjk2=SUM(rjk*rjk)
               if(rjk2.gt.thr)cycle

               rijk2=(rij2*rjk2*rik2)

               rijk=sqrt(rijk2)
               fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
               rijk3=rijk*rijk2
               ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                  *(-rij2+rjk2+rik2)/(rijk3*rijk2) &
                  +1.0_wp/(rijk3)


               r9ijk=ang*fdmp/2.0_wp   !factor 1/2 for doublecounting
               r3 = r3 + r9ijk

               !              iat=jat
               !dfdmp=2._wp*alp9*(0.75_wp*r0av)**(alp9)*fdmp*fdmp
               dfdmp = -(oth*six*par%alp*((cijk/rijk)**oth)**par%alp)*fdmp**2

               !start calculating the derivatives of each part w.r.t. r_ij
               r=sqrt(rij2)

               dang=-0.375_wp*(rij2**3+rij2**2*(rjk2+rik2) &
                  +rij2*(3.0_wp*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2) &
                  -5.0*(rjk2-rik2)**2*(rjk2+rik2)) &
                  /(r*rijk3*rijk2)

               tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)*0.5_wp/r
               sigma = sigma + outer_prod_3x3(rij,rij)*tmp1

               !start calculating the derivatives of each part w.r.t. r_ik
               r=sqrt(rik2)


               dang=-0.375_wp*(rik2**3+rik2**2*(rjk2+rij2) &
                  +rik2*(3.0_wp*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2) &
                  -5.0*(rjk2-rij2)**2*(rjk2+rij2)) &
                  /(r*rijk3*rijk2)

               tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)*0.5_wp/r
               g3(:,1) = g3(:,1) - tmp1*rik
               g3(:,3) = g3(:,3) + tmp1*rik
               sigma = sigma + outer_prod_3x3(rik,rik)*tmp1
               !
               !start calculating the derivatives of each part w.r.t. r_ik
               r=sqrt(rjk2)

               dang=-0.375_wp*(rjk2**3+rjk2**2*(rik2+rij2) &
                  +rjk2*(3.0_wp*rik2**2+2.0*rik2*rij2+3.0*rij2**2) &
                  -5.0*(rik2-rij2)**2*(rik2+rij2)) &
                  /(r*rijk3*rijk2)

               tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)*0.5_wp/r

               g3(:,2) = g3(:,2) - tmp1*rjk
               g3(:,3) = g3(:,3) + tmp1*rjk
               sigma = sigma + outer_prod_3x3(rjk,rjk)*tmp1

            enddo !sz
            enddo !sy
            enddo !sx

         enddo !tz
         enddo !ty
         enddo !tx
         eabc = eabc - r3*c9
         g(:,iat) = g(:,iat) + g3(:,1)
         g(:,jat) = g(:,jat) + g3(:,2)
         g(:,kat) = g(:,kat) + g3(:,3)
         !calculating the CN derivative dE_disp(ijk)/dCN(i)
         dc9=0.5_wp*c9*(dc6ab(iat,jat)/c6ij+dc6ab(iat,kat)/c6ik)
         dc6dcn(iat) = dc6dcn(iat) + r3*dc9
         dc9=0.5_wp*c9*(dc6ab(jat,iat)/c6ij+dc6ab(jat,kat)/c6jk)
         dc6dcn(jat) = dc6dcn(jat) + r3*dc9
         dc9=0.5_wp*c9*(dc6ab(kat,iat)/c6ik+dc6ab(kat,jat)/c6jk)
         dc6dcn(kat) = dc6dcn(kat) + r3*dc9
      enddo kAt_iik
   enddo iAt_iik
!$omp enddo
   ! And now kat=jat, but cycling throug all imagecells without t=s. and jat>iat going though all cells    (iat,jat,jat)
   ! But this counts only 1/2

!$omp do schedule(dynamic)
   iAt_ijj: do iat=2,nat
      jAt_ijj: do jat=1,iat-1
         kat=jat
         ij=iat*(iat-1)/2+jat
         jk=jat*(jat-1)/2+kat
         ik=ij

         c6ij=c6ab(ij)
         c6ik=c6ij

         c6jk=c6ab(jk)
         ikvec=xyz(:,kat)-xyz(:,iat)
         ijvec=ikvec
         jkvec=0.0_wp
         cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
         cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
         cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
         cijk = cij*cjk*cik

         c9=-par%s9*sqrt(c6ij*c6ik*c6jk)
         g3 = 0.0_wp
         r3 = 0.0_wp
         do tx=-rep(1), rep(1)
         do ty=-rep(2), rep(2)
         do tz=-rep(3), rep(3)
            repmin(1)=max(-rep(1),tx-rep(1))
            repmax(1)=min(+rep(1),tx+rep(1))
            repmin(2)=max(-rep(2),ty-rep(2))
            repmax(2)=min(+rep(2),ty+rep(2))
            repmin(3)=max(-rep(3),tz-rep(3))
            repmax(3)=min(+rep(3),tz+rep(3))

            t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
            rij=ijvec+t
            rij2=SUM(rij*rij)
            if(rij2.gt.thr)cycle

            do sx=repmin(1), repmax(1)
            do sy=repmin(2), repmax(2)
            do sz=repmin(3), repmax(3)
               ! every result * 0.5
               if (tx.eq.sx .and. ty.eq.sy .and. tz.eq.sz) cycle
               s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
               rik=ikvec+s
               rik2=SUM(rik*rik)
               if(rik2.gt.thr)cycle

               rjk=jkvec+s-t
               rjk2=SUM(rjk*rjk)
               if(rjk2.gt.thr)cycle

               rijk2=(rij2*rjk2*rik2)

               rijk=sqrt(rijk2)
               fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
               rijk3=rijk*rijk2
               ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                  *(-rij2+rjk2+rik2)/(rijk2*rijk3) &
                  +1.0_wp/(rijk3)
               r9ijk=ang*fdmp/2.0_wp   !factor 1/2 for doublecounting
               r3 = r3 + r9ijk

               !              jat=kat
               !dfdmp=2._wp*alp9*(0.75_wp*r0av)**(alp9)*fdmp*fdmp
               dfdmp = -(oth*six*par%alp*((cijk/rijk)**oth)**par%alp)*fdmp**2
               !start calculating the derivatives of each part w.r.t. r_ij
               r=sqrt(rij2)

               dang=-0.375_wp*(rij2**3+rij2**2*(rjk2+rik2) &
                  +rij2*(3.0_wp*rjk2**2+2.0_wp*rjk2*rik2+3.0_wp*rik2**2) &
                  -5.0_wp*(rjk2-rik2)**2*(rjk2+rik2)) &
                  /(r*rijk3*rijk2)

               tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)*0.5_wp/r
               g3(:,1) = g3(:,1) - tmp1*rij
               g3(:,2) = g3(:,2) + tmp1*rij
               sigma = sigma + outer_prod_3x3(rij,rij)*tmp1

               !start calculating the derivatives of each part w.r.t. r_ik
               r=sqrt(rik2)


               dang=-0.375_wp*(rik2**3+rik2**2*(rjk2+rij2) &
                  +rik2*(3.0_wp*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2) &
                  -5.0*(rjk2-rij2)**2*(rjk2+rij2)) &
                  /(r*rijk3*rijk2)

               tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)*0.5_wp/r
               g3(:,1) = g3(:,1) - tmp1*rik
               g3(:,3) = g3(:,3) + tmp1*rik
               sigma = sigma + outer_prod_3x3(rik,rik)*tmp1
               !
               !start calculating the derivatives of each part w.r.t. r_jk
               r=sqrt(rjk2)

               dang=-0.375_wp*(rjk2**3+rjk2**2*(rik2+rij2) &
                  +rjk2*(3.0_wp*rik2**2+2.0*rik2*rij2+3.0*rij2**2) &
                  -5.0_wp*(rik2-rij2)**2*(rik2+rij2)) &
                  /(r*rijk3*rijk2)

               tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)*0.5_wp/r
               sigma = sigma + outer_prod_3x3(rjk,rjk)*tmp1

            enddo !sz
            enddo !sy
            enddo !sx
         enddo !tz
         enddo !ty
         enddo !tx
         eabc = eabc - r3*c9
         g(:,iat) = g(:,iat) + g3(:,1)
         g(:,jat) = g(:,jat) + g3(:,2)
         g(:,kat) = g(:,kat) + g3(:,3)
         !calculating the CN derivative dE_disp(ijk)/dCN(i)
         dc9=0.5_wp*c9*(dc6ab(iat,jat)/c6ij+dc6ab(iat,kat)/c6ik)
         dc6dcn(iat) = dc6dcn(iat) + r3*dc9
         dc9=0.5_wp*c9*(dc6ab(jat,iat)/c6ij+dc6ab(jat,kat)/c6jk)
         dc6dcn(jat) = dc6dcn(jat) + r3*dc9
         dc9=0.5_wp*c9*(dc6ab(kat,iat)/c6ik+dc6ab(kat,jat)/c6jk)
         dc6dcn(kat) = dc6dcn(kat) + r3*dc9
      enddo jAt_ijj
   enddo iAt_ijj
!$omp enddo

   ! And finally the self interaction iat=jat=kat all

!$omp do schedule(dynamic)
   iAt_iii: do iat=1,nat
      jat=iat
      kat=iat
      ijvec=0.0_wp
      ij=iat*(iat-1)/2+jat
      ik=iat*(iat-1)/2+kat
      jk=jat*(jat-1)/2+kat
      ikvec=ijvec
      jkvec=ikvec
      c6ij=c6ab(ij)
      c6ik=c6ij
      c6jk=c6ij
      c9=-par%s9*sqrt(c6ij*c6ij*c6ij)
      cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
      cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
      cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
      cijk = cij*cjk*cik

      r3 = 0.0_wp
      do tx=-rep(1), rep(1)
      do ty=-rep(2), rep(2)
      do tz=-rep(3), rep(3)
         repmin(1)=max(-rep(1),tx-rep(1))
         repmax(1)=min(+rep(1),tx+rep(1))
         repmin(2)=max(-rep(2),ty-rep(2))
         repmax(2)=min(+rep(2),ty+rep(2))
         repmin(3)=max(-rep(3),tz-rep(3))
         repmax(3)=min(+rep(3),tz+rep(3))
         ! IF iat and jat are the same then cycle
         if ((tx.eq.0) .and.(ty.eq.0) .and.(tz.eq.0))cycle
         t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
         rij=t
         rij2=SUM(rij*rij)
         if(rij2.gt.thr)cycle

         do sx=repmin(1), repmax(1)
         do sy=repmin(2), repmax(2)
         do sz=repmin(3), repmax(3)
            ! if iat and kat are the same then cycle
            if ((sx.eq.0) .and.( sy.eq.0) .and.( sz.eq.0))cycle
            ! If kat and jat are the same then cycle
            if ((sx.eq.tx) .and. (sy.eq.ty) .and. (sz.eq.tz)) cycle

            ! every result * 1/6 becaues every triple is counted twice due
            ! to the two loops t and s going from -rep to rep -> *1/2
            !
            !plus 1/3 becaues every triple is three times in each unitcell
            s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
            rik=s
            rik2=SUM(rik*rik)
            if(rik2.gt.thr)cycle

            rjk=jkvec+s-t
            rjk2=SUM(rjk*rjk)
            if(rjk2.gt.thr)cycle

            rijk2=(rij2*rjk2*rik2)

            rijk=sqrt(rijk2)
            fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
            rijk3=rijk*rijk2
            ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
               *(-rij2+rjk2+rik2)/(rijk2*rijk3) &
               +1.0_wp/(rijk3)
            r9ijk=ang*fdmp/6.0_wp
            r3 = r3 + r9ijk

            !                          iat=jat=kat
            dfdmp = -(oth*six*par%alp*((cijk/rijk)**oth)**par%alp)*fdmp**2
            !start calculating the derivatives of each part w.r.t. r_ij

            r=sqrt(rij2)
            dang=-0.375_wp*(rij2**3+rij2**2*(rjk2+rik2) &
               +rij2*(3.0_wp*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2) &
               -5.0*(rjk2-rik2)**2*(rjk2+rik2)) &
               /(r*rijk3*rijk2)


            tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)/(r*6.0_wp)
            sigma = sigma + outer_prod_3x3(rij,rij)*tmp1

            !start calculating the derivatives of each part w.r.t. r_ik

            r=sqrt(rik2)

            dang=-0.375_wp*(rik2**3+rik2**2*(rjk2+rij2) &
               +rik2*(3.0_wp*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2) &
               -5.0*(rjk2-rij2)**2*(rjk2+rij2)) &
               /(r*rijk3*rijk2)

            tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)/(r*6.0_wp)
            sigma = sigma + outer_prod_3x3(rik,rik)*tmp1
            !
            !start calculating the derivatives of each part w.r.t. r_jk

            r=sqrt(rjk2)
            dang=-0.375_wp*(rjk2**3+rjk2**2*(rik2+rij2) &
               +rjk2*(3.0_wp*rik2**2+2.0*rik2*rij2+3.0*rij2**2) &
               -5.0*(rik2-rij2)**2*(rik2+rij2)) &
               /(r*rijk3*rijk2)

            tmp1=(-dang*c9*fdmp+dfdmp/r*c9*ang)/(r*6.0_wp)
            sigma = sigma + outer_prod_3x3(rjk,rjk)*tmp1

         enddo !sz
         enddo !sy
         enddo !sx
      enddo !tz
      enddo !ty
      enddo !tx
      eabc = eabc - r3*c9
      !calculating the CN derivative dE_disp(ijk)/dCN(i)
      dc9=0.5_wp*c9*(dc6ab(iat,jat)/c6ij+dc6ab(iat,kat)/c6ik)
      dc6dcn(iat) = dc6dcn(iat) + r3*dc9
      dc9=0.5_wp*c9*(dc6ab(jat,iat)/c6ij+dc6ab(jat,kat)/c6jk)
      dc6dcn(jat) = dc6dcn(jat) + r3*dc9
      dc9=0.5_wp*c9*(dc6ab(kat,iat)/c6ik+dc6ab(kat,jat)/c6jk)
      dc6dcn(kat) = dc6dcn(kat) + r3*dc9
   enddo iAt_iii
!$omp enddo
!$omp end parallel

end subroutine dabcappr_3d_dftd3_like_style

end module dftd4
