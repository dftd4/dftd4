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

!> provides different kinds of coordination numbers used in this program
!
!  Implemented is the original DFT-D3 coordination number, a modified version
!  using a better counting function and the DFT-D4 coordination number.
!  The derivative is given as dCNi/dRj -> dcndr(:,j,i), so matrix-vector
!  operations can be used to obtain the molecular gradient
module coordination_number
   use iso_fortran_env, only : wp => real64
   use mctc_param
   implicit none

   real(wp),private,parameter :: cnthr = 1600.0_wp

   real(wp),parameter :: k1 = 16.0_wp !< steepness of counting function
   real(wp),parameter :: k2 = 4.0_wp/3.0_wp

   real(wp),parameter :: ka=10.0_wp !< steepness of first counting function
   real(wp),parameter :: kb=20.0_wp !< steepness of second counting function
   real(wp),parameter :: r_shift=2.0_wp !< offset for second counting function

   real(wp),parameter :: k4=4.10451_wp  !< EN scaling parameter in covCN
   real(wp),parameter :: k5=19.08857_wp !< EN scaling parameter in covCN
   real(wp),parameter :: k6=2*11.28174_wp**2 !< EN scaling parameter in covCN

   real(wp),parameter :: kn=7.50_wp !< steepness of counting function
   real(wp),parameter :: kr=0.25_wp
   real(wp),parameter :: ke=0.05_wp

   integer,private,parameter :: max_elem = 118

   !> pauling EN's
   real(wp),parameter :: en(max_elem) = (/ &
   & 2.20_wp,3.00_wp, & ! H,He
   & 0.98_wp,1.57_wp,2.04_wp,2.55_wp,3.04_wp,3.44_wp,3.98_wp,4.50_wp, & ! Li-Ne
   & 0.93_wp,1.31_wp,1.61_wp,1.90_wp,2.19_wp,2.58_wp,3.16_wp,3.50_wp, & ! Na-Ar
   & 0.82_wp,1.00_wp, & ! K,Ca
   &           1.36_wp,1.54_wp,1.63_wp,1.66_wp,1.55_wp, &
   &           1.83_wp,1.88_wp,1.91_wp,1.90_wp,1.65_wp, & ! Sc-Zn
   &           1.81_wp,2.01_wp,2.18_wp,2.55_wp,2.96_wp,3.00_wp, & ! Ga-Kr
   & 0.82_wp,0.95_wp, & ! Rb,Sr
   &           1.22_wp,1.33_wp,1.60_wp,2.16_wp,1.90_wp, &
   &           2.20_wp,2.28_wp,2.20_wp,1.93_wp,1.69_wp, & ! Y-Cd
   &           1.78_wp,1.96_wp,2.05_wp,2.10_wp,2.66_wp,2.60_wp, & ! In-Xe
   & 0.79_wp,0.89_wp, & ! Cs,Ba
   &      1.10_wp,1.12_wp,1.13_wp,1.14_wp,1.15_wp,1.17_wp,1.18_wp, & ! La-Eu
   &      1.20_wp,1.21_wp,1.22_wp,1.23_wp,1.24_wp,1.25_wp,1.26_wp, & ! Gd-Yb
   &           1.27_wp,1.30_wp,1.50_wp,2.36_wp,1.90_wp, &
   &           2.20_wp,2.20_wp,2.28_wp,2.54_wp,2.00_wp, & ! Lu-Hg
   &           1.62_wp,2.33_wp,2.02_wp,2.00_wp,2.20_wp,2.20_wp, & ! Tl-Rn
   ! only dummies below
   & 1.50_wp,1.50_wp, & ! Fr,Ra
   &      1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Ac-Am
   &      1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Cm-No
   &           1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, &
   &           1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Rf-Cn
   &           1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp /) ! Nh-Og


contains

! ========================================================================
!> original D3 type coordination number from 2010
pure subroutine ncoord_d3(mol,cn,thr)
   use class_molecule

   implicit none

   !> molecular structure information
   type(molecule),intent(in) :: mol
   !> coordination number
   real(wp),intent(out) :: cn(mol%n)
   !> cutoff threshold for neglecting CN count
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j
   real(wp) :: rij(3), r, rco, den, rr, r2

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn = 0.0_wp

   do i = 1, mol%n
      do j = 1, i-1
         rij = mol%xyz(:,j) - mol%xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r=sqrt(r2)
!        covalent distance in bohr
         rco=k2*(covalent_radius(mol%at(j)) + covalent_radius(mol%at(i)))
         rr=rco/r
!        counting function exponential has a better long-range
!        behavior than MHGs inverse damping
         cn(i)=cn(i)+1.0_wp/(1.0_wp+exp(-k1*(rr-1.0_wp)))
         cn(j)=cn(j)+1.0_wp/(1.0_wp+exp(-k1*(rr-1.0_wp)))
      enddo
   enddo

end subroutine ncoord_d3

! ========================================================================
!> original D3 type coordination number from 2010
pure subroutine dncoord_d3(mol,cn,dcndr,thr)
   use class_molecule

   implicit none

   !> molecular structure information
   type(molecule),intent(in) :: mol
   !> coordination number
   real(wp),intent(out) :: cn(mol%n)
   !> derivative of the coordination number w.r.t. nuclear positions
   real(wp),intent(out) :: dcndr(3,mol%n,mol%n)
   !> cutoff threshold for neglecting CN count
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i, j
   real(wp) :: r, r2, rij(3)
   real(wp) :: rcovij
   real(wp) :: expterm
   real(wp) :: dtmp, tmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn  = 0._wp
   dcndr = 0._wp

   do i = 1, mol%n
      do j = 1, i-1
         rij = mol%xyz(:,j) - mol%xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r = sqrt(r2)
         rcovij=k2*(covalent_radius(mol%at(i))+covalent_radius(mol%at(j)))
         expterm=exp(-k1*(rcovij/r-1._wp))
         tmp = 1._wp/(1._wp+expterm)
         dtmp = (-k1*rcovij*expterm)/(r2*((expterm+1._wp)**2))
         cn(i) = cn(i) + tmp
         cn(j) = cn(j) + tmp
         dcndr(:,i,i)=-dtmp*rij/r + dcndr(:,i,i)
         dcndr(:,j,j)= dtmp*rij/r + dcndr(:,j,j)
         dcndr(:,i,j)= dtmp*rij/r
         dcndr(:,j,i)=-dtmp*rij/r
      enddo
   enddo

end subroutine dncoord_d3

!> gradients for pbc coordination number with error function
subroutine pbc_ncoord_erf(mol,cn,thr)
   use iso_fortran_env, wp => real64
   use class_molecule
   use pbc_tools, only : get_realspace_cutoff, outer_prod_3x3

   implicit none

   !> molecular structure information
   type(molecule),intent(in) :: mol
   !> coordination number
   real(wp),intent(out) :: cn(mol%n)
   !> cutoff threshold for neglecting CN count
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j,tx,ty,tz
   real(wp) :: rij(3), r, rco, den, rr, r2, t(3)
   integer  :: rep_cn(3)
   real(wp) :: tmp,dtmp
   real(wp),parameter :: sqrtpi = 1.77245385091_wp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = 1600.0_wp
   endif

   if (mol%npbc > 0) &
   call get_realspace_cutoff(mol%lattice,cn_thr,rep_cn)
   where(.not.mol%pbc) rep_cn = 0

   cn = 0.0_wp

   do i = 1, mol%n
      do j = 1, i-1 ! loop over all atoms for PBC case
         do concurrent(tx = -rep_cn(1):rep_cn(1),&
               &       ty = -rep_cn(2):rep_cn(2),&
               &       tz = -rep_cn(3):rep_cn(3))
            t = [tx,ty,tz]
            rij = mol%xyz(:,j) - mol%xyz(:,i) + matmul(mol%lattice,t)
            r2  = sum(rij**2)
            if (r2.gt.cn_thr) cycle
            r=sqrt(r2)
!           covalent distance in bohr
            rco=k2*(covalent_radius(mol%at(j)) + covalent_radius(mol%at(i)))
            tmp=0.5_wp*(1.0_wp+erf(-kn*(r-rco)/rco))
            dtmp = -kn/sqrtpi/rco*exp(-kn**2*(r-rco)**2/rco**2)
            cn(i)=cn(i) + tmp
            cn(j)=cn(j) + tmp
         enddo
      enddo
      do concurrent(tx = -rep_cn(1):rep_cn(1),&
            &       ty = -rep_cn(2):rep_cn(2),&
            &       tz = -rep_cn(3):rep_cn(3))
         ! avoid self interaction
         if ((tx.eq.0).and.(ty.eq.0).and.(tz.eq.0)) cycle
         t = [tx,ty,tz]
         rij = matmul(mol%lattice,t)
         r2  = sum(rij**2)
         if (r2.gt.cn_thr) cycle
         r=sqrt(r2)
         ! covalent distance in bohr
         rco=k2*2*covalent_radius(mol%at(i))
         tmp=0.5_wp*(1.0_wp+erf(-kn*(r-rco)/rco))
         dtmp = -kn/sqrtpi/rco*exp(-kn**2*(r-rco)**2/rco**2)
         cn(i) = cn(i) + tmp
      enddo
   enddo

end subroutine pbc_ncoord_erf

!> gradients for pbc coordination number with error function
subroutine pbc_dncoord_erf(mol,cn,dcndr,dcndL,thr)
   use iso_fortran_env, wp => real64
   use class_molecule
   use pbc_tools, only : get_realspace_cutoff, outer_prod_3x3

   implicit none

   !> molecular structure information
   type(molecule),intent(in) :: mol
   !> coordination number
   real(wp),intent(out) :: cn(mol%n)
   !> derivative of the coordination number w.r.t. nuclear positions
   real(wp),intent(out) :: dcndr(3,mol%n,mol%n)
   !> derivative of the coordination number w.r.t. lattice parameters
   real(wp),intent(out) :: dcndL(3,3,mol%n)
   !> cutoff threshold for neglecting CN count
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j,tx,ty,tz
   real(wp) :: rij(3), r, rco, den, rr, r2, t(3)
   integer  :: rep_cn(3)
   real(wp) :: tmp,dtmp
   real(wp),parameter :: sqrtpi = 1.77245385091_wp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = 1600.0_wp
   endif

   if (mol%npbc > 0) &
   call get_realspace_cutoff(mol%lattice,cn_thr,rep_cn)
   where(.not.mol%pbc) rep_cn = 0

   cn = 0.0_wp
   dcndr = 0.0_wp
   dcndL = 0.0_wp

   do i = 1, mol%n
      do j = 1, i-1 ! loop over all atoms for PBC case
         do concurrent(tx = -rep_cn(1):rep_cn(1),&
               &       ty = -rep_cn(2):rep_cn(2),&
               &       tz = -rep_cn(3):rep_cn(3))
            t = [tx,ty,tz]
            rij = mol%xyz(:,j) - mol%xyz(:,i) + matmul(mol%lattice,t)
            r2  = sum(rij**2)
            if (r2.gt.cn_thr) cycle
            r=sqrt(r2)
!           covalent distance in bohr
            rco=k2*(covalent_radius(mol%at(j)) + covalent_radius(mol%at(i)))
            tmp=0.5_wp*(1.0_wp+erf(-kn*(r-rco)/rco))
            dtmp = -kn/sqrtpi/rco*exp(-kn**2*(r-rco)**2/rco**2)
            cn(i)=cn(i) + tmp
            cn(j)=cn(j) + tmp
            dcndr(:,i,i)= dtmp*rij/r + dcndr(:,i,i)
            dcndr(:,j,j)=-dtmp*rij/r + dcndr(:,j,j)
            dcndr(:,i,j)= dtmp*rij/r + dcndr(:,i,j)
            dcndr(:,j,i)=-dtmp*rij/r + dcndr(:,j,i)
            dcndL(:,:,j)= dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,j)
            dcndL(:,:,i)= dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,i)
         enddo
      enddo
      do concurrent(tx = -rep_cn(1):rep_cn(1),&
            &       ty = -rep_cn(2):rep_cn(2),&
            &       tz = -rep_cn(3):rep_cn(3))
         ! avoid self interaction
         if ((tx.eq.0).and.(ty.eq.0).and.(tz.eq.0)) cycle
         t = [tx,ty,tz]
         rij = matmul(mol%lattice,t)
         r2  = sum(rij**2)
         if (r2.gt.cn_thr) cycle
         r=sqrt(r2)
         ! covalent distance in bohr
         rco=k2*2*covalent_radius(mol%at(i))
         tmp=0.5_wp*(1.0_wp+erf(-kn*(r-rco)/rco))
         dtmp = -kn/sqrtpi/rco*exp(-kn**2*(r-rco)**2/rco**2)
         cn(i) = cn(i) + tmp
         dcndL(:,:,i) = dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,i)
      enddo
   enddo

end subroutine pbc_dncoord_erf

!> gradients for periodic covalent coordination number
pure subroutine pbc_ncoord_d4(mol,cn,thr)
   use mctc_constants
   use class_molecule
   use pbc_tools, only : outer_prod_3x3, get_realspace_cutoff
   implicit none

   !> molecular structure information
   type(molecule),intent(in) :: mol
   !> number of images to consider for lattice translations
   integer  :: rep_cn(3)
   integer  :: tx,ty,tz
   real(wp) :: t(3)
   !> coordination number
   real(wp),intent(out) :: cn(mol%n)
   !> cutoff threshold for neglecting CN count
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j
   real(wp) :: rij(3), r, rco, den, r2, xn, tmp, dtmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   if (mol%npbc > 0) &
   call get_realspace_cutoff(mol%lattice,cn_thr,rep_cn)
   where(.not.mol%pbc) rep_cn = 0

   cn  = 0.0_wp

   do i=1,mol%n
      do j=1,i-1
         do concurrent(tx = -rep_cn(1):rep_cn(1),&
               &       ty = -rep_cn(2):rep_cn(2),&
               &       tz = -rep_cn(3):rep_cn(3))
            t = [tx,ty,tz]
            rij = mol%xyz(:,j) - mol%xyz(:,i) + matmul(mol%lattice,t)
            r2  = sum(rij**2)
            if (r2.gt.cn_thr) cycle
            r=sqrt(r2)
            ! covalent distance in bohr
            rco=k2*(covalent_radius(mol%at(j)) + covalent_radius(mol%at(i)))
            den=k4*exp(-(abs((en(mol%at(i))-en(mol%at(j))))+ k5)**2/k6 )
            tmp = den * 0.5_wp * (1.0_wp + erf(-kn*(r-rco)/rco))
            dtmp = -den*kn/sqrtpi/rco*exp(-kn**2*(r-rco)**2/rco**2)
            cn(i)=cn(i)+tmp
            cn(j)=cn(j)+tmp
         enddo
      enddo
      do concurrent(tx = -rep_cn(1):rep_cn(1),&
            &       ty = -rep_cn(2):rep_cn(2),&
            &       tz = -rep_cn(3):rep_cn(3))
         ! avoid self interaction
         if ((tx.eq.0).and.(ty.eq.0).and.(tz.eq.0)) cycle
         t = [tx,ty,tz]
         rij = matmul(mol%lattice,t)
         r2  = sum(rij**2)
         if (r2.gt.cn_thr) cycle
         r=sqrt(r2)
         ! covalent distance in bohr
         rco=k2*2*covalent_radius(mol%at(i))
         den=k4*exp(-k5**2/k6)
         tmp = den * 0.5_wp * (1.0_wp + erf(-kn*(r-rco)/rco))
         dtmp = -den*kn/sqrtpi/rco*exp(-kn**2*(r-rco)**2/rco**2)
         cn(i)=cn(i)+tmp
      enddo
   enddo

end subroutine pbc_ncoord_d4

!> gradients for periodic covalent coordination number
pure subroutine pbc_dncoord_d4(mol,cn,dcndr,dcndL,thr)
   use mctc_constants
   use class_molecule
   use pbc_tools, only : outer_prod_3x3, get_realspace_cutoff
   implicit none

   !> molecular structure information
   type(molecule),intent(in) :: mol
   !> number of images to consider for lattice translations
   integer  :: rep_cn(3)
   integer  :: tx,ty,tz
   real(wp) :: t(3)
   !> coordination number
   real(wp),intent(out) :: cn(mol%n)
   !> derivative of the coordination number w.r.t. nuclear positions
   real(wp),intent(out) :: dcndr(3,mol%n,mol%n)
   !> derivative of the coordination number w.r.t. lattice parameters
   real(wp),intent(out) :: dcndL(3,3,mol%n)
   !> cutoff threshold for neglecting CN count
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j
   real(wp) :: rij(3), r, rco, den, r2, xn, tmp, dtmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   if (mol%npbc > 0) &
   call get_realspace_cutoff(mol%lattice,cn_thr,rep_cn)
   where(.not.mol%pbc) rep_cn = 0

   cn  = 0.0_wp
   dcndr = 0.0_wp
   dcndL = 0.0_wp

   do i=1,mol%n
      do j=1,i-1
         do concurrent(tx = -rep_cn(1):rep_cn(1),&
               &       ty = -rep_cn(2):rep_cn(2),&
               &       tz = -rep_cn(3):rep_cn(3))
            t = [tx,ty,tz]
            rij = mol%xyz(:,j) - mol%xyz(:,i) + matmul(mol%lattice,t)
            r2  = sum(rij**2)
            if (r2.gt.cn_thr) cycle
            r=sqrt(r2)
            ! covalent distance in bohr
            rco=k2*(covalent_radius(mol%at(j)) + covalent_radius(mol%at(i)))
            den=k4*exp(-(abs((en(mol%at(i))-en(mol%at(j))))+ k5)**2/k6 )
            tmp = den * 0.5_wp * (1.0_wp + erf(-kn*(r-rco)/rco))
            dtmp = -den*kn/sqrtpi/rco*exp(-kn**2*(r-rco)**2/rco**2)
            cn(i)=cn(i)+tmp
            cn(j)=cn(j)+tmp
            dcndr(:,i,i)=-dtmp*rij/r + dcndr(:,i,i)
            dcndr(:,j,j)= dtmp*rij/r + dcndr(:,j,j)
            dcndr(:,i,j)=-dtmp*rij/r + dcndr(:,i,j)
            dcndr(:,j,i)= dtmp*rij/r + dcndr(:,j,i)
            dcndL(:,:,j)= dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,j)
            dcndL(:,:,i)= dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,i)
         enddo
      enddo
      do concurrent(tx = -rep_cn(1):rep_cn(1),&
            &       ty = -rep_cn(2):rep_cn(2),&
            &       tz = -rep_cn(3):rep_cn(3))
         ! avoid self interaction
         if ((tx.eq.0).and.(ty.eq.0).and.(tz.eq.0)) cycle
         t = [tx,ty,tz]
         rij = matmul(mol%lattice,t)
         r2  = sum(rij**2)
         if (r2.gt.cn_thr) cycle
         r=sqrt(r2)
         ! covalent distance in bohr
         rco=k2*2*covalent_radius(mol%at(i))
         den=k4*exp(-k5**2/k6)
         tmp = den * 0.5_wp * (1.0_wp + erf(-kn*(r-rco)/rco))
         dtmp = -den*kn/sqrtpi/rco*exp(-kn**2*(r-rco)**2/rco**2)
         cn(i)=cn(i)+tmp
         dcndL(:,:,i)= dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,i)
      enddo
   enddo

end subroutine pbc_dncoord_d4

!> cutoff function for large coordination numbers
subroutine dncoord_logcn(n,cn,dcndr,dcndL,cn_max)
   implicit none
   !> number of atoms
   integer, intent(in) :: n
   !> on input coordination number, on output modified CN
   real(wp), intent(inout) :: cn(n)
   !> on input derivative of CN w.r.t. cartesian coordinates,
   !  on output derivative of modified CN
   real(wp), intent(inout), optional :: dcndr(3,n,n)
   !> on input derivative of CN w.r.t. strain deformation,
   !  on output derivative of modified CN
   real(wp), intent(inout), optional :: dcndL(3,3,n)
   !> maximum CN (not strictly obeyed)
   real(wp), intent(in), optional :: cn_max
   !  local
   real(wp) :: cnmax
   real(wp) :: dcnpdcn
   integer  :: i

   if (present(cn_max)) then
      cnmax = max(cn_max,0.0_wp)
   else
      cnmax = 4.5_wp
   endif

   if (present(dcndL)) then
      do i = 1, n
         dcnpdcn = dlog_cn_cut(cn(i),cnmax)
         dcndL(:,:,i) = dcnpdcn*dcndL(:,:,i)
      enddo
   endif

   if (present(dcndr)) then
      do i = 1, n
         dcnpdcn = dlog_cn_cut(cn(i),cnmax)
         dcndr(:,:,i) = dcnpdcn*dcndr(:,:,i)
      enddo
   endif

   do i = 1, n
      cn(i) = log_cn_cut(cn(i),cnmax)
   enddo

end subroutine dncoord_logcn

pure elemental function log_cn_cut(cn,cnmax) result(cnp)
   real(wp), intent(in) :: cn
   real(wp), intent(in) :: cnmax
   real(wp) :: cnp
   cnp = log(1.0_wp + exp(cnmax)) - log(1.0_wp + exp(cnmax - cn))
end function log_cn_cut

pure elemental function dlog_cn_cut(cn,cnmax) result(dcnpdcn)
   real(wp), intent(in) :: cn
   real(wp), intent(in) :: cnmax
   real(wp) :: dcnpdcn
   dcnpdcn = exp(cnmax)/(exp(cnmax) + exp(cn))
end function dlog_cn_cut

pure elemental function erf_count(k,r,r0) result(count)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   count = 0.5_wp * (1.0_wp + erf(-k*(r-r0)/r0))
end function erf_count

pure elemental function derf_count(k,r,r0) result(count)
   use mctc_constants
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   count = -k/sqrtpi/r0*exp(-k**2*(r-r0)**2/r0**2)
end function derf_count

pure elemental function exp_count(k,r,r0) result(count)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   count =1.0_wp/(1.0_wp+exp(-k*(r0/r-1.0_wp)))
end function exp_count

pure elemental function dexp_count(k,r,r0) result(count)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   real(wp) :: expterm
   expterm=exp(-k*(r0/r-1._wp))
   count = (-k*r0*expterm)/(r**2*((expterm+1._wp)**2))
end function dexp_count

end module coordination_number
