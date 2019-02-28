module coordination_number
   use iso_fortran_env, only : wp => real64
   implicit none

   real(wp),private,parameter :: cnthr = 1600.0_wp

   real(wp),parameter :: k1 = 16.0_wp

   real(wp),parameter :: ka=10.0_wp
   real(wp),parameter :: kb=20.0_wp
   real(wp),parameter :: r_shift=2.0_wp

   real(wp),parameter :: k4=4.10451_wp
   real(wp),parameter :: k5=19.08857_wp
   real(wp),parameter :: k6=2*11.28174_wp**2

   real(wp),parameter :: kn=7.50_wp
   real(wp),parameter :: kr=0.25_wp
   real(wp),parameter :: ke=0.05_wp

   integer,private,parameter :: max_elem = 118
!  D3 radii
!   real(wp),parameter :: rcov(max_elem) = (/ &
!  & 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865, &
!  & 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527, &
!  & 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820, &
!  & 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730, &
!  & 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923, &
!  & 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188, &
!  & 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349, &
!  & 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216, &
!  & 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717, &
!  & 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967, &
!  & 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625, &
!  & 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657, &
!  & 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833, &
!  & 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098, &
!  & 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878, &
!  & 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790, &
!  & 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584, &
!  & 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289, &
!  & 3.82984466, 3.85504098, 3.88023730, 3.90543362 /)

!  covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
!  188-197), values for metals decreased by 10 %
   real(wp),private,parameter :: rad(max_elem) = (/  &
   & 0.32,0.46, & ! H,He
   & 1.20,0.94,0.77,0.75,0.71,0.63,0.64,0.67, & ! Li-Ne
   & 1.40,1.25,1.13,1.04,1.10,1.02,0.99,0.96, & ! Na-Ar
   & 1.76,1.54, & ! K,Ca
   &           1.33,1.22,1.21,1.10,1.07,1.04,1.00,0.99,1.01,1.09, & ! Sc-Zn
   &           1.12,1.09,1.15,1.10,1.14,1.17, & ! Ga-Kr
   & 1.89,1.67, & ! Rb,Sr
   &           1.47,1.39,1.32,1.24,1.15,1.13,1.13,1.08,1.15,1.23, & ! Y-Cd
   &           1.28,1.26,1.26,1.23,1.32,1.31, & ! In-Xe
   & 2.09,1.76, & ! Cs,Ba
   &      1.62,1.47,1.58,1.57,1.56,1.55,1.51, & ! La-Eu
   &      1.52,1.51,1.50,1.49,1.49,1.48,1.53, & ! Gd-Yb
   &           1.46,1.37,1.31,1.23,1.18,1.16,1.11,1.12,1.13,1.32, & ! Lu-Hg
   &           1.30,1.30,1.36,1.31,1.38,1.42, & ! Tl-Rn
   & 2.01,1.81, & ! Fr,Ra
   &      1.67,1.58,1.52,1.53,1.54,1.55,1.49, & ! Ac-Am
   &      1.49,1.51,1.51,1.48,1.50,1.56,1.58, & ! Cm-No
   &           1.45,1.41,1.34,1.29,1.27,1.21,1.16,1.15,1.09,1.22, & ! Lr-Cn
   &           1.22,1.29,1.46,1.58,1.48,1.41 /) ! Nh-Og
   real(wp),parameter :: rcov(max_elem) = 4.0_wp/3.0_wp*rad/0.52917726_wp

!  pauling EN's 
   real(wp),parameter :: en(max_elem) = (/ &
   & 2.20,3.00, & ! H,He
   & 0.98,1.57,2.04,2.55,3.04,3.44,3.98,4.50, & ! Li-Ne
   & 0.93,1.31,1.61,1.90,2.19,2.58,3.16,3.50, & ! Na-Ar
   & 0.82,1.00, & ! K,Ca
   &           1.36,1.54,1.63,1.66,1.55,1.83,1.88,1.91,1.90,1.65, & ! Sc-Zn
   &           1.81,2.01,2.18,2.55,2.96,3.00, & ! Ga-Kr
   & 0.82,0.95, & ! Rb,Sr
   &           1.22,1.33,1.60,2.16,1.90,2.20,2.28,2.20,1.93,1.69, & ! Y-Cd
   &           1.78,1.96,2.05,2.10,2.66,2.60, & ! In-Xe
   & 0.79,0.89, & ! Cs,Ba
   &      1.10,1.12,1.13,1.14,1.15,1.17,1.18, & ! La-Eu
   &      1.20,1.21,1.22,1.23,1.24,1.25,1.26, & ! Gd-Yb
   &           1.27,1.30,1.50,2.36,1.90,2.20,2.20,2.28,2.54,2.00, & ! Lu-Hg
   &           1.62,2.33,2.02,2.00,2.20,2.20, & ! Tl-Rn
   ! only dummies below
   & 1.50,1.50, & ! Fr,Ra
   &      1.50,1.50,1.50,1.50,1.50,1.50,1.50, & ! Ac-Am
   &      1.50,1.50,1.50,1.50,1.50,1.50,1.50, & ! Cm-No
   &           1.50,1.50,1.50,1.50,1.50,1.50,1.50,1.50,1.50,1.50, & ! Rf-Cn
   &           1.50,1.50,1.50,1.50,1.50,1.50 /) ! Nh-Og


!  test for PBC case: divergence of GFN2-xTB CN [TODO]

contains

!! ========================================================================
!  modified D3 type coordination number from 2018
!  INPUT
!  nat  :: number of atoms
!  at   :: ordinal number of atoms
!  xyz  :: coordinates in Bohr
!  OUTPUT
!  cn   :: coordination number
!  PARAMETERS: kn,k2
!  NOTE: k2 is already included in rcov
pure subroutine ncoord_erf(mol,cn,thr)
   use class_molecule

   implicit none

   type(molecule),intent(in) :: mol
   real(wp),intent(out) :: cn(mol%nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j
   real(wp) :: rij(3), r, rco, den, rr, r2, tmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn = 0.0_wp

   do i = 1, mol%nat
      do j = 1, i-1
         rij = mol%xyz(:,j) - mol%xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle 
         r=sqrt(r2)
!        covalent distance in bohr
         rco=rcov(mol%at(j)) + rcov(mol%at(i))
!        error function has an even better long range behavior
         tmp = 0.5_wp * (1.0_wp + erf(-kn*(r-rco)/rco))
         cn(i)=cn(i)+tmp
         cn(j)=cn(j)+tmp
      enddo
   enddo

end subroutine ncoord_erf

!! ========================================================================
!  modified D3 type coordination number from 2018
!  INPUT
!  nat  :: number of atoms
!  at   :: ordinal number of atoms
!  xyz  :: coordinates in Bohr
!  OUTPUT
!  cn   :: coordination number
!  dcn  :: derivative of coordination number w.r.t. atom position
!  PARAMETERS: kn,k2
!  NOTE: k2 is already included in rcov
pure subroutine dncoord_erf(mol,cn,dcn,thr)
   use class_molecule

   implicit none

   type(molecule),intent(in) :: mol
   real(wp),intent(out) :: cn(mol%nat)
   real(wp),intent(out) :: dcn(3,mol%nat,mol%nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i, j
   real(wp) :: r, r2, rij(3)
   real(wp) :: rcovij
   real(wp) :: dtmp, tmp
   real(wp),parameter :: hlfosqrtpi = 1.0_wp/1.77245385091_wp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn  = 0._wp
   dcn = 0._wp

   do i = 2, mol%nat
      do j = 1, i-1
         rij = mol%xyz(:,j) - mol%xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r = sqrt(r2)
         rcovij=(rcov(mol%at(i))+rcov(mol%at(j)))
         tmp = 0.5_wp * (1.0_wp + erf(-kn*(r-rcovij)/rcovij))
         dtmp =-hlfosqrtpi*kn*exp(-kn**2*(r-rcovij)**2/rcovij**2)/rcovij
         cn(i) = cn(i) + tmp
         cn(j) = cn(j) + tmp
         dcn(:,j,j)= dtmp*rij/r + dcn(:,j,j)
         dcn(:,i,j)= dtmp*rij/r
         dcn(:,j,i)=-dtmp*rij/r
         dcn(:,i,i)=-dtmp*rij/r + dcn(:,i,i)
      enddo
   enddo

end subroutine dncoord_erf

!! ========================================================================
!  original D3 type coordination number from 2010
!  INPUT
!  mol  :: molecular stucture
!  OUTPUT
!  cn   :: coordination number
!  PARAMETERS: k1,k2
!  NOTE: k2 is already included in rcov
pure subroutine ncoord_d3(mol,cn,thr)
   use class_molecule

   implicit none

   type(molecule),intent(in) :: mol
   real(wp),intent(out) :: cn(mol%nat)
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

   do i = 1, mol%nat
      do j = 1, i-1
         rij = mol%xyz(:,j) - mol%xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle 
         r=sqrt(r2)
!        covalent distance in bohr
         rco=rcov(mol%at(j)) + rcov(mol%at(i))
         rr=rco/r
!        counting function exponential has a better long-range 
!        behavior than MHGs inverse damping
         cn(i)=cn(i)+1.0_wp/(1.0_wp+exp(-k1*(rr-1.0_wp)))
         cn(j)=cn(j)+1.0_wp/(1.0_wp+exp(-k1*(rr-1.0_wp)))
      enddo
   enddo

end subroutine ncoord_d3

!! ========================================================================
!  original D3 type coordination number from 2010
!  INPUT
!  mol  :: molecular stucture
!  OUTPUT
!  cn   :: coordination number
!  dcn  :: derivative of coordination number w.r.t. atom position
!  PARAMETERS: k1,k2
!  NOTE: k2 is already included in rcov
pure subroutine dncoord_d3(mol,cn,dcn,thr)
   use class_molecule

   implicit none

   type(molecule),intent(in) :: mol
   real(wp),intent(out) :: cn(mol%nat)
   real(wp),intent(out) :: dcn(3,mol%nat,mol%nat)
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
   dcn = 0._wp

   do i = 1, mol%nat
      do j = 1, i-1
         rij = mol%xyz(:,j) - mol%xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r = sqrt(r2)
         rcovij=(rcov(mol%at(i))+rcov(mol%at(j)))
         expterm=exp(-k1*(rcovij/r-1._wp))
         tmp = 1._wp/(1._wp+expterm)
         dtmp = (-k1*rcovij*expterm)/(r2*((expterm+1._wp)**2))
         cn(i) = cn(i) + tmp
         cn(j) = cn(j) + tmp
         dcn(:,i,i)=-dtmp*rij/r + dcn(:,i,i)
         dcn(:,j,j)= dtmp*rij/r + dcn(:,j,j)
         dcn(:,i,j)= dtmp*rij/r
         dcn(:,j,i)=-dtmp*rij/r
      enddo
   enddo

end subroutine dncoord_d3

!! ========================================================================
!  covalent coordination number of the DFT-D4 method
!  INPUT
!  mol  :: molecular stucture
!  OUTPUT
!  cn   :: coordination number
!  PARAMETERS: k1,k2,k4,k5,k6
!  NOTE: k2 is already included in rcov
pure subroutine ncoord_d4(mol,cn,thr)
   use class_molecule

   implicit none

   type(molecule),intent(in) :: mol
   real(wp),intent(out) :: cn(mol%nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j
   real(wp) :: rij(3), r, rco, den, rr, r2, xn, tmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn = 0._wp

   do i=1,mol%nat
      do j=1,i-1
         rij = mol%xyz(:,j) - mol%xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle 
         r=sqrt(r2)
!        covalent distance in bohr
         rco=rcov(mol%at(j)) + rcov(mol%at(i))
         rr=rco/r
         den=k4*exp(-(abs((en(mol%at(i))-en(mol%at(j))))+ k5)**2/k6 )
!        counting function exponential has a better long-range 
!        behavior than MHGs inverse damping
         !tmp = den/(1.d0+exp(-k1*(rr-1.0d0)))
!        error function has an even better long range behavior
         tmp = den * 0.5_wp * (1.0_wp + erf(-kn*(r-rco)/rco))
         cn(i)=cn(i)+tmp
         cn(j)=cn(j)+tmp
      enddo
   enddo

end subroutine ncoord_d4

!! ========================================================================
!  derivative of the covalent coordination number of the DFT-D4 method
!  NOTE: the derivative is inlined in the dispgrad in dftd4 by hand
!  INPUT
!  mol  :: molecular stucture
!  OUTPUT
!  cn   :: coordination number
!  dcn  :: derivative of coordination number w.r.t. atom position
!  PARAMETERS: k1,k2,k4,k5,k6
!  NOTE: k2 is already included in rcov
pure subroutine dncoord_d4(mol,cn,dcn,thr)
   use class_molecule

   implicit none

   type(molecule),intent(in) :: mol
   real(wp),intent(out) :: cn(mol%nat)
   real(wp),intent(out) :: dcn(3,mol%nat,mol%nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i, j, ia, ja
   real(wp) :: r, r2, rij(3)
   real(wp) :: rcovij,den
   real(wp) :: expterm
   real(wp) :: dtmp, tmp
   real(wp),parameter :: sqrtpi = 1.77245385091_wp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn  = 0._wp
   dcn = 0._wp

   do i = 1, mol%nat
      ia = mol%at(i)
      do j = 1, i-1
         ja = mol%at(j)
         rij = mol%xyz(:,j) - mol%xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r = sqrt(r2)
         rcovij=(rcov(ia)+rcov(ja))
         den=k4*exp(-(abs((en(ia)-en(ja)))+ k5)**2/k6 )
         !expterm=exp(-k1*(rcovij/r-1._wp))
         !tmp = den/(1._wp+expterm)
         tmp = den * 0.5_wp * (1.0_wp + erf(-kn*(r-rcovij)/rcovij))
         !dtmp = (-k1*rcovij*expterm*den)/(r2*((expterm+1._wp)**2))
         dtmp = -den*kn/sqrtpi/rcovij*exp(-kn**2*(r-rcovij)**2/rcovij**2)
         cn(i) = cn(i) + tmp
         cn(j) = cn(j) + tmp
         dcn(:,i,i)=-dtmp*rij/r + dcn(:,i,i)
         dcn(:,j,j)= dtmp*rij/r + dcn(:,j,j)
         dcn(:,i,j)=-dtmp*rij/r
         dcn(:,j,i)= dtmp*rij/r
      enddo
   enddo

end subroutine dncoord_d4

end module coordination_number
