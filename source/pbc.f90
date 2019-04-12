!> generate a Wigner--Seitz cell from a given structure
!
subroutine generate_wsc(mol,wsc,rep)
   use iso_fortran_env, wp => real64
   use class_molecule
   use class_wsc
   implicit none
   !> molecular structure informtion
   type(molecule),intent(inout) :: mol
   !> Wigner--Seitz cell data type (might be contained in mol)
   type(ws_cell), intent(inout) :: wsc
   !> images of the unit cell to consider
   integer,       intent(in)    :: rep(3)
! ------------------------------------------------------------------------
!  Variables
! ------------------------------------------------------------------------
   integer  :: ii,jj,ich
   integer  :: aa,bb,cc
   integer  :: c,wc
   integer  :: minpos
   integer  :: nminpos
   integer  :: t1,t2,t3
   real(wp) :: mindist
   real(wp) :: nmindist
   real(wp),parameter :: tol = 0.01_wp !< overall WSC tolerance to consider atoms as WSC-images
   real(wp),allocatable,dimension(:,:,:) :: txyz
   real(wp),allocatable,dimension(:)     :: dist
   logical, allocatable,dimension(:)     :: trans

! ------------------------------------------------------------------------
!  allocate space for the WSC first, otherwise associate will
!  reference some dangling pointers which leads to segfaults...
   call wsc%allocate(mol%nat,rep,mol%lattice)

! ------------------------------------------------------------------------
! initialize
! ------------------------------------------------------------------------
   allocate( txyz(3,wsc%cells,mol%nat) ); txyz = 0.0_wp
   allocate( dist(wsc%cells) );     dist = 0.0_wp
   allocate( trans(wsc%cells) );    trans = .true.
! ------------------------------------------------------------------------
! Create the Wigner-Seitz Cell (WSC)
! ------------------------------------------------------------------------
   wsc%at  = 0
   wsc%itbl= 0
!$omp parallel default(none) &
!$omp private(ii,jj,wc,c,dist,trans) &
!$omp shared(mol,wsc,rep,txyz) &
!$omp shared(mindist,minpos,nmindist,nminpos)
!$omp do schedule(dynamic)
   ! Each WSC of one atom consists of n atoms
   do ii=1,mol%nat
      do jj=1,mol%nat
         if (ii.eq.jj) cycle
         ! find according neighbours
         c=1
         dist = 0.0_wp
         do aa=-rep(1),rep(1),1
            do bb=-rep(2),rep(2),1
               do cc=-rep(3),rep(3),1
                  call get_translation(aa,bb,cc,mol%xyz(:,jj),wsc%lattice,txyz(:,c,ii))
                  dist(c)=norm2(mol%xyz(:,ii)-txyz(:,c,ii))
                  c=c+1
               end do
            end do
         end do
         ! get first image with same dist
         ! find minimum in dist-array and assign it to minpos = minimum position
         trans=.true.
         minpos=minloc(dist(:),dim=1)
         ! minimum distance saved in mindist
         mindist=dist(minpos)
         trans(minpos)=.false.
         wc=1
         wsc%xyz(:,wc,jj,ii)=txyz(:,minpos,ii)
         ! get other images with same distance
         find_images : do
            nminpos=minloc(dist(:),dim=1,mask=trans(:))
            nmindist=dist(nminpos)
            if(abs(mindist-nmindist).lt.tol)then
               trans(nminpos)=.false.
               wc=wc+1
               wsc%xyz(:,wc,jj,ii)=txyz(:,nminpos,ii)
            else
               wsc%w(jj,ii)=1.0_wp/real(wc,wp)
               wsc%itbl(jj,ii)=wc
               wsc%at(jj,ii)=jj
               exit find_images
            end if
         end do find_images
      end do
   end do
!$omp enddo
!$omp endparallel

contains

   pure subroutine get_translation(t1,t2,t3,xyz,abc,out)
      implicit none
      integer,  intent(in)  :: t1,t2,t3
      real(wp), intent(in)  :: xyz(3)
      real(wp), intent(in)  :: abc(3,3)
      real(wp), intent(out) :: out(3)
      out=xyz+t1*abc(:,1)+t2*abc(:,2)+t3*abc(:,3)
   end subroutine get_translation

end subroutine generate_wsc

