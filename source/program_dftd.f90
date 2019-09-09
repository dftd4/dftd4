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

program dftd
   use iso_fortran_env, istdout => output_unit, kdp => real64
!$ use omp_lib
! ------------------------------------------------------------------------
!  general purpose library
! ------------------------------------------------------------------------
   use mctc_global
   use mctc_timings
   use mctc_econv

! ------------------------------------------------------------------------
!  class definitions
! ------------------------------------------------------------------------
   use class_molecule
   use class_set
   use class_results

! ------------------------------------------------------------------------
!  interfaces
! ------------------------------------------------------------------------
   use geometry_reader
   use coordination_number
   use eeq_model
   use dftd4
   use dfuncpar
   use dispersion_calculator

   implicit none

! ------------------------------------------------------------------------
!  class declarations
! ------------------------------------------------------------------------
   type(molecule)       :: mol
   type(options)        :: set
   type(dftd_options)   :: dopt
   type(dftd_parameter) :: dparam
   type(chrg_parameter) :: chrgeq
   type(dftd_results)   :: dresults

! ------------------------------------------------------------------------
!  local variables
! ------------------------------------------------------------------------
   integer  :: ndim                      ! matrix dimension
   integer  :: i,j,k,l,ii,jj
   integer  :: err
   real(wp) :: memory
   real(wp) :: etmp,etwo,emany,er,el,es
   real(wp) :: molpol,molc6,molc8        ! molecular Polarizibility
   real(wp) :: lattice_grad(3,3)
   real(wp),allocatable :: q(:)          ! partial charges
   real(wp),allocatable :: dqdr(:,:,:)   ! partial charges
   real(wp),allocatable :: covcn(:)      ! covalent coordination number
   real(wp),allocatable :: dcovcndr(:,:,:) ! covalent coordination number derivative
   real(wp),allocatable :: cn(:)         ! erf coordination number
   real(wp),allocatable :: dcndr(:,:,:)  ! erf coordination number derivative
   real(wp),allocatable :: gweights(:)   ! gaussian weights
   real(wp),allocatable :: refc6(:,:)    ! reference C6 coeffients
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: aw(:,:)
   real(wp),allocatable :: ges(:,:)
   real(wp),allocatable :: gr(:,:)
   real(wp),allocatable :: gl(:,:)
   real(wp),parameter   :: step = 1.0e-5_wp, step2 = 0.5_wp/step

   integer, parameter   :: rep_wsc(3) = [1,1,1]

!! ------------------------------------------------------------------------
!!  signal processing
!! ------------------------------------------------------------------------
!   external :: wSIGTERM
!   external :: wSIGINT
!!  here two important signal handlers are installed, it seems that
!!  FORTRAN by itself cannot handle signals in the way I expected it
!!  to do, but this will force it to die when I hit CTRL^C.
!   call signal(2,wSIGINT)
!   call signal(15,wSIGTERM)

! ------------------------------------------------------------------------
!  general setup
! ------------------------------------------------------------------------
!  initialize the timing system
   call start_timing_run
   call init_timing(10,verb=.true.) ! verbosity allows printing of cputime
   call start_timing(1)

!  initialize the messagebuffer for the error handler
   call init_errorbuffer

!  set this for mctc_global
   name = 'dftd4'

! ------------------------------------------------------------------------
!  command line arguments
! ------------------------------------------------------------------------
   call read_commandline_arguments(set)

! ------------------------------------------------------------------------
!  get molecular geometry
! ------------------------------------------------------------------------
   if (set%lperiodic) then
      call read_geometry(set%fname,mol)
   else
      call get_geometry(mol,set%fname)
   endif
   call generate_wsc(mol,mol%wsc)
   if (set%inchrg) mol%chrg = set%chrg

! ------------------------------------------------------------------------
!  Output:
!  Header, Citation and Licence
! ------------------------------------------------------------------------
   call dftd4_header(set%verbose)
   if (.not.set%silent) then
   call dftd4_citation
   call prdate('S')
   write(istdout,'(a)')
   endif

! ------------------------------------------------------------------------
!  Set defaults
! ------------------------------------------------------------------------
   if (set%inparam) then
      dparam = set%dparam
      if (.not.(set%lgradient.or.set%lhessian)) set%lenergy = .true.
   else if (allocated(set%func)) then
      call d4par(set%func,dparam,set%lmbd)
      if (.not.(set%lgradient.or.set%lhessian)) set%lenergy = .true.
   else
      if (set%lenergy) &
         call raise('E','Dispersion energy requested but no parameters given')
      if (set%lgradient) &
         call raise('E','Dispersion gradient requested but no parameters given')
      if (set%lhessian) &
         call raise('E','Dispersion Hessian requested but no parameters given')
   endif
   if (.not.(set%lenergy.or.set%lgradient.or.set%lhessian)) set%lmolpol = .true.
   if (set%lhessian .and. set%lperiodic) then
      call raise('W',"Cannot calculate Hessian under periodic boundary conditions")
      set%lhessian = .false.
   endif
   if (set%lorca .and. set%lperiodic) then
      call raise('W',"To my knowledge there are no PBC in ORCA!")
      set%lorca = .false.
   endif

   dopt = set%export()

   call d4_calculation(istdout,dopt,mol,dparam,dresults)

! ------------------------------------------------------------------------
!  Output:
! ------------------------------------------------------------------------
if (set%lenergy.or.set%lgradient.or.set%lhessian) &
   call generic_header(istdout,'Results',49,10)
   if (set%lenergy) then
      write(istdout,'('// &
      &      '1x,"Edisp  /kcal,au:",f11.4,1x,f12.8)') &
      &       dresults%energy*autokcal,dresults%energy
      write(istdout,'(a)')
   endif
   if (set%ltmer.or.(.not.set%lorca.and.set%lenergy)) &
      call out_tmer('.EDISP',dresults%energy)

   if (set%lgradient) then
      if (set%verbose) then
         write(istdout,'(1x,a)') &
            "Dispersion gradient"
         call write_gradient(mol,istdout,dresults%gradient)
         write(istdout,'(a)')
      endif
      write(istdout,'(1x,"|G| =",1x,f34.10)') &
      &       norm2(dresults%gradient)
      write(istdout,'(a)')

      if (set%lorca) then
         call orca_gradient(mol,istdout,dresults%gradient)
      else
         call out_gradient(mol,'gradient',dresults%energy,dresults%gradient)
         if (mol%npbc > 0) then
            call out_gradlatt(mol,'gradlatt',dresults%energy,dresults%lattice_gradient)
         endif
      endif
   endif

   if (set%lhessian) then
      call orca_hessian(mol,istdout,dresults%hessian)
   endif

   call raise('F','Some non-fatal runtime exceptions occurred, please check:')

! ------------------------------------------------------------------------
!  Print timings
! ------------------------------------------------------------------------
   if (.not.set%silent) then
   write(istdout,'(a)')
   call stop_timing_run
   call stop_timing(1)
   call prdate('E')
   call prtiming(1,'total')
   endif

   write(istdout,'(a)')
   call terminate(0)

end program dftd
