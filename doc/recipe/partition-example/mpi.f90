program mpi_partition
   use, intrinsic :: iso_fortran_env, only : r8 => real64
   use dftd4, only : d4_model, damping_param, get_dispersion, &
      & get_rational_damping, new_d4_model, new_work_partition, &
      & realspace_cutoff, work_partition
   use mctc_env, only : error_type
   use mctc_io, only : structure_type, new
   use mpi_f08, only : MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, &
      & MPI_SUM, MPI_Abort, MPI_Allreduce, MPI_Comm_rank, MPI_Comm_size, &
      & MPI_Finalize, MPI_Init
   implicit none

   type(structure_type) :: mol
   type(error_type), allocatable :: error
   integer, allocatable :: num(:)
   real(r8), allocatable :: xyz(:, :), gradient(:, :)
   real(r8) :: energy
   type(d4_model) :: disp
   class(damping_param), allocatable :: param
   type(work_partition) :: partition
   integer :: rank, nranks

   call MPI_Init()
   call MPI_Comm_rank(MPI_COMM_WORLD, rank)
   call MPI_Comm_size(MPI_COMM_WORLD, nranks)

   num = [6, 1, 1, 1, 1]
   xyz = reshape([ &  ! coordinates in Bohr
     &  0.0000000_r8, -0.0000000_r8,  0.0000000_r8, &
     & -1.1922080_r8,  1.1922080_r8,  1.1922080_r8, &
     &  1.1922080_r8, -1.1922080_r8,  1.1922080_r8, &
     & -1.1922080_r8, -1.1922080_r8, -1.1922080_r8, &
     &  1.1922080_r8,  1.1922080_r8, -1.1922080_r8], &
     & [3, size(num)])
   call new(mol, num, xyz, charge=0.0_r8, uhf=0)

   call get_rational_damping("pbe0", param, s9=1.0_r8)
   if (.not.allocated(param)) call fatal("No parameters for PBE0 available")
   call new_d4_model(error, disp, mol)
   if (allocated(error)) call fatal(error%message)

   call new_work_partition(error, partition, rank, nranks)
   if (allocated(error)) call fatal(error%message)

   allocate(gradient(3, mol%nat))
   call get_dispersion(mol, disp, param, realspace_cutoff(), energy, gradient, &
      & partition=partition)

   call MPI_Allreduce(MPI_IN_PLACE, energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
      & MPI_COMM_WORLD)
   call MPI_Allreduce(MPI_IN_PLACE, gradient, size(gradient), &
      & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD)

   if (rank == 0) then
      print "(a, f13.10, a)", "PBE0-D4 dispersion energy: ", energy, " Hartree"
   end if

   call MPI_Finalize()

contains

   subroutine fatal(message)
      character(len=*), intent(in) :: message
      print "(2a)", "Error: ", message
      call MPI_Abort(MPI_COMM_WORLD, 1)
   end subroutine fatal

end program mpi_partition
