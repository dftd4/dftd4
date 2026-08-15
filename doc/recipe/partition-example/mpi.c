#include <stdbool.h>
#include <stdio.h>

#include <mpi.h>

#include "dftd4.h"

int main(int argc, char **argv)
{
  dftd4_error error = dftd4_new_error();
  dftd4_structure mol = NULL;
  dftd4_model d4 = NULL;
  dftd4_param param = NULL;
  int rank, nranks;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  int nat = 5;
  int num[5] = {6, 1, 1, 1, 1};
  double xyz[15] = {
     0.00000000, -0.00000000,  0.00000000,
    -1.19220800,  1.19220800,  1.19220800,
     1.19220800, -1.19220800,  1.19220800,
    -1.19220800, -1.19220800, -1.19220800,
     1.19220800,  1.19220800, -1.19220800};

  mol = dftd4_new_structure(error, nat, num, xyz, NULL, NULL, NULL);
  if (dftd4_check_error(error)) goto handle_error;

  param = dftd4_load_rational_damping(error, "pbe0", true);
  if (dftd4_check_error(error)) goto handle_error;

  d4 = dftd4_new_d4_model(error, mol);
  if (dftd4_check_error(error)) goto handle_error;

  dftd4_set_model_work_partition(error, d4, rank, nranks);
  if (dftd4_check_error(error)) goto handle_error;

  double energy;
  double gradient[15];
  dftd4_get_dispersion(error, mol, d4, param, &energy, gradient, NULL);
  if (dftd4_check_error(error)) goto handle_error;

  MPI_Allreduce(MPI_IN_PLACE, &energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, gradient, 3 * nat, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  if (rank == 0) {
    printf("PBE0-D4 dispersion energy: %13.10lf Hartree\n", energy);
  }

  dftd4_delete(error);
  dftd4_delete(mol);
  dftd4_delete(d4);
  dftd4_delete(param);
  MPI_Finalize();
  return 0;

handle_error: {
  char msg[512];
  dftd4_get_error(error, msg, NULL);
  printf("Error: %s\n", msg);

  dftd4_delete(error);
  dftd4_delete(mol);
  dftd4_delete(d4);
  dftd4_delete(param);
  MPI_Abort(MPI_COMM_WORLD, 1);
  return 1;
}
}
