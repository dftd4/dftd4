#include <stdlib.h>
#include <assert.h>

#include "dftd4.h"

int
main (int argc, char **argv)
{
// definitions analogous to DFT-D4 test for API
// must yield exactly the same results as Fortran code!
   double const thr = 1.0e-10;
   int    const natoms = 3;
   int    const attyp[3] = {8,1,1};
   double const charge = 0.0;
   double const coord[3*3] =
      {0.00000000000000,  0.00000000000000, -0.73578586109551,
       1.44183152868459,  0.00000000000000,  0.36789293054775,
      -1.44183152868459,  0.00000000000000,  0.36789293054775};

   DFTD_parameter dparam_b2plyp = (DFTD_parameter){
      0.6400, 1.16888646, 0.0, 0.44154604, 4.73114642, 1.0, 16, 1.0};
   DFTD_parameter dparam_tpss = (DFTD_parameter){
      1.0000, 1.76596355, 0.0, 0.42822303, 4.54257102, 1.0, 16, 1.0};
   DFTD_options opt_1 = (DFTD_options){
      p_mbd_approx_atm, p_refq_goedecker,
      6.0, 3.0, 2.0, false, true, false, true, 0 };
   DFTD_options opt_2 = (DFTD_options){
      p_mbd_approx_atm, p_refq_goedecker,
      6.0, 3.0, 2.0, false, false, true, false, 0 };

   double energy = 0.0;
   double grad[3*3] = {0.0};
   double hess[3*3*3*3] = {0.0};

   D4_calculation(&natoms, attyp, &charge, coord, "-", &dparam_tpss, &opt_1,
         &energy, grad, hess);
   assert(abs(-0.26682682254336E-03 - energy) < thr);

   assert(abs(-0.97182441530696E-05 - hess[3]             ) < thr);
   assert(abs(-0.66109938971051E-05 - hess[7*(3*natoms)+4]) < thr);
   assert(abs( 7.59401653431350E-06 - hess[3*3*3*3-1]) < thr);

   D4_calculation(&natoms, attyp, &charge, coord, "-", &dparam_b2plyp, &opt_2,
         &energy, grad, hess);
   assert(abs(-0.13368190339570E-03 - energy) < thr);

   assert(abs( 0.00000000000000E+00 - grad[0]) < thr);
   assert(abs( 0.39779053826285E-04 - grad[2]) < thr);
   assert(abs(-0.19889526913143E-04 - grad[5]) < thr);
   assert(abs(grad[3]+grad[6])                 < thr);

   return EXIT_SUCCESS;
}

