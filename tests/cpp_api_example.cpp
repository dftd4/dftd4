#include <iostream>
#include <cmath>
#include <cassert>

#include "dftd4.h"

using namespace std;

int
main (int argc, char **argv)
{
// definitions analogous to DFT-D4 test for API
// must yield exactly the same results as Fortran code!
   const double thr = 1.0e-10;
   const int    natoms = 3;
   const int    attyp[natoms] = {8,1,1};
   const double charge = 0.0;
   const double coord[3*natoms] =
      {0.00000000000000,  0.00000000000000, -0.73578586109551,
       1.44183152868459,  0.00000000000000,  0.36789293054775,
      -1.44183152868459,  0.00000000000000,  0.36789293054775};

   dftd::DFTD_parameter dparam_b2plyp = {
      0.6400, 1.16888646, 0.0, 0.44154604, 4.73114642, 1.0, 16, 1.0};
   dftd::DFTD_parameter dparam_tpss = {
      1.0000, 1.76596355, 0.0, 0.42822303, 4.54257102, 1.0, 16, 1.0};
   dftd::DFTD_options opt_1 = {
      dftd::p_mbd_approx_atm, dftd::p_refq_goedecker,
      6.0, 3.0, 2.0, false, true, false, true, 0 };
   dftd::DFTD_options opt_2 = {
      dftd::p_mbd_approx_atm, dftd::p_refq_goedecker,
      6.0, 3.0, 2.0, false, false, true, false, 0 };

   double energy = 0.0;
   double grad[3*natoms] = {0.0};
   double hess[3*natoms*3*natoms] = {0.0};

   dftd::D4_calculation(&natoms, attyp, &charge, coord, "-", &dparam_tpss, &opt_1,
         &energy, nullptr, hess, nullptr, nullptr, nullptr);
   assert(fabs(-0.26678857789318E-03 - energy) < thr);

   assert(fabs(-0.97182441530696E-05 - hess[3]             ) < thr);
   assert(fabs(-0.66109938971051E-05 - hess[7*(3*natoms)+4]) < thr);
   assert(fabs( 7.59401653431350E-06 - hess[3*3*3*3-1]) < thr);

   dftd::D4_calculation(&natoms, attyp, &charge, coord, "-", &dparam_b2plyp, &opt_2,
         &energy, grad, nullptr, nullptr, nullptr, nullptr);
   assert(fabs(-0.13366273625493E-03 - energy) < thr);

   assert(fabs( 0.00000000000000E+00 - grad[0]) < thr);
   assert(fabs( 0.39779053826285E-04 - grad[2]) < thr);
   assert(fabs(-0.19889526913143E-04 - grad[5]) < thr);
   assert(fabs(grad[3]+grad[6])                 < thr);

   return EXIT_SUCCESS;
}

