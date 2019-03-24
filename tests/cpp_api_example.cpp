#include <iostream>
#include <cassert>

#include "dftd_api.h"

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

   const dftd::DFTD_parameter dparam_b2plyp {
      s6 : 0.6400, s8 : 1.16888646, a1 : 0.44154604, a2 : 4.73114642};
   const dftd::DFTD_parameter dparam_tpss {
      s6 : 1.0000, s8 : 1.76596355, a1 : 0.42822303, a2 : 4.54257102};
   const dftd::DFTD_options opt_1 {
      lmbd : dftd::p_mbd_approx_atm, refq : dftd::p_refq_goedecker,
      wf : 6.0, g_a : 3.0, g_c : 2.0,
      lmolpol : false, lenergy : true, lgradient : false, lhessian : true,
      verbose : false, veryverbose : false, silent : false };
   const dftd::DFTD_options opt_2 = {
      lmbd : dftd::p_mbd_approx_atm, refq : dftd::p_refq_goedecker,
      wf : 6.0, g_a : 3.0, g_c : 2.0,
      lmolpol : false, lenergy : false, lgradient : true, lhessian : false,
      verbose : false, veryverbose : false, silent : true };

   double energy {0.0};
   double grad[3*natoms] {0.0};
   double hess[3*natoms*3*natoms] {0.0};

   dftd::D4_calculation(natoms, attyp, charge, coord, dparam_tpss, opt_1,
         energy, grad, hess);
   assert(abs(-0.26682682254336E-03 - energy) < thr);

   assert(abs( 7.9334628241320 - hess[2]             ) < thr);
   assert(abs(-3.2756224894310 - hess[7*(3*natoms)+3]) < thr);
   assert(abs( 0.0000000000000 - hess[2*(3*natoms)+4]) < thr);

   D4_calculation(natoms, attyp, charge, coord, dparam_b2plyp, opt_2,
         energy, grad, hess);
   assert(abs(-0.13368190339570E-03 - energy) < thr);

   assert(abs( 0.00000000000000E+00 - grad[0]) < thr);
   assert(abs( 0.39778648945254E-04 - grad[2]) < thr);
   assert(abs(-0.19889324472627E-04 - grad[5]) < thr);
   assert(abs(grad[3]+grad[6])                 < thr);

   return EXIT_SUCCESS;
}
