/* This file is part of dftd4.
 * SPDX-Identifier: LGPL-3.0-or-later
 *
 * dftd4 is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * dftd4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with dftd4.  If not, see <https://www.gnu.org/licenses/>.
**/

#include <stdlib.h>
#include <assert.h>

#include "dftd4.h"

int
main (void)
{
   int const natoms = 7;
   int const attyp[7] = {6,6,6,1,1,1,1};
   double const coord[21] =
      {0.00000000000000, 0.00000000000000,-1.79755622305860,
       0.00000000000000, 0.00000000000000, 0.95338756106749,
       0.00000000000000, 0.00000000000000, 3.22281255790261,
      -0.96412815539807,-1.66991895015711,-2.53624948351102,
      -0.96412815539807, 1.66991895015711,-2.53624948351102,
       1.92825631079613, 0.00000000000000,-2.53624948351102,
       0.00000000000000, 0.00000000000000, 5.23010455462158};
   double energy;
   double pair_disp2[49];
   double pair_disp3[49];
   double gradient[21];
   double sigma[9];
   double c6[49];

   assert(dftd4_get_version() > 0);

   dftd4_error error;
   dftd4_structure mol;
   dftd4_model disp;
   dftd4_param param;

   error = dftd4_new_error();
   assert(!!error);

   mol = dftd4_new_structure(error, natoms, attyp, coord, NULL, NULL, NULL);
   if (dftd4_check_error(error)) {return 1;};
   assert(!!mol);

   disp = dftd4_new_d4_model(error, mol);
   if (dftd4_check_error(error)) {return 1;}
   assert(!!disp);

   // C6 coefficients
   dftd4_get_properties(error, mol, disp, NULL, NULL, c6, NULL);
   if (dftd4_check_error(error)) {return 1;}

   // PBE-D4
   param = dftd4_new_rational_damping(error, 1.0, 0.95948085, 0.0, 0.38574991, 4.80688534, 16.0);
   if (dftd4_check_error(error)) {return 1;}
   assert(!!param);
   dftd4_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd4_check_error(error)) {return 1;}
   dftd4_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd4_check_error(error)) {return 1;}
   dftd4_get_pairwise_dispersion(error, mol, disp, param, pair_disp2, pair_disp3);
   if (dftd4_check_error(error)) {return 1;}
   dftd4_delete(param);

   // DSD-BLYP-D4-ATM
   param = dftd4_load_rational_damping(error, "dsdblyp", true);
   if (dftd4_check_error(error)) {return 1;}
   assert(!!param);
   dftd4_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd4_check_error(error)) {return 1;}
   dftd4_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd4_check_error(error)) {return 1;}
   dftd4_delete(param);

   dftd4_delete(disp);
   dftd4_delete(mol);
   dftd4_delete(error);

   assert(!param);
   assert(!disp);
   assert(!mol);
   assert(!error);

   return 0;
}
