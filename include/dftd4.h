/* This file is part of dftd4.
 *
 * Copyright (C) 2017-2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher
 *
 * xtb is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * xtb is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with xtb.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#ifdef __cplusplus
namespace dftd {
extern "C" {
#else
#include <stdbool.h>
#endif

// since there is no enumerator in Fortran, we define
static const int p_refq_goedecker        = 5;

static const int p_mbd_none       = 0; // just pair-wise dispersion
static const int p_mbd_rpalike    = 1; // RPA-like (=MBD) non-additivity
static const int p_mbd_approx_atm = 3; // approximate C9 from C6

typedef struct {
   double s6;
   double s8;
   double s10;
   double a1;
   double a2;
   double s9;
   int alp;
   double beta;
} DFTD_parameter;

typedef struct {
   int lmbd;
   int refq;
   double wf;
   double g_a;
   double g_c;
   bool lmolpol;
   bool lenergy;
   bool lgradient;
   bool lhessian;
   bool print_level;
} DFTD_options;

extern int
D4_calculation(const int* natoms, const int* attyp, const double* charge,
      const double* coord, const char* outfile,
      const DFTD_parameter* dparam, const DFTD_options* dopt,
      double* energy, double* grad, double* hess,
      double* polarizibilities, double* c6_coefficients, double* charges);

extern int
D4_PBC_calculation(const int* natoms, const int* attyp, const double* charge,
      const double* coord, const double* lattice, const bool* pbc,
      const char* outfile, const DFTD_parameter* dparam, const DFTD_options* dopt,
      double* energy, double* grad, double* latgrad, double* stress, double* hess,
      double* polarizibilities, double* c6_coefficients, double* charges);

extern int
D4_damping_parameters(
      const char* name, const DFTD_parameter* dparam, const int* lmbd);

extern int
D3_calculation(const int* natoms, const int* attyp,
      const double* coord, const char* outfile,
      const DFTD_parameter* dparam, const DFTD_options* dopt,
      double* energy, double* grad, double* hess,
      double* polarizibilities, double* c6_coefficients);

extern int
D3_PBC_calculation(const int* natoms, const int* attyp,
      const double* coord, const double* lattice, const bool* pbc,
      const char* outfile, const DFTD_parameter* dparam, const DFTD_options* dopt,
      double* energy, double* grad, double* latgrad, double* stress, double* hess,
      double* polarizibilities, double* c6_coefficients);

#ifdef __cplusplus
}
}
#endif
