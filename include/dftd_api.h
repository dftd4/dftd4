#pragma once

#ifdef __cplusplus
namespace dftd {
extern "C" {
#endif

// since there is no enumerator in Fortran, we define
static const int p_refq_goedecker        = 5;

static const int p_mbd_none       = 0; // just pair-wise dispersion
static const int p_mbd_rpalike    = 1; // RPA-like (=MBD) non-additivity
static const int p_mbd_approx_atm = 3; // approximate C9 from C6

struct DFTD_parameter {
   double s6;
   double s8;
   double s10;
   double a1;
   double a2;
   double s9;
   int alp;
   double beta;
};

struct DFTD_options {
   int lmbd;
   int refq;
   double wf;
   double g_a;
   double g_c;
   bool lmolpol;
   bool lenergy;
   bool lgradient;
   bool lhessian;
   bool verbose;
   bool veryverbose;
   bool silent;
};

extern void D4_calculation(const int& natoms, const int* attyp,
      const double& charge, const double* coord,
      const DFTD_parameter& dparam, const DFTD_options& dopt,
      double& energy, double* grad, double* hess);

#ifdef __cplusplus
}
}
#endif
