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
   double s6  {-1.0};
   double s8  {-1.0};
   double s10 { 0.0};
   double a1  {-1.0};
   double a2  {-1.0};
   double s9  { 1.0};
   int alp {16};
   double beta {1.0};
};

struct DFTD_options {
   int lmbd = -1;
   int refq = -1;
   double wf = 0.0;
   double g_a = 0.0;
   double g_c = 0.0;
   bool lmolpol = false;
   bool lenergy = false;
   bool lgradient = false;
   bool lhessian = false;
   bool verbose = false;
   bool veryverbose = false;
   bool silent = false;
};

extern void D4_calculation(const int& natoms, const int* attyp,
      const double& charge, const double* coord,
      const DFTD_parameter& dparam, const DFTD_options& dopt,
      double& energy, double* grad, double* hess);

#ifdef __cplusplus
}
}
#endif
