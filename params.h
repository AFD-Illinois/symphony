#include "dec.h"

//parameters of calculation
//we use Gaussian CGS units
const double mass_electron = 9.1093826e-28;
const double speed_light = 2.99792458e10;
const double theta_e = 10.;
const double electron_charge = 4.80320680e-10;
const double B_field = 30.;
const double n_e = 1.;
const double obs_angle = (M_PI  / 3.);
const int C = 10;
const double n_max = 30.;

//power law parameters
const double power_law_p = 3.;
const double gamma_min = 1.;
const double gamma_max = 1000.;
const double n_e_NT = 1.;
//double gamma_cutoff = 1000.; //also a kappa distribution parameter

//kappa distribution parameters
const double kappa = 3.5;
const double kappa_width = 10.; //width of core of kappa dist.
const double gamma_cutoff = 1000000000000;

//struct to pass parameters to the integrand
struct parameters{
  double n;
  double nu;
};

//need 2 separate n integrations to numerically resolve STOKES_V
static int stokes_v_switch = 0.;

