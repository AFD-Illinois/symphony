#include <math.h>

//parameters of calculation
//we use Gaussian CGS units
#define C_PI (3.1415926535897932384)
const double mass_electron = 9.1093826e-28;
const double h = 6.6260693e-27;
const double speed_light = 2.99792458e10;
const double electron_charge = 4.80320680e-10;
double theta_e = 10.;
int C = 10;
double n_max = 30.;

//power law parameters
double power_law_p = 3.;
double gamma_min = 1.;
double gamma_max = 1000.;
double n_e_NT = 1.;
//double gamma_cutoff = 1000.; //also a kappa distribution parameter

//kappa distribution parameters
double kappa = 3.5;
double kappa_width = 10.;
double gamma_cutoff = 1000000000000;

//need 2 separate n integrations to numerically resolve STOKES_V
static int stokes_v_switch = 0.;
