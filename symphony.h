//GSL libraries
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>

//other C header files
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

//choose distribution function
#define THERMAL (0)
#define POWER_LAW (1)
#define KAPPA_DIST (2)
#define DISTRIBUTION_FUNCTION (POWER_LAW)

//choose absorptivity or emissivity
#define ABSORP (10)
#define EMISS  (11)
#define MODE   (EMISS)

//choose polarization mode
#define STOKES_I (15)
#define STOKES_Q (16)
#define STOKES_U (17)
#define STOKES_V (18)
#define EXTRAORDINARY_MODE (19)
#define ORDINARY_MODE (20)
#define POL_MODE (STOKES_I)

//function declarations
double n_peak(double nu);
double polarization_term(double gamma, double n, double nu);
double my_Bessel_J(double n, double x);
double my_Bessel_dJ(double n, double x);
double maxwell_juttner_f(double gamma);
double integrand_without_extra_factor(double gamma, double n, double nu);
double gamma_integrand(double gamma, void * params);
double gamma_integration_result(double n, void * params);
double n_summation(double nu);
double n_integration(double n_minus, double nu);
double gsl_integrate(double min, double max, double n, double nu);
double normalize_f();
double power_law_to_be_normalized(double gamma, void * params);
double power_law_f(double gamma);
double kappa_to_be_normalized(double gamma, void * params);
double kappa_f(double gamma);
double derivative(double n_start, double nu);
double differential_of_f(double gamma, double nu);
double j_nu(double nu, double B_temp, double n_e_temp, double obs_angl_temp);
double alpha_nu(double nu, double B_temp, double n_e_temp, double obs_angl_temp);
double j_nu_fit(double nu, double B_temp, double n_e_temp, double obs_angl_temp, int pol);
double alpha_nu_fit(double nu, double B_temp, double n_e_temp, double obs_angl_temp, int pol);

//fitting formulae
double thermal_I(double nu);
double thermal_Q(double nu);
double thermal_V(double nu);
double planck_func(double nu);
double thermal_I_abs(double nu);
double thermal_Q_abs(double nu);
double thermal_V_abs(double nu);
double kappa_I(double nu);
double kappa_Q(double nu);
double kappa_V(double nu);
double kappa_I_abs(double nu);
double kappa_Q_abs(double nu);
double kappa_V_abs(double nu);
double power_law_I(double nu);
double power_law_Q(double nu);
double power_law_V(double nu);
double power_law_I_abs(double nu);
double power_law_Q_abs(double nu);
double power_law_V_abs(double nu);

//variable declarations
extern const double mass_electron;
extern const double h;
extern const double speed_light;
extern double theta_e;
extern const double electron_charge;
extern double B;
extern double n_e;
extern double obs_angle;
extern int C;
extern double n_max;
//power law parameters
extern double power_law_p;
extern double gamma_min;
extern double gamma_max;
extern double n_e_NT;
//kappa distribution parameters
extern double kappa;
extern double kappa_width;
extern double gamma_cutoff;

//struct to pass parameters to the integrand
struct parameters{
  double n;
  double nu;
};
