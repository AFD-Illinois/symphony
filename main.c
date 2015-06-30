/* Currently Unnamed Synchrotron Calculator
 * by Alex Pandya, Zhaowei Zhang
 * 6/30/15
 *
 *References:
 * 1) Leung, Gammie, and Noble (2011)
 * 2) Xiao (2006)
*/

//GSL libraries
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>

//other C header files
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

//parameters of calculation
double mass_electron = 9.1093826e-28;
double speed_light = 2.99792458e10;
double theta_e = 10.;
double electron_charge = 4.80320680e-10;
double B_field = 30.;
double n_e = 1.;
double observer_angle = (M_PI  / 3.);
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
double w = 10.; //width of core of kappa dist.
double gamma_cutoff = 1000000000000;

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

//struct to pass parameters to integrand
struct parameters{
  double n;
  double nu;
};

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

/* main: defines nu_c (cyclotron frequency) and
 * loops through values of nu, to give output
 * absorptivity or emissivity vs. nu/nu_c
 */
int main(int argc, char *argv[]) {
  double nu_c = (electron_charge * B_field)
 	       /(2. * M_PI * mass_electron * speed_light);
  int index = 0;
  for (index; index < 31; index++) {
    double nu = pow(10., index/5.) * nu_c;
    printf("\n%e	%e", nu/nu_c, n_summation(nu));
  }
  printf("\n");
  return 0;
}

/*n_peak: gives the location of the peak of the n integrand for 
 *the THERMAL distribution; uses Eq. 68 in [1]
 *
 *@param nu: Input, frequency of emission/absorption
 *@returns n_peak: Output, location of integrand's peak for the 
 * n-integral for the THERMAL distribution 
 */
double n_peak(double nu) {
  double nu_c = (electron_charge * B_field)
	      / (2. * M_PI * mass_electron * speed_light);
  double beta = 0.;
  if (nu <= nu_c * theta_e*theta_e || theta_e < 1.) {
    beta = sqrt((1. - 1./pow((1. + theta_e),2.))); 
    }
  else {
    beta = sqrt(1. - pow((2. * theta_e * nu / nu_c), -2./3.)); //beta68
  }
  double n_peak =  (theta_e + 1. + pow((2. * theta_e * nu / nu_c),1./3.))
                 * (nu/nu_c) * (1. - beta*beta * pow(cos(observer_angle),2.));
  return n_peak;
}

/*polarization_term: term in the gamma integrand that varies based upon the 
 *polarization mode; uses eq. 13-19, 22-25, 27-30, in [1]
 *
 *@param gamma: Input, Lorentz factor 
 *@param n: Input, index n of sum (synchrotron harmonic number)
 *@param nu: Input, frequency of emission/absorption
 *@returns: piece of gamma integrand that determines polarization mode 
 */
double polarization_term(double gamma, double n, double nu) {
  double nu_c = (electron_charge * B_field)
	       /(2. * M_PI * mass_electron * speed_light);
  double beta = sqrt(1. - 1./(gamma*gamma));
  double cos_xi = (gamma * nu - n * nu_c)
		 /(gamma * nu * beta * cos(observer_angle));
  double M = (cos(observer_angle) - beta * cos_xi)/sin(observer_angle);
  double N = beta * sqrt(1 - (cos_xi*cos_xi));
  double z = (nu * gamma * beta * sin(observer_angle) 
       * sqrt(1. - cos_xi*cos_xi))/nu_c;
  double K_xx = M*M * pow(my_Bessel_J(n, z), 2.);
  double K_yy = N*N * pow(my_Bessel_dJ(n, z), 2.);
  double T_x = (2.*nu*cos(observer_angle))/(nu_c*pow(sin(observer_angle), 2.)
          +sqrt(nu_c*nu_c*pow(sin(observer_angle), 4.)
          +4.*nu*nu*pow(cos(observer_angle), 2.)));
  #if POL_MODE == STOKES_I
    double ans = K_xx + K_yy;
  #elif POL_MODE == STOKES_Q
    double ans = K_xx - K_yy;
  #elif POL_MODE == STOKES_U
  double ans = 0.;
  #elif POL_MODE == STOKES_V
    double ans = -2.*M*N*my_Bessel_J(n, z)*my_Bessel_dJ(n, z);
  #elif POL_MODE == EXTRAORDINARY_MODE //factor of 2 comes from eq. 43 vs 42
    double ans = pow(M*T_x*my_Bessel_J(n, z)+N*my_Bessel_dJ(n, z), 2.)
	           /(1. + T_x*T_x);
  #elif POL_MODE == ORDINARY_MODE //factor of 2 again from eq. 43 vs 42
    double ans = pow(M*T_x*my_Bessel_J(n, z)-N*my_Bessel_dJ(n, z), 2.)
	           /(1. + T_x*T_x);
  #endif

  //we need to account for factor of 2 difference between eq. 42 and 43 of [1]
  if(POL_MODE == EXTRAORDINARY_MODE || ORDINARY_MODE && MODE == ABSORP) {
    ans = 2. * ans;
  } 

  return ans;
}

/*differential_of_f: term in gamma integrand only for absorptivity calculation; 
 *it is the differential Df = 2\pi\nu (1/(mc)*d/dgamma + (beta cos(theta)
 * -cos(xi))/(p*beta*c) * d/d(cos(xi))) f 
 *this is eq. 41 of [1]
 *below it is applied for the THERMAL, POWER_LAW, and KAPPA_DIST distributions
 *@param gamma: Input, Lorentz factor
 *@param nu: Input, frequency of emission/absorption
 *@returns: Output, Df term in gamma integrand; depends on distribution function
 */
double differential_of_f(double gamma, double nu) {
  #if DISTRIBUTION_FUNCTION == THERMAL
    double prefactor = (M_PI * nu / (mass_electron*speed_light*speed_light)) 
		     * (n_e/(theta_e * gsl_sf_bessel_Kn(2, 1./theta_e)));
    double body = (-1./theta_e) * exp(-gamma/theta_e);
    double f = prefactor * body;
  #elif DISTRIBUTION_FUNCTION == POWER_LAW
    double pl_norm = 4.* M_PI/(normalize_f());
    double prefactor = (M_PI * nu / (mass_electron*speed_light*speed_light)) 
                     * (n_e_NT*(power_law_p-1.))
                     /((pow(gamma_min, 1.-power_law_p) 
                     - pow(gamma_max, 1.-power_law_p)));
    double term1 = ((-power_law_p-1.)*exp(-gamma/gamma_cutoff)
               *pow(gamma,-power_law_p-2.)/(sqrt(gamma*gamma - 1.)));
    double term2 = (exp(-gamma/gamma_cutoff) * pow(gamma,(-power_law_p-1.))
                  /(gamma_cutoff * sqrt(gamma*gamma - 1.)));
    double term3 = (exp(-gamma/gamma_cutoff) * pow(gamma,-power_law_p))
               /pow((gamma*gamma - 1.), (3./2.));
    double f = pl_norm * prefactor * (term1 - term2 - term3);
  #elif DISTRIBUTION_FUNCTION == KAPPA_DIST
    double prefactor = n_e * (1./normalize_f()) * 4. * M_PI*M_PI * nu 
                     * mass_electron*mass_electron * speed_light;
    double term1 = ((- kappa - 1.) / (kappa * w)) * pow((1. + (gamma - 1.)
                  /(kappa * w)), -kappa-2.);
    double term2 = pow((1. + (gamma - 1.)/(kappa * w)), (- kappa - 1.)) 
                    * (- 1./gamma_cutoff);
    double f = prefactor * (term1 + term2) * exp(-gamma/gamma_cutoff);
  #endif
    return f;
}

/*maxwell_juttner_f: Maxwell-Juttner distribution function in terms of Lorentz 
 * factor gamma; uses eq. 47, 49, 50 of [1]
 *@param gamma: Input, Lorentz factor
 *@returns: THERMAL distribution function, which goes into the gamma integrand 
 */
double maxwell_juttner_f(double gamma) {
  double beta = sqrt(1. - 1./(gamma*gamma));
  double d = (n_e * gamma * sqrt(gamma*gamma-1.) * exp(-gamma/theta_e))
            /(4. * M_PI * theta_e * gsl_sf_bessel_Kn(2, 1./theta_e));
  double ans = 1./(pow(mass_electron, 3.) * pow(speed_light, 3.) * gamma*gamma 
               * beta) * d;
  return ans;
}

/*power_law_to_be_normalized: the power-law distribution is normalized as-is, 
 * but we have added an exponential cutoff e^(-gamma/gamma_cutoff), so it must 
 * be normalized again.  The normalization constant is given by 1 over
 * the integral of power_law_to_be_normalized from 1 to infinity.
 * uses eq. 50 of [1]
 *@param gamma: Input, Lorentz factor
 *@param *params: void pointer to struct parameters, which contains
 * n (harmonic number) and nu (frequency of emission/absorption)
 *@returns: Output, power-law distribution to be normalized via normalize_f()
 */
double power_law_to_be_normalized(double gamma, void * params) {
  double norm_term = 4. * M_PI;
  double prefactor = (power_law_p - 1.) / (pow(gamma_min, 1. - power_law_p) 
                    - pow(gamma_max, 1. - power_law_p));
  //double body = pow(gamma, -p) * exp(- gamma / gamma_cutoff);
  double body = pow(gamma, -power_law_p);
  double ans = norm_term * prefactor * body;
  return ans;
}

/*power_law_f: power-law distribution function, normalized via call to the 
 * normalize_f() function. Uses eq. 50 of [1].
 *@param gamma: Input, Lorentz factor
 *@returns: Ouput, a normalized power-law distribution for the gamma integrand
 */
double power_law_f(double gamma) {
  double beta = sqrt(1. - 1./(gamma*gamma));
  double prefactor = n_e_NT * (power_law_p - 1.) / (pow(gamma_min, 1. 
                     - power_law_p) - pow(gamma_max, 1. - power_law_p));
  //double body = pow(gamma, -p) * exp(- gamma / gamma_cutoff);
  double body = pow(gamma, -power_law_p);
  double ans = 1./normalize_f() * prefactor * body * 1./(pow(mass_electron, 3.) 
               * pow(speed_light, 3.) * gamma*gamma * beta);
  return ans;
}

/*kappa_to_be_normalized: the kappa distribution is not normalized, so we must 
 * find the normalization constant by taking 1 over the integral of the function
 * kappa_to_be_normalized from 1 to infinity.  Uses eq. 42 of [2].
 *@param gamma: Input, Lorentz factor
 *@param *params: void pointer to struct parameters, which contains
 * n (harmonic number) and nu (frequency of emission/absorption)
 *@returns: kappa distribution function to be normalized via normalize_f()
 */
double kappa_to_be_normalized(double gamma, void * params)
{
  double kappa_body = pow((1. + (gamma - 1.)/(kappa * w)), -kappa-1);
  double cutoff = exp(-gamma/gamma_cutoff);
  double norm_term = 4. * M_PI * pow(mass_electron, 3.) * pow(speed_light, 3.) 
                   * gamma * sqrt(gamma*gamma-1.);
  //double ans = kappa_body * cutoff * norm_term;
  double ans = kappa_body * norm_term;
  return ans;
}

/*kappa_f: kappa distribution function, numerically normalized
 *uses eq. 42 of [2]
 *@param gamma: Input, Lorentz factor
 *@returns: normalized kappa distribution function to go into gamma integrand 
 */
double kappa_f(double gamma)
{
  double norm = 1./normalize_f();
  double kappa_body = n_e * pow((1. + (gamma - 1.)/(kappa * w)), -kappa-1);
  double cutoff = exp(-gamma/gamma_cutoff);
  //double ans = norm * kappa_body * cutoff;
  double ans = norm * kappa_body;
  return ans;
}

/*integrand_without_extra_factor: the function I(n, xi, gamma) from eq.
 * 60 of [1].
 *@param gamma: Input, Lorentz factor
 *@param n: Input, harmonic number; index of sum in gamma integrand
 *@param nu: Input, frequency of emission/absorption
 *@returns: Ouput, gamma integrand without terms pulled out of integral in
 * eq. 60, 62 of [1]
 */
double integrand_without_extra_factor(double gamma, double n, double nu)
{
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);
  double beta = sqrt(1. - 1./(gamma*gamma));
  double cos_xi = (gamma * nu - n * nu_c)
                 /(gamma * nu * beta * cos(observer_angle));
  #if DISTRIBUTION_FUNCTION == THERMAL
    double ans = (2. * M_PI * electron_charge*electron_charge * nu*nu)
                / speed_light * (pow(mass_electron, 3.) * pow(speed_light, 3.) 
                * gamma*gamma * beta * 2. * M_PI) * maxwell_juttner_f(gamma) 
                * polarization_term(gamma, n, nu);
  #elif DISTRIBUTION_FUNCTION == POWER_LAW
    double ans = (2. * M_PI * electron_charge*electron_charge * nu*nu)
                / speed_light * (pow(mass_electron, 3.) * pow(speed_light, 3.) 
                * gamma*gamma * beta * 2. * M_PI) * power_law_f(gamma) 
                * polarization_term(gamma, n, nu);
  #elif DISTRIBUTION_FUNCTION == KAPPA_DIST
    double ans = (2. * M_PI * electron_charge*electron_charge * nu*nu)
                / speed_light * (pow(mass_electron, 3.) * pow(speed_light, 3.) 
                * gamma*gamma * beta * 2. * M_PI) * kappa_f(gamma) 
                * polarization_term(gamma, n, nu);
  #endif
    return ans;
}

/*gamma_integrand: full gamma integrand, to be integrated from gamma_minus
 * to gamma_plus, set by eq. 64 of [1]
 *@param gamma: Lorentz factor
 *@param *params: void pointer to struct parameters, which contains
 * n (harmonic number) and nu (frequency of emission/absorption)
 *@returns: Output, integrand to be integrated by gamma_integration_result()
 */
double gamma_integrand(double gamma, void * params)
{
  struct parameters n_and_nu = *(struct parameters*) params;
  double n = n_and_nu.n;
  double nu = n_and_nu.nu;
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);
  double beta = sqrt(1. - 1./(gamma*gamma));
  #if MODE == EMISS
    double cos_xi = (gamma * nu - n * nu_c)
                   /(gamma * nu * beta * cos(observer_angle));
    double prefactor = 1./(nu * beta * fabs(cos(observer_angle)));
    double ans = prefactor * integrand_without_extra_factor(gamma, n, nu);
  #elif MODE == ABSORP
    double prefactor = -speed_light*electron_charge*electron_charge / (2. * nu);
    double ans = prefactor*gamma*gamma*beta*differential_of_f(gamma, nu)
                *polarization_term(gamma, n, nu)
                *(1./(nu*beta*fabs(cos(observer_angle))));
  #endif
    return ans;
}

/*gamma_integration_result: performs gamma integral from gamma_minus to 
 * gamma_plus; calls gsl_integrate(), which is a wrapper for the GSL integrator
 * GSL QAG.
 *@param n: Input, harmonic number and index of sum n
 *@param *params: void pointer to struct parameters, which contains
 * n (harmonic number) and nu (frequency of emission/absorption)
 *@returns: Output, result of integral over gamma
 */
double gamma_integration_result(double n, void * params)
{
  double nu = *(double *) params;
  double nu_c = (electron_charge * B_field)
                /(2. * M_PI * mass_electron * speed_light);
  double gamma_minus = ((n*nu_c)/nu - fabs(cos(observer_angle))
                  *sqrt((pow((n*nu_c)/nu, 2.)) - pow(sin(observer_angle), 2.)))
                  /(pow(sin(observer_angle), 2));
  double gamma_plus  = ((n*nu_c)/nu + fabs(cos(observer_angle))
                  *sqrt((pow((n*nu_c)/nu, 2.)) - pow(sin(observer_angle), 2.)))
                  /(pow(sin(observer_angle), 2));
  double result = 0.;

  //needs help resolving peak for nu/nu_c > 1e6
  //gamma_peak is an accurate numerical fit to peak location
  double gamma_peak = 1.33322780e-06 * n / ((nu/nu_c) / 1.e6);
  double width = 0.;
  if (nu/nu_c < 3.e8) width = 10.;
  else               width = 1000.;

  double gamma_minus_high = gamma_peak - (gamma_peak-gamma_minus)/width;
  double gamma_plus_high = gamma_peak + (gamma_plus-gamma_peak)/width;

  if (nu/nu_c > 1.e6) {
    result = gsl_integrate(gamma_minus_high, gamma_plus_high, n, nu);
  }
  else {
  result = gsl_integrate(gamma_minus, gamma_plus, n, nu);
  }

  if(isnan(result) != 0) result = 0.;

  return result;
}

/*n_integration: j_nu and alpha_nu are given by an integral over gamma of
 * an integrand that contains a sum over n; we do the integral over gamma
 * and then the sum over n.  For numerical accuracy and speed, we only sum up
 * to n = n_max = 30, and integrate from n_max to the point where our integral
 * stops giving appreciable contributions.  n_integration() does the n integral.
 *@params n_minus: Input, minimum n for which the integrand is real
 *@params nu: Input, frequency of absorption/emission
 *@returns: Output, result of integration over n
 */
double n_integration(double n_minus, double nu)
{
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);
  if (DISTRIBUTION_FUNCTION == THERMAL && nu/nu_c < 1e6) {
    n_max = (int)(n_max + n_minus + 1.);
    double ans = gsl_integrate(n_max, C * n_peak(nu), -1, nu);
    return ans;
    }
  else {
    double n_start = (int)(n_max + n_minus + 1.);
    double ans = 0.;
    double contrib = 0.;
    double delta_n = 1.e2;
    double deriv_tol = 1.e-5;
    double tolerance = 1.e10;

    while (fabs(contrib) >= fabs(ans/tolerance)) {
      double deriv = derivative(n_start, nu);

      if(fabs(deriv) < deriv_tol) delta_n = 100. * delta_n;

      contrib = gsl_integrate(n_start, (n_start + delta_n), -1, nu);
      ans = ans + contrib;
      n_start = n_start + delta_n;

    }

  return ans;
  }
}

/*n_summation: performs the sum over n (harmonic number) from n = 1 to value 
 * n_max = 30, set heuristically in [1].  Past n_max we integrate over n using
 * n_integration().
 *@param nu: Input, frequency of absorption/emission
 *@returns: j_nu or alpha_nu, the emissivity or absorptivity, respectively
 */
double n_summation(double nu)
{
  #if POL_MODE == STOKES_U
    return 0.;
  #endif
  double ans = 0.;
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);
  double n_minus = (nu/nu_c) * fabs(sin(observer_angle));
  int x = (int)(n_minus+1.);
  for (x; x <= n_max + (int)n_minus ; x++) {
    ans += gamma_integration_result(x, &nu);
  }

  ans = ans + n_integration(n_minus, nu);
  return ans;
}

/*derivative: wrapper for the GSL derivative gsl_deriv_central(); used to
 *help speed up the integration over harmonic number n
 *@param n_start: Input, value of n at which the derivative is to be performed
 *@param nu: Input, frequency of absorption/emission
 */
double derivative(double n_start, double nu)
{
  gsl_function F;
  double result;
  double abserr;
  F.function = gamma_integration_result;
  F.params = &nu;
  gsl_deriv_central(&F, n_start, 1e-8, &result, &abserr);
  return result;
}

/*normalize_f: normalizes the distribution function (power-law with cutoff
 * or kappa) using GSL's QAGIU integrator.
 *@returns: 1 over the normalization constant for the chosen distribution
 */
double normalize_f()
{
  static double ans = 0;
  if (ans != 0) return ans;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
  double result, error;
  gsl_function F;
  #if DISTRIBUTION_FUNCTION == POWER_LAW
    F.function = &power_law_to_be_normalized;
  #elif DISTRIBUTION_FUNCTION == KAPPA_DIST
    F.function = &kappa_to_be_normalized;
  #else
    return 0;
  #endif

  double unused = 0.;
  F.params = &unused;
  gsl_integration_qagiu(&F, 1, 0, 1e-8, 1000,
                         w, &result, &error);
  gsl_integration_workspace_free (w);
  ans = result;
  return result;
}

/*gsl_integrate: wrapper for the GSL QAG integration routine
 *@param min: Input, lower bound of integral
 *@param max: Input, upper bound of integral
 *@param n: Input, harmonic number and index of sum over n
 *@param nu: Input, frequency of emission/absorption
 *@returns: Output, result of either gamma integral or n integral
 */
double gsl_integrate(double min, double max, double n, double nu)
{
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);
  if (nu/nu_c > 1.e6 || POL_MODE == STOKES_V) {
    gsl_set_error_handler_off();
  }
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
  double result, error;
  gsl_function F;
  if (n < 0) { //do n integration
    F.function = &gamma_integration_result;
    F.params = &nu;
    gsl_integration_qag(&F, min, max, 0., 1.e-3, 1000,
                         3,  w, &result, &error);
  }
  else {  //do gamma integration
  struct parameters n_and_nu;
  n_and_nu.n = n;
  n_and_nu.nu = nu;
  F.function = &gamma_integrand;
  F.params = &n_and_nu;
  gsl_integration_qag(&F, min, max, 0., 1.e-3, 1000,
                       3,  w, &result, &error);
  }

  gsl_integration_workspace_free (w);
  return result;
}
