#include "symphony.h"

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
  double nu_c = (electron_charge * B)
               /(2. * M_PI * mass_electron * speed_light);
  if (nu/nu_c > 1.e6) 
  {
    gsl_set_error_handler_off();
  }

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
  double result, error;
  gsl_function F;
  if (n < 0) //do n integration
  {
    F.function = &gamma_integration_result;
    F.params = &nu;
    gsl_integration_qag(&F, min, max, 0., 1.e-3, 1000,
                         3,  w, &result, &error);
  }
  else      //do gamma integration
  {
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
