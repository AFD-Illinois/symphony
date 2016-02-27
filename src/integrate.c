#include "integrate.h"

/*derivative_of_n: wrapper for the GSL derivative gsl_deriv_central(); used to
 *                 help speed up the integration over harmonic number n.
 *
 *@params: n_start is value of n at which the derivative is to be performed,
 *         struct of parameters params (to pass to gamma_integration_result)
 *@returns: the derivative in n-space at a point in the n integral, used
 *          to speed up n integration.
 */
double derivative_of_n(double n_start, struct parameters * params)
{
  gsl_function F;
  double result;
  double abserr;

  F.function = gamma_integration_result;
  F.params = params;

  gsl_deriv_central(&F, n_start, 1e-8, &result, &abserr);
  return result;
}

/*normalize_f: normalizes the distribution function using GSL's 
 *             QAGIU integrator.  TODO: modify this to take a function 
 *             pointer to the supplied distribution function.
 *
 *@params: struct of parameters to pass to distribution function
 *@returns: 1 over the normalization constant for the chosen distribution
 */
double normalize_f(struct parameters * params)
{
  static double ans = 0;
  if (ans != 0) return ans;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
  double result, error;
  gsl_function F;
  
  if(params->distribution == params->POWER_LAW)
  {
    F.function = &power_law_to_be_normalized;
  }
  else if(params->distribution == params->KAPPA_DIST)
  {
    F.function = &kappa_to_be_normalized;
  }
  else
  {
    return 0;
  }

  F.params = params;
  // TODO: Remove hard coded limits
  gsl_integration_qagiu(&F, 1, 0, 1e-8, 1000,
                         w, &result, &error
                       );
  gsl_integration_workspace_free (w);
  ans = result;
  //TODO: Debug mode 
  //printf("\nRESULT of NORM = %e\n", result);
  return result;
}

/*gamma_integral: evaluates an integral from gamma_minus to gamma_plus
 *                of the gamma integrand using GSL's QAG integrator.
 *
 *@params: min (lower bound of integral), max (upper bound of integral),
 *         n (harmonic number), struct of parameters params to be passed
 *         to gamma_integrand
 *@returns: result of the gamma integral.
 */
double gamma_integral(double min, 
                      double max, 
                      double n,
                      struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  struct parametersGSL paramsGSL;
  paramsGSL.params = *params;
  paramsGSL.n      = n;

  //TODO: describe why this is necessary
  if(params->nu/nu_c >= 1.e6)
  {
     gsl_set_error_handler_off();
  }

  double result, error;

  gsl_function F;
  F.function = &gamma_integrand;
  F.params = &paramsGSL;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);

  // TODO: Remove hard coded limits
  gsl_integration_qag(&F, min, max, 0., 1.e-3, 1000,
                         3,  w, &result, &error); 

  gsl_integration_workspace_free (w);

  return result;
}

/*n_integral: evaluates an integral over the n integrand between two
 *            limits provided by the adaptive routine in n_integration().
 *            This is evaluated using GSL's QAG integrator.
 *
 *@params: min (lower bound of integral), max (upper bound of integral),
 *         n (harmonic number), struct of parameters params to be passed
 *         to gamma_integration_result
 *@returns: result of the n integral.
 */

double n_integral(double min,
                  double max, 
                  struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  //TODO: describe why this is necessary
  if(params->nu/nu_c >= 1.e6)
  {
  gsl_set_error_handler_off();
  }

  double result, error;

  gsl_function F;
  F.function = &gamma_integration_result;
  F.params = params;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);

  // TODO: Remove hard coded limits
  gsl_integration_qag(&F, min, max, 0., 1.e-3, 1000,
                         3,  w, &result, &error);

  gsl_integration_workspace_free (w);

  return result;
}
