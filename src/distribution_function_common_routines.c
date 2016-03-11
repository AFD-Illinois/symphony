#include "distribution_function_common_routines.h"

/*normalize_f: normalizes the distribution function using GSL's 
 *             QAGIU integrator.   
 *
 *@params: struct of parameters to pass to distribution function
 *@returns: 1 over the normalization constant for the chosen distribution
 */
double normalize_f(double (*distribution)(double, void *),
                   struct parameters * params
                  )
{
  static double ans = 0;
  if (ans != 0) return ans;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
  double result, error;
  gsl_function F;

  F.function = distribution;

  F.params = params;
  // TODO: Remove hard coded limits
  gsl_integration_qagiu(&F, 1, 0, 1e-8, 1000,
                         w, &result, &error
                       );
  gsl_integration_workspace_free (w);
  ans = result;
  return result;
}


/*differential_of_f: The absorptivity integrand ([1] eq. 12) has a term 
 *                   dependent on a differential operator ([1] eq. 13) 
 *                   applied to the distribution function. This function
 *                   calls each individual differential, evaluated 
 *                   analytically, for the 3 distributions studied. 
 *                   These are evaluated analytically for speed, but
 *                   a calculator to evaluate this differential without
 *                   the analytic result is also provided. This calculator
 *                   works with any gyrotropic distribution function.
 * 
 *@params: Lorentz factor gamma, struct of parameters params
 *@returns: differential of the distribution function term in the gamma
 *          integrand.  
 */
double differential_of_f(double gamma, struct parameters * params) 
{
  /*described in Section 2 of [1] */

  double Df = 0.;

  if(params->distribution == params->MAXWELL_JUETTNER)
  {
    Df = differential_of_maxwell_juettner(gamma, params);
  }
  else if(params->distribution == params->POWER_LAW)
  {
    Df = differential_of_power_law(gamma, params);
  }
  else if(params->distribution == params->KAPPA_DIST)
  {
    Df = differential_of_kappa(gamma, params);
  }

  return Df;
}

//double num_differential_of_f(double gamma, struct parameters * params)
//{
//  double Df = 0.;
//  double epsilon = 0.0000001;
//
//  double prefactor = (2. * params->pi * params->nu
//                      / (params->mass_electron
//                         *params->speed_light*params->speed_light))
//                     * params->pi * pow(params->mass_electron, 3.) 
//                     * pow(params->speed_light, 3.);
//
//  Df =  (params->distribution_function(gamma+epsilon, params)
//         - params->distribution_function(gamma-epsilon, params))
//          / (2. * epsilon);
//
//  return prefactor * Df;
//}

