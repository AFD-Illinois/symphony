#include "symphony.h"

///*derivative: wrapper for the GSL derivative gsl_deriv_central(); used to
// *help speed up the integration over harmonic number n
// *@param n_start: Input, value of n at which the derivative is to be performed
// *@param nu: Input, frequency of absorption/emission
// */
double derivative(double n_start, double nu)
{
  gsl_function F;
  double result;
  double abserr;
  F.function = gamma_integration_result;
  F.params = &nu;
  gsl_deriv_central(&F, n_start, 1e-8, &result, &abserr);
  return result;

  return 0.;
}
//
///*normalize_f: normalizes the distribution function (power-law with cutoff
// * or kappa) using GSL's QAGIU integrator.
// *@returns: 1 over the normalization constant for the chosen distribution
// */
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

  //double unused = 0.;
  //F.params = &unused;
  F.params = &params;
  gsl_integration_qagiu(&F, 1, 0, 1e-8, 1000,
                         w, &result, &error);
  gsl_integration_workspace_free (w);
  ans = result;

  return result;
}

///*gsl_integrate: wrapper for the GSL QAG integration routine
// *@param min: Input, lower bound of integral
// *@param max: Input, upper bound of integral
// *@param n: Input, harmonic number and index of sum over n
// *@param nu: Input, frequency of emission/absorption
// *@returns: Output, result of either gamma integral or n integral
// */
double gsl_integrate(double min, double max, double n, 
                     struct parameters * params
                    )
{
  double nu_c = get_nu_c(*params);
  if (params->nu/nu_c > 1.e6) 
  {
    gsl_set_error_handler_off();
  }

  struct parametersGSL paramsGSL;
  paramsGSL.params = *params;


  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
  double result, error;
  gsl_function F;
  if (n < 0) //do n integration
  {
    F.function = &gamma_integration_result;
    //F.params = &nu;
    F.params = params;
    gsl_integration_qag(&F, min, max, 0., 1.e-3, 1000,
                         3,  w, &result, &error);
  }
  else      //do gamma integration
  {
    F.function = &gamma_integrand;
    //struct parametersGSL paramsGSL;
    //paramsGSL.params = params;
    paramsGSL.n      = n;
    F.params   = &paramsGSL;
    gsl_integration_qag(&F, min, max, 0., 1.e-3, 1000,
                         3,  w, &result, &error
                       );
  }

  gsl_integration_workspace_free (w);

  return result;
}

//trapezoidal rule integrator, useful while restructuring
double gamma_integral(double min, double max, double n,
                 struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  struct parametersGSL paramsGSL;
  paramsGSL.params = *params;
  paramsGSL.n      = n;
//
//  printf("\nTRAPEZOID: %e\n", gamma_integrand((max+min)/2., &paramsGSL));
//
//  double result = 0.;
//
//  result = (max - min) * (gamma_integrand(max, &paramsGSL) 
//                        + gamma_integrand(min, &paramsGSL))/2.;
//  
//  printf("\n trapezoid result = %e", result);
//  
//  return result;

  double result, error;

  gsl_function F;
  F.function = &gamma_integrand;
  F.params = &paramsGSL;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);

  gsl_integration_qag(&F, min, max, 0., 1.e-7, 1000,
                         3,  w, &result, &error); 

  //printf ("result          = %e\n", result);

  gsl_integration_workspace_free (w);

  return result;


}

