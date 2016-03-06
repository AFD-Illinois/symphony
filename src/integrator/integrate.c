#include "integrate.h"

/*gamma_integration_result: integrates gamma_integrand() from gamma_minus to
 *                          gamma_plus (described in section 4.1 of [1]).
 *                          Calls the function gamma_integral() in file
 *                          integrate.c, which is a wrapper for the GSL
 *                          integrator QAG.
 *
 *@params: harmonic number n, 
 *         void pointer to struct of parameters paramsInput
 *@returns: result of integrating the gamma integrand over gamma.  Note that
 *          this still remains to be summed over n.
 */
double gamma_integration_result(double n, void * paramsInput)
{
  struct parameters * params = (struct parameters*) paramsInput;

  double nu_c = get_nu_c(*params);

  double gamma_minus =  
    (   (n*nu_c)/params->nu - fabs(cos(params->observer_angle))
      * sqrt(  pow((n*nu_c)/params->nu, 2.)
             - pow(sin(params->observer_angle), 2.)
            )
    )
    / (pow(sin(params->observer_angle), 2));

  double gamma_plus =  
    (   (n*nu_c)/params->nu + fabs(cos(params->observer_angle))
      * sqrt(  pow((n*nu_c)/params->nu, 2.)
             - pow(sin(params->observer_angle), 2.)
            )
    )
    / (pow(sin(params->observer_angle), 2));

  double result = 0.;

  /*integrator needs help resolving peak for nu/nu_c > 1e6*/

  double gamma_peak = (gamma_plus+gamma_minus)/2.;

  double width = 1.;

  /* TODO: Remove hard coded limits here and put in params */
  if (params->nu/nu_c > 1.e6 && params->nu/nu_c <= 3.e8) width = 10.;
  else if (params->nu/nu_c > 3.e8)               width = 1000.;

  double gamma_minus_high = gamma_peak - (gamma_peak - gamma_minus)/width;
  double gamma_plus_high  = gamma_peak - (gamma_peak - gamma_plus) /width;

  /*Stokes V is hard to resolve; do 2 separate integrations to make it easier*/
  if(params->polarization == params->STOKES_V && params->stokes_v_switch >= 0) 
  {
    double neg_result = gamma_integral(gamma_minus_high, gamma_peak, n, params);
    double pos_result = gamma_integral(gamma_peak, gamma_plus_high,  n, params);

    if (params->stokes_v_switch == 0)
    {
      result = pos_result;
    }
    else
    {
      result = neg_result;
    }
  }

  if (params->polarization != params->STOKES_V || params->stokes_v_switch < 0) 
  {
    result = gamma_integral(gamma_minus_high, gamma_plus_high, n, params);
  }

  /*GSL QAG sometimes erroneously gives NaN instead of small values; 
    return 0 instead */
  if(isnan(result) != 0) result = 0.;

  return result;
}

/*n_integration: j_nu() and alpha_nu() are given by an integral over gamma of
 *               an integrand that contains a sum over n; we do the integral 
 *               over gamma and then the sum over n.  For numerical accuracy 
 *               and speed, we only sum up to n = n_max = 30, and integrate 
 *               from n_max to the point where our integral stops giving 
 *               appreciable contributions.  n_integration() performs the
 *               n integral by calling n_integral(), in integrate.c, which
 *               is a wrapper for GSL's QAG integrator.
 *
 *@params: n_minus is minimum n for which the integrand is real,
 *         struct of parameters params
 *@returns: result of integral over n, evaluated adaptively from n_minus until the integral
 *          contributions become approximately negligible.
 */
double n_integration(double n_minus, 
                     double (*n_peak)(struct parameters * params),
                     struct parameters * params
                    )
{
  double nu_c = get_nu_c(*params);

  double n_start = (int)(params->n_max + n_minus + 1.);

  /*For the MAXWELL_JUETTNER distribution, the n-space peak location is known
    analytically; this speeds up evaluation of MAXWELL_JUETTNER until
    nu/nu_c = 1e6, where the approximation breaks down and the adaptive 
    procedure works better.*/
  if (params->use_n_peak == 1 && params->nu/nu_c < 1e6) 
  {
    double ans = n_integral(n_start, params->C * n_peak(params), params);
    return ans;
  }
  else /*For other distributions, the n-space peak is found adaptively. */
  {
    double ans = 0.;
    double contrib = 0.;

    /*set parameters on adaptive n integration routine*/
    double delta_n = 1.e5;
    double deriv_tol = 1.e-5;
    double tolerance = 1.e5;

    // TODO: Put a debug mode which activates the diagnostic printf statements.
    //printf("GOT INTO N INTEGRATION");

    /*keep taking steps and integrating in n until the integral stops giving
      contributions greater than tolerance */
    while (fabs(contrib) >= fabs(ans/tolerance)) 
    {
      double deriv = derivative_of_n(n_start, params);
      if(fabs(deriv) < deriv_tol) delta_n = 100. * delta_n;

      contrib = n_integral(n_start, (n_start + delta_n), params);
      ans = ans + contrib;

      n_start = n_start + delta_n;
    }

    return ans;
  }
}

/*n_summation: performs the sum over n (harmonic number) from n = 1 to value 
 *             n_max = 30, set heuristically in [2].  Past n_max we integrate 
 *             over n using n_integration().
 *
 *@params: struct of parameters params
 *@returns: j_nu or alpha_nu, the emissivity or absorptivity, respectively,
 *          depending if the user has called j_nu() or alpha_nu().
 */
double n_summation(struct parameters *params)
{

  if(params->polarization == params->STOKES_U)
  {
    return 0.;
  }

  double ans = 0.;

  double nu_c    = get_nu_c(*params);
  double n_minus = (params->nu/nu_c) * fabs(sin(params->observer_angle)); 

  //need 2 separate n integrations to numerically resolve STOKES_V
  params->stokes_v_switch = -1; // TODO: describe: -1: , 1: 

  /*perform n summation by summing the result of the gamma integral for 
    each value of n from 1 to n_max*/
  for (int n=(int)(n_minus+1.); n <= params->n_max + (int)n_minus ; n++) 
  {
     ans += gamma_integration_result(n, params);
  }

  params->stokes_v_switch = 0;

  /*add result of n sum from 1 to n_max to an integral over n from n_max to
    the point where the integral no longer gives appreciable contributions*/
  ans += n_integration(n_minus, params->n_peak, params);

  /*if doing Stokes V, must perform two n integrals: one over the positive
    part of gamma integrand and one over the negative part; this helps GSL QAG
    resolve the n integrals*/
  if(params->polarization == params->STOKES_V)
  {
    params->stokes_v_switch = 1;
    ans += n_integration(n_minus, params->n_peak, params);
  }
  
  return ans;
}

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
                      struct parameters * params
                     )
{
  double nu_c = get_nu_c(*params);

  struct parametersGSL paramsGSL;
  paramsGSL.params = *params;
  paramsGSL.n      = n;

//  set_distribution_function(&paramsGSL);  //TODO: move this?

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
                  struct parameters * params
                 )
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
