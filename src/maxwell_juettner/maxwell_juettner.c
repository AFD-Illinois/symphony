#include "maxwell_juettner.h"

/*n_peak: Calculates and returns the location of the peak of the n integrand for 
 *        the MAXWELL_JUETTNER distribution; uses Eq. 68 in [2] to do this. 
 *        This analytic estimate of peak location in n-space speeds up the 
 *        evaluation of j_nu() and alpha_nu() for the MAXWELL_JUETTNER 
 *        distribution.  Analytic locations of the peak in n-space for
 *        other distributions are not known, so peak locations are found
 *        adaptively in n_integration().
 *
 *@params: struct of parameters params
 *@returns: location of the peak of the n integrand for the 
 *          MAXWELL_JUETTNER distribution 
 */
double maxwell_juettner_n_peak(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double beta = 0.;

  if (   params->nu <= nu_c * params->theta_e * params->theta_e
      || params->theta_e < 1.
     ) /*beta low nu limit*/
  {
    beta = sqrt((1. - 1./pow((1. + params->theta_e),2.)));
  }
  else /*beta high nu limit */
  {
    beta = sqrt(1. - pow((2. * params->theta_e * params->nu / nu_c), -2./3.));
  }

  double n_peak =  (params->theta_e + 1. + pow((2. * params->theta_e * params->nu / nu_c),1./3.))
                 * (params->nu/nu_c) * (1. - beta*beta * pow(cos(params->observer_angle),2.));

  return n_peak;
}


/*maxwell_juettner_f: Relativistic thermal (Maxwell-Juettner) distribution
 *                    function (eq. 14 and 15 of [1]).
 *
 *@params: Lorentz factor gamma, struct of parameters params
 *@returns: the Maxwell-Juttner distribution function evaluated at the
 *          input Lorentz factor gamma and for provided parameters in
 *          struct params (such as dimensionless electron temperature
 *          theta_e). 
 */
double maxwell_juettner_f(double gamma, struct parameters * params) 
{
  /*we have an analytic estimate for where the n-space peak is*/
  if(params->use_n_peak != 1)
  {
    params->use_n_peak = 1;
    double (*MJ_n_peak)(struct parameters *) = &maxwell_juettner_n_peak;
    params->n_peak     = MJ_n_peak;
  }

  double beta = sqrt(1. - 1./(gamma*gamma));

  double term1 = (  params->electron_density * gamma * sqrt(gamma*gamma-1.) 
                  * exp(-gamma/params->theta_e)
                 ) 
                 / (4. * params->pi 
                       * params->theta_e 
                       * gsl_sf_bessel_Kn(2, 1./params->theta_e)
                   );

  double ans = 1./(  pow(params->mass_electron, 3.)
                   * pow(params->speed_light, 3.) 
                   * gamma*gamma * beta
                  ) * term1;

  return ans;
}

/*differential_of_maxwell_juettner: The integrand for the absorptivity
 *                                  calculation ([1] eq. 12) depends on
 *                                  a differential of the distribution 
 *                                  function ([1] eq. 13).  For the Maxwell-
 *                                  Juettner distribution, this is evaluated
 *                                  analytically for speed and accuracy. 
 *                                  TODO: numerical derivative calculator
 *                                  for any gyrotropic distribution function.
 *
 *@params: Lorentz factor gamma, struct of parameters params
 *@returns: the differential of the Maxwell-Juettner distribution function
 *          ([1] eq. 12 and 13) for the alpha_nu() calculation.
 *
 */
//double differential_of_maxwell_juettner(double gamma, struct parameters * params)
//{
//  double Df = 0.;
//
//  double prefactor = (params->pi * params->nu 
//                   / (params->mass_electron
//                      *params->speed_light*params->speed_light)) 
//	     * (params->electron_density
//                /(params->theta_e * gsl_sf_bessel_Kn(2, 1./params->theta_e)));
//
//  double body = (-1./params->theta_e) * exp(-gamma/params->theta_e);
//
//  Df = prefactor * body;
//
//  return Df;
//}

