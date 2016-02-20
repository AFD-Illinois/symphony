#include "maxwell_juettner.h"

/*maxwell_juttner_f: Maxwell-Juttner distribution function in terms of Lorentz 
 * factor gamma; uses eq. 47, 49, 50 of [1]
 *@param gamma: Input, Lorentz factor
 *@returns: THERMAL distribution function, which goes into the gamma integrand 
 */
double maxwell_juttner_f(double gamma, struct parameters * params) 
{
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

double differential_of_maxwell_juttner(double gamma, struct parameters * params)
{
  double Df = 0.;

  // TODO: Follow equation formatting
  double prefactor = (params->pi * params->nu 
                   / (params->mass_electron*params->speed_light*params->speed_light)) 
	     * (params->electron_density/(params->theta_e * gsl_sf_bessel_Kn(2, 1./params->theta_e)));

  double body = (-1./params->theta_e) * exp(-gamma/params->theta_e);

  Df = prefactor * body;

  return Df;
}
