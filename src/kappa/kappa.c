#include "kappa.h"

/*kappa_to_be_normalized: the kappa distribution is not normalized, so we must 
 * find the normalization constant by taking 1 over the integral of the function
 * kappa_to_be_normalized from 1 to infinity.  Uses eq. 42 of [2].
 *@param gamma: Input, Lorentz factor
 *@param *params: void pointer to struct parameters, which contains
 * n (harmonic number) and nu (frequency of emission/absorption)
 *@returns: kappa distribution function to be normalized via normalize_f()
 */
double kappa_to_be_normalized(double gamma, void * paramsInput)
{
  struct parameters params = *(struct parameters*) paramsInput;

  double kappa_body = pow((1. + (gamma - 1.)
                           /(params.kappa * params.kappa_width)), 
                          -params.kappa-1);

  double cutoff = exp(-gamma/params.gamma_cutoff);

  double norm_term = 4. * params.pi * pow(params.mass_electron, 3.) 
                     * pow(params.speed_light, 3.) 
                     * gamma * sqrt(gamma*gamma-1.);

  double ans = kappa_body * cutoff * norm_term;

  return ans;
}

/*kappa_f: kappa distribution function, numerically normalized
 *uses eq. 42 of [2]
 *@param gamma: Input, Lorentz factor
 *@returns: normalized kappa distribution function to go into gamma integrand 
 */
double kappa_f(double gamma, struct parameters * params)
{
  double norm = 1./normalize_f(params);

  double kappa_body = params->electron_density * pow((1. + (gamma - 1.)
                     /(params->kappa * params->kappa_width)), -params->kappa-1);

  double cutoff = exp(-gamma/params->gamma_cutoff);

  double ans = norm * kappa_body * cutoff;

  return ans;
}

double differential_of_kappa(double gamma, struct parameters * params) 
{
  double prefactor = params->electron_density * (1./normalize_f(params)) 
                     * 4. * params->pi*params->pi * params->nu 
                     * params->mass_electron * params->mass_electron 
                     * params->speed_light;

  double term1 = ((- params->kappa - 1.) 
                  / (params->kappa * params->kappa_width)) 
                * pow((1. + (gamma - 1.) 
                      / (params->kappa * params->kappa_width)), 
                      -params->kappa-2.);

  double term2 = pow((1. + (gamma - 1.)
                     /(params->kappa * params->kappa_width)), 
                     (- params->kappa - 1.)) * (- 1./params->gamma_cutoff);

  double Df = prefactor * (term1 + term2) * exp(-gamma/params->gamma_cutoff);

  return Df;
}
