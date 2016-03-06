#include "kappa.h"

/*kappa_law_to_be_normalized: Unnormalized kappa distribution function 
 *                            (eq. 18 of [1]), along with an exponential
 *                            cutoff because the kappa distribution has a
 *                            power-law tail at high energy.  The kappa
 *                            distribution is normalized numerically by 
 *                            the functionnormalize_f(), in the file
 *                            integrate.c.  One can "turn off" the exponential 
 *                            cutoff by setting the parameter gamma_cutoff to 
 *                            be very large, in which case the cutoff term is 
 *                            approximately equal to 1.
 *
 *@params: Lorentz factor gamma, void pointer to struct of parameters
 *@returns: unnormalized kappa distribution function with exponential
 *          cutoff.
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

/*kappa_f: normalized kappa distribution function with exponential cutoff.  
 *         This is normalized by calling normalize_f(), in the file 
 *         integrate.c, which uses GSL's QAGIU integrator.
 *
 *@params: Lorentz factor gamma, struct of parameters params
 *@returns: normalized kappa distribution function with exponential cutoff.
 */

double kappa_f(double gamma, struct parameters * params)
{

  /*no analytic estimate for where n-space peak is; must find it adaptively*/
//  if(params->use_n_peak != 0) params->use_n_peak = 0;

  double norm = 1./normalize_f(&kappa_to_be_normalized, params);

  double kappa_body = params->electron_density * pow((1. + (gamma - 1.)
                     /(params->kappa * params->kappa_width)), -params->kappa-1);

  double cutoff = exp(-gamma/params->gamma_cutoff);

  double ans = norm * kappa_body * cutoff;

  return ans;
}

/*differential_of_kappa: The integrand for the absorptivity calculation 
 *                       ([1] eq. 12) depends on a differential of the 
 *                       distribution function ([1] eq. 13).  For the 
 *                       kappa distribution, this is evaluated analytically 
 *                       for speed and accuracy. 
 *
 *@params: Lorentz factor gamma, struct of parameters params
 *@returns: the differential of the kappa distribution function
 *          ([1] eq. 12, 13, and 18) for the alpha_nu() calculation.
 *
 */
double differential_of_kappa(double gamma, struct parameters * params) 
{
  double prefactor = params->electron_density * (1./normalize_f(&kappa_to_be_normalized, params)) 
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
