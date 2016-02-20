#include "power_law.h"

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
double power_law_to_be_normalized(double gamma, void * paramsInput) 
{
  struct parameters * params = (struct parameters*) paramsInput;

  double norm_term = 4. * params->pi;

  double prefactor = (params->power_law_p - 1.) / 
                     (pow(params->gamma_min, 1. - params->power_law_p) 
                      - pow(params->gamma_max, 1. 
                            - params->power_law_p));

  double body = pow(gamma, -params->power_law_p) 
                * exp(- gamma / params->gamma_cutoff);

  double ans = norm_term * prefactor * body;

  return ans;
}

/*power_law_f: power-law distribution function, normalized via call to the 
 * normalize_f() function. Uses eq. 50 of [1].
 *@param gamma: Input, Lorentz factor
 *@returns: Ouput, a normalized power-law distribution for the gamma integrand
 */
double power_law_f(double gamma, struct parameters * params) 
{
  double beta = sqrt(1. - 1./(gamma*gamma));

  double prefactor = params->electron_density * (params->power_law_p - 1.) 
                     / (pow(params->gamma_min, 1. - params->power_law_p) 
                        - pow(params->gamma_max, 1. - params->power_law_p));

  double body = pow(gamma, -params->power_law_p) 
                * exp(- gamma / params->gamma_cutoff);

  double ans = 1./normalize_f(params) * prefactor * body 
                              * 1./(pow(params->mass_electron, 3.) 
                              * pow(params->speed_light, 3.) 
                              * gamma*gamma * beta);

  return ans;

  /* TODO: Can this be done using the following?
   * \int f d^3p = m_e n_e c^2
   * double power_law_unnormalized() 
   * { return pow(gamma, -params->power_law_p) };
   *
   * return power_law_unnormalized/normalize_f(power_law_to_be_normalized);
   */
}

double differential_of_power_law(double gamma, struct parameters * params)
{
  double pl_norm = 4.* params->pi/(normalize_f(params));

  double prefactor = (  params->pi * params->nu 
                      / (params->mass_electron * pow(params->speed_light, 2.))
                     )
                     * (params->electron_density*(params->power_law_p-1.))
                     / ( (  pow(params->gamma_min, 1.-params->power_law_p) 
                          - pow(params->gamma_max, 1.-params->power_law_p)
                         )
                       );

  double term1 = ((-params->power_law_p-1.)*exp(-gamma/params->gamma_cutoff)
             *pow(gamma,-params->power_law_p-2.)/(sqrt(gamma*gamma - 1.)));

  double term2 = (exp(-gamma/params->gamma_cutoff) * pow(gamma,(-params->power_law_p-1.))
                /(params->gamma_cutoff * sqrt(gamma*gamma - 1.)));

  double term3 = (exp(-gamma/params->gamma_cutoff) * pow(gamma,-params->power_law_p))
             /pow((gamma*gamma - 1.), (3./2.));

  double Df = pl_norm * prefactor * (term1 - term2 - term3);

  return Df;
}
