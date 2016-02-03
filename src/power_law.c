#include "symphony.h"

/*differential_of_f: term in gamma integrand only for absorptivity calculation; 
 *it is the differential Df = 2\pi\nu (1/(mc)*d/dgamma + (beta cos(theta)
 * -cos(xi))/(p*beta*c) * d/d(cos(xi))) f 
 *this is eq. 41 of [1]
 *below it is applied for the THERMAL, POWER_LAW, and KAPPA_DIST distributions
 *@param gamma: Input, Lorentz factor
 *@param nu: Input, frequency of emission/absorption
 *@returns: Output, Df term in gamma integrand; depends on distribution function
 */
double differential_of_f(double gamma, double nu) 
{
  /*described in Section 3 of [1] */

  #if DISTRIBUTION_FUNCTION == THERMAL
    double prefactor = (M_PI * nu / (electron_mass*speed_light*speed_light)) 
		     * (n_e/(theta_e * gsl_sf_bessel_Kn(2, 1./theta_e)));

    double body = (-1./theta_e) * exp(-gamma/theta_e);

    double f = prefactor * body;

  #elif DISTRIBUTION_FUNCTION == POWER_LAW
    double pl_norm = 4.* M_PI/(normalize_f(power_law_to_be_normalized));

    double prefactor = (M_PI * nu / (electron_mass*speed_light*speed_light)) 
                     * (n_e_NT*(power_law_p-1.))
                     /((pow(gamma_min, 1.-power_law_p) 
                     - pow(gamma_max, 1.-power_law_p)));

    double term1 = ((-power_law_p-1.)*exp(-gamma/gamma_cutoff)
               *pow(gamma,-power_law_p-2.)/(sqrt(gamma*gamma - 1.)));

    double term2 = (exp(-gamma/gamma_cutoff) * pow(gamma,(-power_law_p-1.))
                  /(gamma_cutoff * sqrt(gamma*gamma - 1.)));

    double term3 = (exp(-gamma/gamma_cutoff) * pow(gamma,-power_law_p))
               /pow((gamma*gamma - 1.), (3./2.));

    double f = pl_norm * prefactor * (term1 - term2 - term3);

  #elif DISTRIBUTION_FUNCTION == KAPPA_DIST
    double prefactor = n_e * (1./normalize_f(kappa_to_be_normalized)) * 4. * M_PI*M_PI * nu
                     * electron_mass*electron_mass * speed_light;

    double term1 = ((- kappa - 1.) / (kappa * kappa_width)) 
		  * pow((1. + (gamma - 1.)/(kappa * kappa_width)), -kappa-2.);

    double term2 = pow((1. + (gamma - 1.)/(kappa * kappa_width)),(- kappa - 1.)) 
                    * (- 1./gamma_cutoff);

    double f = prefactor * (term1 + term2) * exp(-gamma/gamma_cutoff);

  #endif

    return f;
}


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
double power_law_to_be_normalized(double gamma, void * params) 
{
  double norm_term = 4. * M_PI;

  double prefactor = (power_law_p - 1.) / (pow(gamma_min, 1. - power_law_p) 
                    - pow(gamma_max, 1. - power_law_p));

  double body = pow(gamma, -power_law_p) * exp(- gamma / gamma_cutoff);
  //double body = pow(gamma, -power_law_p); for PL w/o cutoff

  double ans = norm_term * prefactor * body;

  return ans;
}

/*power_law_f: power-law distribution function, normalized via call to the 
 * normalize_f() function. Uses eq. 50 of [1].
 *@param gamma: Input, Lorentz factor
 *@returns: Ouput, a normalized power-law distribution for the gamma integrand
 */
double power_law_f(double gamma) 
{
  double beta = sqrt(1. - 1./(gamma*gamma));

  double prefactor = n_e_NT * (power_law_p - 1.) / (pow(gamma_min, 1. 
                     - power_law_p) - pow(gamma_max, 1. - power_law_p));

  double body = pow(gamma, -power_law_p) * exp(- gamma / gamma_cutoff);
  //double body = pow(gamma, -power_law_p); //for PL w/o cutoff

  double ans = 1./normalize_f(power_law_to_be_normalized) 
               * prefactor * body * 1./(pow(electron_mass, 3.) 
               * pow(speed_light, 3.) * gamma*gamma * beta);

  return ans;
}

/*kappa_to_be_normalized: the kappa distribution is not normalized, so we must 
 * find the normalization constant by taking 1 over the integral of the function
 * kappa_to_be_normalized from 1 to infinity.  Uses eq. 42 of [2].
 *@param gamma: Input, Lorentz factor
 *@param *params: void pointer to struct parameters, which contains
 * n (harmonic number) and nu (frequency of emission/absorption)
 *@returns: kappa distribution function to be normalized via normalize_f()
 */
double kappa_to_be_normalized(double gamma, void * params)
{
  double kappa_body = pow((1. + (gamma - 1.)/(kappa * kappa_width)), -kappa-1);

  double cutoff = exp(-gamma/gamma_cutoff);

  double norm_term = 4. * M_PI * pow(electron_mass, 3.) * pow(speed_light, 3.) 
                   * gamma * sqrt(gamma*gamma-1.);

  double ans = kappa_body * cutoff * norm_term;
  //double ans = kappa_body * norm_term; //for kappa w/o cutoff

  return ans;
}

/*kappa_f: kappa distribution function, numerically normalized
 *uses eq. 42 of [2]
 *@param gamma: Input, Lorentz factor
 *@returns: normalized kappa distribution function to go into gamma integrand 
 */
double kappa_f(double gamma)
{
  double norm = 1./normalize_f(kappa_to_be_normalized);

  double kappa_body = n_e * pow((1. + (gamma - 1.)
                     /(kappa * kappa_width)), -kappa-1);

  double cutoff = exp(-gamma/gamma_cutoff);

  double ans = norm * kappa_body * cutoff;
  //double ans = norm * kappa_body; //for kappa w/o cutoff

  return ans;
}

/*gamma_integrand: full gamma integrand, to be integrated from gamma_minus
 * to gamma_plus, set by eq. 64 of [1]
 *@param gamma: Lorentz factor
 *@param *params: void pointer to struct parameters, which contains
 * n (harmonic number) and nu (frequency of emission/absorption)
 *@returns: Output, integrand to be integrated by gamma_integration_result()
 */
double gamma_integrand(double gamma, void * params)
{
  struct parameters paramsLocal = *(struct parameters*) params;

  double n                = paramsLocal.n;
  double nu               = paramsLocal.nu;
  double stokes_v_switch  = paramsLocal.stokes_v_switch;
  double magnetic_field   = paramsLocal.magnetic_field;
  double electron_charge  = paramsLocal.electron_charge;
  double electron_mass    = paramsLocal.electron_mass;
  double electron_density = paramsLocal.electron_density;
  double observer_angle   = paramsLocal.observer_angle;
  double mode             = paramsLocal.mode; /* Emissivities or absorptivities
                                                 TODO: Rotativities */
  double polarization_mode= paramsLocal.polarization_mode;
  double *function_pointer= paramsLocal.function_pointer;



  double nu_c = (electron_charge * B)
               /(2. * M_PI * electron_mass * speed_light);

  double beta = sqrt(1. - 1./(gamma*gamma));

  if (mode == EMISSIVITY)
  {
    double cos_xi = (gamma * nu - n * nu_c)
                   /(gamma * nu * beta * cos(observer_angle));

    double prefactor = 1./(nu * beta * fabs(cos(observer_angle)));

    double ans = prefactor * integrand(gamma, n, nu, function_pointer);
  }
  else if (MODE == ABSORPTIVITY)
  {
    double prefactor = -speed_light*electron_charge*electron_charge / (2. * nu);

    double ans = prefactor*gamma*gamma*beta*differential_of_f(gamma, nu, function_pointer)
                *polarization_term(gamma, n, nu)
                *(1./(nu*beta*fabs(cos(observer_angle))));
  }
  
  return ans;
}

