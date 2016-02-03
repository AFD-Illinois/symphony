#include "symphony.h"

/*maxwell_juttner_f: Maxwell-Juttner distribution function in terms of Lorentz 
 * factor gamma; uses eq. 47, 49, 50 of [1]
 *@param gamma: Input, Lorentz factor
 *@returns: THERMAL distribution function, which goes into the gamma integrand 
 */
double maxwell_juttner(double gamma) 
{
  double beta = sqrt(1. - 1./(gamma*gamma));

  double d = (n_e * gamma * sqrt(gamma*gamma-1.) * exp(-gamma/theta_e))
            /(4. * M_PI * theta_e * gsl_sf_bessel_Kn(2, 1./theta_e));

  double ans = 1./(pow(electron_mass, 3.) * pow(speed_light, 3.) * gamma*gamma 
               * beta) * d;

  return ans;
}

/*differential_of_f: term in gamma integrand only for absorptivity calculation; 
 *it is the differential Df = 2\pi\nu (1/(mc)*d/dgamma + (beta cos(theta)
 * -cos(xi))/(p*beta*c) * d/d(cos(xi))) f 
 *this is eq. 41 of [1]
 *below it is applied for the THERMAL, POWER_LAW, and KAPPA_DIST distributions
 *@param gamma: Input, Lorentz factor
 *@param nu: Input, frequency of emission/absorption
 *@returns: Output, Df term in gamma integrand; depends on distribution function
 */
double differential_of_maxwell_juttner(double gamma, double nu) 
{
  /*described in Section 3 of [1] */

  double prefactor = (M_PI * nu / (electron_mass*speed_light*speed_light)) 
	     * (n_e/(theta_e * gsl_sf_bessel_Kn(2, 1./theta_e)));

  double body = (-1./theta_e) * exp(-gamma/theta_e);

  double f = prefactor * body;

  return f;
}

/*integrand_without_extra_factor: the function I(n, xi, gamma) from eq.
 * 60 of [1].
 *@param gamma: Input, Lorentz factor
 *@param n: Input, harmonic number; index of sum in gamma integrand
 *@param nu: Input, frequency of emission/absorption
 *@returns: Ouput, gamma integrand without terms pulled out of integral in
 * eq. 60, 62 of [1]
 */
double integrand(double gamma, double n, double nu, function_pointer)
{
  /*This is the function I(n, xi, gamma) in Eq. 60 of [1]*/

  double nu_c = (electron_charge * B)
               /(2. * M_PI * electron_mass * speed_light);

  double beta = sqrt(1. - 1./(gamma*gamma));

  double cos_xi = (gamma * nu - n * nu_c)
                 /(gamma * nu * beta * cos(observer_angle));

  double ans = (2. * M_PI * electron_charge*electron_charge * nu*nu)
              / speed_light * (pow(electron_mass, 3.) * pow(speed_light, 3.) 
              * gamma*gamma * beta * 2. * M_PI) * function_pointer(gamma) 
              * polarization_term(gamma, n, nu);

  return ans;
}

