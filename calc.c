#include "symphony.h"
#include "params.h"

double B;
double n_e;
double obs_angle;
double j_nu(double nu, double B_temp, double n_e_temp, double obs_angl_temp)
{
  B = B_temp;
  n_e = n_e_temp;
  obs_angle = obs_angl_temp;
  return n_summation(nu);
}

double alpha_nu(double nu, double B_temp, double n_e_temp, double obs_angl_temp)
{
  B = B_temp;
  n_e = n_e_temp;
  obs_angle = obs_angl_temp;
  return n_summation(nu);
}

/*n_peak: gives the location of the peak of the n integrand for 
 *the THERMAL distribution; uses Eq. 68 in [1]
 *
 *@param nu: Input, frequency of emission/absorption
 *@returns n_peak: Output, location of integrand's peak for the 
 * n-integral for the THERMAL distribution 
 */
double n_peak(double nu)
{
  double nu_c = (electron_charge * B)
	      / (2. * M_PI * mass_electron * speed_light);

  double beta = 0.;

  if (nu <= nu_c * theta_e*theta_e || theta_e < 1.) 
  {
    beta = sqrt((1. - 1./pow((1. + theta_e),2.))); //beta in low nu limit
  }
  else 
  {
    beta = sqrt(1. - pow((2. * theta_e * nu / nu_c), -2./3.));//beta for high nu
  }
  
  double n_peak =  (theta_e + 1. + pow((2. * theta_e * nu / nu_c),1./3.))
                 * (nu/nu_c) * (1. - beta*beta * pow(cos(obs_angle),2.));

  return n_peak;
}

/*polarization_term: term in the gamma integrand that varies based upon the 
 *polarization mode; uses eq. 13-19, 22-25, 27-30, in [1]
 *
 *@param gamma: Input, Lorentz factor 
 *@param n: Input, index n of sum (synchrotron harmonic number)
 *@param nu: Input, frequency of emission/absorption
 *@returns: piece of gamma integrand that determines polarization mode 
 */
double polarization_term(double gamma, double n, double nu) 
{
  /*below calculation is described in Section 3.1 of [1]*/
  
  double nu_c = (electron_charge * B)
	       /(2. * M_PI * mass_electron * speed_light);

  double beta = sqrt(1. - 1./(gamma*gamma));

  double cos_xi = (gamma * nu - n * nu_c)
		 /(gamma * nu * beta * cos(obs_angle));

  double M = (cos(obs_angle) - beta * cos_xi)/sin(obs_angle);

  double N = beta * sqrt(1 - (cos_xi*cos_xi));

  double z = (nu * gamma * beta * sin(obs_angle) 
       * sqrt(1. - cos_xi*cos_xi))/nu_c;

  double K_xx = M*M * pow(my_Bessel_J(n, z), 2.);

  double K_yy = N*N * pow(my_Bessel_dJ(n, z), 2.);

  double T_x = (2.*nu*cos(obs_angle))/(nu_c*pow(sin(obs_angle), 2.)
          +sqrt(nu_c*nu_c*pow(sin(obs_angle), 4.)
          +4.*nu*nu*pow(cos(obs_angle), 2.)));

  #if POL_MODE == STOKES_I
    double ans = K_xx + K_yy;

  #elif POL_MODE == STOKES_Q
    double ans = K_xx - K_yy;

  #elif POL_MODE == STOKES_U
  double ans = 0.;

  #elif POL_MODE == STOKES_V
    double ans = -2.*M*N*my_Bessel_J(n, z)*my_Bessel_dJ(n, z);
  #endif

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
double differential_of_f(double gamma, double nu) 
{
  /*described in Section 3 of [1] */

  #if DISTRIBUTION_FUNCTION == THERMAL
    double prefactor = (M_PI * nu / (mass_electron*speed_light*speed_light)) 
		     * (n_e/(theta_e * gsl_sf_bessel_Kn(2, 1./theta_e)));

    double body = (-1./theta_e) * exp(-gamma/theta_e);

    double f = prefactor * body;

  #elif DISTRIBUTION_FUNCTION == POWER_LAW
    double pl_norm = 4.* M_PI/(normalize_f());

    double prefactor = (M_PI * nu / (mass_electron*speed_light*speed_light)) 
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
    double prefactor = n_e * (1./normalize_f()) * 4. * M_PI*M_PI * nu
                     * mass_electron*mass_electron * speed_light;

    double term1 = ((- kappa - 1.) / (kappa * kappa_width)) 
		  * pow((1. + (gamma - 1.)/(kappa * kappa_width)), -kappa-2.);

    double term2 = pow((1. + (gamma - 1.)/(kappa * kappa_width)),(- kappa - 1.)) 
                    * (- 1./gamma_cutoff);

    double f = prefactor * (term1 + term2) * exp(-gamma/gamma_cutoff);

  #endif

    return f;
}

/*maxwell_juttner_f: Maxwell-Juttner distribution function in terms of Lorentz 
 * factor gamma; uses eq. 47, 49, 50 of [1]
 *@param gamma: Input, Lorentz factor
 *@returns: THERMAL distribution function, which goes into the gamma integrand 
 */
double maxwell_juttner_f(double gamma) 
{
  double beta = sqrt(1. - 1./(gamma*gamma));

  double d = (n_e * gamma * sqrt(gamma*gamma-1.) * exp(-gamma/theta_e))
            /(4. * M_PI * theta_e * gsl_sf_bessel_Kn(2, 1./theta_e));

  double ans = 1./(pow(mass_electron, 3.) * pow(speed_light, 3.) * gamma*gamma 
               * beta) * d;

  return ans;
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

  double ans = 1./normalize_f() * prefactor * body * 1./(pow(mass_electron, 3.) 
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

  double norm_term = 4. * M_PI * pow(mass_electron, 3.) * pow(speed_light, 3.) 
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
  double norm = 1./normalize_f();

  double kappa_body = n_e * pow((1. + (gamma - 1.)
                     /(kappa * kappa_width)), -kappa-1);

  double cutoff = exp(-gamma/gamma_cutoff);

  double ans = norm * kappa_body * cutoff;
  //double ans = norm * kappa_body; //for kappa w/o cutoff

  return ans;
}

/*integrand_without_extra_factor: the function I(n, xi, gamma) from eq.
 * 60 of [1].
 *@param gamma: Input, Lorentz factor
 *@param n: Input, harmonic number; index of sum in gamma integrand
 *@param nu: Input, frequency of emission/absorption
 *@returns: Ouput, gamma integrand without terms pulled out of integral in
 * eq. 60, 62 of [1]
 */
double integrand_without_extra_factor(double gamma, double n, double nu)
{
  /*This is the function I(n, xi, gamma) in Eq. 60 of [1]*/

  double nu_c = (electron_charge * B)
               /(2. * M_PI * mass_electron * speed_light);

  double beta = sqrt(1. - 1./(gamma*gamma));

  double cos_xi = (gamma * nu - n * nu_c)
                 /(gamma * nu * beta * cos(obs_angle));

  #if DISTRIBUTION_FUNCTION == THERMAL
    double ans = (2. * M_PI * electron_charge*electron_charge * nu*nu)
                / speed_light * (pow(mass_electron, 3.) * pow(speed_light, 3.) 
                * gamma*gamma * beta * 2. * M_PI) * maxwell_juttner_f(gamma) 
                * polarization_term(gamma, n, nu);

  #elif DISTRIBUTION_FUNCTION == POWER_LAW
    double ans = (2. * M_PI * electron_charge*electron_charge * nu*nu)
                / speed_light * (pow(mass_electron, 3.) * pow(speed_light, 3.) 
                * gamma*gamma * beta * 2. * M_PI) * power_law_f(gamma) 
                * polarization_term(gamma, n, nu);

  #elif DISTRIBUTION_FUNCTION == KAPPA_DIST
    double ans = (2. * M_PI * electron_charge*electron_charge * nu*nu)
                / speed_light * (pow(mass_electron, 3.) * pow(speed_light, 3.) 
                * gamma*gamma * beta * 2. * M_PI) * kappa_f(gamma) 
                * polarization_term(gamma, n, nu);
  #endif

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
  struct parameters n_and_nu = *(struct parameters*) params;
  double n = n_and_nu.n;
  double nu = n_and_nu.nu;

  double nu_c = (electron_charge * B)
               /(2. * M_PI * mass_electron * speed_light);

  double beta = sqrt(1. - 1./(gamma*gamma));

  #if MODE == EMISS
    double cos_xi = (gamma * nu - n * nu_c)
                   /(gamma * nu * beta * cos(obs_angle));

    double prefactor = 1./(nu * beta * fabs(cos(obs_angle)));

    double ans = prefactor * integrand_without_extra_factor(gamma, n, nu);

  #elif MODE == ABSORP
    double prefactor = -speed_light*electron_charge*electron_charge / (2. * nu);

    double ans = prefactor*gamma*gamma*beta*differential_of_f(gamma, nu)
                *polarization_term(gamma, n, nu)
                *(1./(nu*beta*fabs(cos(obs_angle))));
  #endif

    return ans;
}

/*gamma_integration_result: performs gamma integral from gamma_minus to 
 * gamma_plus; calls gsl_integrate(), which is a wrapper for the GSL integrator
 * GSL QAG.
 *@param n: Input, harmonic number and index of sum n
 *@param *params: void pointer to struct parameters, which contains
 * n (harmonic number) and nu (frequency of emission/absorption)
 *@returns: Output, result of integral over gamma
 */
double gamma_integration_result(double n, void * params)
{
  /*This evaluates the gamma integral*/
  double nu = *(double *) params;

  double nu_c = (electron_charge * B)
                /(2. * M_PI * mass_electron * speed_light);

  double gamma_minus = ((n*nu_c)/nu - fabs(cos(obs_angle))
                  *sqrt((pow((n*nu_c)/nu, 2.)) - pow(sin(obs_angle), 2.)))
                  /(pow(sin(obs_angle), 2));

  double gamma_plus  = ((n*nu_c)/nu + fabs(cos(obs_angle))
                  *sqrt((pow((n*nu_c)/nu, 2.)) - pow(sin(obs_angle), 2.)))
                  /(pow(sin(obs_angle), 2));

  double result = 0.;

  /*integrator needs help resolving peak for nu/nu_c > 1e6*/

  double gamma_peak = (gamma_plus+gamma_minus)/2.;

  double width = 1.;

  if (nu/nu_c > 1.e6 && nu/nu_c <= 3.e8) width = 10.;
  else if (nu/nu_c > 3.e8)               width = 1000.;

  double gamma_minus_high = gamma_peak - (gamma_peak-gamma_minus)/width;
  double gamma_plus_high = gamma_peak + (gamma_plus-gamma_peak)/width;


  /*Stokes V is hard to resolve; do 2 separate integrations to make it easier*/
  if(POL_MODE == STOKES_V && stokes_v_switch >= 0) 
  {
    double neg_result = gsl_integrate(gamma_minus_high, gamma_peak, n, nu);
    double pos_result = gsl_integrate(gamma_peak, gamma_plus_high, n, nu);

    if(stokes_v_switch == 0) result = pos_result;
    else                     result = neg_result;
  }

  if (POL_MODE != STOKES_V || stokes_v_switch < 0) 
  {
    result = gsl_integrate(gamma_minus_high, gamma_plus_high, n, nu);
  }

  /*GSL QAG sometimes erroneously gives NaN instead of small values; 
    return 0 instead */
  if(isnan(result) != 0) result = 0.;

  return result;
}

/*n_integration: j_nu and alpha_nu are given by an integral over gamma of
 * an integrand that contains a sum over n; we do the integral over gamma
 * and then the sum over n.  For numerical accuracy and speed, we only sum up
 * to n = n_max = 30, and integrate from n_max to the point where our integral
 * stops giving appreciable contributions.  n_integration() does the n integral.
 *@params n_minus: Input, minimum n for which the integrand is real
 *@params nu: Input, frequency of absorption/emission
 *@returns: Output, result of integration over n
 */
double n_integration(double n_minus, double nu)
{
  double nu_c = (electron_charge * B)
               /(2. * M_PI * mass_electron * speed_light);

  if (DISTRIBUTION_FUNCTION == THERMAL && nu/nu_c < 1e6) 
  {
    double n_start = (int)(n_max + n_minus + 1.);
    double ans = gsl_integrate(n_start, C * n_peak(nu), -1, nu);
    return ans;
  }

  else 
  {
    double n_start = (int)(n_max + n_minus + 1.);
    double ans = 0.;
    double contrib = 0.;

    /*set parameters on adaptive n integration routine*/
    double delta_n = 1.e5;
    double deriv_tol = 1.e-5;
    double tolerance = 1.e5;

    /*keep taking steps and integrating in n until the integral stops giving
      contributions greater than tolerance */
    while (fabs(contrib) >= fabs(ans/tolerance)) 
    {
      double deriv = derivative(n_start, nu);
      if(fabs(deriv) < deriv_tol) delta_n = 100. * delta_n;

      contrib = gsl_integrate(n_start, (n_start + delta_n), -1, nu);
      ans = ans + contrib;

      n_start = n_start + delta_n;
    }

  return ans;
  }
}

/*n_summation: performs the sum over n (harmonic number) from n = 1 to value 
 * n_max = 30, set heuristically in [1].  Past n_max we integrate over n using
 * n_integration().
 *@param nu: Input, frequency of absorption/emission
 *@returns: j_nu or alpha_nu, the emissivity or absorptivity, respectively
 */
double n_summation(double nu)
{

  #if POL_MODE == STOKES_U
    return 0.;
  #endif

  double ans = 0.;

  double nu_c = (electron_charge * B)
               /(2. * M_PI * mass_electron * speed_light);

  double n_minus = (nu/nu_c) * fabs(sin(obs_angle));

  stokes_v_switch = -1;

  /*perform n summation by summing the result of the gamma integral for 
    each value of n from 1 to n_max*/
  for (int x=(int)(n_minus+1.); x <= n_max + (int)n_minus ; x++) 
  {
    ans += gamma_integration_result(x, &nu);
  }

  stokes_v_switch = 0;

  /*add result of n sum from 1 to n_max to an integral over n from n_max to
    the point where the integral no longer gives appreciable contributions*/
  ans += n_integration(n_minus, nu);

  /*if doing Stokes V, must perform two n integrals: one over the positive
    part of gamma integrand and one over the negative part; this helps GSL QAG
    resolve the n integrals*/
  #if POL_MODE == STOKES_V
    stokes_v_switch = 1;
    ans += n_integration(n_minus, nu);
  #endif

  return ans;
}
