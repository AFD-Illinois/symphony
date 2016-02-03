#include "symphony.h"
#include "params.h"

/* Wrapper for n_summation */
double j_nu(double nu, 
            double magnetic_field,
            double electron_density,
            double observer_angle
           )
{
  return n_summation(nu, 
                     magnetic_field,
                     electron_density,
                     observer_angle,
                     0 /* Compute j_nu */
                     );
}

/* Wrapper for n_summation */
double alpha_nu(double nu, 
                double magnetic_field,
                double electron_density,
                double observer_angle
               )
{
  return n_summation(nu, 
                     magnetic_field,
                     electron_density,
                     observer_angle,
                     1 /* Compute alpha_nu */
                     );
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
	      / (2. * M_PI * electron_mass * speed_light);

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
                 * (nu/nu_c) * (1. - beta*beta * pow(cos(observer_angle),2.));

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
	       /(2. * M_PI * electron_mass * speed_light);

  double beta = sqrt(1. - 1./(gamma*gamma));

  double cos_xi = (gamma * nu - n * nu_c)
		 /(gamma * nu * beta * cos(observer_angle));

  double M = (cos(observer_angle) - beta * cos_xi)/sin(observer_angle);

  double N = beta * sqrt(1 - (cos_xi*cos_xi));

  double z = (nu * gamma * beta * sin(observer_angle) 
       * sqrt(1. - cos_xi*cos_xi))/nu_c;

  double K_xx = M*M * pow(my_Bessel_J(n, z), 2.);

  double K_yy = N*N * pow(my_Bessel_dJ(n, z), 2.);

  double T_x = (2.*nu*cos(observer_angle))/(nu_c*pow(sin(observer_angle), 2.)
          +sqrt(nu_c*nu_c*pow(sin(observer_angle), 4.)
          +4.*nu*nu*pow(cos(observer_angle), 2.)));

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
                /(2. * M_PI * electron_mass * speed_light);

  double gamma_minus = ((n*nu_c)/nu - fabs(cos(observer_angle))
                  *sqrt((pow((n*nu_c)/nu, 2.)) - pow(sin(observer_angle), 2.)))
                  /(pow(sin(observer_angle), 2));

  double gamma_plus  = ((n*nu_c)/nu + fabs(cos(observer_angle))
                  *sqrt((pow((n*nu_c)/nu, 2.)) - pow(sin(observer_angle), 2.)))
                  /(pow(sin(observer_angle), 2));

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
               /(2. * M_PI * electron_mass * speed_light);

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
double n_summation(double nu, 
                   double magnetic_field,
                   double electron_density,
                   double observer_angle,
                   double electron_charge,
                   double electron_mass,
                   double speed_light,
                   int pol_mode,
                   int j_nu_or_alpha_nu
                  )
{
  if (pol_mode == STOKES_U)
  {
    return 0.;
  }

  double ans = 0.;

  double nu_c = (electron_charge * magnetic_field)
               /(2. * M_PI * electron_mass * speed_light);

  double n_minus = (nu/nu_c) * fabs(sin(observer_angle));

  /* TODO: Describe why the stokes_v_switch is changing */
  int stokes_v_switch = -1;

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
