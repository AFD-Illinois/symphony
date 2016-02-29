#include "symphony.h"

/*j_nu: wrapper for the emissivity calculation; takes in values of all
 *      necessary paramters and sets a struct of parameters using the input
 *      values.  It then passes this struct to n_summation(), which begins
 *      the emissivity calculation.
 *
 *@params: nu, magnetic_field, electron_density, observer_angle,
 *         distribution, polarization, theta_e, power_law_p,
 *         gamma_min, gamma_max, gamma_cutoff, kappa, 
 *         kappa_width
 *
 *@returns: n_summation(&params), which takes the struct of 
 *          parameters (now populated with values) and 
 *          performs the integration to evaluate j_nu().
 */
double j_nu(double nu, 
            double magnetic_field, 
            double electron_density,
            double observer_angle,
            int distribution,
            int polarization,
            double theta_e,
            double power_law_p,
            double gamma_min,
            double gamma_max,
            double gamma_cutoff,
            double kappa,
            double kappa_width
           )
{
/*fill the struct with values*/
  struct parameters params;
  setConstParams(&params);
  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.polarization       = polarization;
  params.mode               = params.EMISSIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

  return n_summation(&params);
}

/*alpha_nu: wrapper for the absorptivity calculation; takes in values of all
 *          necessary paramters and sets a struct of parameters using the input
 *          values.  It then passes this struct to n_summation(), which begins
 *          the absorptivity calculation.
 *
 *@params: nu, magnetic_field, electron_density, observer_angle,
 *         distribution, polarization, theta_e, power_law_p,
 *         gamma_min, gamma_max, gamma_cutoff, kappa, 
 *         kappa_width
 *
 *@returns: n_summation(&params), which takes the struct of 
 *          parameters (now populated with values) and 
 *          performs the integration to evaluate alpha_nu().
 */

double alpha_nu(double nu,
                double magnetic_field,
                double electron_density,
                double observer_angle,
                int distribution,
                int polarization,
                double theta_e,
                double power_law_p,
                double gamma_min,
                double gamma_max,
                double gamma_cutoff,
                double kappa,
                double kappa_width
               )
{
/*fill the struct with values*/
  struct parameters params;
  setConstParams(&params);
  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.polarization       = polarization;
  params.mode               = params.ABSORPTIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

  return n_summation(&params);
}

/*get_nu_c: takes in values of electron_charge, magnetic_field, mass_electron,
 *          and speed_light, and returns the cyclotron frequency nu_c.  
 * 
 *@params:  struct parameters params, contains parameters mentioned above
 *@returns: cyclotron frequency, nu_c, for the provided parameters
 */
double get_nu_c(struct parameters params)
{
  return  (params.electron_charge * params.magnetic_field)
        / (2. * M_PI * params.mass_electron * params.speed_light);
}

/*distribution_function: wrapper for the electron momentum-space distribution
 *                       function.  The value of params->distribution 
 *                       determines which distribution function is called.
 *                       The procedure to add a new gyrotropic distribution
 *                       function is outlined in the README.
 *
 *@params: electron Lorentz factor gamma, struct parameters * params
 *         described above
 *@returns: the selected distribution function, evaluated at the 
 *          input Lorentz factor gamma, for the parameters set in
 *          params.
 */
double distribution_function(double gamma, struct parameters * params)
{
  if(params->distribution == params->MAXWELL_JUETTNER)
  {
    return maxwell_juettner_f(gamma, params);
  }
  else if(params->distribution == params->POWER_LAW)
  {
    return power_law_f(gamma, params);
  }
  else if(params->distribution == params->KAPPA_DIST)
  {
    return kappa_f(gamma, params);
  }

  return 0.;
}


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
double n_peak(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double beta = 0.;

  if (params->nu <= nu_c * params->theta_e * params->theta_e 
      || params->theta_e < 1.) /*beta low nu limit*/
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

/*polarization_term: term in the gamma integrand that varies based upon the 
 *                   Stokes parameter; this term is denoted K_S in [1],
 *                   eq. 3 and 12.  The functional form K_S takes is
 *                   described in eq. 4-7 of [1].
 *
 *@params: Lorentz factor gamma, harmonic number n,
 *         struct of parameters params 
 *@returns: result of evaluating term K_S (the term dependent on Stokes
 *          parameter) in the gamma integrand. 
 */
double polarization_term(double gamma, double n,
                         struct parameters * params
                        ) 
{
  /*below calculation is described in Section 2 of [1]*/
  
  double nu_c = get_nu_c(*params);

  double beta = sqrt(1. - 1./(gamma*gamma));

  /*xi = the angle between the electron velocity vector and the magnetic
         field vector */
  double cos_xi =   (gamma * params->nu - n * nu_c)
		  / ( gamma * params->nu * beta 
                     * cos(params->observer_angle)
                    );

  double M = (cos(params->observer_angle) - beta * cos_xi)
             / sin(params->observer_angle);

  double N = beta * sqrt(1 - (cos_xi*cos_xi));

  double z = (params->nu * gamma * beta * sin(params->observer_angle) 
              * sqrt(1. - cos_xi*cos_xi))/nu_c;

  double K_xx = M*M * pow(my_Bessel_J(n, z), 2.);

  double K_yy = N*N * pow(my_Bessel_dJ(n, z), 2.);

  double ans = 0.;

  if(params->polarization == params->STOKES_I)
  {
    ans = K_xx + K_yy;
  }

  else if(params->polarization == params->STOKES_Q)
  {
    ans = K_xx - K_yy;
  }

  else if(params->polarization == params->STOKES_U)
  {
    ans = 0.;
  }

  else if(params->polarization == params->STOKES_V)
  {
    ans = -2.*M*N*my_Bessel_J(n, z)*my_Bessel_dJ(n, z);
  }

  return ans;
}

/*differential_of_f: The absorptivity integrand ([1] eq. 12) has a term 
 *                   dependent on a differential operator ([1] eq. 13) 
 *                   applied to the distribution function. This function
 *                   calls each individual differential, evaluated 
 *                   analytically, for the 3 distributions studied. 
 *                   These are evaluated analytically for speed, but
 *                   a calculator to evaluate this differential without
 *                   the analytic result is also provided. This calculator
 *                   works with any gyrotropic distribution function.
 * 
 *@params: Lorentz factor gamma, struct of parameters params
 *@returns: differential of the distribution function term in the gamma
 *          integrand.  
 */
//double differential_of_f(double gamma, struct parameters * params) 
//{
//  /*described in Section 2 of [1] */
//
//  double Df = 0.;
//
//  if(params->distribution == params->MAXWELL_JUETTNER)
//  {
//    Df = differential_of_maxwell_juettner(gamma, params);
//  }
//  else if(params->distribution == params->POWER_LAW)
//  {
//    Df = differential_of_power_law(gamma, params);
//  }
//  else if(params->distribution == params->KAPPA_DIST)
//  {
//    Df = differential_of_kappa(gamma, params);
//  }
//
//  return Df;
//}

double num_differential_of_f(double gamma, struct parameters * params)
{
  double Df = 0.;
  double epsilon = 0.0000001;
  double prefactor = (2. * params->pi * params->nu
                      / (params->mass_electron
                         *params->speed_light*params->speed_light))
                     * params->pi * pow(params->mass_electron, 3.) 
                     * pow(params->speed_light, 3);

  if(params->distribution == params->MAXWELL_JUETTNER)
  {
    Df =  (maxwell_juettner_f(gamma+epsilon, params)
           - maxwell_juettner_f(gamma-epsilon, params))
          / epsilon;
  }
  else if(params->distribution == params->POWER_LAW)
  {
    Df =  (power_law_f(gamma+epsilon, params)
           - power_law_f(gamma-epsilon, params))
          / epsilon;
  }
  else if(params->distribution == params->KAPPA_DIST)
  {
    Df =  (kappa_f(gamma+epsilon, params)
           - kappa_f(gamma-epsilon, params))
          / epsilon;
  }

  return prefactor * Df;
}

/*gamma_integrand: full gamma integrand (eq. 3, 12 of [1]) with the summation
 *                 moved outside the integral and the delta function evaluated
 *                 by setting eq. 10 of [1] to zero. Also, d^3p is converted
 *                 to coordinates gamma (Lorentz factor), xi (angle between 
 *                 electron velocity vector and magnetic field vector), 
 *                 and phi (gyrophase).
 *
 *@params: Lorentz factor gamma, 
 *         void pointer to paramsGSLInput (struct of params similar to
 *         parameters, except also contains harmonic number n.  This is 
 *         necessary because the gamma integral is moved inside the n
 *         summation in order to get a smooth integrand.)
 *returns: gamma integrand, which still needs to be integrated over
 *         gamma and then summed over n.
 */
double gamma_integrand(double gamma, void * paramsGSLInput)
{
  struct parametersGSL * paramsGSL = (struct parametersGSL*) paramsGSLInput;
  struct parameters * params       = &(paramsGSL->params);

  double beta = sqrt(1. - 1./(gamma*gamma));

  double func_I = 
      (2. * params->pi * pow(params->electron_charge * params->nu, 2.) )
    / params->speed_light * (  pow(params->mass_electron * params->speed_light, 3.) 
                             * gamma*gamma * beta * 2. * params->pi
                            )
    * distribution_function(gamma, params) 
    * polarization_term(gamma, paramsGSL->n, params);

  double ans = 0.;

  if (params->mode == params->EMISSIVITY)
  {
    double prefactor = 1./(params->nu * beta * fabs(cos(params->observer_angle)));

    ans = prefactor * func_I;
  }
  else if (params->mode == params->ABSORPTIVITY)
  {
    double prefactor = -  params->speed_light
                        * params->electron_charge
                        * params->electron_charge 
                        / (2. * params->nu);

    ans = prefactor*gamma*gamma*beta*num_differential_of_f(gamma, params)
                   *polarization_term(gamma, paramsGSL->n, params)
                   *(1./(params->nu*beta*fabs(cos(params->observer_angle))));
  }

  return ans;
}

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
double n_integration(double n_minus, struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  /*For the MAXWELL_JUETTNER distribution, the n-space peak location is known
    analytically; this speeds up evaluation of MAXWELL_JUETTNER until
    nu/nu_c = 1e6, where the approximation breaks down and the adaptive 
    procedure works better.*/
  if (params->distribution == params->MAXWELL_JUETTNER && params->nu/nu_c < 1e6) 
  {
    double n_start = (int)(params->n_max + n_minus + 1.);
    double ans = n_integral(n_start, params->C * n_peak(params), params);
    return ans;
  }
  else /*For other distributions, the n-space peak is found adaptively. */
  {
    double n_start = (int)(params->n_max + n_minus + 1.);
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
  ans += n_integration(n_minus, params);

  /*if doing Stokes V, must perform two n integrals: one over the positive
    part of gamma integrand and one over the negative part; this helps GSL QAG
    resolve the n integrals*/
  if(params->polarization == params->STOKES_V)
  {
    params->stokes_v_switch = 1;
    ans += n_integration(n_minus, params);
  }
  
  return ans;
}
