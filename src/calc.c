#include "symphony.h"

double j_nu(double nu, 
            double magnetic_field, 
            double electron_density,
            double observer_angle,
            int distribution,
            int polarization
           )
{
//fill the struct with values
  struct parameters params;
  setConstParams(&params);
  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.polarization       = polarization;
  params.mode               = params.EMISSIVITY;
  params.theta_e            = 10.; //TODO: PASS THESE IN AS PARAMETERS
  params.power_law_p        = 3.5;
  params.gamma_min          = 1.;
  params.gamma_max          = 100.;
  params.gamma_cutoff       = 100000000000.;

  //printf("\n POWER_LAW:  %e \n", power_law_f(1.5, &params));

  return n_summation(&params);
}

double alpha_nu(double nu, 
                double magnetic_field, 
                double electron_density,
                double observer_angle,
                int distribution,
                int polarization
               )
{
//fill the struct with values
  struct parameters params;
  setConstParams(&params);
  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.polarization       = polarization;
  params.mode               = params.ABSORPTIVITY;

  return n_summation(&params);
}

double get_nu_c(struct parameters params)
{
  return  (params.electron_charge * params.magnetic_field)
        / (2. * M_PI * params.mass_electron * params.speed_light);
}


double distribution_function(double gamma, struct parameters * params)
{
  if(params->distribution == params->THERMAL)
  {
    return maxwell_juttner_f(gamma, params);
  }
//  else if(params->distribution == params->POWER_LAW)
//  {
//    return power_law_f(gamma, params);
//  }
//  else if(params->distribution == params->KAPPA_DIST)
//  {
//    return kappa_f(gamma, params);
//  }

  return 0.;

}


/////*n_peak: gives the location of the peak of the n integrand for 
//// *the THERMAL distribution; uses Eq. 68 in [1]
//// *
//// *@param nu: Input, frequency of emission/absorption
//// *@returns n_peak: Output, location of integrand's peak for the 
//// * n-integral for the THERMAL distribution 
//// */
double n_peak(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double beta = 0.;

  if (params->nu <= nu_c * params->theta_e*params->theta_e || params->theta_e < 1.) 
  {
    beta = sqrt((1. - 1./pow((1. + params->theta_e),2.))); //beta in low nu limit
  }
  else 
  {
    beta = sqrt(1. - pow((2. * params->theta_e * params->nu / nu_c), -2./3.));//beta for high nu
  }
  
  double n_peak =  (params->theta_e + 1. + pow((2. * params->theta_e * params->nu / nu_c),1./3.))
                 * (params->nu/nu_c) * (1. - beta*beta * pow(cos(params->observer_angle),2.));

  return n_peak;
}

/////*polarization_term: term in the gamma integrand that varies based upon the 
//// *polarization mode; uses eq. 13-19, 22-25, 27-30, in [1]
//// *
//// *@param gamma: Input, Lorentz factor 
//// *@param n: Input, index n of sum (synchrotron harmonic number)
//// *@param nu: Input, frequency of emission/absorption
//// *@returns: piece of gamma integrand that determines polarization mode 
//// */
double polarization_term(double gamma, double n,
                         struct parameters * params
                        ) 
{
  /*below calculation is described in Section 3.1 of [1]*/
  
  double nu_c = get_nu_c(*params);

  double beta = sqrt(1. - 1./(gamma*gamma));

  double cos_xi =   (gamma * params->nu - n * nu_c)
		     / (gamma * params->nu * beta 
                        * cos(params->observer_angle));

  double M = (cos(params->observer_angle) - beta * cos_xi)
             /sin(params->observer_angle);

  double N = beta * sqrt(1 - (cos_xi*cos_xi));

  double z = (params->nu * gamma * beta * sin(params->observer_angle) 
       * sqrt(1. - cos_xi*cos_xi))/nu_c;

  double K_xx = M*M * pow(my_Bessel_J(n, z), 2.);

  double K_yy = N*N * pow(my_Bessel_dJ(n, z), 2.);

//  double T_x = (2.*params->nu*cos(params->observer_angle))
//               /(nu_c*pow(sin(params->observer_angle), 2.)
//               +sqrt(nu_c*nu_c*pow(sin(params->observer_angle), 4.)
//               +4.*params->nu*params->nu*pow(cos(params->observer_angle), 2.)));

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

/////*differential_of_f: term in gamma integrand only for absorptivity calculation; 
//// *it is the differential Df = 2\pi\nu (1/(mc)*d/dgamma + (beta cos(theta)
//// * -cos(xi))/(p*beta*c) * d/d(cos(xi))) f 
//// *this is eq. 41 of [1]
//// *below it is applied for the THERMAL, POWER_LAW, and KAPPA_DIST distributions
//// *@param gamma: Input, Lorentz factor
//// *@param nu: Input, frequency of emission/absorption
//// *@returns: Output, Df term in gamma integrand; depends on distribution function
//// */
double differential_of_f(double gamma, struct parameters params) 
{
  /*described in Section 3 of [1] */

  double f = 0.;

  if(params.distribution == params.THERMAL)
  {
    double prefactor = (params.pi * params.nu 
                     / (params.mass_electron*params.speed_light*params.speed_light)) 
		     * (params.electron_density/(params.theta_e * gsl_sf_bessel_Kn(2, 1./params.theta_e)));

    double body = (-1./params.theta_e) * exp(-gamma/params.theta_e);

    f = prefactor * body;
  }

  else if(params.distribution == params.POWER_LAW)
  {
    double pl_norm = 4.* M_PI/(normalize_f(&params));

    double prefactor = (M_PI * params.nu / (params.mass_electron*params.speed_light*params.speed_light)) 
                     * (params.electron_density*(params.power_law_p-1.))
                     /((pow(params.gamma_min, 1.-params.power_law_p) 
                     - pow(params.gamma_max, 1.-params.power_law_p)));

    double term1 = ((-params.power_law_p-1.)*exp(-gamma/params.gamma_cutoff)
               *pow(gamma,-params.power_law_p-2.)/(sqrt(gamma*gamma - 1.)));

    double term2 = (exp(-gamma/params.gamma_cutoff) * pow(gamma,(-params.power_law_p-1.))
                  /(params.gamma_cutoff * sqrt(gamma*gamma - 1.)));

    double term3 = (exp(-gamma/params.gamma_cutoff) * pow(gamma,-params.power_law_p))
               /pow((gamma*gamma - 1.), (3./2.));

    f = pl_norm * prefactor * (term1 - term2 - term3);
  }

  else if(params.distribution == params.KAPPA_DIST)
  {
    double prefactor = params.electron_density * (1./normalize_f(&params)) 
                       * 4. * M_PI*M_PI * params.nu * params.mass_electron
                       * params.mass_electron * params.speed_light;

    double term1 = ((- params.kappa - 1.) / (params.kappa * params.kappa_width)) 
		  * pow((1. + (gamma - 1.)/(params.kappa * params.kappa_width)), -params.kappa-2.);

    double term2 = pow((1. + (gamma - 1.)/(params.kappa * params.kappa_width)),(- params.kappa - 1.)) 
                    * (- 1./params.gamma_cutoff);

    f = prefactor * (term1 + term2) * exp(-gamma/params.gamma_cutoff);
  }

    return f;
}

/////*maxwell_juttner_f: Maxwell-Juttner distribution function in terms of Lorentz 
//// * factor gamma; uses eq. 47, 49, 50 of [1]
//// *@param gamma: Input, Lorentz factor
//// *@returns: THERMAL distribution function, which goes into the gamma integrand 
//// */
double maxwell_juttner_f(double gamma, struct parameters * params) 
{
  double beta = sqrt(1. - 1./(gamma*gamma));

  double d = (params->electron_density * gamma * sqrt(gamma*gamma-1.) 
              * exp(-gamma/params->theta_e)) / (4. * params->pi 
                    * params->theta_e 
                    * gsl_sf_bessel_Kn(2, 1./params->theta_e));

  double ans = 1./(pow(params->mass_electron, 3.) * pow(params->speed_light, 3.) 
                   * gamma*gamma * beta) * d;

  return ans;
}

/////*power_law_to_be_normalized: the power-law distribution is normalized as-is, 
//// * but we have added an exponential cutoff e^(-gamma/gamma_cutoff), so it must 
//// * be normalized again.  The normalization constant is given by 1 over
//// * the integral of power_law_to_be_normalized from 1 to infinity.
//// * uses eq. 50 of [1]
//// *@param gamma: Input, Lorentz factor
//// *@param *params: void pointer to struct parameters, which contains
//// * n (harmonic number) and nu (frequency of emission/absorption)
//// *@returns: Output, power-law distribution to be normalized via normalize_f()
//// */
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
//double body = pow(gamma, -power_law_p); for PL w/o cutoff

  double ans = norm_term * prefactor * body;

  return ans;
}

////*power_law_f: power-law distribution function, normalized via call to the 
//// * normalize_f() function. Uses eq. 50 of [1].
//// *@param gamma: Input, Lorentz factor
//// *@returns: Ouput, a normalized power-law distribution for the gamma integrand
//// */
double power_law_f(double gamma, struct parameters * params) 
{
  double beta = sqrt(1. - 1./(gamma*gamma));

  double prefactor = params->electron_density * (params->power_law_p - 1.) 
                     / (pow(params->gamma_min, 1. - params->power_law_p) 
                        - pow(params->gamma_max, 1. - params->power_law_p));

  double body = pow(gamma, -params->power_law_p) 
                * exp(- gamma / params->gamma_cutoff);
//double body = pow(gamma, -params->power_law_p); //for PL w/o cutoff

  double ans = 1./normalize_f(params) * prefactor * body 
                              * 1./(pow(params->mass_electron, 3.) 
                              * pow(params->speed_light, 3.) 
                              * gamma*gamma * beta);

  return ans;
}

/////*kappa_to_be_normalized: the kappa distribution is not normalized, so we must 
//// * find the normalization constant by taking 1 over the integral of the function
//// * kappa_to_be_normalized from 1 to infinity.  Uses eq. 42 of [2].
//// *@param gamma: Input, Lorentz factor
//// *@param *params: void pointer to struct parameters, which contains
//// * n (harmonic number) and nu (frequency of emission/absorption)
//// *@returns: kappa distribution function to be normalized via normalize_f()
//// */
double kappa_to_be_normalized(double gamma, void * paramsInput)
{
  struct parameters params = *(struct parameters*) paramsInput;

  double kappa_body = pow((1. + (gamma - 1.)/(params.kappa * params.kappa_width)), -params.kappa-1);

  double cutoff = exp(-gamma/params.gamma_cutoff);

  double norm_term = 4. * M_PI * pow(params.mass_electron, 3.) * pow(params.speed_light, 3.) 
                   * gamma * sqrt(gamma*gamma-1.);

  double ans = kappa_body * cutoff * norm_term;
  //double ans = kappa_body * norm_term; //for kappa w/o cutoff

  return ans;
}

/////*kappa_f: kappa distribution function, numerically normalized
//// *uses eq. 42 of [2]
//// *@param gamma: Input, Lorentz factor
//// *@returns: normalized kappa distribution function to go into gamma integrand 
//// */
double kappa_f(double gamma, struct parameters * params)
{
  double norm = 1./normalize_f(params);

  double kappa_body = params->electron_density * pow((1. + (gamma - 1.)
                     /(params->kappa * params->kappa_width)), -params->kappa-1);

  double cutoff = exp(-gamma/params->gamma_cutoff);

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
double gamma_integrand(double gamma, void * paramsGSLInput)
{
  struct parametersGSL * paramsGSL = (struct parametersGSL*) paramsGSLInput;
  struct parameters * params       = &(paramsGSL->params);

  double nu_c = get_nu_c(*params);

  double beta = sqrt(1. - 1./(gamma*gamma));

//inserting integrand_without_extra_factor here

  double cos_xi = (gamma * params->nu - paramsGSL->n * nu_c)
                 /(gamma * params->nu * beta * cos(params->observer_angle));

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

    ans = prefactor*gamma*gamma*beta*differential_of_f(gamma, *params)
                *polarization_term(gamma, paramsGSL->n, params)
                *(1./(params->nu*beta*fabs(cos(params->observer_angle))));
  }

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
double gamma_integration_result(double n, void * paramsInput)
{
  /*This evaluates the gamma integral*/
  struct parameters * params = (struct parameters*) paramsInput;

  double nu_c = get_nu_c(*params);

  /* Eqn () in Leung et. al. */
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
   // double neg_result = gsl_integrate(gamma_minus_high, gamma_peak, n, params);
   // double pos_result = gsl_integrate(gamma_peak, gamma_plus_high,  n, params);
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
    //result = gsl_integrate(gamma_minus_high, gamma_plus_high, n, params);
    result = gamma_integral(gamma_minus_high, gamma_plus_high, n, params);
  }

  /*GSL QAG sometimes erroneously gives NaN instead of small values; 
    return 0 instead */
  if(isnan(result) != 0) result = 0.;

  return result;
}

/////*n_integration: j_nu and alpha_nu are given by an integral over gamma of
//// * an integrand that contains a sum over n; we do the integral over gamma
//// * and then the sum over n.  For numerical accuracy and speed, we only sum up
//// * to n = n_max = 30, and integrate from n_max to the point where our integral
//// * stops giving appreciable contributions.  n_integration() does the n integral.
//// *@params n_minus: Input, minimum n for which the integrand is real
//// *@params nu: Input, frequency of absorption/emission
//// *@returns: Output, result of integration over n
//// */
double n_integration(double n_minus, struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  if (params->distribution == params->THERMAL && params->nu/nu_c < 1e6) 
  {
    double n_start = (int)(params->n_max + n_minus + 1.);
    //double ans = gsl_integrate(n_start, params->C * n_peak(params), -1, params);
    double ans = n_integral(n_start, params->C * n_peak(params), -1, params);
    return ans;
  }

  else 
  {
    double n_start = (int)(params->n_max + n_minus + 1.);
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
      double deriv = derivative(n_start, params->nu);
      if(fabs(deriv) < deriv_tol) delta_n = 100. * delta_n;

      //contrib = gsl_integrate(n_start, (n_start + delta_n), -1, params);
      contrib = n_integral(n_start, params->C * n_peak(params), -1, params);
      ans = ans + contrib;

      n_start = n_start + delta_n;
    }

  return ans;
  }
}

///*n_summation: performs the sum over n (harmonic number) from n = 1 to value 
// * n_max = 30, set heuristically in [1].  Past n_max we integrate over n using
// * n_integration().
// *@param nu: Input, frequency of absorption/emission
// *@returns: j_nu or alpha_nu, the emissivity or absorptivity, respectively
// */
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
  params->stokes_v_switch = -1; // TODO: describe: 0: , 1: 

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
