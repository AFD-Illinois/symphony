#include "integrands.h"

/*polarization_term: term in the gamma integrand that varies based upon the 
 *                   Stokes parameter; this term is denoted K_S in [1],
 *                   eq. 3 and 12.  The functional form K_S takes is
 *                   described in eq. 4-7 of [1].  This term is a part
 *                   of the gamma integrand, and thus is called by
 *                   gamma_integrand() below. 
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

  double ans = 0.;

  if (params->mode == params->EMISSIVITY)
  {
    double func_I =
      (2. * params->pi * pow(params->electron_charge * params->nu, 2.) )
    / params->speed_light * (  pow(params->mass_electron * params->speed_light, 3.)
                             * gamma*gamma * beta * 2. * params->pi
                            )
    * params->distribution_function(gamma, params)
    * polarization_term(gamma, paramsGSL->n, params);

    double prefactor = 1./(params->nu * beta * fabs(cos(params->observer_angle)));

    ans = prefactor * func_I;
  }
  else if (params->mode == params->ABSORPTIVITY)
  {
    double prefactor = -  params->speed_light
                        * params->electron_charge
                        * params->electron_charge 
                        / (2. * params->nu);

    ans = prefactor * gamma * gamma * beta
         * numerical_differential_of_f(gamma,  params)
         * polarization_term(gamma, paramsGSL->n, params)
         * (1./(params->nu*beta*fabs(cos(params->observer_angle))));

  }

  return ans;

}

void set_distribution_function(struct parameters * params)
{
  if(params->distribution == params->MAXWELL_JUETTNER)
  {
    params->distribution_function = &maxwell_juettner_f;
    params->use_n_peak            = 1;
    params->n_peak                = &maxwell_juettner_n_peak;
    params->analytic_differential = &differential_of_maxwell_juettner;
  }
  else if(params->distribution == params->POWER_LAW)
  {
    params->distribution_function = &power_law_f;
    params->use_n_peak            = 0;
    params->analytic_differential = &differential_of_power_law;
  }
  else if(params->distribution == params->KAPPA_DIST)
  {
    params->distribution_function = &kappa_f;
    params->use_n_peak            = 0;
    params->analytic_differential = differential_of_kappa;
  }
}
