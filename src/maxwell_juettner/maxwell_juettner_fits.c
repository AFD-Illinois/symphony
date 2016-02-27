#include "maxwell_juettner.h"

/*maxwell_juettner_I: fitting formula for the emissivity (polarized in Stokes I)
 *                    produced by a Maxwell-Juettner (relativistic thermal) 
 *                    distribution of electrons. (Eq. 29, 31 of [1])
 *
 *@params: struct of parameters params
 *@returns: fit to the emissivity, polarized in Stokes I, for the given 
 *          parameters for a Maxwell-Juettner distribution.
 */
double maxwell_juettner_I(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_s = (2./9.)*nu_c*sin(params->observer_angle)*params->theta_e
                *params->theta_e;

  double X = params->nu/nu_s;

  double prefactor = (params->electron_density 
                      * pow(params->electron_charge, 2.) 
                      * nu_c)/params->speed_light;

  double term1 = sqrt(2.)*params->pi/27. * sin(params->observer_angle);

  double term2 = pow(pow(X, 0.5)+pow(2., 11./12.)*pow(X, 1./6.), 2.);
  
  double term3 = exp(-pow(X, 1./3.));
  
  double ans = prefactor * term1 * term2 * term3;

  return ans;
}

/*maxwell_juettner_Q: fitting formula for the emissivity (polarized in Stokes Q)
 *                    produced by a Maxwell-Juettner (relativistic thermal) 
 *                    distribution of electrons. (Eq. 29, 31 of [1])
 *
 *@params: struct of parameters params
 *@returns: fit to the emissivity, polarized in Stokes Q, for the given 
 *          parameters for a Maxwell-Juettner distribution.
 */
double maxwell_juettner_Q(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_s = (2./9.)*nu_c*sin(params->observer_angle)
                *params->theta_e*params->theta_e;

  double X = params->nu/nu_s;

  double prefactor = (params->electron_density 
                      * pow(params->electron_charge, 2.) 
                      * nu_c)/params->speed_light;

  double term1 = sqrt(2.)*params->pi/27. * sin(params->observer_angle);

  double term2 = (7.*pow(params->theta_e, 24./25.)+35.)
		/(10.*pow(params->theta_e, 24./25.)+75.);

  double term3 = pow(pow(X, 0.5)+term2*pow(2., 11./12.)*pow(X, 1./6.), 2.);

  double ans = prefactor*term1*term3*exp(-pow(X, 1./3.));

  return -ans;
}

/*maxwell_juettner_V: fitting formula for the emissivity (polarized in Stokes V)
 *                    produced by a Maxwell-Juettner (relativistic thermal) 
 *                    distribution of electrons. (Eq. 29, 31 of [1])
 *
 *@params: struct of parameters params
 *@returns: fit to the emissivity, polarized in Stokes V, for the given 
 *          parameters for a Maxwell-Juettner distribution.
 */
double maxwell_juettner_V(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_s = (2./9.)*nu_c*sin(params->observer_angle)*params->theta_e
                *params->theta_e;

  double X = params->nu/nu_s;

  double prefactor = (params->electron_density 
                      * pow(params->electron_charge, 2.) 
                      * nu_c)/params->speed_light;

  double term1 = (37.-87.*sin(params->observer_angle-28./25.))
                /(100.*(params->theta_e+1.));

  double term2 = pow(1.+(pow(params->theta_e, 3./5.)/25.+7./10.)
		*pow(X, 9./25.), 5./3.);

  double ans = prefactor*term1*term2*exp(-pow(X, 1./3.));

  return -ans;
}

/*planck_func: The Planck function (used in eq. 25 of [1]) can be used to
 *             obtain alpha_nu() fitting formulae from the j_nu() fitting 
 *             formulae for the Maxwell-Juettner (relativistic thermal)
 *             distribution.
 *
 *@params: struct of parameters params
 *@returns: Planck function evaluated for the supplied parameters
 */
double planck_func(struct parameters * params)
{
  double term1 = (2.*params->plancks_constant*pow(params->nu, 3.))
                /pow(params->speed_light, 2.);

  double term2 = (exp(params->plancks_constant*params->nu
                  /(params->theta_e*params->mass_electron
                    *pow(params->speed_light, 2.)))-1.);

  double ans = term1 / term2;

  return ans;
}

/*maxwell_juettner_I_abs: Fitting formula for the absorptivity, polarized in
 *                        Stokes I, for a Maxwell-Juettner electron momentum 
 *                        distribution.  Uses eq. 30, 31, 32 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to the absorptivity, in Stokes I, for a Maxwell-
 *          Juettner distribution of electrons.
 */
double maxwell_juettner_I_abs(struct parameters * params)
{
  double ans = maxwell_juettner_I(params)/planck_func(params);
  return ans;
}

/*maxwell_juettner_Q_abs: Fitting formula for the absorptivity, polarized in
 *                        Stokes Q, for a Maxwell-Juettner electron momentum 
 *                        distribution.  Uses eq. 30, 31, 32 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to the absorptivity, in Stokes Q, for a Maxwell-
 *          Juettner distribution of electrons.
 */
double maxwell_juettner_Q_abs(struct parameters * params)
{
  double ans = maxwell_juettner_Q(params)/planck_func(params);
  return ans;
}

/*maxwell_juettner_V_abs: Fitting formula for the absorptivity, polarized in
 *                        Stokes V, for a Maxwell-Juettner electron momentum 
 *                        distribution.  Uses eq. 30, 31, 32 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to the absorptivity, in Stokes V, for a Maxwell-
 *          Juettner distribution of electrons.
 */
double maxwell_juettner_V_abs(struct parameters * params)
{
  double ans = maxwell_juettner_V(params)/planck_func(params);
  return ans;
}
