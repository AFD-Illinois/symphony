#include "maxwell_juettner.h"

// TODO: Change names thermal -> maxwell_juettner
double thermal_I(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_s = (2./9.)*nu_c*sin(params->observer_angle)*params->theta_e*params->theta_e;
  double X = params->nu/nu_s;
  double prefactor = (params->electron_density * pow(params->electron_charge, 2.) * nu_c)/params->speed_light;
  double term1 = sqrt(2.)*params->pi/27. * sin(params->observer_angle);
  double term2 = pow(pow(X, 0.5)+pow(2., 11./12.)*pow(X, 1./6.), 2.);
  double term3 = exp(-pow(X, 1./3.));
  double ans = prefactor * term1 * term2 * term3;

  return ans;
}

double thermal_Q(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_s = (2./9.)*nu_c*sin(params->observer_angle)*params->theta_e*params->theta_e;
  double X = params->nu/nu_s;
  double prefactor = (params->electron_density * pow(params->electron_charge, 2.) * nu_c)/params->speed_light;
  double term1 = sqrt(2.)*params->pi/27. * sin(params->observer_angle);
  double term2 = (7.*pow(params->theta_e, 24./25.)+35.)
		/(10.*pow(params->theta_e, 24./25.)+75.);
  double term3 = pow(pow(X, 0.5)+term2*pow(2., 11./12.)*pow(X, 1./6.), 2.);
  double ans = prefactor*term1*term3*exp(-pow(X, 1./3.));

  return -ans;
}

double thermal_V(struct parameters * params)
{
  double nu_c = get_nu_c(*params);

  double nu_s = (2./9.)*nu_c*sin(params->observer_angle)*params->theta_e*params->theta_e;
  double X = params->nu/nu_s;
  double prefactor = (params->electron_density * pow(params->electron_charge, 2.) * nu_c)/params->speed_light;
  double term1 = (37.-87.*sin(params->observer_angle-28./25.))/(100.*(params->theta_e+1.));
  double term2 = pow(1.+(pow(params->theta_e, 3./5.)/25.+7./10.)
		*pow(X, 9./25.), 5./3.);
  double ans = prefactor*term1*term2*exp(-pow(X, 1./3.));

  return -ans;
}

double planck_func(struct parameters * params)
{
  double term1 = (2.*params->plancks_constant*pow(params->nu, 3.))/pow(params->speed_light, 2.);
  double term2 = (exp(params->plancks_constant*params->nu
                  /(params->theta_e*params->mass_electron*pow(params->speed_light, 2.)))-1.);
  double ans = term1 / term2;
  return ans;
}

double thermal_I_abs(struct parameters * params)
{
  double ans = thermal_I(params)/planck_func(params);
  return ans;
}

double thermal_Q_abs(struct parameters * params)
{
  double ans = thermal_Q(params)/planck_func(params);
  return ans;
}

double thermal_V_abs(struct parameters * params)
{
  double ans = thermal_V(params)/planck_func(params);
  return ans;
}
