#include "dec.h"

double thermal_I(double nu)
{
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);

  double nu_s = (2./9.)*nu_c*sin(obs_angle)*theta_e*theta_e;
  double X = nu/nu_s;
  double prefactor = (n_e * pow(electron_charge, 2.) * nu_c)/speed_light;
  double term1 = sqrt(2.)*M_PI/27. * sin(obs_angle);
  double term2 = pow(pow(X, 0.5)+pow(2., 11./12.)*pow(X, 1./6.), 2.);
  double term3 = exp(-pow(X, 1./3.));
  double ans = prefactor * term1 * term2 * term3;

  return ans;
}

double thermal_Q(double nu)
{
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);

  double nu_s = (2./9.)*nu_c*sin(obs_angle)*theta_e*theta_e;
  double X = nu/nu_s;
  double prefactor = (n_e * pow(electron_charge, 2.) * nu_c)/speed_light;
  double term1 = sqrt(2.)*M_PI/27. * sin(obs_angle);
  double term2 = (7.*pow(theta_e, 24./25.)+35.)
		/(10.*pow(theta_e, 24./25.)+75.);
  double term3 = pow(pow(X, 0.5)+term2*pow(2., 11./12.)*pow(X, 1./6.), 2.);
  double ans = prefactor*term1*term3*exp(-pow(X, 1./3.));

  return -ans;
}

double thermal_V(double nu)
{
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);

  double nu_s = (2./9.)*nu_c*sin(obs_angle)*theta_e*theta_e;
  double X = nu/nu_s;
  double prefactor = (n_e * pow(electron_charge, 2.) * nu_c)/speed_light;
  double term1 = (37.-87.*sin(obs_angle-28./25.))/(100.*(theta_e+1.));
  double term2 = pow(1.+(pow(theta_e, 3./5.)/25.+7./10.)
		*pow(X, 9./25.), 5./3.);
  double ans = prefactor*term1*term2*exp(-pow(X, 1./3.));

  return -ans;
}

double planck_func(double nu)
{
  double term1 = (2.*h*pow(nu, 3.))/pow(speed_light, 2.);
  double term2 = (exp(h*nu/(theta_e*mass_electron*pow(speed_light, 2.)))-1.);
  double ans = term1 / term2;
  return ans;
}

double thermal_I_abs(double nu)
{
  double ans = thermal_I(nu)/planck_func(nu);
  return ans;
}

double thermal_Q_abs(double nu)
{
  double ans = thermal_Q(nu)/planck_func(nu);
  return ans;
}

double thermal_V_abs(double nu)
{
  double ans = thermal_V(nu)/planck_func(nu);
  return ans;
}

double kappa_I(double nu)
{
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);

  double nu_w = pow(kappa_width*kappa, 2.)*nu_c*sin(obs_angle);
  double X_k = nu/nu_w;
  double prefactor = (n_e*pow(electron_charge, 2.)*nu_c*sin(obs_angle))
	            /speed_light;
  double Nlow = 4.*M_PI*tgamma(kappa-4./3.)/(pow(3., 7./3.)*tgamma(kappa-2.));
  double Nhigh = (1./4.)*pow(3., (kappa-1.)/2.)*(kappa-2.)*(kappa-1.)
		*tgamma(kappa/4.-1./3.)*tgamma(kappa/4.+4./3.);
  double x = 3.*pow(kappa, -3./2.);
  double ans = prefactor*Nlow*pow(X_k, 1./3.)*pow(1.+pow(X_k, x*(3.*kappa-4.)
	      /6.)*pow(Nlow/Nhigh, x), -1./x);

  return ans;

}

double kappa_Q(double nu)
{
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);

  double nu_w = pow(kappa_width*kappa, 2.)*nu_c*sin(obs_angle);
  double X_k = nu/nu_w;
  double prefactor = (n_e*pow(electron_charge, 2.)*nu_c*sin(obs_angle))
	            /speed_light;
  double Nlow = -(1./2.)*4.*M_PI*tgamma(kappa-4./3.)/(pow(3., 7./3.)
	        *tgamma(kappa-2.));
  double Nhigh = -(pow(4./5., 2)+kappa/50.)*(1./4.)*pow(3., (kappa-1.)/2.)*
                  (kappa-2.)*(kappa-1.)*tgamma(kappa/4.-1./3.)
                 *tgamma(kappa/4.+4./3.);
  double x = (37./10.)*pow(kappa, -8./5.);
  double ans = prefactor*Nlow*pow(X_k, 1./3.)
              *pow(1.+pow(X_k, x*(3.*kappa-4.)/6.)*pow(Nlow/Nhigh, x), -1./x);

  return ans;

}

double kappa_V(double nu)
{
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);

  double nu_w = pow(kappa_width*kappa, 2.)*nu_c*sin(obs_angle);
  double X_k = nu/nu_w;
  double prefactor = (n_e*pow(electron_charge, 2.)*nu_c*sin(obs_angle))
                    /speed_light;
  double Nlow = -pow(3./4., 2.)*pow(pow(sin(obs_angle), -12./5.)-1., 12./25.)
                *(pow(kappa, -66./125.)/kappa_width)*pow(X_k, -7./20.)*4.*M_PI
                *tgamma(kappa-4./3.)/(pow(3., 7./3.)*tgamma(kappa-2.));
  double Nhigh = -pow(7./8., 2.)*pow(pow(sin(obs_angle), -5./2.)-1., 11./25.)
                 *(pow(kappa, -11./25.)/kappa_width)*pow(X_k, -1./2.)
                 *(1./4.)*pow(3., (kappa-1.)/2.)*(kappa-2.)*(kappa-1.)
                 *tgamma(kappa/4.-1./3.)*tgamma(kappa/4.+4./3.);
  double x = 3.*pow(kappa, -3./2.);
  double ans = prefactor*Nlow*pow(X_k, 1./3.)
              *pow(1.+pow(X_k, x*(3.*kappa-4.)/6.)*pow(Nlow/Nhigh, x), -1./x);

  return ans;

}

double kappa_I_abs(double nu)
{
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);

  double nu_w = pow(kappa_width*kappa, 2.)*nu_c*sin(obs_angle);
  double X_k = nu/nu_w;
  double prefactor = n_e*electron_charge/(B_field*sin(obs_angle));
  double a = kappa - 1./3.; 
  double b = kappa + 1.;
  double c = kappa + 2./3.;
  double z = -kappa*kappa_width;
  double hyp2f1 = pow(1.-z, -a)*tgamma(c)*tgamma(b-a)/(tgamma(b)*tgamma(c-a))
                 *gsl_sf_hyperg_2F1(a, c-b, a-b+1., 1./(1.-z))+pow(1.-z, -b)
                 *tgamma(c)*tgamma(a-b)/(tgamma(a)*tgamma(c-b))
                 *gsl_sf_hyperg_2F1(b, c-a, b-a+1., 1./(1.-z));
  double Nlow = pow(3., 1./6.)*(10./41.)*pow(2.*M_PI, 2.)
               /pow(kappa_width*kappa, 16./3.-kappa)*(kappa-2.)*(kappa-1.)
               *kappa/(3.*kappa-1.)*tgamma(5./3.)*hyp2f1;
  double Nhigh = 2.*pow(M_PI, 5./2.)/3.*(kappa-2.)*(kappa-1.)*kappa
                 /pow(kappa_width*kappa, 5.)*(2*tgamma(2.+kappa/2.)
                 /(2.+kappa)-1.)*(pow(3./kappa, 19./4.)+3./5.);
  double x = pow(-7./4. + 8.*kappa/5., -43./50.);
  double ans = prefactor*Nlow*pow(X_k, -5./3.)
              *pow(1.+pow(X_k, x*(3.*kappa-1.)/6.)*pow(Nlow/Nhigh, x), -1./x);

  return ans;
}

double kappa_Q_abs(double nu)
{
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);

  double nu_w = pow(kappa_width*kappa, 2.)*nu_c*sin(obs_angle);
  double X_k = nu/nu_w;
  double prefactor = n_e*electron_charge/(B_field*sin(obs_angle));
  double a = kappa - 1./3.; 
  double b = kappa + 1.;
  double c = kappa + 2./3.;
  double z = -kappa*kappa_width;
  double hyp2f1 = pow(1.-z, -a)*tgamma(c)*tgamma(b-a)/(tgamma(b)*tgamma(c-a))
                 *gsl_sf_hyperg_2F1(a, c-b, a-b+1., 1./(1.-z))+pow(1.-z, -b)
                 *tgamma(c)*tgamma(a-b)/(tgamma(a)*tgamma(c-b))
                 *gsl_sf_hyperg_2F1(b, c-a, b-a+1., 1./(1.-z));
  double Nlow = -(25./48.)*pow(3., 1./6.)*(10./41.)*pow(2.*M_PI, 2.)
               /pow(kappa_width*kappa, 16./3.-kappa)*(kappa-2.)*(kappa-1.)
               *kappa/(3.*kappa-1.)*tgamma(5./3.)*hyp2f1;
  double Nhigh = -(pow(21., 2.)*pow(kappa, -144./25.)+11./20.)*2.*pow(M_PI, 5./2.)/3.*(kappa-2.)*(kappa-1.)*kappa
                 /pow(kappa_width*kappa, 5.)*(2*tgamma(2.+kappa/2.)
                 /(2.+kappa)-1.);
  double x = (7./5.)*pow(kappa, -23./20.);
  double ans = prefactor*Nlow*pow(X_k, -5./3.)
              *pow(1.+pow(X_k, x*(3.*kappa-1.)/6.)*pow(Nlow/Nhigh, x), -1./x);

  return ans;
}

double kappa_V_abs(double nu)
{
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);

  double nu_w = pow(kappa_width*kappa, 2.)*nu_c*sin(obs_angle);
  double X_k = nu/nu_w;
  double prefactor = n_e*electron_charge/(B_field*sin(obs_angle));
  double a = kappa - 1./3.; 
  double b = kappa + 1.;
  double c = kappa + 2./3.;
  double z = -kappa*kappa_width;
  double hyp2f1 = pow(1.-z, -a)*tgamma(c)*tgamma(b-a)/(tgamma(b)*tgamma(c-a))
                 *gsl_sf_hyperg_2F1(a, c-b, a-b+1., 1./(1.-z))+pow(1.-z, -b)
                 *tgamma(c)*tgamma(a-b)/(tgamma(a)*tgamma(c-b))
                 *gsl_sf_hyperg_2F1(b, c-a, b-a+1., 1./(1.-z));
  double Nlow = -(77./(100.*kappa_width))*pow(pow(sin(obs_angle), -114./50.)-1., 223./500.)*pow(X_k, -7./20.)*pow(kappa, -7./10)*pow(3., 1./6.)*(10./41.)*pow(2.*M_PI, 2.)
               /pow(kappa_width*kappa, 16./3.-kappa)*(kappa-2.)*(kappa-1.)
               *kappa/(3.*kappa-1.)*tgamma(5./3.)*hyp2f1;
  double Nhigh = -(143./10. * pow(kappa_width, -116./125.))*pow(pow(sin(obs_angle), -41./20.)-1., 1./2.)*(13.*13.*pow(kappa, -8.)+13./(2500.)*kappa - 263./5000.+47./(200.*kappa))*pow(X_k, -1./2.)*2.*pow(M_PI, 5./2.)/3.*(kappa-2.)*(kappa-1.)*kappa
                 /pow(kappa_width*kappa, 5.)*(2*tgamma(2.+kappa/2.)
                 /(2.+kappa)-1.);
  double x = (61./50.)*pow(kappa, -142./125.)+7./1000.;
  double ans = prefactor*Nlow*pow(X_k, -5./3.)
              *pow(1.+pow(X_k, x*(3.*kappa-1.)/6.)*pow(Nlow/Nhigh, x), -1./x);

  return ans;
}

double power_law_I(double nu)
{
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);

  double prefactor = (n_e_NT*pow(electron_charge,2.)*nu_c)/speed_light;
  double term1 = pow(3., power_law_p/2.)*(power_law_p-1.)*sin(obs_angle);
  double term2 = 2.*(power_law_p+1.)*(pow(gamma_min, 1.-power_law_p)-pow(gamma_max, 1.-power_law_p));
  double term3 = tgamma((3.*power_law_p-1.)/12.)*tgamma((3.*power_law_p+19.)/12.);
  double term4 = pow(nu/(nu_c*sin(obs_angle)), -(power_law_p-1.)/2.);
  double ans = prefactor*term1/term2*term3*term4;

 return ans;
}

double power_law_I_abs(double nu)
{
  double nu_c = (electron_charge * B_field)
               /(2. * M_PI * mass_electron * speed_light);

  double prefactor = (n_e_NT*pow(electron_charge,2.))/(nu*mass_electron*speed_light);
  double term1 = pow(3., (power_law_p+1.)/2.)*(power_law_p-1.);
  double term2 = 4.*(pow(gamma_min, 1.-power_law_p)-pow(gamma_max, 1.-power_law_p));
  double term3 = tgamma((3.*power_law_p+2.)/12.)*tgamma((3.*power_law_p+22.)/12.);
  double term4 = pow(nu/(nu_c*sin(obs_angle)), -(power_law_p+2.)/2.);
  double ans = prefactor*term1/term2*term3*term4;

 return ans;
}
