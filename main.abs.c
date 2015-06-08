#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>

/* other C header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

//global variables
double m = 9.1093826e-28;
double c = 2.99792458e10;
double theta_e = 10.;
double e = 4.80320680e-10;
double B = 30.;
double n_e = 1.;
double theta = (M_PI  / 3.);
int C = 10;
double n_max = 30.;

//power law parameters
double p = 3.;
double gamma_min = 1.;
double gamma_max = 1000.;
double n_e_NT = 1.;
double gamma_cutoff = 1000.;

//function declarations
double n_peak(double nu);
double K_s(double gamma, double n, double nu);
double K_q(double gamma, double n, double nu);
double K_u(double gamma, double n, double nu);
double K_v(double gamma, double n, double nu);
double D_thermal(double gamma, double nu);
double D_pl(double gamma, double nu);
double my_Bessel_J(double n, double x);
double my_Bessel_dJ(double n, double x);
double MJ_f(double gamma);
double I(double gamma, double n, double nu);
double absorptivity_integrand(double gamma, double n, double nu);
double gamma_integrand(double gamma, void * params);
double gamma_integration_result(double n, void * params);
double n_summation(double nu);
double n_integration(double n_minus, double nu);
double integrate(double min, double max, double n, double nu);
double gsl_integrate(double min, double max, double n, double nu);
double normalize_f();
double power_law_to_be_normalized(double gamma, void * params);
double power_law_f(double gamma);
double n_integration_adaptive(double n_max, double n_minus);
double derivative(double n_start, double nu);

//struct to pass parameters to integrand
struct parameters
{
	double n;
	double nu;
};

int main(int argc, char *argv[])
{
	//define parameters of calculation
	double nu_c = (e * B)/(2. * M_PI * m * c);
	int index = 0;
	double nu = 1. * nu_c;
	//for(index; index < 7; index++)
	//{
	//	double nu = pow(10., index) * nu_c;
	//	n_summation(nu);
	//	printf("\n%e	%e", nu/nu_c, n_summation(nu));
	//}
	//printf("\n");
	//n_summation(nu);
	//struct parameters p;
	//p.n = 10.;
	//p.nu = nu;
	printf("\n%e\n", D_pl(10., nu));
	//printf("\n%e\n", absorptivity_integrand(10., 10., nu));
	return 0;
}

double n_peak(double nu)
{
	double nu_c = (e * B)/(2. * M_PI * m * c);
	if(nu <= nu_c * theta_e*theta_e){
		double beta_67 = sqrt((1. - 1./pow((1. + theta_e),2.)));
		double n_peak = (theta_e + 1. + pow((2. * theta_e * nu / nu_c),1./3.)) * (nu/nu_c) * (1. - beta_67*beta_67 * pow(cos(theta),2.));
		return n_peak;
	}
	else{
		double beta_68 = sqrt(1. - pow((2. * theta_e * nu / nu_c), -2./3.));
		double n_peak = (theta_e + 1. + pow((2. * theta_e * nu / nu_c),1./3.)) * (nu/nu_c) * (1. - beta_68*beta_68 * pow(cos(theta),2.));
		return n_peak;
	}
}

double K_s(double gamma, double n, double nu)
{
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double beta = sqrt(1. - 1./(gamma*gamma));
	double cos_xi = (gamma * nu - n * nu_c)/(gamma * nu * beta * cos(theta));
	double M = (cos(theta) - beta * cos_xi)/sin(theta);
	double N = beta * sqrt(1 - (cos_xi*cos_xi));
	double z = (nu * gamma * beta * sin(theta) * sqrt(1. - cos_xi*cos_xi))/nu_c;
	double K_xx = M*M * pow(my_Bessel_J(n, z), 2.);
	double K_yy = N*N * pow(my_Bessel_dJ(n, z), 2.);
	double ans = K_xx + K_yy;
	return ans;
}

double K_q(double gamma, double n, double nu)
{
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double beta = sqrt(1. - 1./(gamma*gamma));
	double cos_xi = (gamma * nu - n * nu_c)/(gamma * nu * beta * cos(theta));
	double M = (cos(theta) - beta * cos_xi)/sin(theta);
	double N = beta * sqrt(1 - (cos_xi*cos_xi));
	double z = (nu * gamma * beta * sin(theta) * sqrt(1. - cos_xi*cos_xi))/nu_c;
	double K_xx = M*M * pow(my_Bessel_J(n, z), 2.);
	double K_yy = N*N * pow(my_Bessel_dJ(n, z), 2.);
	double ans = K_xx - K_yy;
	return ans;
}

double K_u(double gamma, double n, double nu)
{
	double ans = 0.;
	return ans;
}

double K_v(double gamma, double n, double nu)
{
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double beta = sqrt(1. - 1./(gamma*gamma));
	double cos_xi = (gamma * nu - n * nu_c)/(gamma * nu * beta * cos(theta));
	double M = (cos(theta) - beta * cos_xi)/sin(theta);
	double N = beta * sqrt(1 - (cos_xi*cos_xi));
	double z = (nu * gamma * beta * sin(theta) * sqrt(1. - cos_xi*cos_xi))/nu_c;
	double ans = -2.*M*N*my_Bessel_J(n, z)*my_Bessel_dJ(n, z);
	return ans;
}

//this is specific to the thermal (MJ) distribution
double D_thermal(double gamma, double nu)
{
	double prefactor = (M_PI * nu / (m*c*c)) * (n_e/(theta_e * gsl_sf_bessel_Kn(2, 1./theta_e)));
	double body = (-1./theta_e) * exp(-gamma/theta_e);
	double f = prefactor * body;
	return f;
}

//this is specific to the power-law distribution
double D_pl(double gamma, double nu)
{
	double pl_norm = 1./normalize_f();
	double prefactor = (M_PI * nu / (m*c*c)) * (n_e_NT*(p-1.))/((pow(gamma_min, 1.-p) - pow(gamma_max, 1.-p)));
	double term1 = ((-p-1.)*exp(-gamma/gamma_cutoff)*pow(gamma,-p-2.)/(sqrt(gamma*gamma - 1.)));
	double term2 = (exp(-gamma/gamma_cutoff) * pow(gamma,(-p-1.))/(gamma_cutoff * sqrt(gamma*gamma - 1.)));
	double term3 = (exp(-gamma/gamma_cutoff) * pow(gamma,-p))/pow((gamma*gamma - 1.), (3./2.));
	double f = pl_norm * prefactor * (term1 - term2 - term3);
	return f;
}

double MJ_f(double gamma)
{
	double beta = sqrt(1. - 1./(gamma*gamma));
	double d = (n_e * gamma * sqrt(gamma*gamma-1.) * exp(-gamma/theta_e))/(4. * M_PI * theta_e * gsl_sf_bessel_Kn(2, 1./theta_e));
	double ans = 1./(pow(m, 3.) * pow(c, 3.) * gamma*gamma * beta) * d;
	return ans;
}

//power law distribution function
double power_law_to_be_normalized(double gamma, void * params)
{
	double norm_term = 4. * M_PI;
	double prefactor = n_e_NT * (p - 1.) / (pow(gamma_min, 1. - p) - pow(gamma_max, 1. - p));
	double body = pow(gamma, -p) * exp(- gamma / gamma_cutoff);
	double ans = norm_term * prefactor * body;
	return ans;
}

double power_law_f(double gamma)
{
	double beta = sqrt(1. - 1./(gamma*gamma));
	double prefactor = n_e_NT * (p - 1.) / (pow(gamma_min, 1. - p) - pow(gamma_max, 1. - p));
	double body = pow(gamma, -p) * exp(- gamma / gamma_cutoff);
	double ans = 1./normalize_f() * prefactor * body * 1./(pow(m, 3.) * pow(c, 3.) * gamma*gamma * beta);
	return ans;
}

double I(double gamma, double n, double nu)
{
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double beta = sqrt(1. - 1./(gamma*gamma));
	double cos_xi = (gamma * nu - n * nu_c)/(gamma * nu * beta * cos(theta));
	//distribution function goes in here
	double ans = (2. * M_PI * e*e * nu*nu)/c * (pow(m, 3.) * pow(c, 3.) * gamma*gamma * beta * 2. * M_PI) * MJ_f(gamma) * K_s(gamma, n, nu);
	return ans;
}

double absorptivity_integrand(double gamma, double n, double nu)
{
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double beta = sqrt(1. - 1./(gamma*gamma));
	double cos_xi = (gamma * nu - n * nu_c)/(gamma * nu * beta * cos(theta));
	double prefactor = - c * e*e / (2. * nu);
	//polarization mode goes in below
	double ans = prefactor*gamma*gamma*beta*D_thermal(gamma, nu)*K_s(gamma, n, nu)*(1./(nu*beta*fabs(cos(theta))));
	return ans;
}


//modified for the absorptivity
double gamma_integrand(double gamma, void * params)
{
	struct parameters n_and_nu = *(struct parameters*) params;
	double n = n_and_nu.n;
	double nu = n_and_nu.nu;
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double beta = sqrt(1. - 1./(gamma*gamma));
	double prefactor = -c*e*e / (2. * nu);
	//polarization mode goes in below
	double ans = prefactor*gamma*gamma*beta*D_thermal(gamma, nu)*K_s(gamma, n, nu)*(1./(nu*beta*fabs(cos(theta))));
	return ans;
}

//double gamma_integration_result(double n, double nu)
double gamma_integration_result(double n, void * params)
{
	double nu = *(double *) params;
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double gamma_minus = ((n*nu_c)/nu - fabs(cos(theta))*sqrt((pow((n*nu_c)/nu, 2.)) - pow(sin(theta), 2.)))/(pow(sin(theta), 2));
	double gamma_plus  = ((n*nu_c)/nu + fabs(cos(theta))*sqrt((pow((n*nu_c)/nu, 2.)) - pow(sin(theta), 2.)))/(pow(sin(theta), 2));
	double result = gsl_integrate(gamma_minus, gamma_plus, n, nu);

	return result;
}

double n_integration(double n_minus, double nu)
{
	if(n_max < n_minus)
	{
		n_max = (int) (n_minus+1);
	}
	
	double ans = gsl_integrate(n_max, C * n_peak(nu), -1, nu);
	return ans;
}

double n_integration_adaptive(double n_minus, double nu)
{
	//if(n_max < n_minus)
	//{
	//	n_max = (int)(n_minus + 1.)
	//}
	double n_start = (int)(n_max + n_minus + 1.);
	double ans = 0.;
	double contrib = 0.;
	int i = 0;
	double delta_n = 1.e5;
	double deriv_tol = 1.e-10;
	double tolerance = 1.e13;

	while(contrib >= ans/tolerance)
	{
		double deriv = derivative(n_start, nu);
		if(fabs(deriv) < deriv_tol)
		{
			delta_n = 100. * delta_n;
			//delta_n = 1. * delta_n;
		}

		contrib = gsl_integrate(n_start, (n_start + delta_n), -1, nu);
		ans = ans + contrib;
		n_start = n_start + delta_n;
		i++;
		//printf("\n%d\n", i);
	}

	return ans;
}

double n_summation(double nu)
{
	double j_nu = 0.;
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double n_minus = (nu/nu_c) * fabs(sin(theta));
	int x = (int)(n_minus+1.);
	for(x; x <= n_max + (int)n_minus ; x++)
	{
		j_nu += gamma_integration_result(x, &nu);
	}

	j_nu = j_nu + n_integration(n_minus, nu);
	//printf("\n%e\n", j_nu);
	return j_nu;
}

double derivative(double n_start, double nu)
{
	double dx = 0.01;
	double yhigh = gamma_integration_result(n_start + dx/2., &nu);
	double ylow  = gamma_integration_result(n_start - dx/2., &nu);
	double ans = (yhigh - ylow)/dx;
	return ans;
}
