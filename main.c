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


double n_peak(double nu);
double K_s(double gamma, double n, double nu);
double my_Bessel_J(double n, double x);
double my_Bessel_dJ(double n, double x);
double MJ_f(double gamma);
double I(double gamma, double n, double nu);
double trapez_gamma(double min, double max, double n, double nu);
double trapez_n(double min, double max, double nu);
double gamma_integrand(double gamma, double n, double nu);
double gamma_integration_result(double n, double nu);
double n_summation(double nu);
double n_integration(double n_minus, double nu);
//double integrand(double gamma);

int main(int argc, char *argv[])
{
	//define parameters of calculation
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double nu = 100. * nu_c;
	n_summation(nu);
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

double MJ_f(double gamma)
{
	double beta = sqrt(1. - 1./(gamma*gamma));
	double d = (n_e * gamma * sqrt(gamma*gamma-1.) * exp(-gamma/theta_e))/(4. * M_PI * theta_e * gsl_sf_bessel_Kn(2, 1./theta_e));
	double ans = 1./(pow(m, 3.) * pow(c, 3.) * gamma*gamma * beta) * d;
	return ans;
}

double I(double gamma, double n, double nu)
{
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double beta = sqrt(1. - 1./(gamma*gamma));
	double cos_xi = (gamma * nu - n * nu_c)/(gamma * nu * beta * cos(theta));
	double ans = (2. * M_PI * e*e * nu*nu)/c * (pow(m, 3.) * pow(c, 3.) * gamma*gamma * beta * 2. * M_PI) * MJ_f(gamma) * K_s(gamma, n, nu);
	return ans;
}

double gamma_integrand(double gamma, double n, double nu)
{
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double beta = sqrt(1. - 1./(gamma*gamma));
	double cos_xi = (gamma * nu - n * nu_c)/(gamma * nu * beta * cos(theta));
	double prefactor = 1./(nu * beta * fabs(cos(theta)));
	double ans = prefactor * I(gamma, n, nu);
	return ans;
}

double gamma_integration_result(double n, double nu)
{
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double gamma_minus = ((n*nu_c)/nu - fabs(cos(theta))*sqrt((pow((n*nu_c)/nu, 2.)) - pow(sin(theta), 2.)))/(pow(sin(theta), 2));
	double gamma_plus  = ((n*nu_c)/nu + fabs(cos(theta))*sqrt((pow((n*nu_c)/nu, 2.)) - pow(sin(theta), 2.)))/(pow(sin(theta), 2));
	double result = trapez_gamma(gamma_minus, gamma_plus, n, nu);
	return result;
}

double trapez_gamma (double min, double max, double n, double nu)
{
	int i;
	float interval, sum=0., x;
	int divisions = 1000;

	interval = ((max-min) / (divisions-1));

	for (i=2; i<divisions; i++)
   	{
      		x    = min + interval * (i-1);
      		sum += gamma_integrand(x, n, nu)*interval;
   	}

	//sum += 0.5 *(gamma_integrand(min, n, nu) + gamma_integrand(max, n, nu)) * interval;
   	return (sum);
}

double trapez_n(double min, double max, double nu)
{
	int i;
	float interval, sum=0., x;
	int divisions = 1000;

	interval = ((max-min) / (divisions-1));

	for (i=2; i<divisions; i++)
   	{
      		x    = min + interval * (i-1);
      		sum += gamma_integration_result(x, nu)*interval;
   	}

	//sum += 0.5 *(gamma_integration_result(min, nu) + gamma_integration_result(max, nu)) * interval;
   	return (sum);
}

double n_integration(double n_minus, double nu)
{
	if(n_max < n_minus)
	{
		n_max = n_minus;
	}

	double ans = trapez_n(n_max, C * n_peak(nu), nu);
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
		j_nu += gamma_integration_result(x, nu);
	}

	printf("\n%e\n", j_nu);
	j_nu = j_nu + n_integration(n_minus, nu);
	printf("\n%e\n", j_nu);
	return j_nu;
}
