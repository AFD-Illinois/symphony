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


double n_peak(double nu);
double K_s(double gamma, int n, double nu);
double my_Bessel_J(double n, double x);
double my_Bessel_dJ(double n, double x);
double MJ_f(double gamma);
double I(double gamma, double n, double nu);

int main(int argc, char *argv[])
{
	//define parameters of calculation
	double nu_c = (e * B)/(2. * M_PI * m * c);
	//printf("\n%f\n", nu_c);
	double nu = 1. * nu_c;
	double * jn;
	//printf("\n%f\n", K_s(10, 10, nu));
	printf("\n%f\n", I(5, 4., nu)*1e18 );
	//printf("\n%f\n", my_Bessel_dJ(2000., 3000.));
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

double K_s(double gamma, int n, double nu)
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
	//double ans = N;
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
