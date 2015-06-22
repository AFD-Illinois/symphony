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
//double gamma_cutoff = 1000.; also a kappa distribution parameter

//kappa distribution parameters
double kappa = 3.5;
double gamma_cutoff = 1000;

//function declarations
double n_peak(double nu);
double K_s(double gamma, double n, double nu);
double my_Bessel_J(double n, double x);
double my_Bessel_dJ(double n, double x);
double MJ_f(double gamma);
double I(double gamma, double n, double nu);
double gamma_integrand(double gamma, void * params);
double gamma_integration_result(double n, void * params);
double n_summation(double nu);
double n_integration(double n_minus, double nu);
double integrate(double min, double max, double n, double nu);
double gsl_integrate(double min, double max, double n, double nu);
double s_integrate(double min, double max, double n, double nu);
double normalize_f();
double power_law_to_be_normalized(double gamma, void * params);
double power_law_f(double gamma);
double kappa_to_be_normalized(double gamma, void * params);
double kappa_f(double gamma);
double n_integration_adaptive(double n_max, double n_minus);
double derivative(double n_start, double nu);
double D_thermal(double gamma, double nu);
double D_pl(double gamma, double nu);
double D_kappa(double gamma, double nu);

//struct to pass parameters to integrand
struct parameters
{
	double n;
	double nu;
};

//choose distribution function
#define MJ (0)
#define POWER_LAW (1)
#define KAPPA_DIST (2)
#define DISTRIBUTION_FUNCTION (MJ)

//choose absorptivity or emissivity
#define ABSORP (10)
#define EMISS  (11)
#define MODE   (ABSORP)

//choose polarization mode
#define K_I (15)
#define K_Q (16)
#define K_U (17)
#define K_V (18)
#define K_XX (19)
#define K_OO (20)
#define POL_MODE (K_OO)

int main(int argc, char *argv[])
{
	//define parameters of calculation
	double nu_c = (e * B)/(2. * M_PI * m * c);
	int index = 0;
	//double nu = 1. * nu_c;
	for(index; index < 151; index++)
	{
		double nu = pow(10., index/100.) * nu_c;
		//double nu = nu_c * index/5.;
		printf("\n%e	%e", nu/nu_c, n_summation(nu));
	}
	printf("\n");
	return 0;
}

double n_peak(double nu)
{
	double nu_c = (e * B)/(2. * M_PI * m * c);
	if(nu <= nu_c * theta_e*theta_e || theta_e < 1.){
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
	double T_x = (2.*nu*cos(theta))/(nu_c*pow(sin(theta), 2.)+sqrt(nu_c*nu_c*pow(sin(theta), 4.)+4.*nu*nu*pow(cos(theta), 2.)));
#if POL_MODE == K_I
	double ans = K_xx + K_yy;
#elif POL_MODE == K_Q
	double ans = K_xx - K_yy;
#elif POL_MODE == K_U
	double ans = 0.;
#elif POL_MODE == K_V
	double ans = -2.*M*N*my_Bessel_J(n, z)*my_Bessel_dJ(n, z);
#elif POL_MODE == K_XX //factor of 2 out front comes from eq. 43 as opposed to 42
	double ans = 2.*pow(M*T_x*my_Bessel_J(n, z)+N*my_Bessel_dJ(n, z), 2.)/(1. + T_x*T_x);
#elif POL_MODE == K_OO //factor of 2 again from eq. 43
	double ans = 2.*pow(M*T_x*my_Bessel_J(n, z)-N*my_Bessel_dJ(n, z), 2.)/(1. + T_x*T_x);
#endif
	return ans;
}

double D_thermal(double gamma, double nu)
{
	double prefactor = (M_PI * nu / (m*c*c)) * (n_e/(theta_e * gsl_sf_bessel_Kn(2, 1./theta_e)));
	double body = (-1./theta_e) * exp(-gamma/theta_e);
	double f = prefactor * body;
	return f;
}

double D_pl(double gamma, double nu)
{
	double pl_norm = 1./(normalize_f());
	double prefactor = (M_PI * nu / (m*c*c)) * (n_e_NT*(p-1.))/((pow(gamma_min, 1.-p) - pow(gamma_max, 1.-p)));
	double term1 = ((-p-1.)*exp(-gamma/gamma_cutoff)*pow(gamma,-p-2.)/(sqrt(gamma*gamma - 1.)));
	double term2 = (exp(-gamma/gamma_cutoff) * pow(gamma,(-p-1.))/(gamma_cutoff * sqrt(gamma*gamma - 1.)));
	double term3 = (exp(-gamma/gamma_cutoff) * pow(gamma,-p))/pow((gamma*gamma - 1.), (3./2.));
	double f = pl_norm * prefactor * (term1 - term2 - term3);
	return f;
}

double D_kappa(double gamma, double nu)
{
	//double prefactor = (1./normalize_f()) * 4. * M_PI*M_PI * nu * m*m * c;
	//double term1 = ((- kappa - 1.) / (kappa * theta_e)) * pow((1. + (gamma - 1.)/(kappa * theta_e)), -kappa-2.);
	//double term2 = pow((1. + (gamma - 1.)/(kappa * theta_e)), (- kappa - 1.)) * (- 1./gamma_cutoff);
	//double f = prefactor * (term1 + term2) * exp(-gamma/gamma_cutoff);
	//below is differential for kappa WITHOUT cutoff
	double term1 = (-1.-kappa)*(-2.+kappa)*(-1.+kappa) * nu * M_PI * n_e;
	double term2 = pow((1. + (gamma-1.)/(kappa*theta_e)), -2.-kappa);
	double term3 = 2. * c*c * pow(kappa, 3.) * pow(theta_e, 3.) * m;
	double f = (term1 * term2)/term3;
	return f;
}


double MJ_f(double gamma)
{
	double beta = sqrt(1. - 1./(gamma*gamma));
	double d = (n_e * gamma * sqrt(gamma*gamma-1.) * exp(-gamma/theta_e))/(4. * M_PI * theta_e * gsl_sf_bessel_Kn(2, 1./theta_e));
	double ans = 1./(pow(m, 3.) * pow(c, 3.) * gamma*gamma * beta) * d;
	return ans;
}

double power_law_to_be_normalized(double gamma, void * params)
{
	double norm_term = 4. * M_PI;
	double prefactor = n_e_NT * (p - 1.) / (pow(gamma_min, 1. - p) - pow(gamma_max, 1. - p));
	//double body = pow(gamma, -p) * exp(- gamma / gamma_cutoff);
	double body = pow(gamma, -p);
	double ans = norm_term * prefactor * body;
	return ans;
}

double power_law_f(double gamma)
{
	double beta = sqrt(1. - 1./(gamma*gamma));
	double prefactor = n_e_NT * (p - 1.) / (pow(gamma_min, 1. - p) - pow(gamma_max, 1. - p));
	//double body = pow(gamma, -p) * exp(- gamma / gamma_cutoff);
	double body = pow(gamma, -p);
	double ans = 1./normalize_f() * prefactor * body * 1./(pow(m, 3.) * pow(c, 3.) * gamma*gamma * beta);
	return ans;
}

double kappa_to_be_normalized(double gamma, void * params)
{
	double kappa_body = pow((1. + (gamma - 1.)/(kappa * theta_e)), -kappa-1);
	//double kappa_body = pow((1. + gamma/(kappa*theta_e)), -kappa-1.);
	double cutoff = exp(-gamma/gamma_cutoff);
	double norm_term = 4. * M_PI * pow(m, 3.) * pow(c, 3.) * gamma * sqrt(gamma*gamma-1.);
	//double ans = kappa_body * cutoff * norm_term;
	double ans = kappa_body * norm_term;
	return ans;
}

double kappa_f(double gamma)
{
	double norm = 1./normalize_f();
	double kappa_body = pow((1. + (gamma - 1.)/(kappa * theta_e)), -kappa-1);
	//double kappa_body = pow((1.+ gamma/(kappa*theta_e)), -kappa-1.);
	double cutoff = exp(-gamma/gamma_cutoff);
	//double ans = norm * kappa_body * cutoff;
	double ans = norm * kappa_body;
	return ans;
}

double I(double gamma, double n, double nu)
{
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double beta = sqrt(1. - 1./(gamma*gamma));
	double cos_xi = (gamma * nu - n * nu_c)/(gamma * nu * beta * cos(theta));
//I(gamma, n, nu) depends on the distribution function
#if DISTRIBUTION_FUNCTION == MJ
	double ans = (2. * M_PI * e*e * nu*nu)/c * (pow(m, 3.) * pow(c, 3.) * gamma*gamma * beta * 2. * M_PI) * MJ_f(gamma) * K_s(gamma, n, nu);
#elif DISTRIBUTION_FUNCTION == POWER_LAW
	double ans = (2. * M_PI * e*e * nu*nu)/c * (pow(m, 3.) * pow(c, 3.) * gamma*gamma * beta * 2. * M_PI) * power_law_f(gamma) * K_s(gamma, n, nu);
#elif DISTRIBUTION_FUNCTION == KAPPA_DIST
	double ans = (2. * M_PI * e*e * nu*nu)/c * (pow(m, 3.) * pow(c, 3.) * gamma*gamma * beta * 2. * M_PI) * kappa_f(gamma) * K_s(gamma, n, nu);
#else
	double ans = 0;
#endif
	return ans;
}

#if MODE == EMISS
double gamma_integrand(double gamma, void * params)
{
	struct parameters n_and_nu = *(struct parameters*) params;
	double n = n_and_nu.n;
	double nu = n_and_nu.nu;
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double beta = sqrt(1. - 1./(gamma*gamma));
	double cos_xi = (gamma * nu - n * nu_c)/(gamma * nu * beta * cos(theta));
	double prefactor = 1./(nu * beta * fabs(cos(theta)));
	double ans = prefactor * I(gamma, n, nu);
	return ans;
}
#endif

#if MODE == ABSORP
double gamma_integrand(double gamma, void * params)
{
	struct parameters n_and_nu = *(struct parameters*) params;
	double n = n_and_nu.n;
	double nu = n_and_nu.nu;
	double nu_c = (e * B)/(2. * M_PI * m * c);
	double beta = sqrt(1. - 1./(gamma*gamma));
	double prefactor = -c*e*e / (2. * nu);
#if DISTRIBUTION_FUNCTION == MJ
	double ans = prefactor*gamma*gamma*beta*D_thermal(gamma, nu)*K_s(gamma, n, nu)*(1./(nu*beta*fabs(cos(theta))));
#elif DISTRIBUTION_FUNCTION == POWER_LAW
	double ans = prefactor*gamma*gamma*beta*D_pl(gamma, nu)*K_s(gamma, n, nu)*(1./(nu*beta*fabs(cos(theta))));
#elif DISTRIBUTION_FUNCTION == KAPPA_DIST
	double ans = prefactor*gamma*gamma*beta*D_kappa(gamma, nu)*K_s(gamma, n, nu)*(1./(nu*beta*fabs(cos(theta))));
#else
	double ans = 0.;
#endif
	return ans;
}
#endif

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
	double n_start = (int)(n_max + n_minus + 1.);
	double ans = 0.;
	double contrib = 0.;
	int i = 0;
	double delta_n = 1.e5;// if using Simpson's rule
	//double delta_n = 1.e5;
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

#if DISTRIBUTION_FUNCTION == MJ //we know where peak is analytically
	j_nu = j_nu + n_integration(n_minus, nu);
#elif DISTRIBUTION_FUNCTION != MJ //need to use adaptive integration
	j_nu = j_nu + n_integration_adaptive(n_minus, nu);
#else
	j_nu = 0.;
#endif
	return j_nu;
}

double derivative(double n_start, double nu)
{
	gsl_function F;
	double result, abserr;
	F.function = gamma_integration_result;
	F.params = &nu;

	gsl_deriv_central(&F, n_start, 1e-8, &result, &abserr);
	return result;
}

double normalize_f()
{
	static double ans = 0;
	if(ans != 0)
	{
		return ans;
	}

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	double result, error;

	gsl_function F;
#if DISTRIBUTION_FUNCTION == POWER_LAW
	F.function = &power_law_to_be_normalized;
#elif DISTRIBUTION_FUNCTION == KAPPA_DIST
	F.function = &kappa_to_be_normalized;
#else
	return 0;
#endif
	double unused = 0.;
	F.params = &unused;

	gsl_integration_qagiu(&F, 1, 0, 1e-8, 1000,
	                       w, &result, &error);
	gsl_integration_workspace_free (w);
	ans = result;
	return result;
}


