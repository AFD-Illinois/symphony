#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

struct parameters
{
	double n;
	double nu;
};

double gamma_integrand(double gamma, void * params);
double gamma_integration_result(double n, void * params);
double power_law_to_be_normalized(double gamma, void * params);
double kappa_to_be_normalized(double gamma, void * params);

double integrate(double min, double max, double n, double nu)
{
	//n integration
	if(n < 0)
	{
		int i;
		float interval, sum=0., x;
		int divisions = 1200;
		interval = ((max-min) / (divisions-1));

		for (i=2; i<divisions; i++)
   		{
      			x    = min + interval * (i-1);
      			sum += gamma_integration_result(x, &nu)*interval;
   		}
		return sum;
	}
	//gamma integration
	else
	{
		int i;
		float interval, sum=0., x;
		int divisions = 1200;
		struct parameters n_and_nu;
		n_and_nu.n = n;
		n_and_nu.nu = nu;

		interval = ((max-min) / (divisions-1));

		for (i=2; i<divisions; i++)
   		{
      			x    = min + interval * (i-1);
      			sum += gamma_integrand(x, &n_and_nu)*interval;
   		}
		//sum += 0.5 *(gamma_integrand(min, n, nu) + gamma_integrand(max, n, nu)) * interval;
   		return (sum);
	}
}

double gsl_integrate(double min, double max, double n, double nu)
{
	//n integration
	if(n < 0)
	{
		double nu_c = (e * B)/(2. * M_PI * m * c);
		if(nu/nu_c > 1.e7)
		{
			gsl_set_error_handler_off();
		}

		gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
 		double result, error;
		gsl_function F;
  		F.function = &gamma_integration_result;
		F.params = &nu;
		gsl_integration_qag(&F, min, max, 0., 1.e-3, 1000,
                      3,  w, &result, &error);
		gsl_integration_workspace_free (w);
		return result;
	}
	else
	{
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
		double result, error;
		struct parameters n_and_nu;
		n_and_nu.n = n;
		n_and_nu.nu = nu;
		gsl_function F;
		F.function = &gamma_integrand;
		F.params = &n_and_nu;
		gsl_integration_qag(&F, min, max, 0., 1.e-3, 1000,
                      3,  w, &result, &error);
		gsl_integration_workspace_free (w);
		return result; 
	}
}
double s_integrate(double min, double max, double n, double nu)
{
	if(n < 0)
	{
		int i = 1500;
		if(i % 2 != 0)
		{
			i = i + 1;
		}
		double h = (max - min)/i;
		double s = gamma_integration_result(min, &nu) + gamma_integration_result(max, &nu);

		int index = 1;

		for(index; index < i; index+=2)
		{
			s += 4. * gamma_integration_result(min+index*h, &nu);
		}
		index = 0;
		for(index; index < i-1; index+=2)
		{
			s += 2 * gamma_integration_result(min+index*h, &nu);
		}
		return s * h / 3.;
	}

	else
	{
		struct parameters n_and_nu;
		n_and_nu.n = n;
		n_and_nu.nu = nu;

		int i = 1500;
		if(i % 2 != 0)
		{
			i = i + 1;
		}
		double h = (max - min)/i;
		double s = 0.;
		//double s = gamma_integrand(min+0.1, &n_and_nu) + gamma_integrand(max-0.1, &n_and_nu);
		//printf("\n%e\n", s);
		int index= 1;

		for(index; index < i; index+=2)
		{
			s += 4. * gamma_integrand(min+index*h, &n_and_nu);
		}
		index = 2;
		for(index; index < i-1; index+=2)
		{
			s += 2 * gamma_integrand(min+index*h, &n_and_nu);
		}
		//printf("\n%e\n", s);
		return s * h / 3.;
	}
}
