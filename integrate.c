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

double integrate(double min, double max, double n, double nu)
{
	//n integration
	if(n < 0)
	{
		int i;
		float interval, sum=0., x;
		int divisions = 1000;
		//printf("\n%e\n", test());
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
		int divisions = 1000;
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
		//return 5;
	}
}

double gsl_integrate(double min, double max, double n, double nu)
{
	//n integration
	if(n < 0)
	{
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
 		double result, error;

		gsl_function F;
  		F.function = &gamma_integration_result;
		F.params = &nu;

		gsl_integration_qag(&F, min, max, 0.0, 1e-3, 1000,
                      1,  w, &result, &error);
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
		//printf("\n\n%e\n\n\n", n_and_nu.nu);

		gsl_function F;
		F.function = &gamma_integrand;
		F.params = &n_and_nu;

		gsl_integration_qag(&F, min, max, 0.0, 1e-3, 1000,
                      1,  w, &result, &error);
		  gsl_integration_workspace_free (w);
		
		return result; 
	}
}

double normalize_f()
{
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	double result, error;

	gsl_function F;
	F.function = &power_law_to_be_normalized;
	double unused = 0.;
	F.params = &unused;

	gsl_integration_qagiu(&F, 1, 0.0, 1e-3, 1000,
	                       w, &result, &error);
	gsl_integration_workspace_free (w);
	return result;

}
