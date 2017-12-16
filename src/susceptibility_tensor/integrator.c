#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "susceptibility_tensor.h"

/*tau_integrator: integrator for the first integral (the tau integral) for the
 *                components of the susceptibility tensor.  The chi_33 integral
 *                is faster when we use a fixed-order Gaussian quadrature
 *                method GSL QNG, rather than the GSL QAWO adaptive integrator
 *                used for the other components.
 *
 *@params: double gamma (the second integration is carried out over this
 *         variable), void * parameters (a pointer to the struct of parameters)
 *
 *@returns: numerically evaluated tau integral at a given gamma of a given
 *          component of the susceptibility tensor
 */
double tau_integrator(double gamma, void * parameters)
{
	struct parameters * params = (struct parameters*) parameters;

	if(gamma == 1.)
	{
		return 0.;
	}

        double ans_tot  = 0.;
	double ans_step = 0.;
	double error    = 0.;
	double step;
	double start    = 0.;
        size_t n        = 50;
        size_t limit    = 5000;
        double epsabs   = 0.;
        double epsrel   = 1e-8;
	double sign_correction;
	enum gsl_integration_qawo_enum gsl_weight;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);

	if(params->real == 1)
        {
		gsl_weight      = GSL_INTEG_SINE;
		sign_correction = -1.;
        }
	else
        {
		gsl_weight      = GSL_INTEG_COSINE;
		sign_correction = 1.;
        }

	if(params->tau_integrand == &chi_33_integrand)
        {
                step = 2. * params->pi / gamma;
                sign_correction = 1.;
        }
        else
        {
                step     = params->pi/gamma;
        }

	gsl_integration_qawo_table * table =
                               gsl_integration_qawo_table_alloc(gamma, step, gsl_weight, n);


	//need to update value of gamma
	params-> gamma = gamma;

	gsl_set_error_handler_off();
	gsl_function F;
	F.function = params->tau_integrand;
	F.params   = params;

	int i            = 0;
	int max_counter  = 2500;
	double tolerance = 1e-7;
	int counts       = 0;

	/*TODO: explain this */
	int small_counter = 0;
	double small_tol  = 1e-14;  //TODO: want to set this adaptively
	int max_small_counter = 1000;

	while(i == 0 || counts < max_counter)
	{

		if(params->tau_integrand == &chi_33_integrand)
		{
			gsl_integration_qng(&F, i*step, (i+1)*step, epsabs, epsrel, &ans_step, &error, &limit);
		}
		else
		{
			gsl_integration_qawo(&F, i*step, epsabs, epsrel, limit, w, table, &ans_step, &error);
		}
		ans_tot += ans_step;
		i       += 1;


		if(fabs(ans_step / ans_tot) < tolerance)
		{
			counts += 1;
		}

		if(fabs(ans_tot) < small_tol)
		{
			small_counter++;
		}
		if(small_counter >= max_small_counter)
		{
			return 0.;
		}

	}

	gsl_integration_qawo_table_free(table);
	gsl_integration_workspace_free(w);
	ans_tot = ans_tot * sign_correction;

	return ans_tot;
}

/*trapezoidal: trapezoidal rule integrator for gamma integral, which can be
 *             faster than adaptive approaches because the integrand is
 *             smooth; one can also closely monitor the number of function
 *             calls to tau_integrator, which is important because these calls
 *             are rather expensive.
 *
 *@params: pointer to struct of parameters *p, start (starting point for gamma
 *         integral), end (ending point for gamma integral), samples (number of
 *         calls to tau_integrator allowed)
 *
 *@returns: gamma integral for a given component of chi_ij, evaluated using
 *          a trapezoidal rule integrator
 */
double trapezoidal(struct parameters *p, double start, double end, int samples)
{

        int i = 0;

	double step      = (end - start)/samples;
        double tolerance = 1e-6;
        double ans_step  = 0.;
        double ans_tot   = 0.;

	double p1 = p->gamma_integrand(start+i*step, p);
        double p2 = p->gamma_integrand(start+(i+1)*step, p);

        while(start + i * step <= end)
        {

                ans_step = step * (p1 + p2)/2.;

		printf("\nSTART: %e     END: %e", start+i*step, start + (i+1)*step);

                ans_tot += ans_step;
                i++;

                p1 = p2;
                p2 = p->gamma_integrand(start+(i+1)*step, p);
        }

	return ans_tot;
}

/*trapezoidal_adaptive: trapezoidal rule integrator for gamma integral, similar
 *             to the above function, except it determines the ending point
 *             adaptively.
 *
 *@params: pointer to struct of parameters *p, start (starting point for gamma
 *         integral), step (step size between calls to tau_integrator)
 *
 *@returns: gamma integral for a given component of chi_ij, evaluated using
 *          an adaptive trapezoidal rule integrator
 */
double trapezoidal_adaptive(struct parameters *p, double start, double step)
{

        int i = 0;
        double tolerance = 1e-6;
        double ans_step  = 0.;
        double ans_tot   = 0.;

        double p1 = p->gamma_integrand(start+i*step, p);
        double p2 = p->gamma_integrand(start+(i+1)*step, p);

	while(ans_tot == 0 || fabs(ans_step/ans_tot) > tolerance)
        {

                ans_step = step * (p1 + p2)/2.;

                printf("\nSTART: %e     END: %e", start+i*step, start + (i+1)*step);

                ans_tot += ans_step;
                i++;

                p1 = p2;
                p2 = p->gamma_integrand(start+(i+1)*step, p);
        }

        return ans_tot;
}

/*end_approx: approximate ending point for the gamma integral, beyond which
 *            contributions to the integral are assumed to be negligible.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: approximate end location for gamma integral
 */
double end_approx(struct parameters *p)
{
	double end;

        double MJ_max        = 0.5 * (3. * p->theta_e + sqrt(4. + 9. * p->theta_e * p->theta_e));
        double PL_max_real   = sqrt((1. + p->power_law_p)/p->power_law_p);
	double PL_max_moving = 50./sqrt(p->omega/p->omega_c) + 9. * pow(p->omega/p->omega_c, 1./3.);
	double kappa_max     = (-3. + 3. * p->kappa_width * p->kappa 
			      + sqrt(1. - 4. * p->kappa 
			  	     - 18. * p->kappa_width * p->kappa 
                                     + 4. * pow(p->kappa, 2.) 
                                     + 9. * pow(p->kappa_width * p->kappa, 2.))) 
			   / (2. * (p->kappa - 2.));

        if(p->distribution == p->MAXWELL_JUETTNER)
        {
                end = 7. * MJ_max;
        }
        else if(p->distribution == p->POWER_LAW)
        {
		end = PL_max_moving;
        }
        else if(p->distribution == p->KAPPA_DIST)
	{
		end = 7. * kappa_max;
	}
	else
        {
                printf("\ndistribution or real/imag is set incorrectly");
                return 1.;
        }

	return end;
}

/*gsl_integrator: wrapper for a GSL integrator, to be used for the gamma
 *                integral.  This function can be very expensive, because
 *                it makes a large number of function calls to a region
 *                of parameter space that has a slowly convergent tau integral.
 *
 *@params: pointer to struct of parameters *p, start (starting point for gamma
 *         integral), end (ending point for gamma integral)
 *
 *@returns: gamma integral for a given component of chi_ij, evaluated using
 *          a GSL adaptive integrator
 */
double gsl_integrator(struct parameters *p, double start, double end)
{
	gsl_function F;
        F.function = p->gamma_integrand;
        F.params   = p;

        gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000);

        double ans    = 0.;
        double error  = 0.;
        size_t limit  = 5000;
        double epsabs = 0.;
        double epsrel = 1e-5;
        int gsl_key   = 1;

        printf("\n%e\n", end);


	gsl_integration_qags(&F, start, end, epsabs, epsrel, limit, w, &ans, &error);

//      gsl_integration_qng(&F, start, end, epsabs, epsrel, &ans, &error, &limit);
//	gsl_integration_qag(&F, start, end, epsabs, epsrel, limit, gsl_key, w, &ans, &error);
//	gsl_integration_qagiu(&F, start, epsabs, epsrel, limit, w, &ans, &error);
	
	gsl_integration_workspace_free(w);

	return ans;
}

/*gamma_integrator: function that evaluates the gamma integral, using one of
 *                  the above integration algorithms.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: gamma integral for a given component of chi_ij
 */
double gamma_integrator(struct parameters *p)
{
        double prefactor = 2. * p->pi * p->omega_p*p->omega_p / (p->omega * p->omega);

        double start  = 1.;
        double end = end_approx(p);
	double ans_tot;

	/*for power-law, there is a singularity at low frequency at gamma = 1,
          which cannot be resolved by a trapezoidal integrator.  We are forced
          to use the (much slower) GSL integrator QAGS in these cases. */
	if(p->distribution == p->POWER_LAW && p->omega/p->omega_c < 5.)
	{
		ans_tot = gsl_integrator(p, start, end);
	}
	else
	{
		ans_tot = trapezoidal(p, start, end, 100);
	}

//	double ans_tot = trapezoidal(p, start, end, 100);
//	double ans_tot = trapezoidal_adaptive(p, start, 1.);
//	double ans_tot = gsl_integrator(p, start, end);

        return prefactor * ans_tot;
}
