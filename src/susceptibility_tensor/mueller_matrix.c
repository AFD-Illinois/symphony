#include <stdio.h>
#include <math.h>
#include "susceptibility_tensor.h"
#include <time.h>

/*set_params: sets values for the constants (permittivity of free space,
 *            electron charge, electron mass, etc.) as well as the 
 *            particle distribution function and distribution function
 *            parameters. 
 *
 *@params: pointer to struct of parameters *p
 * 
 *@returns: 1 to indicate success //TODO: should we just make this return nothing? 
 */
int set_params(struct parameters *p)
{
	p->epsilon0  = 1./(4. * M_PI); //permittivity of free space, CGS units
	p->e         = 4.80320680e-10; //electron charge
	p->m         = 9.1093826e-28;  //electron mass
	p->c         = 2.99792458e10;  //speed of light
	p->epsilon   = -1.;            //sign of electron charge
	
	//parameters
	p->B       = 1.;          //background B strength
	p->n_e     = 1.;          //electron number density cm^-3
	p->theta   = M_PI/3.;     //observer angle

	//derived quantities
	p->omega_p = sqrt(p->n_e * p->e*p->e / (p->m * p->epsilon0));//plasma frequency    
        p->omega_c = p->e * p->B / (p->m * p->c);                 //cyclotron frequency

	//integrator parameters
	p->gamma             = 1.5; //will get reset later in integration
	p->real              = 1;   //real part = 1, imag part = 0

	//distribution function
	p->dist              = 0; //MJ=1, PL=2, kappa=3

	//distribution function parameters
	p->theta_e     = 10.;         //dimensionless electron temp
	p->pl_p        = 3.;          //power-law index, p
	p->gamma_min   = 1.;          //power-law gamma_min
	p->gamma_max   = 1000.;       //power-law gamma_max
	p->kappa       = 3.5;         //kappa index
	p->kappa_width = 10.;         //kappa width, like theta_e
	p->gamma_cutoff = 1e10;       //currently unused

	return 1;
}

double chi_11_symphony(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
            int polarization,
            double theta_e,
            double power_law_p,
            double gamma_min,
            double gamma_max,
            double gamma_cutoff,
            double kappa,
            double kappa_width,
            char **error_message
           )
{
//  gsl_error_handler_t *prev_handler;
  double retval;

/*fill the struct with values*/
  struct parameters params;
  setConstParams(&params);

  set_params(&params);
  params.omega = 2. * M_PI * nu;
  params.real  = 1;

  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.polarization       = polarization;
  params.mode               = params.EMISSIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

//  if (error_message != NULL)
//    *error_message = NULL; /* Initialize the user's error message. */
//
//  global_gsl_error_message = &params.error_message;
//  prev_handler = gsl_set_error_handler (_handle_gsl_error);
//  set_distribution_function(&params);

//  retval = n_summation(&params);
  retval = chi_11(&params);

//  gsl_set_error_handler (prev_handler);
//  global_gsl_error_message = NULL;

  /* Success? */

  if (params.error_message == NULL)
    return retval;

  /* Something went wrong. Give the caller the error message if they
   * provided us with a place to save it. */

  if (error_message != NULL)
    *error_message = params.error_message;

  return NAN;
}

/*alpha_I: returns the absorption coefficient alpha_I, for the total intensity
 *         of light along the ray in question, for the given values of 
 *         parameters within the struct p.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: absorption coefficient for total intensity (Stokes I) 
 */
double alpha_I(struct parameters *p)
{
	p->real          = 0;
        double prefactor = 2. * M_PI * p->epsilon0 * p->omega / p->c;
        double term11    = (chi_11(p) * pow(cos(p->theta), 2.)  
			  + chi_33(p) * pow(sin(p->theta), 2.)
			  - 2. * chi_13(p) * sin(p->theta) * cos(p->theta));
	double term22    = chi_22(p);
        double ans       = prefactor * (term11 + term22);
        return ans;
}

/*alpha_Q: returns the absorption coefficient alpha_Q, for linearly polarized
 *         light along the ray in question, for the given values of parameters 
 *         within the struct p.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: absorption coefficient for linearly polarized light (Stokes Q) 
 */
double alpha_Q(struct parameters *p)
{
        p->real          = 0;
        double prefactor = 2. * M_PI * p->epsilon0 * p->omega / p->c;
        double term11    = (chi_11(p) * pow(cos(p->theta), 2.)
                          + chi_33(p) * pow(sin(p->theta), 2.)
                          - 2. * chi_13(p) * sin(p->theta) * cos(p->theta));
        double term22    = chi_22(p);
        double ans       = prefactor * (term11 - term22);
        return ans;
}

/*rho_Q: returns the Faraday conversion coefficient rho_Q, which corresponds
 *       to the conversion between linearly polarized and circularly
 *       polarized light by the medium.  The coefficient is calculated for
 *       the given values of parameters within the struct p.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: Faraday conversion coefficient rho_Q 
 */
double rho_Q(struct parameters *p)
{
        p->real          = 1;
        double prefactor = 2. * M_PI * p->epsilon0 * p->omega / p->c;
        double term11    = (chi_11(p) * pow(cos(p->theta), 2.)
                          + chi_33(p) * pow(sin(p->theta), 2.)
                          - 2. * chi_13(p) * sin(p->theta) * cos(p->theta));
        double term22    = chi_22(p);
        double ans       = prefactor * (term22 - term11);
        return ans;
}

/*alpha_V: returns the absorption coefficient alpha_V, for the circularly
 *         polarized light along the ray in question, for the given values of 
 *         parameters within the struct p.  Uses the IEEE/IAU convention for
 *         the sign of Stokes V.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: absorption coefficient for circularly polarized light (Stokes V) 
 */
double alpha_V(struct parameters *p)
{
	p->real            = 1;
	double prefactor   = 4. * M_PI * p->epsilon0 * p->omega / p->c;
	double term1     = (chi_12(p) * cos(p->theta) - chi_32(p) * sin(p->theta));
	double ans       = prefactor * term1;
	
        printf("\nprefactor: %e term1:  %e ans: %e", prefactor, term1, ans);

	return ans;
}

/*rho_V: returns the Faraday rotation coefficient rho_V, which rotates the
 *       plane of polarization (EVPA) for linearly polarized light,
 *       for the given values of parameters within the struct p.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: Faraday rotation coefficient rho_V 
 */
double rho_V(struct parameters *p)
{
        p->real          = 0;
        double prefactor = 4. * M_PI * p->epsilon0 * p->omega / p->c;
        double term1     = (chi_12(p) * cos(p->theta) - chi_32(p) * sin(p->theta));
        double ans       = prefactor * term1;
        return ans;
}

/*plotter: prints the values of the gamma integrand for the component of chi_ij
 *         determined by p.tau_integrand, from gamma=start to gamma=end, in
 *         increments of step.  These values are printed to a file called
 *         output.txt, and can be plotted easily by an external plotting
 *         software to determine if the gamma integrand is being properly
 *         resolved.
 *
 *@params: struct of parameters p
 *
 *@returns: 0 when completed and prints the gamma integrand to a file for
 *          plotting //TODO: make this function return nothing? 
 */
double plotter(struct parameters p)
{
	FILE *fp;
	fp = fopen("output.txt", "w");

	double start = 1.;
	double end   = 1.01;
	double i     = start;
	double step  = 0.00001;

	p.tau_integrand = &chi_12_integrand;

        while(i < end)
        {

		fprintf(fp, "\n%e    %e", i, tau_integrator(i, &p));
		printf("\n%e", i);
                i = i + step;
        }
        printf("\n");

	return 0.;
}

/*main: sets parameters, runs some calculation, and prints the CPU time elapsed
 *
 *@params: none
 *
 *@returns: nothing
 */
int main(void)
{
	/*start timer*/
	clock_t start = clock(), diff;
        struct parameters p;

	/*set parameters*/
	set_params(&p);
	p.omega = 1. * p.omega_c;
	p.real  = 1;

	/*print gamma	gamma_integrand(gamma) with the function plotter(params)*/
//	plotter(p);

	/*print omega/omega_c	alpha_I(params)*/
	printf("\n%e    %e\n", p.omega/p.omega_c, alpha_V(&p));

	/*calculate and print elapsed time*/
	diff = clock() - start;
	int msec = diff * 1000 / CLOCKS_PER_SEC;
	printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
}
