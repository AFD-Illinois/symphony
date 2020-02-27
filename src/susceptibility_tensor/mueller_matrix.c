#include <stdio.h>
#include <math.h>
#include "susceptibility_tensor.h"

double chi_11_symphony(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
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

//  set_params(&params);
  params.omega = 2. * params.pi * nu;
  params.real  = real_part;

  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.mode               = params.EMISSIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

  double nu_c = get_nu_c(params);
  double omega_p = get_omega_p(params);

  params.omega_c = 2. * params.pi * nu_c;
  params.omega_p = omega_p;

  if (error_message != NULL)
    *error_message = NULL; /* Initialize the user's error message. */

//  global_gsl_error_message = &params.error_message;
//  prev_handler = gsl_set_error_handler (_handle_gsl_error);
  
  retval = chi_11(&params);

//  gsl_set_error_handler (prev_handler);
//  global_gsl_error_message = NULL;

  /* Success? */

  if (params.error_message == NULL)
    return retval;

  /* Something went wrong. Give the caller the error message if they
   * provided us with a place to save it. */

//  if (error_message != NULL)
//    *error_message = params.error_message;

  return NAN;
}

double chi_12_symphony(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
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

//  set_params(&params);
  params.omega = 2. * params.pi * nu;
  params.real  = real_part;

  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.mode               = params.EMISSIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

  double nu_c = get_nu_c(params);
  double omega_p = get_omega_p(params);

  params.omega_c = 2. * params.pi * nu_c;
  params.omega_p = omega_p;

  if (error_message != NULL)
    *error_message = NULL; /* Initialize the user's error message. */

//  global_gsl_error_message = &params.error_message;
//  prev_handler = gsl_set_error_handler (_handle_gsl_error);
  
  retval = chi_12(&params);

//  gsl_set_error_handler (prev_handler);
//  global_gsl_error_message = NULL;

  /* Success? */

  if (params.error_message == NULL)
    return retval;

  /* Something went wrong. Give the caller the error message if they
   * provided us with a place to save it. */

//  if (error_message != NULL)
//    *error_message = params.error_message;

  return NAN;
}

double chi_13_symphony(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
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

//  set_params(&params);
  params.omega = 2. * params.pi * nu;
  params.real  = real_part;

  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.mode               = params.EMISSIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

  double nu_c = get_nu_c(params);
  double omega_p = get_omega_p(params);

  params.omega_c = 2. * params.pi * nu_c;
  params.omega_p = omega_p;

  if (error_message != NULL)
    *error_message = NULL; /* Initialize the user's error message. */

//  global_gsl_error_message = &params.error_message;
//  prev_handler = gsl_set_error_handler (_handle_gsl_error);
  
  retval = chi_13(&params);

//  gsl_set_error_handler (prev_handler);
//  global_gsl_error_message = NULL;

  /* Success? */

  if (params.error_message == NULL)
    return retval;

  /* Something went wrong. Give the caller the error message if they
   * provided us with a place to save it. */

//  if (error_message != NULL)
//    *error_message = params.error_message;

  return NAN;
}

double chi_22_symphony(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
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

//  set_params(&params);
  params.omega = 2. * params.pi * nu;
  params.real  = real_part;

  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.mode               = params.EMISSIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

  double nu_c = get_nu_c(params);
  double omega_p = get_omega_p(params);

  params.omega_c = 2. * params.pi * nu_c;
  params.omega_p = omega_p;

  if (error_message != NULL)
    *error_message = NULL; /* Initialize the user's error message. */

//  global_gsl_error_message = &params.error_message;
//  prev_handler = gsl_set_error_handler (_handle_gsl_error);
  
  retval = chi_22(&params);

//  gsl_set_error_handler (prev_handler);
//  global_gsl_error_message = NULL;

  /* Success? */

  if (params.error_message == NULL)
    return retval;

  /* Something went wrong. Give the caller the error message if they
   * provided us with a place to save it. */

//  if (error_message != NULL)
//    *error_message = params.error_message;

  return NAN;
}

double chi_32_symphony(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
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

//  set_params(&params);
  params.omega = 2. * params.pi * nu;
  params.real  = real_part;

  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.mode               = params.EMISSIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

  double nu_c = get_nu_c(params);
  double omega_p = get_omega_p(params);

  params.omega_c = 2. * params.pi * nu_c;
  params.omega_p = omega_p;

  if (error_message != NULL)
    *error_message = NULL; /* Initialize the user's error message. */

//  global_gsl_error_message = &params.error_message;
//  prev_handler = gsl_set_error_handler (_handle_gsl_error);
  
  retval = chi_32(&params);

//  gsl_set_error_handler (prev_handler);
//  global_gsl_error_message = NULL;

  /* Success? */

  if (params.error_message == NULL)
    return retval;

  /* Something went wrong. Give the caller the error message if they
   * provided us with a place to save it. */

//  if (error_message != NULL)
//    *error_message = params.error_message;

  return NAN;
}

double chi_33_symphony(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
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

//  set_params(&params);
  params.omega = 2. * params.pi * nu;
  params.real  = real_part;

  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.mode               = params.EMISSIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

  double nu_c = get_nu_c(params);
  double omega_p = get_omega_p(params);

  params.omega_c = 2. * params.pi * nu_c;
  params.omega_p = omega_p;

  if (error_message != NULL)
    *error_message = NULL; /* Initialize the user's error message. */

//  global_gsl_error_message = &params.error_message;
//  prev_handler = gsl_set_error_handler (_handle_gsl_error);
  
  retval = chi_33(&params);

//  gsl_set_error_handler (prev_handler);
//  global_gsl_error_message = NULL;

  /* Success? */

  if (params.error_message == NULL)
    return retval;

  /* Something went wrong. Give the caller the error message if they
   * provided us with a place to save it. */

//  if (error_message != NULL)
//    *error_message = params.error_message;

  return NAN;
}

double chi_21_symphony(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
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

//  set_params(&params);
  params.omega = 2. * params.pi * nu;
  params.real  = real_part;

  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.mode               = params.EMISSIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

  double nu_c = get_nu_c(params);
  double omega_p = get_omega_p(params);

  params.omega_c = 2. * params.pi * nu_c;
  params.omega_p = omega_p;

  if (error_message != NULL)
    *error_message = NULL; /* Initialize the user's error message. */

//  global_gsl_error_message = &params.error_message;
//  prev_handler = gsl_set_error_handler (_handle_gsl_error);
  
  retval = chi_21(&params);

//  gsl_set_error_handler (prev_handler);
//  global_gsl_error_message = NULL;

  /* Success? */

  if (params.error_message == NULL)
    return retval;

  /* Something went wrong. Give the caller the error message if they
   * provided us with a place to save it. */

//  if (error_message != NULL)
//    *error_message = params.error_message;

  return NAN;
}

double chi_23_symphony(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
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

//  set_params(&params);
  params.omega = 2. * params.pi * nu;
  params.real  = real_part;

  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.mode               = params.EMISSIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

  double nu_c = get_nu_c(params);
  double omega_p = get_omega_p(params);

  params.omega_c = 2. * params.pi * nu_c;
  params.omega_p = omega_p;

  if (error_message != NULL)
    *error_message = NULL; /* Initialize the user's error message. */

//  global_gsl_error_message = &params.error_message;
//  prev_handler = gsl_set_error_handler (_handle_gsl_error);
  
  retval = chi_21(&params);

//  gsl_set_error_handler (prev_handler);
//  global_gsl_error_message = NULL;

  /* Success? */

  if (params.error_message == NULL)
    return retval;

  /* Something went wrong. Give the caller the error message if they
   * provided us with a place to save it. */

//  if (error_message != NULL)
//    *error_message = params.error_message;

  return NAN;
}

double chi_31_symphony(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
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

//  set_params(&params);
  params.omega = 2. * params.pi * nu;
  params.real  = real_part;

  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.mode               = params.EMISSIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

  double nu_c = get_nu_c(params);
  double omega_p = get_omega_p(params);

  params.omega_c = 2. * params.pi * nu_c;
  params.omega_p = omega_p;

  if (error_message != NULL)
    *error_message = NULL; /* Initialize the user's error message. */

//  global_gsl_error_message = &params.error_message;
//  prev_handler = gsl_set_error_handler (_handle_gsl_error);
  
  retval = chi_31(&params);

//  gsl_set_error_handler (prev_handler);
//  global_gsl_error_message = NULL;

  /* Success? */

  if (params.error_message == NULL)
    return retval;

  /* Something went wrong. Give the caller the error message if they
   * provided us with a place to save it. */

//  if (error_message != NULL)
//    *error_message = params.error_message;

  return NAN;
}

//TODO: write documentation for this function
double alpha_nu_suscept(struct parameters * params)
{
  /*this method has only been tested for 1 <= nu/nu_c <= 1000*/
  double nu_c = get_nu_c(*params);
  if(params->nu/nu_c >= 1000)
  {
    printf("\nFrequency out of range for this method");
    return 0.;
  }

  if(params->polarization == params->STOKES_I)
  {
    return alpha_I(params);
  }
  else if(params->polarization == params->STOKES_Q)
  {
    return alpha_Q(params);
  }
  else if(params->polarization == params->STOKES_U)
  {
    return 0.; /*alpha_U is exactly zero in this basis; see paper*/
  }
  else
  {
    return alpha_V(params);
  }

  printf("\nERROR IN alpha_nu_suscept()\n");
  return 0.;
}

//TODO: write documentation for this function
double rho_nu_suscept(struct parameters * params)
{
  /*this method has only been tested for 1 <= nu/nu_c <= 1000*/
  double nu_c = get_nu_c(*params);
  if(params->nu/nu_c >= 20000)
  {
    printf("\nFrequency out of range for this method");
    return 0.;
  }

  if(params->polarization == params->STOKES_I)
  {
    printf("The coefficient rho_I does not exist.");
    return 0.;
  }
  else if(params->polarization == params->STOKES_Q)
  {
    return rho_Q(params);
  }
  else if(params->polarization == params->STOKES_U)
  {
    return 0.; /*rho_U is exactly zero in this basis; see paper*/
  }
  else
  {
    return rho_V(params);
  }

  printf("\nERROR IN alpha_nu_suscept()\n");
  return 0.;
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
  double prefactor = 2. * p->pi * p->epsilon0 * p->omega / p->speed_light;
  double term11    = (chi_11(p) * pow(cos(p->observer_angle), 2.)  
  		    + chi_33(p) * pow(sin(p->observer_angle), 2.)
  		  - 2.*chi_13(p)*sin(p->observer_angle)*cos(p->observer_angle));
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
  double prefactor = 2. * p->pi * p->epsilon0 * p->omega / p->speed_light;
  double term11    = (chi_11(p) * pow(cos(p->observer_angle), 2.)
                   + chi_33(p) * pow(sin(p->observer_angle), 2.)
                   -2.*chi_13(p)*sin(p->observer_angle)*cos(p->observer_angle));
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
  double prefactor = 2. * p->pi * p->epsilon0 * p->omega / p->speed_light;
  double term11    =(chi_11(p) * pow(cos(p->observer_angle), 2.)
                   +chi_33(p) * pow(sin(p->observer_angle), 2.)
                   -2.*chi_13(p)*sin(p->observer_angle)*cos(p->observer_angle));
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
  p->real          = 1;
  double prefactor = 4. * p->pi * p->epsilon0 * p->omega / p->speed_light;
  double term1     = (chi_12(p) * cos(p->observer_angle) 
                      - chi_32(p) * sin(p->observer_angle));
  double ans       = prefactor * term1;
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
  double prefactor = 4. * p->pi * p->epsilon0 * p->omega / p->speed_light;
  double term1     = (chi_12(p) * cos(p->observer_angle) 
                      - chi_32(p) * sin(p->observer_angle));
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
 *          plotting
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
