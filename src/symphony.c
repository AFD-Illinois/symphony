#include "symphony.h"

/* GSL error handling. This isn't thread-safe since we have to use a single
 * global variable to track the current best error message destination.
 * Fortunately-ish Python isn't either ... */

void
_symphony_error_trap (char *message)
{
    /* This function does nothing, but you can set a breakpoint on it to enter
     * a debugger when something bad happens in a calculation. */
}

static char **global_gsl_error_message = NULL;

static void
_handle_gsl_error (const char *reason, const char *file, int line, int gsl_errno)
{
    const size_t buf_size = 4096; /* arbitrary */

    if (global_gsl_error_message == NULL) {
       fprintf (stderr, "unhandled GSL error: %s (%d; %s:%d)\n", reason,
                gsl_errno, file, line);
       return;
    }

    *global_gsl_error_message = (char *) calloc (buf_size, 1);
    snprintf (*global_gsl_error_message, buf_size - 1,
	      "GSL error: %s (%d; %s:%d)", reason, gsl_errno, file, line);
    _symphony_error_trap (*global_gsl_error_message);
}


/*j_nu: wrapper for the emissivity calculation; takes in values of all
 *      necessary paramters and sets a struct of parameters using the input
 *      values.  It then passes this struct to n_summation(), which begins
 *      the emissivity calculation.
 *
 *@params: nu, magnetic_field, electron_density, observer_angle,
 *         distribution, polarization, theta_e, power_law_p,
 *         gamma_min, gamma_max, gamma_cutoff, kappa,
 *         kappa_width
 *@returns: n_summation(&params), which takes the struct of
 *          parameters (now populated with values) and
 *          performs the integration to evaluate j_nu(). If an
 *          error occurred and error_message is not NULL,
 *          *error_message will be set to a malloc()ed string
 *          explaining the error.
 */
double j_nu(double nu,
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
  gsl_error_handler_t *prev_handler;
  double retval;

/*fill the struct with values*/
  struct parameters params;
  setConstParams(&params);
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

  if (error_message != NULL)
    *error_message = NULL; /* Initialize the user's error message. */

  global_gsl_error_message = &params.error_message;
  prev_handler = gsl_set_error_handler (_handle_gsl_error);
  set_distribution_function(&params);
  retval = n_summation(&params);
  gsl_set_error_handler (prev_handler);
  global_gsl_error_message = NULL;

  /* Success? */
  if (params.error_message == NULL)
    return retval;

  /* Something went wrong. Give the caller the error message if they
   * provided us with a place to save it. */

  if (error_message != NULL)
    *error_message = params.error_message;

  return NAN;
}

/*alpha_nu: wrapper for the absorptivity calculation; takes in values of all
 *          necessary parameters and sets a struct of parameters using the input
 *          values.  It then passes this struct to n_summation(), which begins
 *          the absorptivity calculation.
 *
 *@params: nu, magnetic_field, electron_density, observer_angle,
 *         distribution, polarization, theta_e, power_law_p,
 *         gamma_min, gamma_max, gamma_cutoff, kappa,
 *         kappa_width
 *@returns: n_summation(&params), which takes the struct of
 *          parameters (now populated with values) and
 *          performs the integration to evaluate alpha_nu().
 */
double alpha_nu(double nu,
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
		int chi_method,
		char **error_message
               )
{
  gsl_error_handler_t *prev_handler;
  double retval;

/*fill the struct with values*/
  struct parameters params;
  setConstParams(&params);
  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.polarization       = polarization;
  params.mode               = params.ABSORPTIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

  if (error_message != NULL)
    *error_message = NULL; /* Initialize the user's error message. */

  global_gsl_error_message = &params.error_message;
  prev_handler = gsl_set_error_handler (_handle_gsl_error);
  set_distribution_function(&params);

  /* Choose method to compute alpha_nu */
  if(chi_method == params.SUSCEPT_METHOD)
  {
    params.omega   = 2. * params.pi * params.nu;
    params.omega_c = 2. * params.pi * get_nu_c(params);
    params.omega_p = get_omega_p(params);

    retval = alpha_nu_suscept(&params);
  }
  else
  {
    retval = n_summation(&params);
  }
  gsl_set_error_handler (prev_handler);
  global_gsl_error_message = NULL;

  /* Success? */

  if (params.error_message == NULL)
    return retval;

  /* Something went wrong. Give the caller the error message if they
   * provided us with a place to save it. */

  if (error_message != NULL)
    *error_message = params.error_message;

  return NAN;
}

/*rho_nu: wrapper for the rotativities calculation; takes in values of all
 *          necessary parameters and sets a struct of parameters using the input
 *          values.  It then passes this struct to n_summation(), which begins
 *          the absorptivity calculation.
 *
 *@params: nu, magnetic_field, electron_density, observer_angle,
 *         distribution, polarization, theta_e, power_law_p,
 *         gamma_min, gamma_max, gamma_cutoff, kappa,
 *         kappa_width
 *@returns: n_summation(&params), which takes the struct of
 *          parameters (now populated with values) and
 *          performs the integration to evaluate alpha_nu().
 */
double rho_nu(double nu,
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
	      char **error_message)
{
  gsl_error_handler_t *prev_handler;
  double retval;

/*fill the struct with values*/
  struct parameters params;
  setConstParams(&params);
  params.nu                 = nu;
  params.magnetic_field     = magnetic_field;
  params.observer_angle     = observer_angle;
  params.electron_density   = electron_density;
  params.distribution       = distribution;
  params.polarization       = polarization;
  params.mode               = params.ABSORPTIVITY;
  params.theta_e            = theta_e;
  params.power_law_p        = power_law_p;
  params.gamma_min          = gamma_min;
  params.gamma_max          = gamma_max;
  params.gamma_cutoff       = gamma_cutoff;
  params.kappa              = kappa;
  params.kappa_width        = kappa_width;

  if (error_message != NULL)
    *error_message = NULL; /* Initialize the user's error message. */

  global_gsl_error_message = &params.error_message;
  prev_handler = gsl_set_error_handler (_handle_gsl_error);
  set_distribution_function(&params);

  params.omega   = 2. * params.pi * params.nu;
  params.omega_c = 2. * params.pi * get_nu_c(params);
  params.omega_p = get_omega_p(params);

  retval = rho_nu_suscept(&params);

  gsl_set_error_handler (prev_handler);
  global_gsl_error_message = NULL;

  /* Success? */

  if (params.error_message == NULL)
    return retval;

  /* Something went wrong. Give the caller the error message if they
   * provided us with a place to save it. */

  if (error_message != NULL)
    *error_message = params.error_message;

  return NAN;
}
