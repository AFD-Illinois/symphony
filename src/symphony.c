#include "symphony.h"

/*j_nu: wrapper for the emissivity calculation; takes in values of all
 *      necessary paramters and sets a struct of parameters using the input
 *      values.  It then passes this struct to n_summation(), which begins
 *      the emissivity calculation.
 *
 *@params: nu, magnetic_field, electron_density, observer_angle,
 *         distribution, polarization, theta_e, power_law_p,
 *         gamma_min, gamma_max, gamma_cutoff, kappa, 
 *         kappa_width
 *
 *@returns: n_summation(&params), which takes the struct of 
 *          parameters (now populated with values) and 
 *          performs the integration to evaluate j_nu().
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
            double kappa_width
           )
{
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
  set_distribution_function(&params);


  return n_summation(&params);
}

/*alpha_nu: wrapper for the absorptivity calculation; takes in values of all
 *          necessary paramters and sets a struct of parameters using the input
 *          values.  It then passes this struct to n_summation(), which begins
 *          the absorptivity calculation.
 *
 *@params: nu, magnetic_field, electron_density, observer_angle,
 *         distribution, polarization, theta_e, power_law_p,
 *         gamma_min, gamma_max, gamma_cutoff, kappa, 
 *         kappa_width
 *
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
                double kappa_width
               )
{
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

  return n_summation(&params);
}

