#include "symphony.h"

/*fitting formulae*/

///*wrapper for emissivity fitting formulae; takes in arguments nu, B, n_e, 
//  and observer angle theta*/
double j_nu_fit(double nu,
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
                double kappa_width)
{
//fill the struct with values
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


//  check_for_errors(params); TODO: fix this
  
 if(params.distribution == params.MAXWELL_JUETTNER)
  {
    if     (params.polarization == params.STOKES_I) return maxwell_juettner_I(&params); 
    else if(params.polarization == params.STOKES_Q) return maxwell_juettner_Q(&params);
    else if(params.polarization == params.STOKES_U) return 0.;
    else if(params.polarization == params.STOKES_V) return maxwell_juettner_V(&params);
  }

  else if(params.distribution == params.POWER_LAW)
  {
    if     (params.polarization == params.STOKES_I) return power_law_I(&params);
    else if(params.polarization == params.STOKES_Q) return power_law_Q(&params);
    else if(params.polarization == params.STOKES_U) return 0.;
    else if(params.polarization == params.STOKES_V) return power_law_V(&params);
  }

  else if(params.distribution == params.KAPPA_DIST)
  {
    if     (params.polarization == params.STOKES_I) return kappa_I(&params);
    else if(params.polarization == params.STOKES_Q) return kappa_Q(&params);
    else if(params.polarization == params.STOKES_U) return 0.;
    else if(params.polarization == params.STOKES_V) return kappa_V(&params);
  }

  return 0.;
}

///*wrapper for the absorptivity fitting formulae; takes in nu, B, n_e, and 
//  observer angle theta*/
double alpha_nu_fit(double nu,
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
                    double kappa_width)
{
//fill the struct with values
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

//  check_for_errors(nu, B, n_e, obs_angle);

  if(params.distribution == params.MAXWELL_JUETTNER)
  {
    if     (params.polarization == params.STOKES_I) return maxwell_juettner_I_abs(&params);
    else if(params.polarization == params.STOKES_Q) return maxwell_juettner_Q_abs(&params);
    else if(params.polarization == params.STOKES_U) return 0.;
    else if(params.polarization == params.STOKES_V) return maxwell_juettner_V_abs(&params);
  }

  else if(params.distribution == params.POWER_LAW)
  {
    if     (params.polarization == params.STOKES_I) return power_law_I_abs(&params);
    else if(params.polarization == params.STOKES_Q) return power_law_Q_abs(&params);
    else if(params.polarization == params.STOKES_U) return 0.;
    else if(params.polarization == params.STOKES_V) return power_law_V_abs(&params);
  }

  else if(params.distribution == params.KAPPA_DIST)
  {
    if     (params.polarization == params.STOKES_I) return kappa_I_abs(&params);
    else if(params.polarization == params.STOKES_Q) return kappa_Q_abs(&params);
    else if(params.polarization == params.STOKES_U) return 0.;
    else if(params.polarization == params.STOKES_V) return kappa_V_abs(&params);
  }

  return 0.;
}

double check_for_errors(struct parameters * params)
{

  double nu_c = get_nu_c(*params);

  /* catch potential errors */
  if(params->nu/nu_c > 3e10)
  {
    printf("\n ERROR: nu out of range\n");
    exit(0);
  }
  if(params->magnetic_field < 0)
  {
    printf("\n ERROR: B out of range\n");
    exit(0);
  }
  if(params->electron_density < 0)
  {
    printf("\n ERROR: electron density out of range\n");
    exit(0);
  }
  if(params->kappa < 2.5 || params->kappa > 7.5)
  {
    printf("\n ERROR: kappa out of range of fitting formula\n");
    exit(0);
  }
  if(params->kappa_width < 3 || params->kappa_width > 200)
  {
    printf("\n WARNING: w out of range; fitting formula may be inaccurate\n");
    exit(0);
  }
  if(params->gamma_min < 1)
  {
    printf("\n ERROR: gamma_min < 1\n");
    exit(0);
  }
  if(params->observer_angle < 5.*(params->pi)/180. 
     || params->observer_angle == 90.*(params->pi)/180.)
  {
    printf("\n ERROR: theta out of range; fitting formula may be inaccurate\n");
    exit(0);
  }

  return 0.;
}
