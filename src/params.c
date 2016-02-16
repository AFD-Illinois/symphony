#include "params.h"

void setConstParams(struct parameters *params)
{
  params->pi               = 3.1415926535897932384;
  params->mass_electron    = 9.1093826e-28;
  params->plancks_constant = 6.6260693e-27;
  params->speed_light      = 2.99792458e10;
  params->electron_charge  = 4.80320680e-10;
  params->n_max            = 30;
  params->C                = 10;
  // Keys for the distributions
  params->THERMAL          = 0;
  params->POWER_LAW        = 1;
  params->KAPPA_DIST       = 2;
  // Keys for the polarization modes
  params->STOKES_I         = 15;
  params->STOKES_Q         = 16;
  params->STOKES_U         = 17;
  params->STOKES_V         = 18;
  // Keys for the mode: absorptivity or emmisivity
  params->ABSORPTIVITY     = 10;
  params->EMISSIVITY       = 11;

}

void setUserParams(struct parameters *params)
{
  // USER PARAMS SET BELOW:
  params->nu               = 230e9; // GHz
  params->magnetic_field   = 100;   // Gauss
  params->electron_density = 1;     // g/cc
  params->observer_angle   = params->pi/3;  // rad  
  params->distribution     = params->THERMAL;
  params->polarization     = params->STOKES_I;
  params->gamma_cutoff     = 1000000000000;

  // Thermal distribution parameters
  params->theta_e = 10.;

  //power law parameters
  params->power_law_p = 3.;
  params->gamma_min   = 1.;
  params->gamma_max   = 1000.;

  //kappa distribution parameters
  params->kappa       = 3.5;
  params->kappa_width = 10.;

  //need 2 separate n integrations to numerically resolve STOKES_V
  params->stokes_v_switch = 0.; // TODO: describe: 0: , 1: 

}
