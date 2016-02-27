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
  /* Keys for the distributions */
  params->MAXWELL_JUETTNER = 0;
  params->POWER_LAW        = 1;
  params->KAPPA_DIST       = 2;
  /* Keys for the polarization parameter */
  params->STOKES_I         = 15;
  params->STOKES_Q         = 16;
  params->STOKES_U         = 17;
  params->STOKES_V         = 18;
  /* Keys for the mode: absorptivity or emissivity */
  params->ABSORPTIVITY     = 10;
  params->EMISSIVITY       = 11;

}

