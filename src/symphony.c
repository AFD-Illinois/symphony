/* Symphony
 * by Alex Pandya, Zhaowei Zhang
 * 6/30/15
 *
 *References:
 * 1) Leung, Gammie, and Noble (2011)
 * 2) Xiao (2006)
*/

#include "symphony.h"

/* main: defines nu_c (cyclotron frequency) and loops through values of nu, 
 * to give output absorptivity or emissivity vs. nu/nu_c; can also be modified
 * to give abs/emiss as a function of its other parameters, like obs. angle
 *
 * The distribution function, polarization mode, and emissivity/absorptivity
 * must all be set in the file symphony.h; parameters can be set in params.h
 */

int main(int argc, char *argv[]) 
{
  struct parameters paramsM;
  setConstParams(&paramsM);

//set some the parameters
  paramsM.magnetic_field     = 30.;
  paramsM.electron_density   = 1.;
  paramsM.observer_angle     = paramsM.pi/3.;
  paramsM.distribution       = paramsM.THERMAL;
  paramsM.polarization       = paramsM.STOKES_I;
  paramsM.theta_e            = 10.;
  paramsM.power_law_p        = 3.5;
  paramsM.gamma_min          = 1.;
  paramsM.gamma_max          = 1000.;
  paramsM.gamma_cutoff       = 10000000000.;
  paramsM.kappa              = 3.5;
  paramsM.kappa_width        = 10.;

  double nu_c = get_nu_c(paramsM);

//  printf("\nDIST: %d, MODE: %d, POL: %d", DISTRIBUTION_FUNCTION,
//		 	                  MODE, POL_MODE);

  double max_nuratio = 1e6;
  int points_per_pow_10 = 1;
  int max_index = (int) log10(max_nuratio)*points_per_pow_10;

  printf("\nnu/nu_c         ans             fit");

  for (int index=0; index <= max_index; index++) 
  {

    paramsM.nu = pow(10., (double)index/(double)points_per_pow_10) * nu_c;

    printf("\n%e	%e	%e", paramsM.nu/nu_c, 
                                     j_nu(paramsM.nu, 
                                          paramsM.magnetic_field, 
                                          paramsM.electron_density, 
                                          paramsM.observer_angle, 
                                          paramsM.distribution, 
                                          paramsM.polarization,
                                          paramsM.theta_e,            
                                          paramsM.power_law_p,        
                                          paramsM.gamma_min,          
                                          paramsM.gamma_max,        
                                          paramsM.gamma_cutoff,      
                                          paramsM.kappa,              
                                          paramsM.kappa_width        
                                          ), 
                                 j_nu_fit(paramsM.nu,
                                          paramsM.magnetic_field,
                                          paramsM.electron_density,
                                          paramsM.observer_angle,
                                          paramsM.distribution,
                                          paramsM.polarization,
                                          paramsM.theta_e,           
                                          paramsM.power_law_p,       
                                          paramsM.gamma_min,         
                                          paramsM.gamma_max,       
                                          paramsM.gamma_cutoff,     
                                          paramsM.kappa,             
                                          paramsM.kappa_width 
                                          ));

  }
  printf("\n");
  return 0;
}

