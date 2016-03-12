/* Symphony
 * by Alex Pandya, Zhaowei Zhang
 * 2/26/2016
 *
 *References:
 * 1) Pandya, Zhang, Chandra, Gammie (Accepted to ApJ 2016)
 * 2) Leung, Gammie, and Noble (2011)
 * 3) Xiao (2006)
 */

#include "symphony.h"
#include "fits.h"
#include "params.h"
#include "distribution_function_common_routines.h"

/* main: runs through a demo of symphony's capabilities; currently
 * sets parameters to sample values and outputs nu/nu_c 
 * (where nu_c is the cyclotron frequency) vs. j_nu(). 
 */

int main(int argc, char *argv[]) 
{
  struct parameters paramsM;
  setConstParams(&paramsM);

//set some the parameters
  paramsM.magnetic_field     = 30.;
  paramsM.electron_density   = 1.;
  paramsM.observer_angle     = paramsM.pi/3.;
  paramsM.distribution       = paramsM.MAXWELL_JUETTNER;
  paramsM.polarization       = paramsM.STOKES_I;
  paramsM.theta_e            = 10.;
  paramsM.power_law_p        = 3.5;
  paramsM.gamma_min          = 1.;
  paramsM.gamma_max          = 1000.;
  paramsM.gamma_cutoff       = 10000000000.;
  paramsM.kappa              = 3.5;
  paramsM.kappa_width        = 10.;

  double nu_c = get_nu_c(paramsM);

  double max_nuratio = 1e6;
  int points_per_pow_10 = 1;
  int max_index = (int) log10(max_nuratio)*points_per_pow_10;

//  printf("\nnu/nu_c         j_nu()          j_nu_fit()");

  for (int index=0; index <= max_index; index++) 
  {

    paramsM.nu = pow(10., (double)index/(double)points_per_pow_10) * nu_c;

    printf("\n%e	%e	%e", paramsM.nu/nu_c, 
                                     alpha_nu(paramsM.nu, 
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
                                 alpha_nu_fit(paramsM.nu,
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

