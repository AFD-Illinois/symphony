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
 * can all be set in the file dec.h
 */

int main(int argc, char *argv[]) 
{
  double B = 30.;
  double n_e = 1.;
  double obs_angle = 60.*M_PI/180.;

  double nu_c = (electron_charge * B)
 	       /(2. * M_PI * mass_electron * speed_light);

  printf("\nDIST: %d, MODE: %d, POL: %d", DISTRIBUTION_FUNCTION,
		 	                  MODE, POL_MODE);

  double max_nuratio = 1e6;
  int points_per_pow_10 = 1;
  int max_index = (int) log10(max_nuratio)*points_per_pow_10;

  printf("\nnu/nu_c         ans             fit");

  for (int index=0; index <= max_index; index++) 
  {

    double nu = pow(10., (double)index/(double)points_per_pow_10) * nu_c;

    #if MODE == EMISS
    printf("\n%e	%e	%e", nu/nu_c, j_nu(nu, B, n_e, obs_angle), 
                                              j_nu_fit(nu, B, n_e, obs_angle, POL_MODE));

    #elif MODE == ABSORP
    printf("\n%e	%e	%e", nu/nu_c, alpha_nu(nu, B, n_e, obs_angle), 
                                              alpha_nu_fit(nu, B, n_e, obs_angle, POL_MODE));
    #endif

  }
  printf("\n");
  return 0;
}
