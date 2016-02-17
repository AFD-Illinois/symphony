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
//  struct parameters params;
//  setConstParams(&params);
//  setUserParams(&params);

//set some the parameters
  double nu = 230e5;
  double magnetic_field = 30.;
  double electron_density = 1.;
  double observer_angle = 1.0472;
  int distribution = 0;
  int polarization = 15;

  //double nu_c = get_nu_c(params);

//  printf("\nDIST: %d, MODE: %d, POL: %d", DISTRIBUTION_FUNCTION,
//		 	                  MODE, POL_MODE);

  double max_nuratio = 1e0;
  int points_per_pow_10 = 1;
  int max_index = (int) log10(max_nuratio)*points_per_pow_10;

//  printf("\nnu/nu_c         ans             fit");

//  for (int index=0; index <= max_index; index++) 
//  {

    //double nu = pow(10., (double)index/(double)points_per_pow_10) * nu_c;
    //nu = 5. * index * nu;

//    printf("\n%e	%e	%e", nu, 
                               j_nu    (nu, magnetic_field, electron_density, 
                                        observer_angle, distribution, 
                                        polarization);//, 
//                               j_nu_fit(nu, magnetic_field, electron_density, 
//                                        observer_angle, distribution, 
//                                        polarization));

//    printf("\n%e\n", j_nu_fit(nu, magnetic_field, electron_density,
//                              observer_angle, distribution, polarization));



//   printf("\n%e	%e	%e", nu/nu_c, alpha_nu(nu, B, n_e, obs_angle, dist, pol), 
//                               alpha_nu_fit(nu, B, n_e, obs_angle, dist, pol));

//  }
  printf("\n");
  return 0;
}
