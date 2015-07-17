/* Symphony
 * by Alex Pandya, Zhaowei Zhang
 * 6/30/15
 *
 *References:
 * 1) Leung, Gammie, and Noble (2011)
 * 2) Xiao (2006)
*/

#include "dec.h"

//need 2 separate n integrations to numerically resolve STOKES_V
//static int stokes_v_switch = 0.;

/* main: defines nu_c (cyclotron frequency) and loops through values of nu, 
 * to give output absorptivity or emissivity vs. nu/nu_c; can also be modified
 * to give abs/emiss as a function of its other parameters, like obs. angle
 */


int main(int argc, char *argv[]) {

  double nu_c = (electron_charge * B_field)
 	       /(2. * M_PI * mass_electron * speed_light);

  int index = 0;
  printf("\nDIST: %d, MODE: %d, POL: %d", DISTRIBUTION_FUNCTION,
		 	                  MODE, POL_MODE);
  int num_points = 31;
  for (index; index <= num_points; index++) {

    double nu = pow(10., index/5.) * nu_c;
    printf("\n%e	%e", nu/nu_c, n_summation(nu));

  }
  printf("\n");
  return 0;

}

