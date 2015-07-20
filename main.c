/* Symphony
 * by Alex Pandya, Zhaowei Zhang
 * 6/30/15
 *
 *References:
 * 1) Leung, Gammie, and Noble (2011)
 * 2) Xiao (2006)
*/

#include "dec.h"

/* main: defines nu_c (cyclotron frequency) and loops through values of nu, 
 * to give output absorptivity or emissivity vs. nu/nu_c; can also be modified
 * to give abs/emiss as a function of its other parameters, like obs. angle
 *
 * The distribution function, polarization mode, and emissivity/absorptivity
 * can all be set in the file dec.h
 */

int main(int argc, char *argv[]) 
{

  double nu_c = (electron_charge * B_field)
 	       /(2. * M_PI * mass_electron * speed_light);

  printf("\nDIST: %d, MODE: %d, POL: %d", DISTRIBUTION_FUNCTION,
		 	                  MODE, POL_MODE);
  double max_nuratio = 1e10;
  int points_per_pow_10 = 1;
  int max_index = (int) log10(max_nuratio)*points_per_pow_10;

  for (int index=0; index <= max_index; index++) 
  {

    double nu = pow(10., (double)index/(double)points_per_pow_10) * nu_c;
    printf("\n%e	%e", nu/nu_c, n_summation(nu));

  }
  printf("\n");
  return 0;
}
