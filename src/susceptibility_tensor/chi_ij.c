#include <stdio.h>
#include <math.h>
#include "susceptibility_tensor.h"

/*chi_11: evaluates the component chi_11 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: the component chi_11, evaluated using the parameters in struct p
 */
double chi_11(struct parameters * p)
{
  p->tau_integrand = &chi_11_integrand;
  p->gamma_integrand = &tau_integrator;
  
  return gamma_integrator(p);
}

/*chi_12: evaluates the component chi_12 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: the component chi_12, evaluated using the parameters in struct p
 */
double chi_12(struct parameters * p)
{
  p->tau_integrand = &chi_12_integrand;
  p->gamma_integrand = &tau_integrator;
  
  return gamma_integrator(p);
}

/*chi_32: evaluates the component chi_32 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: the component chi_32, evaluated using the parameters in struct p
 */
double chi_32(struct parameters * p)
{
  p->tau_integrand = &chi_32_integrand;
  p->gamma_integrand = &tau_integrator;
  
  return gamma_integrator(p);
}

/*chi_13: evaluates the component chi_13 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: the component chi_13, evaluated using the parameters in struct p
 */
double chi_13(struct parameters * p)
{
  p->tau_integrand = &chi_13_integrand;
  p->gamma_integrand = &tau_integrator;
  
  return gamma_integrator(p);
}

/*chi_22: evaluates the component chi_22 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: the component chi_22, evaluated using the parameters in struct p
 */
double chi_22(struct parameters * p)
{
  double ans = 0.;
  p->gamma_integrand = &tau_integrator;
  
  if(p->real == 0)
  {
    p->tau_integrand = &chi_22_integrand_p1;
    ans = gamma_integrator(p);
    p->tau_integrand = &chi_22_integrand_p2;
    ans += gamma_integrator(p);
  }
  else
  {
    p->tau_integrand = &chi_22_integrand_real;
    ans = gamma_integrator(p);
  }
  
  return ans;
}

/*chi_33: evaluates the component chi_33 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: the component chi_33, evaluated using the parameters in struct p
 */
double chi_33(struct parameters * p)
{
  p->tau_integrand = &chi_33_integrand;
  p->gamma_integrand = &tau_integrator;
  
  return gamma_integrator(p);
}

/*chi_21: evaluates the component chi_21 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: the component chi_21, evaluated using the parameters in struct p
 */
double chi_21(struct parameters * p)
{
  return -chi_12(p);
}

/*chi_23: evaluates the component chi_23 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: the component chi_23, evaluated using the parameters in struct p
 */
double chi_23(struct parameters * p)
{
  return -chi_32(p);
}

/*chi_31: evaluates the component chi_31 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: the component chi_31, evaluated using the parameters in struct p
 */
double chi_31(struct parameters * p)
{
  return chi_13(p);
}

/*chi_rho_Q: evaluates the combination of the four chi components of the susceptibility tensor needed to evaluate rho_Q.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: the value of the four combined chi components, evaluated using the parameters in struct p
 */
double chi_rho_Q(struct parameters * p){

  p->tau_integrand = &chi_rho_Q_integrand;
  p->gamma_integrand = &tau_integrator;

  return gamma_integrator(p);
}

