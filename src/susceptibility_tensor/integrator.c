#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "susceptibility_tensor.h"

double tau_trapezoidal(struct parameters *p, double start, double end, int samples)
{
  int k                = 1;
  double tau           = 0;
  double eval          = 0;
  double tau_step      = (end - start)/samples;
  double tau_ans_step  = 0.;
  double osc_ans_tot   = (chi_rho_Q_integrand(start, p) + chi_rho_Q_integrand(end, p)) * tau_step / 2;
        
  while(start + (k * tau_step) < end)
  {
    eval = chi_rho_Q_integrand(start + (k * tau_step), p);
    tau_ans_step = tau_step * eval;
    osc_ans_tot += tau_ans_step;
    tau = start + (k * tau_step);
    k++;
  }
  
  return osc_ans_tot;
}



/*tau_integrator: integrator for the first integral (the tau integral) for the
 *                components of the susceptibility tensor.  The chi_33 integral
 *                is faster when we use a fixed-order Gaussian quadrature
 *                method GSL QNG, rather than the GSL QAWO adaptive integrator
 *                used for the other components.
 *
 *@params: double gamma (the second integration is carried out over this
 *         variable), void * parameters (a pointer to the struct of parameters)
 *
 *@returns: numerically evaluated tau integral at a given gamma of a given
 *          component of the susceptibility tensor
 */
double tau_integrator(double gamma, void * parameters)
{
  struct parameters * params = (struct parameters*) parameters;
  
  if(gamma == 1.)
  {
    return 0.;
  }
  
  double ans_tot  = 0.;
  double ans_step = 0.;
  double error    = 0.;
  double step;
  double start    = 0.;
  size_t n        = 50;
  size_t limit    = 5000;
  double epsabs   = 0.;
  double epsrel   = 1e-8;
  double sign_correction;
  enum gsl_integration_qawo_enum gsl_weight;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
  
  if(params->real == 1)
  {
    gsl_weight      = GSL_INTEG_SINE;
    sign_correction = -1.;
  }
  else
  {
    gsl_weight      = GSL_INTEG_COSINE;
    sign_correction = 1.;
  }
  
  if(params->tau_integrand == &chi_33_integrand)
  {
    step = 2. * params->pi / gamma;
    sign_correction = 1.;
  }
  else
  {
    step     = params->pi/gamma;
  }
  
  gsl_integration_qawo_table * table =
                   gsl_integration_qawo_table_alloc(gamma, step, gsl_weight, n);
  
  
  //need to update value of gamma
  params-> gamma = gamma;
  
  gsl_set_error_handler_off();
  gsl_function F;
  F.function = params->tau_integrand;
  F.params   = params;
  
  int i            = 0;
  int max_counter  = 2500;
  double tolerance = 1e-7;
  int counts       = 0;
  
  /*count the number of contributions smaller than small_tol; if there are more 
    than max_small_counter of them, then stop integrating. */
  int small_counter = 0;
  double small_tol  = 1e-14;  
  int max_small_counter = 1000;
  
  while(i == 0 || counts < max_counter)
  {
  
    if(params->tau_integrand == &chi_33_integrand)
    {
      gsl_integration_qng(&F, i*step, (i+1)*step, epsabs, epsrel, 
                          &ans_step, &error, &limit);
    }
    //Insert chi_rho_Q integrator here
    else
    {
      gsl_integration_qawo(&F, i*step, epsabs, epsrel, limit, w, table, 
                           &ans_step, &error);
    }

    ans_tot += ans_step;
    i       += 1;    
    
    if(fabs(ans_step / ans_tot) < tolerance)
    {
      counts += 1;
    }
    
    if(fabs(ans_tot) < small_tol)
    {
      small_counter++;
    }
    if(small_counter >= max_small_counter)
    {
      return 0.;
    }
  
  }
  
  gsl_integration_qawo_table_free(table);
  gsl_integration_workspace_free(w);
  ans_tot = ans_tot * sign_correction;
  
  return ans_tot;
}

/*trapezoidal: trapezoidal rule integrator for gamma integral, which can be
 *             faster than adaptive approaches because the integrand is
 *             smooth; one can also closely monitor the number of function
 *             calls to tau_integrator, which is important because these calls
 *             are rather expensive.
 *
 *@params: pointer to struct of parameters *p, start (starting point for gamma
 *         integral), end (ending point for gamma integral), samples (number of
 *         calls to tau_integrator allowed)
 *
 *@returns: gamma integral for a given component of chi_ij, evaluated using
 *          a trapezoidal rule integrator
 */
double trapezoidal(struct parameters *p, double start, double end, int samples)
{
  int i = 0;
  
  double step      = (end - start)/samples;
  double tolerance = 1e-6;
  double ans_step  = 0.;
  double ans_tot   = 0.;
  
  double p1 = p->gamma_integrand(start+i*step, p);
  double p2 = p->gamma_integrand(start+(i+1)*step, p);
  
  while(start + i * step <= end)
  {
    ans_step = step * (p1 + p2)/2.;
    
    
    ans_tot += ans_step;
    i++;
    
    p1 = p2;
    p2 = p->gamma_integrand(start+(i+1)*step, p);
  }
  
  return ans_tot;
}

/*trapezoidal_adaptive: trapezoidal rule integrator for gamma integral, similar
 *             to the above function, except it determines the ending point
 *             adaptively.
 *
 *@params: pointer to struct of parameters *p, start (starting point for gamma
 *         integral), step (step size between calls to tau_integrator)
 *
 *@returns: gamma integral for a given component of chi_ij, evaluated using
 *          an adaptive trapezoidal rule integrator
 */
double trapezoidal_adaptive(struct parameters *p, double start, double step)
{
  int i = 0;
  double tolerance = 1e-6;
  double ans_step  = 0.;
  double ans_tot   = 0.;
  
  double p1 = p->gamma_integrand(start+i*step, p);
  double p2 = p->gamma_integrand(start+(i+1)*step, p);
  
  while(ans_tot == 0 || fabs(ans_step/ans_tot) > tolerance)
  {
    ans_step = step * (p1 + p2)/2.;
    
    ans_tot += ans_step;
    i++;
    
    p1 = p2;
    p2 = p->gamma_integrand(start+(i+1)*step, p);
  }
  
  return ans_tot;
}

/*end_approx: approximate ending point for the gamma integral, beyond which
 *            contributions to the integral are assumed to be negligible.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: approximate end location for gamma integral
 */
double end_approx(struct parameters *p)
{
  double end;
  
  double MJ_max        = 0.5 * (3. * p->theta_e 
                                + sqrt(4. + 9. * p->theta_e * p->theta_e));
  double PL_max_real   = sqrt((1. + p->power_law_p)/p->power_law_p);
  double PL_max_moving = 50./sqrt(p->omega/p->omega_c) 
                         + 9. * pow(p->omega/p->omega_c, 1./3.);
  double kappa_max     = (-3. + 3. * p->kappa_width * p->kappa 
  		      + sqrt(1. - 4. * p->kappa 
  		  	     - 18. * p->kappa_width * p->kappa 
                               + 4. * pow(p->kappa, 2.) 
                               + 9. * pow(p->kappa_width * p->kappa, 2.))) 
  		   / (2. * (p->kappa - 2.));
  
  if(p->distribution == p->MAXWELL_JUETTNER)
  {
    end = 7. * MJ_max;
  }
  else if(p->distribution == p->POWER_LAW)
  {
    end = PL_max_moving;
  }
  else if(p->distribution == p->KAPPA_DIST)
  {
    end = 7. * kappa_max;
  }
  else
  {
    printf("\ndistribution or real/imag is set incorrectly");
    return 1.;
  }
  
  return end;
}

/*gsl_integrator: wrapper for a GSL integrator, to be used for the gamma
 *                integral.  This function can be very expensive, because
 *                it makes a large number of function calls to a region
 *                of parameter space that has a slowly convergent tau integral.
 *
 *@params: pointer to struct of parameters *p, start (starting point for gamma
 *         integral), end (ending point for gamma integral)
 *
 *@returns: gamma integral for a given component of chi_ij, evaluated using
 *          a GSL adaptive integrator
 */
double gsl_integrator(struct parameters *p, double start, double end)
{
  gsl_function F;
  F.function = p->gamma_integrand;
  F.params   = p;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000);
  
  double ans    = 0.;
  double error  = 0.;
  size_t limit  = 5000;
  double epsabs = 0.;
  double epsrel = 1e-5;
  int gsl_key   = 1;
  
  gsl_integration_qags(&F, start, end, epsabs, epsrel, limit, w, &ans, &error);
  
  gsl_integration_workspace_free(w);
  
  return ans;
}

#define NUM_QUAD 21
/*gauss_legendre: Fast 21-point Gauss-Legendre integrator for gamma integral.
 *
 *@params: pointer to struct of parameters *p, start (starting point for gamma
 *         integral), end (ending point for gamma integral)
 *
 *@returns: gamma integral for a given component of chi_ij, evaluated using
 *          a 21-point Gauss-Legendre quadrature method
 */
double gauss_legendre(struct parameters * p, double start, double end)
{
  double quadPts[NUM_QUAD] = \
        {-9.93752171e-01,  -9.67226839e-01,  -9.20099334e-01,
         -8.53363365e-01,  -7.68439963e-01,  -6.67138804e-01,
         -5.51618836e-01,  -4.24342120e-01,  -2.88021317e-01,
         -1.45561854e-01,   1.98918497e-16,   1.45561854e-01,
          2.88021317e-01,   4.24342120e-01,   5.51618836e-01,
          6.67138804e-01,   7.68439963e-01,   8.53363365e-01,
          9.20099334e-01,   9.67226839e-01,   9.93752171e-01};
  
  double weights[NUM_QUAD] = \
        {0.01601723,  0.03695379,  0.05713443,  0.07610011,  0.09344442,
         0.1087973 ,  0.12183142,  0.13226894,  0.13988739,  0.1445244 ,
         0.14608113,  0.1445244 ,  0.13988739,  0.13226894,  0.12183142,
         0.1087973 ,  0.09344442,  0.07610011,  0.05713443,  0.03695379,
         0.01601723};
  
  /*first we change integration variables to x = 1/(gamma^2 - 1), where the
    upper integration bound b = 1 and the lower bound a = 0.  We then apply
    the transformation: 
    \int_a^b f(x) dx = (b-a)/2 \int_{-1}^1 f((b-a)x/2 + (a+b)/2) dx
                     =     1/2 \int_{-1}^1 f(     x/2 +     1/2) dx */
  double x      = 0.;
  int i         = 0;
  double weight = 0.;
  double sum    = 0.;
  int n = NUM_QUAD;
  
  for(i = 0; i < n; i++)
  {
    x        = quadPts[i];
    weight   = weights[i];

    sum = sum + (end - start)/2. 
                * p->gamma_integrand((end - start)/2. * x + (end + start)/2., p)
                * weight;
  }

  return sum;
}

/*gamma_integrator: function that evaluates the gamma integral, using one of
 *                  the above integration algorithms.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: gamma integral for a given component of chi_ij
 */
double gamma_integrator(struct parameters *p)
{
  double prefactor = 2. * p->pi * p->omega_p*p->omega_p / (p->omega * p->omega);
  
  double start  = 1.;
  double end = end_approx(p);
  double ans_tot;
  
  /*for power-law, there is a singularity at low frequency at gamma = 1,
    which cannot be resolved by a trapezoidal integrator.  We are forced
    to use the (much slower) GSL integrator QAGS in these cases. */
  if(p->distribution == p->POWER_LAW && p->omega/p->omega_c < 5.)
  {
    ans_tot = gsl_integrator(p, start, end);
  }
  else if(p->distribution == p->POWER_LAW && p->omega/p->omega_c > 5.)
  {
    ans_tot = gauss_legendre(p, start, end);
  }
  else
  {
    ans_tot = trapezoidal(p, start, end, 100);
  }
  
  return prefactor * ans_tot;
}
