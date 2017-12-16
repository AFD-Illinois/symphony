#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "susceptibility_tensor.h"

/*MJ: returns the value of the Maxwell-Juettner (relativistic thermal)
 *    distribution function for the parameters in the struct p.
 *    TODO: I think a minus sign and other parts of Df are absorbed
 *          into the integrands themselves, and should be moved here.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: MJ(params->gamma, params->theta_e)
 */
double MJ(struct parameters * params)
{
	double ans = exp(-params->gamma/params->theta_e) 
		   / (4. * params->pi * params->theta_e*params->theta_e 
	  	      * gsl_sf_bessel_Kn(2, 1./params->theta_e));
	return ans;
}

/*PL: returns the value of the power-law distribution function for 
 *    the parameters in the struct p.
 *    TODO: I think a minus sign and other parts of Df are absorbed
 *          into the integrands themselves, and should be moved here.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: PL(params->gamma, params->power_law_p, params->gamma_min,
 *             params->gamma_max)
 */
double PL(struct parameters * params)
{
	if(params->gamma > params->gamma_max || params->gamma < params->gamma_min)
	{
		return 0.;
	}

	double beta = sqrt(1. - 1./pow(params->gamma, 2.));

	double ans = (params->pl_p - 1.) * (-1 + 2. * params->gamma * params->gamma 
					 + params->pl_p * (params->gamma*params->gamma - 1.))
		    / (4. * params->pi * (pow(params->gamma_min, -1. - params->pl_p) - pow(params->gamma_max, -1. - params->pl_p))
			* beta * (params->gamma*params->gamma - 1.)) * pow(params->gamma, -3. - params->pl_p);
	return ans;	
}

/*kappa_to_be_normalized: returns the value of the unnormalized relativistic
 *                        kappa distribution for the parameters in the struct p.
 *                        This function is integrated over gamma to determine
 *                        the normalization constant in normalize_f(), below.
 *    TODO: I think a minus sign and other parts of Df are absorbed
 *          into the integrands themselves, and should be moved here.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: kappa_to_be_normalized(gamma, void * parameters)
 */
double kappa_to_be_normalized(double gamma, void * parameters)
{
	struct parameters * params = (struct parameters*) parameters;

        double beta = sqrt(1. - 1./pow(gamma, 2.));

        double body = pow((1. + (gamma - 1.)/(params->kappa * params->kappa_width)), -1. - params->kappa);

	double d3p = 4. * params->pi * gamma*gamma * beta;

        double ans = body * d3p;

        return ans;

}

/*normalize_f: normalizes the distribution function using GSL's 
 *             QAGIU integrator.   
 *
 *@params: struct of parameters to pass to distribution function
 *@returns: 1 over the normalization constant for the chosen distribution
 */
double normalize_f(double (*distribution)(double, void *),
                   struct parameters * params
                  )
{
  
	/*set GSL QAGIU integrator parameters */
	double lower_bound = 1.;
	double absolute_error = 0.;
	double relative_error = 1e-8;
	int limit  = 1000;
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	double result, error;
	gsl_function F;
	
	F.function = distribution;
	
	F.params = params;
	
	gsl_integration_qagiu(&F, lower_bound, absolute_error, 
	                      relative_error, limit, w, &result, &error
	                     );
	
	
	gsl_integration_workspace_free(w);
	
	return result;
}

/*kappa: returns the value of the relativistic kappa distribution function for 
 *    the parameters in the struct p, with the normalization handled
 *    numerically by normalize_f().
 *    TODO: I think a minus sign and other parts of Df are absorbed
 *          into the integrands themselves, and should be moved here.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: PL(params->gamma, params->power_law_p, params->gamma_min,
 *             params->gamma_max)
 */
double kappa(struct parameters * params)
{
	static double norm                  = 0.;
	static double previous_kappa        = 0.;
	static double previous_kappa_width  = 0.;
	static double previous_gamma_cutoff = 0.;
	if(norm == 0. || previous_kappa_width != params->kappa_width
	              || previous_kappa       != params->kappa)
	{
	  norm                  = 1./normalize_f(&kappa_to_be_normalized, params);
	  previous_kappa        = params->kappa;
	  previous_kappa_width  = params->kappa_width;
	  previous_gamma_cutoff = params->gamma_cutoff;
	}

	double beta = sqrt(1. - 1./pow(params->gamma, 2.));

	double body = -pow((1. + (params->gamma - 1.)/(params->kappa * params->kappa_width)), -2. - params->kappa)
		     *(-1. - params->kappa) / (params->kappa_width * params->kappa);

	double ans = norm * body;

	return ans;

}

/*Df: chooses the distribution function given by the parameter params->dist,
 *    and returns that distribution function evaluated with the given
 *    parameters in the struct params.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: MJ(params), PL(params), or kappa(params), depending on params->dist
 */
double Df(struct parameters * params)
{
	if(params->dist == 0)
	{
		return MJ(params);
	}
	else if(params->dist == 1)
	{
		return PL(params);
	}
	else if(params->dist == 2)
	{
		return kappa(params);
	}

	return 0.;

}

/*I_1_analytic: analytic solution to the first cos_xi integral with Bessel
 *              function order 0, which we have named I_1(0).  This term 
 *              is used in chi_{11, 12, 21, 22}.
 *
 *@params: double alpha, double delta; these parameters are defined in terms of
 *         relevant quantities in the chi_ij integrands.
 *
 *@returns: analytic solution to cos_xi integral I_1(0)
 */
double I_1_analytic(double alpha, double delta)
{
	double A     = sqrt(alpha*alpha + delta*delta);

	if(alpha == 0. || delta == 0.)
	{
		return 0.;
	}

	double ans = 2. * ( (2. * alpha*alpha + (alpha*alpha - 1.)*delta*delta + pow(delta, 4.))*sin(A) 
                - (2. * alpha*alpha - delta*delta) * A * cos(A)) / pow(A, 5.);
    	return ans;

}

/*I_1_of_2: analytic solution to the first cos_xi integral with Bessel function
 *          order 2, which we have named I_1(2).  This term is used in 
 *          chi_{11, 22}.
 *
 *@params: double alpha, double delta; these parameters are defined in terms of
 *         relevant quantities in the chi_ij integrands.
 *
 *@returns: analytic solution to cos_xi integral I_1(2)
 */
double I_1_of_2(double alpha, double delta)
{
	double A   = sqrt(alpha*alpha + delta*delta);
	double ans = -2. * delta*delta * (3. * A * cos(A) + (-3. + A*A) * sin(A)) / pow(A, 5.);
	return ans;
}

/*I_2_analytic: analytic solution to the second cos_xi integral with Bessel 
 *          function order 1, which we have named I_2(1).  This term is used in 
 *          chi_{13, 23, 31, 32}.
 *
 *@params: double alpha, double delta; these parameters are defined in terms of
 *         relevant quantities in the chi_ij integrands.
 *
 *@returns: analytic solution to cos_xi integral I_2(1)
 */
double I_2_analytic(double alpha, double delta)
{
	double A     = sqrt(alpha*alpha + delta*delta);

	if(alpha == 0. || delta == 0.)
        {
                return 0.;
        }

	double num   = 2. * alpha * delta * (3. * A * cos(A) + (-3. + A*A) * sin(A));
	double denom = pow(A, 5.);
	double ans   = num / denom;
	return ans;
}

/*I_3_analytic: analytic solution to the third cos_xi integral with Bessel 
 *          function order 0, which we have named I_3(1).  This term is used in 
 *          chi_33.
 *
 *@params: double alpha, double delta; these parameters are defined in terms of
 *         relevant quantities in the chi_ij integrands.
 *
 *@returns: analytic solution to cos_xi integral I_3(0)
 */
double I_3_analytic(double alpha, double delta)
{
	if(alpha == 0. || delta == 0.)
	{
		return 0.;
	}	

	double A        = sqrt(alpha*alpha + delta*delta);
	double term1    = 6. * alpha*alpha * cos(A) / pow(A, 4.);
	double term2    = -2. * cos(A) / (A*A);
	double term3    = 6. * delta*delta * sin(A) / pow(A, 5.);
	double term4    = -4. * sin(A) / pow(A, 3.);
	double term5    = 2. * alpha*alpha * sin(A) / pow(A, 3.);
	double ans      = term1 + term2 + term3 + term4 + term5;
	return ans;
}

/*I_3_limit: analytic solution to a limiting case of the I_3(0) cos_xi
 *           integral, namely the one where A = alpha.  We subtract this
 *           limiting form off of the existing tau integral for chi_33,
 *           which makes the tau integral converge more rapidly.  We then
 *           evaluate the limit term analytically and add it back on.
 *           This function is the analytic integral of that limit term.
 *
 *@params: double alpha, double delta; these parameters are defined in terms of
 *         relevant quantities in the chi_ij integrands.
 *
 *@returns: analytic solution to the I_3(0) limit case, used to improve the
 *          rate of convergence of the chi_33 tau integral.
 */
double I_3_limit(double alpha, double delta)
{
	if(alpha == 0.)
	{
		return 0.;
	}

        double A        = fabs(alpha);
        double term5    = 2. * alpha*alpha * sin(A) / pow(A, 3.);
        double ans      = term5;
        return ans;
}

/*chi_11_integrand: integrand for the component chi_11 of the susceptibility
 *                  tensor.  The term e^(i*gamma*tau) determines the real
 *                  and imaginary parts, and is left out of the below function
 *                  because the integrator QAWO includes an implicit factor of
 *                  sin(gamma*tau) or cos(gamma*tau), which is appropriately
 *                  chosen in tau_integrator.
 *
 *@params: double tau_prime (the first integration is carried out over this
 *         variable), void * parameters (a pointer to the struct of parameters)
 *
 *@returns: chi_11_integrand at tau_prime with parameters params
 */
double chi_11_integrand(double tau_prime, void * parameters)
{
	struct parameters * params = (struct parameters*) parameters;

	double prefactor  = 1.; //should be 1j
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(params->epsilon * params->omega_c) 
			   * sin(params->theta) * params->gamma * beta 
			   * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * Df(params);
//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = -sin(tau_prime * params->gamma) 
//			    * sin((epsilon * params->omega_c / params->omega) * tau_prime);
	double tauxi_term = 0.5 * (cos((params->epsilon * params->omega_c / params->omega) * tau_prime) * I_1_analytic(alpha, delta)
				   - I_1_of_2(alpha, delta));
	double ans        = prefactor * gamma_term * tauxi_term * params->gamma*params->gamma * beta;
	
	return ans;
}

/*chi_12_integrand: integrand for the component chi_12 of the susceptibility
 *                  tensor.  The term e^(i*gamma*tau) determines the real
 *                  and imaginary parts, and is left out of the below function
 *                  because the integrator QAWO includes an implicit factor of
 *                  sin(gamma*tau) or cos(gamma*tau), which is appropriately
 *                  chosen in tau_integrator.
 *
 *@params: double tau_prime (the first integration is carried out over this
 *         variable), void * parameters (a pointer to the struct of parameters)
 *
 *@returns: chi_12_integrand at tau_prime with parameters params
 */
double chi_12_integrand(double tau_prime, void * parameters)
{
	struct parameters * params = (struct parameters*) parameters;

	double prefactor  = 1.; //should be 1j
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(params->epsilon * params->omega_c) 
			   * sin(params->theta) * params->gamma * beta 
			   * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * Df(params);
//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = -sin(tau_prime * params->gamma) 
//			    * sin((epsilon * params->omega_c / params->omega) * tau_prime);
	double tau_term   = sin( - (params->epsilon * params->omega_c / params->omega) * tau_prime);
	double xi_term    = -0.5 * I_1_analytic(alpha, delta);
	double ans        = prefactor * gamma_term * xi_term * tau_term * params->gamma*params->gamma * beta;
	
	return ans;
}

/*chi_13_integrand: integrand for the component chi_13 of the susceptibility
 *                  tensor.  The term e^(i*gamma*tau) determines the real
 *                  and imaginary parts, and is left out of the below function
 *                  because the integrator QAWO includes an implicit factor of
 *                  sin(gamma*tau) or cos(gamma*tau), which is appropriately
 *                  chosen in tau_integrator.
 *
 *@params: double tau_prime (the first integration is carried out over this
 *         variable), void * parameters (a pointer to the struct of parameters)
 *
 *@returns: chi_13_integrand at tau_prime with parameters params
 */
double chi_13_integrand(double tau_prime, void * parameters)
{
	struct parameters * params = (struct parameters*) parameters;

	double prefactor  = 1.; //should be 1j
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(params->epsilon * params->omega_c) 
			   * sin(params->theta) * params->gamma * beta 
			   * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * Df(params);//* exp(-params->gamma/params->theta_e);
				//explicit imag part
//	double tau_term   = cos(tau_prime * params->gamma) * cos((params->epsilon * params->omega_c / params->omega) * tau_prime/2.);

//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = sin(tau_prime * params->gamma) 
//			    * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / 2.);

	/*K_13 is K_32 except -S_2 -> C_2*/
	double tau_term   = cos((params->epsilon * params->omega_c / params->omega) * tau_prime / 2.);
	double xi_term    = I_2_analytic(alpha, delta); //NOTE: Graf thm error changes sign of this term TODO: write better explanation 
	double ans        = prefactor * gamma_term * xi_term * tau_term * params->gamma*params->gamma * beta;

	return ans;
}

/*chi_22_integrand: integrand for the component chi_22 of the susceptibility
 *                  tensor.  This term has two components, and it seems that
 *                  integrating them together is faster for the real part,
 *                  and integrating them separately is faster for the imaginary
 *                  part.
 *
 *@params: double tau_prime (the first integration is carried out over this
 *         variable), void * parameters (a pointer to the struct of parameters)
 *
 *@returns: chi_22_integrand at tau_prime with parameters params
 */
double chi_22_integrand(double tau_prime, void * parameters)
{
	struct parameters * params = (struct parameters*) parameters;

	if(params->real == 1)
	{
		return chi_22_integrand_real(tau_prime, parameters);
	}
	else
	{
		return chi_22_integrand_p1(tau_prime, parameters)
		      +chi_22_integrand_p2(tau_prime, parameters);
	}

	return 0.;
}

/*note: for imaginary part of chi_22 splitting the integrand up is faster (40sec vs 170sec)
	but for real part the combined integrand is faster (3sec vs 10sec)*/
/*chi_22_integrand_real: integrand for the component chi_22 of the 
 *                       susceptibility tensor.  Integrating the whole term is
 *                       faster for the real part, so the real part's integrand
 *                       is the whole term.
 *
 *@params: double tau_prime (the first integration is carried out over this
 *         variable), void * parameters (a pointer to the struct of parameters)
 *
 *@returns: chi_22_integrand_real at tau_prime with parameters params
 */
double chi_22_integrand_real(double tau_prime, void * parameters)
{
	struct parameters * params = (struct parameters*) parameters;

	double prefactor  = 1.; //should be 1j
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(params->epsilon * params->omega_c) 
			   * sin(params->theta) * params->gamma * beta 
			   * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * Df(params);
//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = -sin(tau_prime * params->gamma) 
//			    * sin((epsilon * params->omega_c / params->omega) * tau_prime);
	double tauxi_term = 0.5 * (cos((params->epsilon * params->omega_c / params->omega) * tau_prime) * I_1_analytic(alpha, delta)
				   + I_1_of_2(alpha, delta));

	double ans        = prefactor * gamma_term * tauxi_term * params->gamma*params->gamma * beta;
	
	return ans;
}

/*chi_22_integrand_p1: integrand for the first part of the component chi_22 of 
 *                     the susceptibility tensor.  Integrating the two parts
 *                     separately is faster for the imaginary part, so the
 *                     below function corresponds to just the first term of the
 *                     imaginary part integrand.
 *
 *@params: double tau_prime (the first integration is carried out over this
 *         variable), void * parameters (a pointer to the struct of parameters)
 *
 *@returns: chi_22_integrand_p1 at tau_prime with parameters params
 */
double chi_22_integrand_p1(double tau_prime, void * parameters)
{
        struct parameters * params = (struct parameters*) parameters;

        double prefactor  = 1.; //should be 1j
        double beta       = sqrt(1. - pow(params->gamma, -2.));
        double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
        double delta      = 2. * params->omega/(params->epsilon * params->omega_c)
                           * sin(params->theta) * params->gamma * beta 
                           * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

        double gamma_term = beta*beta * params->gamma * Df(params);
//      double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//      double tau_term   = -sin(tau_prime * params->gamma) 
//                          * sin((epsilon * params->omega_c / params->omega) * tau_prime);
        double tauxi_term = 0.5 * (cos((params->epsilon * params->omega_c / params->omega) * tau_prime) * I_1_analytic(alpha, delta)
                                   + 0.);

        double ans        = prefactor * gamma_term * tauxi_term * params->gamma*params->gamma * beta;

        return ans;
}

/*chi_22_integrand_p2: integrand for the second part of the component chi_22 of 
 *                     the susceptibility tensor.  Integrating the two parts
 *                     separately is faster for the imaginary part, so the
 *                     below function corresponds to just the second term of the
 *                     imaginary part integrand.
 *
 *@params: double tau_prime (the first integration is carried out over this
 *         variable), void * parameters (a pointer to the struct of parameters)
 *
 *@returns: chi_22_integrand_p2 at tau_prime with parameters params
 */
double chi_22_integrand_p2(double tau_prime, void * parameters)
{
        struct parameters * params = (struct parameters*) parameters;

        double prefactor  = 1.; //should be 1j
        double beta       = sqrt(1. - pow(params->gamma, -2.));
        double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
        double delta      = 2. * params->omega/(params->epsilon * params->omega_c)
                           * sin(params->theta) * params->gamma * beta 
                           * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

        double gamma_term = beta*beta * params->gamma * Df(params);
//      double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//      double tau_term   = -sin(tau_prime * params->gamma) 
//                          * sin((epsilon * params->omega_c / params->omega) * tau_prime);
        double tauxi_term = 0.5 * (0.
                                   + I_1_of_2(alpha, delta));

        double ans        = prefactor * gamma_term * tauxi_term * params->gamma*params->gamma * beta;

        return ans;
}

/*chi_32_integrand: integrand for the component chi_32 of the susceptibility
 *                  tensor.  The term e^(i*gamma*tau) determines the real
 *                  and imaginary parts, and is left out of the below function
 *                  because the integrator QAWO includes an implicit factor of
 *                  sin(gamma*tau) or cos(gamma*tau), which is appropriately
 *                  chosen in tau_integrator.
 *
 *@params: double tau_prime (the first integration is carried out over this
 *         variable), void * parameters (a pointer to the struct of parameters)
 *
 *@returns: chi_32_integrand at tau_prime with parameters params
 */
double chi_32_integrand(double tau_prime, void * parameters)
{
	struct parameters * params = (struct parameters*) parameters;

	double prefactor  = 1.; // should be 1j 
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(params->epsilon * params->omega_c) 
			   * sin(params->theta) * params->gamma * beta 
			   * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * Df(params);//* exp(-params->gamma/params->theta_e);
//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = sin(tau_prime * params->gamma) 
//			    * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / 2.);
	double tau_term   = sin((params->epsilon * params->omega_c / params->omega) * tau_prime / 2.);
	double xi_term    = I_2_analytic(alpha, delta); //should be times 1j * -1j = +1
	double ans        = prefactor * gamma_term * xi_term * tau_term * params->gamma*params->gamma * beta;

	return ans;
}

/*chi_33_integrand: integrand for the component chi_33 of the susceptibility
 *                  tensor.  The rate of convergence for the imaginary part's
 *                  tau integral can be increased significantly by subtracting
 *                  off a limit term, which is then integrated analytically and
 *                  added onto the result of the numerical tau integral.
 *
 *@params: double tau_prime (the first integration is carried out over this
 *         variable), void * parameters (a pointer to the struct of parameters)
 *
 *@returns: chi_13_integrand at tau_prime with parameters params
 */
double chi_33_integrand(double tau_prime, void * parameters)
{
	struct parameters * params = (struct parameters*) parameters;

	double prefactor  = 1.; //should be 1j
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(params->epsilon * params->omega_c) 
			   * sin(params->theta) * params->gamma * beta 
			   * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * Df(params);
//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = -sin(tau_prime * params->gamma) 
//			    * sin((epsilon * params->omega_c / params->omega) * tau_prime);
	double tau_term   = 1.;
	
	double xi_term;
	/*subtracting off this term on the imaginary part speeds convergence
	  by a factor of 2-3.  The term integrates to zero, so we get this
	  speed boost basically for free.  The real part of this term is
	  nonzero, however, and actually slows convergence by a factor of 10.*/
	if(params->real == 1)
	{
		xi_term = I_3_analytic(alpha, delta);
		xi_term *= -sin(params->gamma * tau_prime);
	}
	else
	{
		xi_term = I_3_analytic(alpha, delta) - I_3_limit(alpha, delta);
		xi_term *= cos(params->gamma * tau_prime);
	}

	double ans        = prefactor * gamma_term * xi_term * tau_term * params->gamma*params->gamma * beta;
	
	return ans;
}
