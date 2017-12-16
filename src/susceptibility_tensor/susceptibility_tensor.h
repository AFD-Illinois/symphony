#ifndef SUSCEPTIBILITY_H_
#define SUSCEPTIBILITY_H_

#include "../params.h"

//#define M_PI 3.1415926535897932384

double chi_11_symphony(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
            int polarization,
            double theta_e,
            double power_law_p,
            double gamma_min,
            double gamma_max,
            double gamma_cutoff,
            double kappa,
            double kappa_width,
            char **error_message
           );

double I_1_of_2(double alpha, double delta);
double I_1_analytic(double alpha, double delta);
double I_2_analytic(double alpha, double delta);
double MJ(struct parameters * params);
double Df(struct parameters * params);

double chi_11_integrand(double tau_prime, void * parameters);
double chi_12_integrand(double tau_prime, void * parameters);
double chi_32_integrand(double tau_prime, void * parameters);
double chi_13_integrand(double tau_prime, void * parameters);
double chi_33_integrand(double tau_prime, void * parameters);
double chi_22_integrand(double tau_prime, void * parameters);
double chi_22_integrand_p1(double tau_prime, void * parameters);
double chi_22_integrand_p2(double tau_prime, void * parameters);
double chi_22_integrand_real(double tau_prime, void * parameters);

double alpha_V_integrand(double tau_prime, void * parameters);

double tau_integrator(double gamma, void * parameters);
double gamma_integrator(struct parameters * p);
double end_approx(struct parameters * p);

int set_params(struct parameters * p);
double alpha_I(struct parameters * p);
double alpha_Q(struct parameters * p);
double alpha_V(struct parameters * p);
double rho_Q(struct parameters * p);
double rho_V(struct parameters * p);
double plotter(struct parameters p);

double chi_11(struct parameters * p);
double chi_22(struct parameters * p);
double chi_33(struct parameters * p);
double chi_12(struct parameters * p);
double chi_21(struct parameters * p);
double chi_23(struct parameters * p);
double chi_32(struct parameters * p);
double chi_13(struct parameters * p);
double chi_31(struct parameters * p);
#endif /* SUSCEPTIBILITY_H_ */
