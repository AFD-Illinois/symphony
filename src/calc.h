#ifndef SYMPHONY_CALC_H_
#define SYMPHONY_CALC_H_

#include "params.h"
#include "maxwell_juettner/maxwell_juettner.h"
#include "power_law/power_law.h"
#include "kappa/kappa.h"

double get_nu_c(struct parameters params);
double n_summation(struct parameters *params);
double gamma_integrand(double gamma, void * paramsInput);
double distribution_function(double gamma, struct parameters * params);
double my_Bessel_J(double n, double z);
double my_Bessel_dJ(double n, double z);
double j_nu(double nu,
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
            double kappa_width);
double alpha_nu(double nu,
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
                double kappa_width);
double gamma_integration_result(double n, void * paramsInput);

#endif /* SYMPHONY_CALC_H_ */
