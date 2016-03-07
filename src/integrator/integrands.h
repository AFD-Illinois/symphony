#ifndef SYMPHONY_INTEGRANDS_H_
#define SYMPHONY_INTEGRANDS_H_

#include "params.h"
#include "maxwell_juettner/maxwell_juettner.h"
#include "power_law/power_law.h"
#include "kappa/kappa.h"

double n_summation(struct parameters *params);
double gamma_integrand(double gamma, void * paramsInput);
void   set_distribution_function(struct parameters * params);
double my_Bessel_J(double n, double z);
double my_Bessel_dJ(double n, double z);
double gamma_integration_result(double n, void * paramsInput);
#endif /* SYMPHONY_INTEGRANDS_H_ */
