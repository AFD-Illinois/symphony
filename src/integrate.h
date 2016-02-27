#ifndef SYMPHONY_INTEGRATE_H_
#define SYMPHONY_INTEGRATE_H_

#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include "params.h"
#include "calc.h"
#include "maxwell_juettner/maxwell_juettner.h"
#include "power_law/power_law.h"
#include "kappa/kappa.h"

double gamma_integral(double min, double max, double n,
                 struct parameters * params);
double n_integral(double min, double max, double n,
                  struct parameters * params);
double normalize_f(struct parameters * params);
double derivative_of_n(double n_start, struct parameters * params);
#endif /* SYMPHONY_INTEGRATE_H_ */
