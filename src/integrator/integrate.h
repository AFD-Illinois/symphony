#ifndef SYMPHONY_INTEGRATE_H_
#define SYMPHONY_INTEGRATE_H_

#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include "integrands.h"

double gamma_integral(double min, double max, double n,
                      struct parameters * params
                     );
double n_integral(double min, double max,
                  struct parameters * params
                 );
double derivative_of_n(double n_start, struct parameters * params);
double n_summation(struct parameters *params);
#endif /* SYMPHONY_INTEGRATE_H_ */
