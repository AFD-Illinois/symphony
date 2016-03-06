//#ifndef SYMPHONY_DISTRIBUTION_FUNCTION_COMMON_ROUTINES_H_
//#define SYMPHONY_DISTRIBUTION_FUNCTION_COMMON_ROUTINES_H_

#include <gsl/gsl_integration.h>
#include "params.h"

double normalize_f(double (*distribution)(double, void *),
                   struct parameters * params
                  );

double num_differential_of_f(double gamma, struct parameters * params);

//#endif /* SYMPHONY_DISTRIBUTION_FUNCTION_COMMON_ROUTINES_H_ */
