#ifndef SYMPHONY_INTEGRATE_H_
#define SYMPHONY_INTEGRATE_H_

double gsl_integrate(double min, double max, double n, 
                     struct parameters params
                    );
double normalize_f(struct parameters params);
double derivative(double n_start, double nu);
#endif /* SYMPHONY_INTEGRATE_H_ */
