#ifndef SYMPHONY_INTEGRATE_H_
#define SYMPHONY_INTEGRATE_H_

double gsl_integrate(double min, double max, double n, 
                     struct parameters * params
                    );
double gamma_integral(double min, double max, double n,
                 struct parameters * params);
double n_integral(double min, double max, double n,
                  struct parameters * params);
double normalize_f(struct parameters * params);
double derivative(double n_start, struct parameters * params);
#endif /* SYMPHONY_INTEGRATE_H_ */
