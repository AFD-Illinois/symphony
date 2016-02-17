#ifndef SYMPHONY_CALC_H_
#define SYMPHONY_CALC_H_

double get_nu_c(struct parameters params);
double n_summation(struct parameters *params);
double gamma_integrand(double gamma, void * paramsInput);
double distribution_function(double gamma, struct parameters * params);
double my_Bessel_J(double n, double z);
double my_Bessel_dJ(double n, double z);
double maxwell_juttner_f(double gamma, struct parameters * params);
double power_law_f(double gamma, struct parameters * params);
double kappa_f(double gamma, struct parameters * params);
double j_nu(double nu, 
            double magnetic_field, 
            double electron_density,
            double observer_angle,
            int distribution,
            int polarization);
double gamma_integration_result(double n, void * paramsInput);
double power_law_to_be_normalized(double gamma, void * paramsInput);
double kappa_to_be_normalized(double gamma, void * paramsInput);

#endif /* SYMPHONY_CALC_H_ */
