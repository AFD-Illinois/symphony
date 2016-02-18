#ifndef SYMPHONY_FITS_H_
#define SYMPHONY_FITS_H_

double thermal_I(struct parameters * params);
double thermal_Q(struct parameters * params);
double thermal_V(struct parameters * params);
double thermal_I_abs(struct parameters * params);
double thermal_Q_abs(struct parameters * params);
double thermal_V_abs(struct parameters * params);
double power_law_I(struct parameters * params);
double power_law_Q(struct parameters * params);
double power_law_V(struct parameters * params);
double power_law_I_abs(struct parameters * params);
double power_law_Q_abs(struct parameters * params);
double power_law_V_abs(struct parameters * params);
double kappa_I(struct parameters * params);
double kappa_Q(struct parameters * params);
double kappa_V(struct parameters * params);
double kappa_I_abs(struct parameters * params);
double kappa_Q_abs(struct parameters * params);
double kappa_V_abs(struct parameters * params);
double check_for_errors(struct parameters * params);
double j_nu_fit(double nu, 
                double magnetic_field, 
                double electron_density, 
                double observer_angle, 
                int distribution, 
                int polarization
                );
double alpha_nu_fit(double nu,
                double magnetic_field,
                double electron_density,
                double observer_angle,
                int distribution,
                int polarization
                );
#endif /* SYMPHONY_FITS_H_ */

