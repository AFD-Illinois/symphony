#ifndef SYMPHONY_MAXWELL_JUETTNER_H_
#define SYMPHONY_MAXWELL_JUETTNER_H_
#include "../params.h"

double maxwell_juttner_f(double gamma, struct parameters * params);
double differential_of_maxwell_juttner(double gamma, struct parameters * params);

/* Fits */

/* Emissivities */
double thermal_I(struct parameters * params);
double thermal_Q(struct parameters * params);
double thermal_V(struct parameters * params);

double planck_func(struct parameters * params);

/* Absorptivities */
double thermal_I_abs(struct parameters * params);
double thermal_Q_abs(struct parameters * params);
double thermal_V_abs(struct parameters * params);

#endif /* SYMPHONY_MAXWELL_JUETTNER_H_ */
