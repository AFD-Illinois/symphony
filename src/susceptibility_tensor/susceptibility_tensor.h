#ifndef SUSCEPTIBILITY_H_
#define SUSCEPTIBILITY_H_

#define M_PI 3.1415926535897932384
struct paramsS
{
        double epsilon0;
	double epsilon;
	double e;
	double m;
	double c;
	double B;
	double n_e;
	double theta;
	double theta_e;
	double pl_p;
	double gamma_min;
	double gamma_max;
	double kappa;
	double kappa_width;
	double gamma_cutoff;
	double omega_c;
	double omega_p;
        double omega;
        double gamma;
	int real;
	int dist;
	double (*tau_integrand)(double, void * parameters);
	double (*gamma_integrand)(double, void * parameters);
};

double I_1_of_2(double alpha, double delta);
double I_1_analytic(double alpha, double delta);
double I_2_analytic(double alpha, double delta);
double MJ(struct paramsS * params);
double Df(struct paramsS * params);

double chi_11_integrand(double tau_prime, void * parameters);
double chi_12_integrand(double tau_prime, void * parameters);
double chi_32_integrand(double tau_prime, void * parameters);
double chi_13_integrand(double tau_prime, void * parameters);
double chi_33_integrand(double tau_prime, void * parameters);
double chi_22_integrand(double tau_prime, void * parameters);
double chi_22_integrand_p1(double tau_prime, void * parameters);
double chi_22_integrand_p2(double tau_prime, void * parameters);
double chi_22_integrand_real(double tau_prime, void * parameters);

double alpha_V_integrand(double tau_prime, void * parameters);

double tau_integrator(double gamma, void * parameters);
double gamma_integrator(struct paramsS * p);
double end_approx(struct paramsS * p);

double chi_11(struct paramsS * p);
double chi_22(struct paramsS * p);
double chi_33(struct paramsS * p);
double chi_12(struct paramsS * p);
double chi_21(struct paramsS * p);
double chi_23(struct paramsS * p);
double chi_32(struct paramsS * p);
double chi_13(struct paramsS * p);
double chi_31(struct paramsS * p);
#endif /* SUSCEPTIBILITY_H_ */
