from symphonyHeaders cimport j_nu, alpha_nu, j_nu_fit, alpha_nu_fit, rho_nu_fit
#from symphonyHeaders cimport differential_of_f, numerical_differential_of_f
#from symphonyHeaders cimport setConstParams

def j_nu_py(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
            int polarization,
            double theta_e,
            double power_law_p,
            double gamma_min,
            double gamma_max,
            double gamma_cutoff,
            double kappa,
            double kappa_width):

  """Returns j_nu(nu, magnetic_field, electron_density, observer_angle, 
                  distribution, polarization, theta_e, power_law_p, 
                  gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width).
     Keys for distribution functions: symphonyPy.MAXWELL_JUETTNER, 
                                      symphonyPy.POWER_LAW, 
                                      symphonyPy.KAPPA_DIST
     Keys for Stokes parameter: symphonyPy.STOKES_I,
                                symphonyPy.STOKES_Q,
                                symphonyPy.STOKES_U,
                                symphonyPy.STOKES_V"""

  return j_nu(nu, magnetic_field, electron_density,
              observer_angle, distribution, polarization,
              theta_e, power_law_p, gamma_min, gamma_max,
              gamma_cutoff, kappa, kappa_width)

def alpha_nu_py(double nu,
                double magnetic_field,
                double electron_density,
                double observer_angle,
                int distribution,
                int polarization,
                double theta_e,
                double power_law_p,
                double gamma_min,
                double gamma_max,
                double gamma_cutoff,
                double kappa,
                double kappa_width):

  """Returns alpha_nu(nu, magnetic_field, electron_density, observer_angle,
                      distribution, polarization, theta_e, power_law_p, 
                      gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width).
     Keys for distribution functions: symphonyPy.MAXWELL_JUETTNER, 
                                      symphonyPy.POWER_LAW, 
                                      symphonyPy.KAPPA_DIST
     Keys for Stokes parameter: symphonyPy.STOKES_I,
                                symphonyPy.STOKES_Q,
                                symphonyPy.STOKES_U,
                                symphonyPy.STOKES_V"""

  return alpha_nu(nu, magnetic_field, electron_density,
              observer_angle, distribution, polarization,
              theta_e, power_law_p, gamma_min, gamma_max,
              gamma_cutoff, kappa, kappa_width)

def j_nu_fit_py(double nu,
                double magnetic_field,
                double electron_density,
                double observer_angle,
                int distribution,
                int polarization,
                double theta_e,
                double power_law_p,
                double gamma_min,
                double gamma_max,
                double gamma_cutoff,
                double kappa,
                double kappa_width):

  """Returns j_nu_fit(nu, magnetic_field, electron_density, observer_angle, 
                    distribution, polarization, theta_e, power_law_p, 
                    gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width).
     Keys for distribution functions: symphonyPy.MAXWELL_JUETTNER, 
                                      symphonyPy.POWER_LAW, 
                                      symphonyPy.KAPPA_DIST
     Keys for Stokes parameter: symphonyPy.STOKES_I,
                                symphonyPy.STOKES_Q,
                                symphonyPy.STOKES_U,
                                symphonyPy.STOKES_V"""

  return j_nu_fit(nu, magnetic_field, electron_density,
                  observer_angle, distribution, polarization,
                  theta_e, power_law_p, gamma_min, gamma_max,
                  gamma_cutoff, kappa, kappa_width)

def alpha_nu_fit_py(double nu,
                    double magnetic_field,
                    double electron_density,
                    double observer_angle,
                    int distribution,
                    int polarization,
                    double theta_e,
                    double power_law_p,
                    double gamma_min,
                    double gamma_max,
                    double gamma_cutoff,
                    double kappa,
                    double kappa_width):

  """Returns alpha_nu_fit(nu, magnetic_field, electron_density, observer_angle, 
                        distribution, polarization, theta_e, power_law_p, 
                        gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width).
     Keys for distribution functions: symphonyPy.MAXWELL_JUETTNER, 
                                      symphonyPy.POWER_LAW, 
                                      symphonyPy.KAPPA_DIST
     Keys for Stokes parameter: symphonyPy.STOKES_I,
                                symphonyPy.STOKES_Q,
                                symphonyPy.STOKES_U,
                                symphonyPy.STOKES_V"""

  return alpha_nu_fit(nu, magnetic_field, electron_density,
                      observer_angle, distribution, polarization,
                      theta_e, power_law_p, gamma_min, gamma_max,
                      gamma_cutoff, kappa, kappa_width)


def rho_nu_fit_py(double nu,
                  double magnetic_field,
                  double electron_density,
                  double observer_angle,
                  int distribution,
                  int polarization,
                  double theta_e,
                  double power_law_p,
                  double gamma_min,
                  double gamma_max,
                  double gamma_cutoff,
                  double kappa,
                  double kappa_width):

  """Returns rho_nu_fit(nu, magnetic_field, electron_density, observer_angle, 
                        distribution, polarization, theta_e, power_law_p, 
                        gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width).
     Keys for distribution functions: symphonyPy.MAXWELL_JUETTNER, 
                                      symphonyPy.POWER_LAW, 
                                      symphonyPy.KAPPA_DIST
     Keys for Stokes parameter: symphonyPy.STOKES_I,
                                symphonyPy.STOKES_Q,
                                symphonyPy.STOKES_U,
                                symphonyPy.STOKES_V"""

  return rho_nu_fit(nu, magnetic_field, electron_density,
                      observer_angle, distribution, polarization,
                      theta_e, power_law_p, gamma_min, gamma_max,
                      gamma_cutoff, kappa, kappa_width)



#def numerical_differential_of_f_py(double gamma,
#                                   double nu,
#                  	  	   double magnetic_field,
#                 	 	   double electron_density,
#                  	 	   double observer_angle,
#                  		   int distribution,
#                   		   int polarization,
#                    		   double theta_e,
#                    		   double power_law_p,
#                    		   double gamma_min,
#                    		   double gamma_max,
#                   		   double gamma_cutoff,
#                    		   double kappa,
#                    		   double kappa_width):
#
#  
#  struct parameters params;
#  setConstParams(&params);
#  params.nu                 = nu;
#  params.magnetic_field     = magnetic_field;
#  params.observer_angle     = observer_angle;
#  params.electron_density   = electron_density;
#  params.distribution       = distribution;
#  params.polarization       = polarization;
#  params.mode               = params.EMISSIVITY;
#  params.theta_e            = theta_e;
#  params.power_law_p        = power_law_p;
#  params.gamma_min          = gamma_min;
#  params.gamma_max          = gamma_max;
#  params.gamma_cutoff       = gamma_cutoff;
#  params.kappa              = kappa;
#  params.kappa_width        = kappa_width;
#  set_distribution_function(&params);
#
#  return numerical_differential_of_f(gamma, &params)



#DEFINE KEYS FOR DISTRIBUTION FUNCTIONS
MAXWELL_JUETTNER = 0
POWER_LAW        = 1
KAPPA_DIST       = 2

#DEFINE KEYS FOR STOKES PARAMETERS
STOKES_I         = 15
STOKES_Q         = 16
STOKES_U         = 17
STOKES_V         = 18
