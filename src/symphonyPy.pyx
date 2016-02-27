from symphonyHeaders cimport j_nu, alpha_nu, j_nu_fit, alpha_nu_fit

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
                  gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width)."""

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
                      gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width)."""

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

  return alpha_nu_fit(nu, magnetic_field, electron_density,
                      observer_angle, distribution, polarization,
                      theta_e, power_law_p, gamma_min, gamma_max,
                      gamma_cutoff, kappa, kappa_width)
