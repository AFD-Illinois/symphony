cdef extern from "symphony.h":
    
    double j_nu(double nu,
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
                double kappa_width)

    double alpha_nu(double nu,
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
                    double kappa_width)

    double j_nu_fit(double nu,
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
                    double kappa_width)

    double alpha_nu_fit(double nu,
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
                        double kappa_width)

    double rho_nu_fit(double nu,
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
                        double kappa_width)

