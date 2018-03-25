from susceptibilityHeaders cimport chi_11_symphony, chi_12_symphony, chi_13_symphony, chi_21_symphony, chi_22_symphony, chi_23_symphony, chi_31_symphony, chi_32_symphony, chi_33_symphony

def chi_11_symphony_py(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
            double theta_e,
            double power_law_p,
            double gamma_min,
            double gamma_max,
            double gamma_cutoff,
            double kappa,
            double kappa_width):

  """Returns chi_11_nu(nu, magnetic_field, electron_density, observer_angle, 
                       distribution, real_part, theta_e, power_law_p, gamma_min, 
                       gamma_max, gamma_cutoff, kappa, kappa_width).
     Keys for distribution functions: symphonyPy.MAXWELL_JUETTNER, 
                                      symphonyPy.POWER_LAW, 
                                      symphonyPy.KAPPA_DIST
  """

  cdef char* error_message = NULL
  result = chi_11_symphony(nu, magnetic_field, electron_density,
                observer_angle, distribution, real_part, theta_e, 
		power_law_p, gamma_min, gamma_max, gamma_cutoff, 
		kappa, kappa_width, &error_message)
  if error_message:
    raise RuntimeError (error_message)
  return result

def chi_12_symphony_py(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
            double theta_e,
            double power_law_p,
            double gamma_min,
            double gamma_max,
            double gamma_cutoff,
            double kappa,
            double kappa_width):

  """Returns chi_12_nu(nu, magnetic_field, electron_density, observer_angle, 
                       distribution, real_part, theta_e, power_law_p, gamma_min, 
                       gamma_max, gamma_cutoff, kappa, kappa_width).
     Keys for distribution functions: symphonyPy.MAXWELL_JUETTNER, 
                                      symphonyPy.POWER_LAW, 
                                      symphonyPy.KAPPA_DIST
  """

  cdef char* error_message = NULL
  result = chi_12_symphony(nu, magnetic_field, electron_density,
                observer_angle, distribution, real_part, theta_e, 
		power_law_p, gamma_min, gamma_max, gamma_cutoff, 
		kappa, kappa_width, &error_message)
  if error_message:
    raise RuntimeError (error_message)
  return result

def chi_13_symphony_py(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
            double theta_e,
            double power_law_p,
            double gamma_min,
            double gamma_max,
            double gamma_cutoff,
            double kappa,
            double kappa_width):

  """Returns chi_13_nu(nu, magnetic_field, electron_density, observer_angle, 
                       distribution, real_part, theta_e, power_law_p, gamma_min, 
                       gamma_max, gamma_cutoff, kappa, kappa_width).
     Keys for distribution functions: symphonyPy.MAXWELL_JUETTNER, 
                                      symphonyPy.POWER_LAW, 
                                      symphonyPy.KAPPA_DIST
  """

  cdef char* error_message = NULL
  result = chi_13_symphony(nu, magnetic_field, electron_density,
                observer_angle, distribution, real_part, theta_e, 
		power_law_p, gamma_min, gamma_max, gamma_cutoff, 
		kappa, kappa_width, &error_message)
  if error_message:
    raise RuntimeError (error_message)
  return result

def chi_21_symphony_py(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
            double theta_e,
            double power_law_p,
            double gamma_min,
            double gamma_max,
            double gamma_cutoff,
            double kappa,
            double kappa_width):

  """Returns chi_21_nu(nu, magnetic_field, electron_density, observer_angle, 
                       distribution, real_part, theta_e, power_law_p, gamma_min, 
                       gamma_max, gamma_cutoff, kappa, kappa_width).
     Keys for distribution functions: symphonyPy.MAXWELL_JUETTNER, 
                                      symphonyPy.POWER_LAW, 
                                      symphonyPy.KAPPA_DIST
  """

  cdef char* error_message = NULL
  result = chi_21_symphony(nu, magnetic_field, electron_density,
                observer_angle, distribution, real_part, theta_e, 
		power_law_p, gamma_min, gamma_max, gamma_cutoff, 
		kappa, kappa_width, &error_message)
  if error_message:
    raise RuntimeError (error_message)
  return result

def chi_22_symphony_py(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
            double theta_e,
            double power_law_p,
            double gamma_min,
            double gamma_max,
            double gamma_cutoff,
            double kappa,
            double kappa_width):

  """Returns chi_22_nu(nu, magnetic_field, electron_density, observer_angle, 
                       distribution, real_part, theta_e, power_law_p, gamma_min, 
                       gamma_max, gamma_cutoff, kappa, kappa_width).
     Keys for distribution functions: symphonyPy.MAXWELL_JUETTNER, 
                                      symphonyPy.POWER_LAW, 
                                      symphonyPy.KAPPA_DIST
  """

  cdef char* error_message = NULL
  result = chi_22_symphony(nu, magnetic_field, electron_density,
                observer_angle, distribution, real_part, theta_e, 
		power_law_p, gamma_min, gamma_max, gamma_cutoff, 
		kappa, kappa_width, &error_message)
  if error_message:
    raise RuntimeError (error_message)
  return result

def chi_23_symphony_py(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
            double theta_e,
            double power_law_p,
            double gamma_min,
            double gamma_max,
            double gamma_cutoff,
            double kappa,
            double kappa_width):

  """Returns chi_23_nu(nu, magnetic_field, electron_density, observer_angle, 
                       distribution, real_part, theta_e, power_law_p, gamma_min, 
                       gamma_max, gamma_cutoff, kappa, kappa_width).
     Keys for distribution functions: symphonyPy.MAXWELL_JUETTNER, 
                                      symphonyPy.POWER_LAW, 
                                      symphonyPy.KAPPA_DIST
  """

  cdef char* error_message = NULL
  result = chi_23_symphony(nu, magnetic_field, electron_density,
                observer_angle, distribution, real_part, theta_e, 
		power_law_p, gamma_min, gamma_max, gamma_cutoff, 
		kappa, kappa_width, &error_message)
  if error_message:
    raise RuntimeError (error_message)
  return result

def chi_31_symphony_py(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
            double theta_e,
            double power_law_p,
            double gamma_min,
            double gamma_max,
            double gamma_cutoff,
            double kappa,
            double kappa_width):

  """Returns chi_31_nu(nu, magnetic_field, electron_density, observer_angle, 
                       distribution, real_part, theta_e, power_law_p, gamma_min, 
                       gamma_max, gamma_cutoff, kappa, kappa_width).
     Keys for distribution functions: symphonyPy.MAXWELL_JUETTNER, 
                                      symphonyPy.POWER_LAW, 
                                      symphonyPy.KAPPA_DIST
  """

  cdef char* error_message = NULL
  result = chi_31_symphony(nu, magnetic_field, electron_density,
                observer_angle, distribution, real_part, theta_e, 
		power_law_p, gamma_min, gamma_max, gamma_cutoff, 
		kappa, kappa_width, &error_message)
  if error_message:
    raise RuntimeError (error_message)
  return result

def chi_32_symphony_py(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
            double theta_e,
            double power_law_p,
            double gamma_min,
            double gamma_max,
            double gamma_cutoff,
            double kappa,
            double kappa_width):

  """Returns chi_32_nu(nu, magnetic_field, electron_density, observer_angle, 
                       distribution, real_part, theta_e, power_law_p, gamma_min, 
                       gamma_max, gamma_cutoff, kappa, kappa_width).
     Keys for distribution functions: symphonyPy.MAXWELL_JUETTNER, 
                                      symphonyPy.POWER_LAW, 
                                      symphonyPy.KAPPA_DIST
  """

  cdef char* error_message = NULL
  result = chi_32_symphony(nu, magnetic_field, electron_density,
                observer_angle, distribution, real_part, theta_e, 
		power_law_p, gamma_min, gamma_max, gamma_cutoff, 
		kappa, kappa_width, &error_message)
  if error_message:
    raise RuntimeError (error_message)
  return result

def chi_33_symphony_py(double nu,
            double magnetic_field,
            double electron_density,
            double observer_angle,
            int distribution,
	    int real_part,
            double theta_e,
            double power_law_p,
            double gamma_min,
            double gamma_max,
            double gamma_cutoff,
            double kappa,
            double kappa_width):

  """Returns chi_33_nu(nu, magnetic_field, electron_density, observer_angle, 
                       distribution, real_part, theta_e, power_law_p, gamma_min, 
                       gamma_max, gamma_cutoff, kappa, kappa_width).
     Keys for distribution functions: symphonyPy.MAXWELL_JUETTNER, 
                                      symphonyPy.POWER_LAW, 
                                      symphonyPy.KAPPA_DIST
  """

  cdef char* error_message = NULL
  result = chi_33_symphony(nu, magnetic_field, electron_density,
                observer_angle, distribution, real_part, theta_e, 
		power_law_p, gamma_min, gamma_max, gamma_cutoff, 
		kappa, kappa_width, &error_message)
  if error_message:
    raise RuntimeError (error_message)
  return result
