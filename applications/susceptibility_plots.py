import numpy as np
import pylab as pl

import sys
sys.path.insert(0, '../src/susceptibility_tensor')
sys.path.insert(0, '../build')
from susceptibility_interpolator import chi_ij
import susceptibility_tensor.susceptibilityPy as susp
import symphonyPy as sp

# Set plot parameters to make beautiful plots
pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 20  
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex']     = True
pl.rcParams['axes.linewidth']  = 1.5
pl.rcParams['axes.titlesize']  = 'medium'
pl.rcParams['axes.labelsize']  = 'large'

pl.rcParams['xtick.major.size'] = 8     
pl.rcParams['xtick.minor.size'] = 4     
pl.rcParams['xtick.major.pad']  = 8     
pl.rcParams['xtick.minor.pad']  = 8     
pl.rcParams['xtick.color']      = 'k'     
pl.rcParams['xtick.labelsize']  = 'large'
pl.rcParams['xtick.direction']  = 'in'    

pl.rcParams['ytick.major.size'] = 8     
pl.rcParams['ytick.minor.size'] = 4     
pl.rcParams['ytick.major.pad']  = 8     
pl.rcParams['ytick.minor.pad']  = 8     
pl.rcParams['ytick.color']      = 'k'     
pl.rcParams['ytick.labelsize']  = 'large'
pl.rcParams['ytick.direction']  = 'in'

#define constants
epsilon0  = 1./(4. * np.pi)
e         = 4.80320680e-10
m         = 9.1093826e-28
c         = 2.99792458e10
epsilon   = -1.

def nu_c(magnetic_field):
    ans = e * magnetic_field / (2. * np.pi * m * c)
    return ans

nuratio = np.logspace(0., 1., 3)

magnetic_field = 30.
nu = nu_c(magnetic_field)
electron_density = 1.
observer_angle = np.pi/3.
distribution = sp.MAXWELL_JUETTNER
polarization = sp.STOKES_I
theta_e = 10.
power_law_p = 3.
gamma_min = 1.
gamma_max = 1000.
gamma_cutoff = 1.e10
kappa = 3.5
kappa_width = 10.
chi_method = sp.SYMPHONY_METHOD

symphony_result = sp.alpha_nu_py(nu, 
				 magnetic_field, 
		 		 electron_density, 
		 		 observer_angle, 
		 		 distribution, 
		 		 polarization, 
		 		 theta_e, 
		 		 power_law_p, 
		 		 gamma_min, 
		 		 gamma_max, 
		 		 gamma_cutoff, 
		 		 kappa, 
		 		 kappa_width,
		 		 chi_method)

def interp_alpha_I(nu,
                   magnetic_field,
                   electron_density,
                   observer_angle,
                   distribution,
                   polarization,
                   theta_e,
                   power_law_p,
                   gamma_min,
                   gamma_max,
                   gamma_cutoff,
                   kappa,
                   kappa_width):
	real_part = 0
	omega = 2. * np.pi * nu
	prefactor = 2. * np.pi * epsilon0 * omega/c
	chi_11 = chi_ij(nu,
                             magnetic_field,
                             electron_density,
                             observer_angle,
                             distribution,
                             real_part,
                             theta_e,
                             power_law_p,
                             gamma_min,
                             gamma_max,
                             gamma_cutoff,
                             kappa,
                             kappa_width,
                             11)
	chi_13 = chi_ij(nu,
                             magnetic_field,
                             electron_density,
                             observer_angle,
                             distribution,
                             real_part,
                             theta_e,
                             power_law_p,
                             gamma_min,
                             gamma_max,
                             gamma_cutoff,
                             kappa,
                             kappa_width,
                             13)
	chi_22 = chi_ij(nu,
                             magnetic_field,
                             electron_density,
                             observer_angle,
                             distribution,
                             real_part,
                             theta_e,
                             power_law_p,
                             gamma_min,
                             gamma_max,
                             gamma_cutoff,
                             kappa,
                             kappa_width,
                             22)
	chi_33 = chi_ij(nu,
                             magnetic_field,
                             electron_density,
                             observer_angle,
                             distribution,
                             real_part,
                             theta_e,
                             power_law_p,
                             gamma_min,
                             gamma_max,
                             gamma_cutoff,
                             kappa,
                             kappa_width,
                             33)

	term11 = (chi_11 * np.cos(observer_angle)**2. 
		  + chi_33*np.sin(observer_angle)**2. 
                  - 2. * chi_13 * np.sin(observer_angle) * np.cos(observer_angle))
	term22 = chi_22
	ans = prefactor * (term11 + term22)
	return ans

interp_result = interp_alpha_I(nu,
                               magnetic_field,
                               electron_density,
                               observer_angle,
                               distribution,
                               polarization,
                               theta_e,
                               power_law_p,
                               gamma_min,
                               gamma_max,
                               gamma_cutoff,
                               kappa,
                               kappa_width)

print symphony_result, interp_result
