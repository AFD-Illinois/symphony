import symphonyPy 
import numpy as np
import pylab as pl

"""This is a short demo on how to use the symphony python interface"""

#Define some sample parameters (in CGS units)
nu = 230e9
magnetic_field = 30.
electron_density = 1.
observer_angle = np.pi/3.
distribution = symphonyPy.MAXWELL_JUETTNER
polarization = symphonyPy.STOKES_I
theta_e = 10.
power_law_p = 3.5
gamma_min = 1.
gamma_max = 1000.
gamma_cutoff = 100000000000.
kappa = 3.5
kappa_width = 10.

#Output
print 'nu        ', 'j_nu_py()           ', 'j_nu_fit_py()'

print nu, 
print '  ',
print symphonyPy.j_nu_py(nu, magnetic_field, electron_density, observer_angle,
                         distribution, polarization, theta_e, power_law_p, 
                         gamma_min, gamma_max, gamma_cutoff, kappa, 
                         kappa_width),
print '  ',
print symphonyPy.j_nu_fit_py(nu, magnetic_field, electron_density, 
                             observer_angle, distribution, polarization,
                             theta_e, power_law_p, gamma_min, gamma_max, 
                             gamma_cutoff, kappa, kappa_width)

x = np.arange(0., np.pi, 0.01)
pl.plot(x, np.sin(x))
pl.show()
