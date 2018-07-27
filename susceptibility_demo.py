import build.symphonyPy as sp
import numpy as np
import build.susceptibility_tensor.susceptibilityPy as susp
import sys
sys.path.insert(0, './src/susceptibility_tensor')
sys.path.insert(0, './build')
from susceptibility_interpolator import chi_ij, spline_selector

#define parameters
e = 4.80320680e-10
m = 9.1093826e-28
c = 2.99792458e10
B = 1. 
epsilon0 = 1./(4.*np.pi)
nu_c = (e * B) / (2. * np.pi * m * c)
nu = 100.*nu_c
real = 1
theta = np.pi/3.
n_e = 1.
theta_e = 10.
p = 3.
gamma_min = 1.
gamma_max = 1000.
gamma_cutoff = 1e10
kappa = 3.5
kappa_width = 10.
dist = sp.KAPPA_DIST

#-----load interpolated component-----#
print '-----compare interpolated component with integrated component-----'
component = 11
print 'component:', component, 'real part (1=real, 0=imag):', real
chi_interp = chi_ij(nu, B, n_e, theta, dist, real, theta_e, p, gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width, component)
chi_int    = susp.chi_11_symphony_py(nu, B, n_e, theta, dist, real, theta_e, p, gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width)

print 'interpolated:', chi_interp, 'integrated:', chi_int, 'error:', (chi_interp-chi_int)/chi_int

if(real == 0):
	print '------------------------------------------------------------------'
	print '-----compare integrated alpha_Q to symphony-----'
	print 'NOTE: this can take up to 15 min.'
	print 'defined variables'
	chi_11 = susp.chi_11_symphony_py(nu, B, n_e, theta, dist, real, theta_e, p, gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width)
	print 'chi_11 done,  25% done with chi_ij components'
	chi_13 = susp.chi_13_symphony_py(nu, B, n_e, theta, dist, real, theta_e, p, gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width)
	print 'chi_13 done,  50% done with chi_ij components'
	chi_33 = susp.chi_33_symphony_py(nu, B, n_e, theta, dist, real, theta_e, p, gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width)
	print 'chi_33 done,  75% done with chi_ij components'
	chi_22 = susp.chi_22_symphony_py(nu, B, n_e, theta, dist, real, theta_e, p, gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width)
	print 'chi_22 done, 100% done with chi_ij components'
	
	alpha_Q_chi_int = 2.*np.pi*epsilon0*(2.*np.pi*nu)/c * (chi_11*np.cos(theta)**2. + chi_33*np.sin(theta)**2. - 2.*chi_13*np.sin(theta)*np.cos(theta) - chi_22)

	print 'computing symphony alpha_Q'
	symphony_alpha_Q = sp.alpha_nu_py(nu, B, n_e, theta, dist, sp.STOKES_Q, theta_e, p, gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width, sp.SYMPHONY_METHOD)

	print 'chi_ij method:', alpha_Q_chi_int, 'symphony method:', symphony_alpha_Q, 'error:', (alpha_Q_chi_int - symphony_alpha_Q)/symphony_alpha_Q


else:
	print '------------------------------------------------------------------'
	print '-----compare integrated alpha_V to symphony-----'
	print 'NOTE: this can take up to 10 min.'
	chi_12 = susp.chi_12_symphony_py(nu, B, n_e, theta, dist, real, theta_e, p, gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width)
        print 'chi_12 done,  50% done with chi_ij components'
        chi_32 = susp.chi_32_symphony_py(nu, B, n_e, theta, dist, real, theta_e, p, gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width)
        print 'chi_32 done, 100% done with chi_ij components'
	alpha_V_chi_int = 4.*np.pi*epsilon0*(2.*np.pi*nu)/c * (chi_12*np.cos(theta) - chi_32*np.sin(theta))
	symphony_alpha_V = sp.alpha_nu_py(nu, B, n_e, theta, dist, sp.STOKES_V, theta_e, p, gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width, sp.SYMPHONY_METHOD)
	print 'chi_ij method:', alpha_V_chi_int, 'symphony method:', symphony_alpha_V, 'error:', (alpha_V_chi_int - symphony_alpha_V)/symphony_alpha_V
