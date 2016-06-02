import sys
#symphony_build_path = '/home/mani/work/symphony/build'
symphony_build_path = '/home/alex/Documents/Spring_2016/symphony/symphony/build'
sys.path.append(symphony_build_path)
import symphonyPy as sp
import numpy as np

#--------------------Testing parameters---------------------------------------#
theta_e      = 10.
gamma_min    = 1.
gamma_max    = 1000.
gamma_cutoff = 1e10
power_law_p  = 2.5
kappa        = 3.5
kappa_width  = 10.
B            = 30.
n_e          = 1.
nu           = 230.e9
obs_angle    = np.pi/3.
#-----------------------------------------------------------------------------#

#----------------------expected values and actual values----------------------#
MJ_I_exp_fit    = 1.31261797658e-22
MJ_Q_exp_fit    = -1.04527441106e-22 
MJ_V_exp_fit    = -2.81598011687e-24
PL_I_exp_fit    = 2.03808334021e-24
PL_Q_exp_fit    = -1.04527441106e-22
PL_V_exp_fit    = -3.88375258782e-26
Kappa_I_exp_fit = 2.77087098467e-22
Kappa_Q_exp_fit = -1.79193573531e-22
Kappa_V_exp_fit = -4.31126410186e-24

MJ_I_exp_fit_abs    = 1.36195875355e-19
MJ_Q_exp_fit_abs    = -1.0845658519e-19
MJ_V_exp_fit_abs    =-2.92183166648e-21
PL_I_exp_fit_abs    = 1.90379277421e-21
PL_Q_exp_fit_abs    = -1.47291138324e-21
PL_V_exp_fit_abs    = -4.24640504782e-23
Kappa_I_exp_fit_abs = 1.01739581553e-19
Kappa_Q_exp_fit_abs = -7.08859966782e-20
Kappa_V_exp_fit_abs = -1.65908484702e-21

MJ_I_exp    = 1.29871898683e-22
MJ_Q_exp    = -9.87577471434e-23
MJ_V_exp    = -2.7138066784e-24
PL_I_exp    = 2.03654269253e-24
PL_Q_exp    = -1.47432480277e-24
PL_V_exp    = -3.86129355636e-26
Kappa_I_exp = 2.78933923611e-22
Kappa_Q_exp = -1.7966105582e-22
Kappa_V_exp = -3.84237974681e-24

MJ_I_exp_abs    = 1.34753662893e-19
MJ_Q_exp_abs    = -1.02469959257e-19
MJ_V_exp_abs    = -2.81581615428e-21
PL_I_exp_abs    = 1.90164946304e-21
PL_Q_exp_abs    = -1.46655689386e-21
PL_V_exp_abs    = -4.27534148958e-23
Kappa_I_exp_abs = 1.03817164764e-19
Kappa_Q_exp_abs = -7.09134633078e-20
Kappa_V_exp_abs = -1.66592280833e-21

#----------------------------tests--------------------------------------------#
print '-------------------------------------------------------------------'
print '                symphony automated testing program               '
print '-------------------------------------------------------------------'
print ''

print 'Testing emissivity fitting formulae'
print '-------------------------------------------------------------------'
print 'Maxwell-Juettner Emissivities'

MJ_I_fit = sp.j_nu_fit_py(nu, B, n_e, obs_angle, sp.MAXWELL_JUETTNER,
                          sp.STOKES_I, theta_e, power_law_p, gamma_min,
                          gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(MJ_I_fit - MJ_I_exp_fit)/MJ_I_exp_fit > 0.01):
	print 'STOKES_I                                     FAIL'
else:
        print 'STOKES_I                                     PASS' 

MJ_Q_fit = sp.j_nu_fit_py(nu, B, n_e, obs_angle, sp.MAXWELL_JUETTNER,
                          sp.STOKES_Q, theta_e, power_law_p, gamma_min,
                          gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(MJ_Q_fit - MJ_Q_exp_fit)/MJ_Q_exp_fit > 0.01):
        print 'STOKES_Q                                     FAIL'
else:
        print 'STOKES_Q                                     PASS'

MJ_V_fit = sp.j_nu_fit_py(nu, B, n_e, obs_angle, sp.MAXWELL_JUETTNER,
                          sp.STOKES_V, theta_e, power_law_p, gamma_min,
                          gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(MJ_V_fit - MJ_V_exp_fit)/MJ_V_exp_fit > 0.01):
        print 'STOKES_V                                     FAIL'
else:
        print 'STOKES_V                                     PASS'

print ''
print 'Power-law Emissivities'

PL_I_fit = sp.j_nu_fit_py(nu, B, n_e, obs_angle, sp.POWER_LAW,
                          sp.STOKES_I, theta_e, power_law_p, gamma_min,
                          gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(PL_I_fit - PL_I_exp_fit)/PL_I_exp_fit > 0.01):
        print 'STOKES_I                                     FAIL'
else:
        print 'STOKES_I                                     PASS'

PL_Q_fit = sp.j_nu_fit_py(nu, B, n_e, obs_angle, sp.POWER_LAW,
                          sp.STOKES_Q, theta_e, power_law_p, gamma_min,
                          gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(PL_Q_fit - PL_Q_exp_fit)/PL_Q_exp_fit > 0.01):
        print 'STOKES_Q                                     FAIL'
else:
        print 'STOKES_Q                                     PASS'

PL_V_fit = sp.j_nu_fit_py(nu, B, n_e, obs_angle, sp.POWER_LAW,
                          sp.STOKES_V, theta_e, power_law_p, gamma_min,
                          gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(PL_V_fit - PL_V_exp_fit)/PL_V_exp_fit > 0.01):
        print 'STOKES_V                                     FAIL'
else:
        print 'STOKES_V                                     PASS'

print ''
print 'Kappa Emissivities'

Kappa_I_fit = sp.j_nu_fit_py(nu, B, n_e, obs_angle, sp.KAPPA_DIST,
                             sp.STOKES_I, theta_e, power_law_p, gamma_min,
                             gamma_max, gamma_cutoff, kappa, kappa_width)

if(np.abs(Kappa_I_fit - Kappa_I_exp_fit)/Kappa_I_exp_fit > 0.01):
        print 'STOKES_I                                     FAIL'
else:
        print 'STOKES_I                                     PASS'

Kappa_Q_fit = sp.j_nu_fit_py(nu, B, n_e, obs_angle, sp.KAPPA_DIST,
                             sp.STOKES_Q, theta_e, power_law_p, gamma_min,
                             gamma_max, gamma_cutoff, kappa, kappa_width)

if(np.abs(Kappa_Q_fit - Kappa_Q_exp_fit)/Kappa_Q_exp_fit > 0.01):
        print 'STOKES_Q                                     FAIL'
else:
        print 'STOKES_Q                                     PASS'

Kappa_V_fit = sp.j_nu_fit_py(nu, B, n_e, obs_angle, sp.KAPPA_DIST,
                             sp.STOKES_V, theta_e, power_law_p, gamma_min,
                             gamma_max, gamma_cutoff, kappa, kappa_width)

if(np.abs(Kappa_V_fit - Kappa_V_exp_fit)/Kappa_V_exp_fit > 0.01):
        print 'STOKES_V                                     FAIL'
else:
        print 'STOKES_V                                     PASS'

print ''
print 'Maxwell-Juettner Absorptivities'

MJ_I_fit_abs = sp.alpha_nu_fit_py(nu, B, n_e, obs_angle, sp.MAXWELL_JUETTNER,
                          sp.STOKES_I, theta_e, power_law_p, gamma_min,
                          gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(MJ_I_fit_abs - MJ_I_exp_fit_abs)/MJ_I_exp_fit_abs > 0.01):
	print 'STOKES_I                                     FAIL'
else:
        print 'STOKES_I                                     PASS' 

MJ_Q_fit_abs = sp.alpha_nu_fit_py(nu, B, n_e, obs_angle, sp.MAXWELL_JUETTNER,
                          sp.STOKES_Q, theta_e, power_law_p, gamma_min,
                          gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(MJ_Q_fit_abs - MJ_Q_exp_fit_abs)/MJ_Q_exp_fit_abs > 0.01):
        print 'STOKES_Q                                     FAIL'
else:
        print 'STOKES_Q                                     PASS'

MJ_V_fit_abs = sp.alpha_nu_fit_py(nu, B, n_e, obs_angle, sp.MAXWELL_JUETTNER,
                          sp.STOKES_V, theta_e, power_law_p, gamma_min,
                          gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(MJ_V_fit_abs - MJ_V_exp_fit_abs)/MJ_V_exp_fit_abs > 0.01):
        print 'STOKES_V                                     FAIL'
else:
        print 'STOKES_V                                     PASS'

print ''
print 'Power-law Absorptivities'

PL_I_fit_abs = sp.alpha_nu_fit_py(nu, B, n_e, obs_angle, sp.POWER_LAW,
                          sp.STOKES_I, theta_e, power_law_p, gamma_min,
                          gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(PL_I_fit_abs - PL_I_exp_fit_abs)/PL_I_exp_fit_abs > 0.01):
        print 'STOKES_I                                     FAIL'
else:
        print 'STOKES_I                                     PASS'

PL_Q_fit_abs = sp.alpha_nu_fit_py(nu, B, n_e, obs_angle, sp.POWER_LAW,
                          sp.STOKES_Q, theta_e, power_law_p, gamma_min,
                          gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(PL_Q_fit_abs - PL_Q_exp_fit_abs)/PL_Q_exp_fit_abs > 0.01):
        print 'STOKES_Q                                     FAIL'
else:
        print 'STOKES_Q                                     PASS'

PL_V_fit_abs = sp.alpha_nu_fit_py(nu, B, n_e, obs_angle, sp.POWER_LAW,
                          sp.STOKES_V, theta_e, power_law_p, gamma_min,
                          gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(PL_V_fit_abs - PL_V_exp_fit_abs)/PL_V_exp_fit_abs > 0.01):
        print 'STOKES_V                                     FAIL'
else:
        print 'STOKES_V                                     PASS'

print ''
print 'Kappa Absorptivities'

Kappa_I_fit_abs = sp.alpha_nu_fit_py(nu, B, n_e, obs_angle, sp.KAPPA_DIST,
                             sp.STOKES_I, theta_e, power_law_p, gamma_min,
                             gamma_max, gamma_cutoff, kappa, kappa_width)

if(np.abs(Kappa_I_fit_abs - Kappa_I_exp_fit_abs)/Kappa_I_exp_fit_abs > 0.01):
        print 'STOKES_I                                     FAIL'
else:
        print 'STOKES_I                                     PASS'

Kappa_Q_fit_abs = sp.alpha_nu_fit_py(nu, B, n_e, obs_angle, sp.KAPPA_DIST,
                             sp.STOKES_Q, theta_e, power_law_p, gamma_min,
                             gamma_max, gamma_cutoff, kappa, kappa_width)

if(np.abs(Kappa_Q_fit_abs - Kappa_Q_exp_fit_abs)/Kappa_Q_exp_fit_abs > 0.01):
        print 'STOKES_Q                                     FAIL'
else:
        print 'STOKES_Q                                     PASS'

Kappa_V_fit_abs = sp.alpha_nu_fit_py(nu, B, n_e, obs_angle, sp.KAPPA_DIST,
                             sp.STOKES_V, theta_e, power_law_p, gamma_min,
                             gamma_max, gamma_cutoff, kappa, kappa_width)

if(np.abs(Kappa_V_fit_abs - Kappa_V_exp_fit_abs)/Kappa_V_exp_fit_abs > 0.01):
        print 'STOKES_V                                     FAIL'
else:
        print 'STOKES_V                                     PASS'


print ''
print 'Testing integrated values'
print '-------------------------------------------------------------------'

print 'Maxwell-Juettner Emissivities'

MJ_I = sp.j_nu_py(nu, B, n_e, obs_angle, sp.MAXWELL_JUETTNER,
                  sp.STOKES_I, theta_e, power_law_p, gamma_min,
                  gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(MJ_I - MJ_I_exp)/MJ_I_exp > 0.01):
	print 'STOKES_I                                     FAIL'
else:
        print 'STOKES_I                                     PASS'

MJ_Q = sp.j_nu_py(nu, B, n_e, obs_angle, sp.MAXWELL_JUETTNER,
                  sp.STOKES_Q, theta_e, power_law_p, gamma_min,
                  gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(MJ_Q - MJ_Q_exp)/MJ_Q_exp > 0.01):
        print 'STOKES_Q                                     FAIL'
else:
        print 'STOKES_Q                                     PASS'

MJ_V = sp.j_nu_py(nu, B, n_e, obs_angle, sp.MAXWELL_JUETTNER,
                  sp.STOKES_V, theta_e, power_law_p, gamma_min,
                  gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(MJ_V - MJ_V_exp)/MJ_V_exp > 0.01):
        print 'STOKES_V                                     FAIL'
else:
        print 'STOKES_V                                     PASS'

print ''
print 'Power-law Emissivities'

PL_I = sp.j_nu_py(nu, B, n_e, obs_angle, sp.POWER_LAW,
                  sp.STOKES_I, theta_e, power_law_p, gamma_min,
                  gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(PL_I - PL_I_exp)/PL_I_exp > 0.01):
        print 'STOKES_I                                     FAIL'
else:
        print 'STOKES_I                                     PASS'

PL_Q = sp.j_nu_py(nu, B, n_e, obs_angle, sp.POWER_LAW,
                  sp.STOKES_Q, theta_e, power_law_p, gamma_min,
                  gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(PL_Q - PL_Q_exp)/PL_Q_exp > 0.01):
        print 'STOKES_Q                                     FAIL'
else:
        print 'STOKES_Q                                     PASS'

PL_V = sp.j_nu_py(nu, B, n_e, obs_angle, sp.POWER_LAW,
                  sp.STOKES_V, theta_e, power_law_p, gamma_min,
                  gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(PL_V - PL_V_exp)/PL_V_exp > 0.01):
        print 'STOKES_V                                     FAIL'
else:
        print 'STOKES_V                                     PASS'

print ''
print 'Kappa Emissivities'

Kappa_I = sp.j_nu_py(nu, B, n_e, obs_angle, sp.KAPPA_DIST,
                     sp.STOKES_I, theta_e, power_law_p, gamma_min,
                     gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(Kappa_I - Kappa_I_exp)/Kappa_I_exp > 0.01):
        print 'STOKES_I                                     FAIL'
else:
        print 'STOKES_I                                     PASS'

Kappa_Q = sp.j_nu_py(nu, B, n_e, obs_angle, sp.KAPPA_DIST,
                     sp.STOKES_Q, theta_e, power_law_p, gamma_min,
                     gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(Kappa_Q - Kappa_Q_exp)/Kappa_Q_exp > 0.01):
        print 'STOKES_Q                                     FAIL'
else:
        print 'STOKES_Q                                     PASS'

Kappa_V = sp.j_nu_py(nu, B, n_e, obs_angle, sp.KAPPA_DIST,
                     sp.STOKES_V, theta_e, power_law_p, gamma_min,
                     gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(Kappa_V - Kappa_V_exp)/Kappa_V_exp > 0.01):
        print 'STOKES_V                                     FAIL'
else:
        print 'STOKES_V                                     PASS'

print ''
print 'Maxwell-Juettner Absorptivities'

MJ_I_abs = sp.alpha_nu_py(nu, B, n_e, obs_angle, sp.MAXWELL_JUETTNER,
                  sp.STOKES_I, theta_e, power_law_p, gamma_min,
                  gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(MJ_I_abs - MJ_I_exp_abs)/MJ_I_exp_abs > 0.01):
	print 'STOKES_I                                     FAIL'
else:
        print 'STOKES_I                                     PASS'

MJ_Q_abs = sp.alpha_nu_py(nu, B, n_e, obs_angle, sp.MAXWELL_JUETTNER,
                  sp.STOKES_Q, theta_e, power_law_p, gamma_min,
                  gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(MJ_Q_abs - MJ_Q_exp_abs)/MJ_Q_exp_abs > 0.01):
        print 'STOKES_Q                                     FAIL'
else:
        print 'STOKES_Q                                     PASS'

MJ_V_abs = sp.alpha_nu_py(nu, B, n_e, obs_angle, sp.MAXWELL_JUETTNER,
                  sp.STOKES_V, theta_e, power_law_p, gamma_min,
                  gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(MJ_V_abs - MJ_V_exp_abs)/MJ_V_exp_abs > 0.01):
        print 'STOKES_V                                     FAIL'
else:
        print 'STOKES_V                                     PASS'

print ''
print 'Power-law Absorptivities'

PL_I_abs = sp.alpha_nu_py(nu, B, n_e, obs_angle, sp.POWER_LAW,
                  sp.STOKES_I, theta_e, power_law_p, gamma_min,
                  gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(PL_I_abs - PL_I_exp_abs)/PL_I_exp_abs > 0.01):
        print 'STOKES_I                                     FAIL'
else:
        print 'STOKES_I                                     PASS'

PL_Q_abs = sp.alpha_nu_py(nu, B, n_e, obs_angle, sp.POWER_LAW,
                  sp.STOKES_Q, theta_e, power_law_p, gamma_min,
                  gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(PL_Q_abs - PL_Q_exp_abs)/PL_Q_exp_abs > 0.01):
        print 'STOKES_Q                                     FAIL'
else:
        print 'STOKES_Q                                     PASS'

PL_V_abs = sp.alpha_nu_py(nu, B, n_e, obs_angle, sp.POWER_LAW,
                  sp.STOKES_V, theta_e, power_law_p, gamma_min,
                  gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(PL_V_abs - PL_V_exp_abs)/PL_V_exp_abs > 0.01):
        print 'STOKES_V                                     FAIL'
else:
        print 'STOKES_V                                     PASS'

print ''
print 'Kappa Absorptivities'

Kappa_I_abs = sp.alpha_nu_py(nu, B, n_e, obs_angle, sp.KAPPA_DIST,
                     sp.STOKES_I, theta_e, power_law_p, gamma_min,
                     gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(Kappa_I_abs - Kappa_I_exp_abs)/Kappa_I_exp_abs > 0.01):
        print 'STOKES_I                                     FAIL'
else:
        print 'STOKES_I                                     PASS'

Kappa_Q_abs = sp.alpha_nu_py(nu, B, n_e, obs_angle, sp.KAPPA_DIST,
                     sp.STOKES_Q, theta_e, power_law_p, gamma_min,
                     gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(Kappa_Q_abs - Kappa_Q_exp_abs)/Kappa_Q_exp_abs > 0.01):
        print 'STOKES_Q                                     FAIL'
else:
        print 'STOKES_Q                                     PASS'

Kappa_V_abs = sp.alpha_nu_py(nu, B, n_e, obs_angle, sp.KAPPA_DIST,
                     sp.STOKES_V, theta_e, power_law_p, gamma_min,
                     gamma_max, gamma_cutoff, kappa, kappa_width)
if(np.abs(Kappa_V_abs - Kappa_V_exp_abs)/Kappa_V_exp_abs > 0.01):
        print 'STOKES_V                                     FAIL'
else:
        print 'STOKES_V                                     PASS'
