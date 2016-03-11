import sys

#symphony_build_path = '/home/mani/work/symphony/build'
symphony_build_path = '/home/alex/Documents/Spring_2016/1symphony/symphony/build'
sys.path.append(symphony_build_path)

import symphonyPy as sp
import numpy as np
import pylab as pl

#set constant parameters for the calculation
m = 9.1093826e-28
c = 2.99792458e10
electron_charge = 4.80320680e-10
theta_e = 10. #Do we need to do anything about electron temp?
e = 4.80320680e-10
h = 6.6260693e-27
gamma_min = 1.
gamma_max = 1000.
gamma_cutoff = 1e10
power_law_p = 3.5
kappa = 3.5
kappa_width = 10.
nuratio = 1.e2
B_scale = 30. #TODO: check if sim data is actually normalized to init. val.

#import data from Dr. Kunz's simulation
N2 = 1152
N1 = 1152
rank = 0
size = 1

datafiles_path = '/home/alex/Documents/Spring_2016/'
B_x = np.loadtxt(datafiles_path + 'mirror_bx.out') * B_scale
B_y = np.loadtxt(datafiles_path + 'mirror_by.out') * B_scale
B_z = np.loadtxt(datafiles_path + 'mirror_bz.out') * B_scale
B_mag = np.sqrt(B_x**2. + B_y**2. + B_z**2.)
n_e = np.loadtxt(datafiles_path + 'mirror_d.out')

#still thinking about how to do observer_angle.  Is there a better way than this?
obs_angle = np.arccos(1./(1.) * (0.*B_x + -1.*B_y + 0.*B_z)/B_mag)

#Generate all averages
B_x_avg        = np.mean(B_x)
B_y_avg        = np.mean(B_y)
B_z_avg        = np.mean(B_z)
n_e_avg        = np.mean(n_e)
obs_angle_avg  = np.mean(obs_angle)
B_mag_avg      = np.mean(B_mag)
nu_c           = electron_charge * (B_mag * B_scale) / (2. * np.pi * m * c)
nu_c_avg       = np.mean(nu_c)
nu_avg         = nuratio * nu_c_avg

#-------------------------------MJ_I-------------------------------------------#
MJ_I_exact_avg = 0
MJ_I_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        MJ_I_exact[j][i] = sp.j_nu_fit_py(nuratio * nu_c[j][i],
                                          B_mag[j][i],
                                          n_e[j][i], 
                                          obs_angle[j][i],
                                          sp.MAXWELL_JUETTNER, 
                                          sp.STOKES_I, theta_e, power_law_p, 
                                          gamma_min, gamma_max, gamma_cutoff, 
                                          kappa, kappa_width
                                         )


#np.savetxt("MJ_I_using_symphony_fits.txt", MJ_I_exact)
#pl.contourf(MJ_I_exact, 100)
#pl.colorbar()
#pl.show()

MJ_I_exact_avg = np.mean(MJ_I_exact)

MJ_I_avgs  = sp.j_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                            obs_angle_avg, sp.MAXWELL_JUETTNER,
                            sp.STOKES_I, theta_e, power_law_p,
                            gamma_min, gamma_max, gamma_cutoff,
                            kappa, kappa_width)

print 'MJ_I   ', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (MJ_I_exact_avg 
                                                  - MJ_I_avgs)/MJ_I_exact_avg

#-------------------------------MJ_Q-------------------------------------------#
MJ_Q_exact_avg = 0
MJ_Q_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        MJ_Q_exact[j][i] = sp.j_nu_fit_py(nuratio * nu_c[j][i],
                                          B_mag[j][i],
                                          n_e[j][i], 
                                          obs_angle[j][i],
                                          sp.MAXWELL_JUETTNER, 
                                          sp.STOKES_Q, theta_e, power_law_p, 
                                          gamma_min, gamma_max, gamma_cutoff, 
                                          kappa, kappa_width
                                         )


#np.savetxt("MJ_Q_using_symphony_fits.txt", MJ_Q_exact)
#pl.contourf(MJ_Q_exact, 100)
#pl.colorbar()
#pl.show()

MJ_Q_exact_avg = np.mean(MJ_Q_exact)

MJ_Q_avgs  = sp.j_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                            obs_angle_avg, sp.MAXWELL_JUETTNER,
                            sp.STOKES_Q, theta_e, power_law_p,
                            gamma_min, gamma_max, gamma_cutoff,
                            kappa, kappa_width)

print 'MJ_Q   ', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (MJ_Q_exact_avg 
                                                  - MJ_Q_avgs)/MJ_Q_exact_avg

#-------------------------------MJ_V-------------------------------------------#
MJ_V_exact_avg = 0
MJ_V_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        MJ_V_exact[j][i] = sp.j_nu_fit_py(nuratio * nu_c[j][i],
                                          B_mag[j][i],
                                          n_e[j][i], 
                                          obs_angle[j][i],
                                          sp.MAXWELL_JUETTNER, 
                                          sp.STOKES_V, theta_e, power_law_p, 
                                          gamma_min, gamma_max, gamma_cutoff, 
                                          kappa, kappa_width
                                         )


#np.savetxt("MJ_V_using_symphony_fits.txt", MJ_V_exact)
#pl.contourf(MJ_V_exact, 100)
#pl.colorbar()
#pl.show()

MJ_V_exact_avg = np.mean(MJ_V_exact)

MJ_V_avgs  = sp.j_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                            obs_angle_avg, sp.MAXWELL_JUETTNER,
                            sp.STOKES_V, theta_e, power_law_p,
                            gamma_min, gamma_max, gamma_cutoff,
                            kappa, kappa_width)

print 'MJ_V   ', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (MJ_V_exact_avg 
                                                  - MJ_V_avgs)/MJ_V_exact_avg


#-------------------------------PL_I-------------------------------------------#
PL_I_exact_avg = 0
PL_I_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        PL_I_exact[j][i] = sp.j_nu_fit_py(nuratio * nu_c[j][i],
                                          B_mag[j][i],
                                          n_e[j][i], 
                                          obs_angle[j][i],
                                          sp.POWER_LAW, 
                                          sp.STOKES_I, theta_e, power_law_p, 
                                          gamma_min, gamma_max, gamma_cutoff, 
                                          kappa, kappa_width
                                         )


#np.savetxt("PL_I_using_symphony_fits.txt", PL_I_exact)
#pl.contourf(PL_I_exact, 100)
#pl.colorbar()
#pl.show()

PL_I_exact_avg = np.mean(PL_I_exact)

PL_I_avgs  = sp.j_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                            obs_angle_avg, sp.POWER_LAW,
                            sp.STOKES_I, theta_e, power_law_p,
                            gamma_min, gamma_max, gamma_cutoff,
                            kappa, kappa_width)

print 'PL_I   ', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (PL_I_exact_avg 
                                                  - PL_I_avgs)/PL_I_exact_avg

#-------------------------------PL_Q-------------------------------------------#
PL_Q_exact_avg = 0
PL_Q_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        PL_Q_exact[j][i] = sp.j_nu_fit_py(nuratio * nu_c[j][i],
                                          B_mag[j][i],
                                          n_e[j][i], 
                                          obs_angle[j][i],
                                          sp.POWER_LAW, 
                                          sp.STOKES_Q, theta_e, power_law_p, 
                                          gamma_min, gamma_max, gamma_cutoff, 
                                          kappa, kappa_width
                                         )


#np.savetxt("PL_Q_using_symphony_fits.txt", PL_Q_exact)
#pl.contourf(PL_Q_exact, 100)
#pl.colorbar()
#pl.show()

PL_Q_exact_avg = np.mean(PL_Q_exact)

PL_Q_avgs  = sp.j_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                            obs_angle_avg, sp.POWER_LAW,
                            sp.STOKES_Q, theta_e, power_law_p,
                            gamma_min, gamma_max, gamma_cutoff,
                            kappa, kappa_width)

print 'PL_Q   ', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (PL_Q_exact_avg 
                                                  - PL_Q_avgs)/PL_Q_exact_avg

#-------------------------------PL_V-------------------------------------------#
PL_V_exact_avg = 0
PL_V_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        PL_V_exact[j][i] = sp.j_nu_fit_py(nuratio * nu_c[j][i],
                                          B_mag[j][i],
                                          n_e[j][i], 
                                          obs_angle[j][i],
                                          sp.POWER_LAW, 
                                          sp.STOKES_V, theta_e, power_law_p, 
                                          gamma_min, gamma_max, gamma_cutoff, 
                                          kappa, kappa_width
                                         )


#np.savetxt("PL_V_using_symphony_fits.txt", PL_V_exact)
#pl.contourf(PL_V_exact, 100)
#pl.colorbar()
#pl.show()

PL_V_exact_avg = np.mean(PL_V_exact)

PL_V_avgs  = sp.j_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                            obs_angle_avg, sp.POWER_LAW,
                            sp.STOKES_V, theta_e, power_law_p,
                            gamma_min, gamma_max, gamma_cutoff,
                            kappa, kappa_width)

print 'PL_V   ', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (PL_V_exact_avg 
                                                  - PL_V_avgs)/PL_V_exact_avg

#-------------------------------kappa_I-------------------------------------------#
kappa_I_exact_avg = 0
kappa_I_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        kappa_I_exact[j][i] = sp.j_nu_fit_py(nuratio * nu_c[j][i],
                                             B_mag[j][i],
                                             n_e[j][i], 
                                             obs_angle[j][i],
                                             sp.KAPPA_DIST, 
                                             sp.STOKES_I, theta_e, power_law_p, 
                                             gamma_min, gamma_max, gamma_cutoff, 
                                             kappa, kappa_width
                                            )


#np.savetxt("kappa_I_using_symphony_fits.txt", kappa_I_exact)
#pl.contourf(kappa_I_exact, 100)
#pl.colorbar()
#pl.show()

kappa_I_exact_avg = np.mean(kappa_I_exact)

kappa_I_avgs  = sp.j_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                               obs_angle_avg, sp.KAPPA_DIST,
                               sp.STOKES_I, theta_e, power_law_p,
                               gamma_min, gamma_max, gamma_cutoff,
                               kappa, kappa_width)

print 'kappa_I', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (kappa_I_exact_avg 
                                                     - kappa_I_avgs)/kappa_I_exact_avg

#-------------------------------kappa_Q-------------------------------------------#
kappa_Q_exact_avg = 0
kappa_Q_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        kappa_Q_exact[j][i] = sp.j_nu_fit_py(nuratio * nu_c[j][i],
                                             B_mag[j][i],
                                             n_e[j][i], 
                                             obs_angle[j][i],
                                             sp.KAPPA_DIST, 
                                             sp.STOKES_Q, theta_e, power_law_p, 
                                             gamma_min, gamma_max, gamma_cutoff, 
                                             kappa, kappa_width
                                            )


#np.savetxt("kappa_Q_using_symphony_fits.txt", kappa_Q_exact)
#pl.contourf(kappa_Q_exact, 100)
#pl.colorbar()
#pl.show()

kappa_Q_exact_avg = np.mean(kappa_Q_exact)

kappa_Q_avgs  = sp.j_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                               obs_angle_avg, sp.KAPPA_DIST,
                               sp.STOKES_Q, theta_e, power_law_p,
                               gamma_min, gamma_max, gamma_cutoff,
                               kappa, kappa_width)

print 'kappa_Q', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (kappa_Q_exact_avg 
                                                     - kappa_Q_avgs)/kappa_Q_exact_avg

#-------------------------------kappa_V-------------------------------------------#
kappa_V_exact_avg = 0
kappa_V_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        kappa_V_exact[j][i] = sp.j_nu_fit_py(nuratio * nu_c[j][i],
                                             B_mag[j][i],
                                             n_e[j][i], 
                                             obs_angle[j][i],
                                             sp.KAPPA_DIST, 
                                             sp.STOKES_V, theta_e, power_law_p, 
                                             gamma_min, gamma_max, gamma_cutoff, 
                                             kappa, kappa_width
                                            )


#np.savetxt("kappa_V_using_symphony_fits.txt", kappa_V_exact)
#pl.contourf(kappa_V_exact, 100)
#pl.colorbar()
#pl.show()

kappa_V_exact_avg = np.mean(kappa_V_exact)

kappa_V_avgs  = sp.j_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                               obs_angle_avg, sp.KAPPA_DIST,
                               sp.STOKES_V, theta_e, power_law_p,
                               gamma_min, gamma_max, gamma_cutoff,
                               kappa, kappa_width)

print 'kappa_V', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (kappa_V_exact_avg 
                                                     - kappa_V_avgs)/kappa_V_exact_avg


#---------------------------ABSORPTIVITY---------------------------------------#

print ''

#-------------------------------MJ_I-------------------------------------------#
MJ_I_exact_avg = 0
MJ_I_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        MJ_I_exact[j][i] = sp.alpha_nu_fit_py(nuratio * nu_c[j][i],
                                          B_mag[j][i],
                                          n_e[j][i], 
                                          obs_angle[j][i],
                                          sp.MAXWELL_JUETTNER, 
                                          sp.STOKES_I, theta_e, power_law_p, 
                                          gamma_min, gamma_max, gamma_cutoff, 
                                          kappa, kappa_width
                                         )


#np.savetxt("MJ_I_using_symphony_fits.txt", MJ_I_exact)
pl.contourf(MJ_I_exact, 100)
pl.colorbar()
pl.show()

MJ_I_exact_avg = np.mean(MJ_I_exact)

MJ_I_avgs  = sp.alpha_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                            obs_angle_avg, sp.MAXWELL_JUETTNER,
                            sp.STOKES_I, theta_e, power_law_p,
                            gamma_min, gamma_max, gamma_cutoff,
                            kappa, kappa_width)

print 'MJ_I_abs   ', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (MJ_I_exact_avg 
                                                  - MJ_I_avgs)/MJ_I_exact_avg

#-------------------------------MJ_Q-------------------------------------------#
MJ_Q_exact_avg = 0
MJ_Q_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        MJ_Q_exact[j][i] = sp.alpha_nu_fit_py(nuratio * nu_c[j][i],
                                          B_mag[j][i],
                                          n_e[j][i], 
                                          obs_angle[j][i],
                                          sp.MAXWELL_JUETTNER, 
                                          sp.STOKES_Q, theta_e, power_law_p, 
                                          gamma_min, gamma_max, gamma_cutoff, 
                                          kappa, kappa_width
                                         )


#np.savetxt("MJ_Q_using_symphony_fits.txt", MJ_Q_exact)
#pl.contourf(MJ_Q_exact, 100)
#pl.colorbar()
#pl.show()

MJ_Q_exact_avg = np.mean(MJ_Q_exact)

MJ_Q_avgs  = sp.alpha_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                            obs_angle_avg, sp.MAXWELL_JUETTNER,
                            sp.STOKES_Q, theta_e, power_law_p,
                            gamma_min, gamma_max, gamma_cutoff,
                            kappa, kappa_width)

print 'MJ_Q_abs   ', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (MJ_Q_exact_avg 
                                                  - MJ_Q_avgs)/MJ_Q_exact_avg

#-------------------------------MJ_V-------------------------------------------#
MJ_V_exact_avg = 0
MJ_V_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        MJ_V_exact[j][i] = sp.alpha_nu_fit_py(nuratio * nu_c[j][i],
                                          B_mag[j][i],
                                          n_e[j][i], 
                                          obs_angle[j][i],
                                          sp.MAXWELL_JUETTNER, 
                                          sp.STOKES_V, theta_e, power_law_p, 
                                          gamma_min, gamma_max, gamma_cutoff, 
                                          kappa, kappa_width
                                         )


#np.savetxt("MJ_V_using_symphony_fits.txt", MJ_V_exact)
#pl.contourf(MJ_V_exact, 100)
#pl.colorbar()
#pl.show()

MJ_V_exact_avg = np.mean(MJ_V_exact)

MJ_V_avgs  = sp.alpha_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                            obs_angle_avg, sp.MAXWELL_JUETTNER,
                            sp.STOKES_V, theta_e, power_law_p,
                            gamma_min, gamma_max, gamma_cutoff,
                            kappa, kappa_width)

print 'MJ_V_abs   ', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (MJ_V_exact_avg 
                                                  - MJ_V_avgs)/MJ_V_exact_avg


#-------------------------------PL_I-------------------------------------------#
PL_I_exact_avg = 0
PL_I_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        PL_I_exact[j][i] = sp.alpha_nu_fit_py(nuratio * nu_c[j][i],
                                          B_mag[j][i],
                                          n_e[j][i], 
                                          obs_angle[j][i],
                                          sp.POWER_LAW, 
                                          sp.STOKES_I, theta_e, power_law_p, 
                                          gamma_min, gamma_max, gamma_cutoff, 
                                          kappa, kappa_width
                                         )


#np.savetxt("PL_I_using_symphony_fits.txt", PL_I_exact)
#pl.contourf(PL_I_exact, 100)
#pl.colorbar()
#pl.show()

PL_I_exact_avg = np.mean(PL_I_exact)

PL_I_avgs  = sp.alpha_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                            obs_angle_avg, sp.POWER_LAW,
                            sp.STOKES_I, theta_e, power_law_p,
                            gamma_min, gamma_max, gamma_cutoff,
                            kappa, kappa_width)

print 'PL_I_abs   ', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (PL_I_exact_avg 
                                                  - PL_I_avgs)/PL_I_exact_avg

#-------------------------------PL_Q-------------------------------------------#
PL_Q_exact_avg = 0
PL_Q_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        PL_Q_exact[j][i] = sp.alpha_nu_fit_py(nuratio * nu_c[j][i],
                                          B_mag[j][i],
                                          n_e[j][i], 
                                          obs_angle[j][i],
                                          sp.POWER_LAW, 
                                          sp.STOKES_Q, theta_e, power_law_p, 
                                          gamma_min, gamma_max, gamma_cutoff, 
                                          kappa, kappa_width
                                         )


#np.savetxt("PL_Q_using_symphony_fits.txt", PL_Q_exact)
#pl.contourf(PL_Q_exact, 100)
#pl.colorbar()
#pl.show()

PL_Q_exact_avg = np.mean(PL_Q_exact)

PL_Q_avgs  = sp.alpha_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                            obs_angle_avg, sp.POWER_LAW,
                            sp.STOKES_Q, theta_e, power_law_p,
                            gamma_min, gamma_max, gamma_cutoff,
                            kappa, kappa_width)

print 'PL_Q_abs   ', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (PL_Q_exact_avg 
                                                  - PL_Q_avgs)/PL_Q_exact_avg

#-------------------------------PL_V-------------------------------------------#
PL_V_exact_avg = 0
PL_V_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        PL_V_exact[j][i] = sp.alpha_nu_fit_py(nuratio * nu_c[j][i],
                                          B_mag[j][i],
                                          n_e[j][i], 
                                          obs_angle[j][i],
                                          sp.POWER_LAW, 
                                          sp.STOKES_V, theta_e, power_law_p, 
                                          gamma_min, gamma_max, gamma_cutoff, 
                                          kappa, kappa_width
                                         )


#np.savetxt("PL_V_using_symphony_fits.txt", PL_V_exact)
#pl.contourf(PL_V_exact, 100)
#pl.colorbar()
#pl.show()

PL_V_exact_avg = np.mean(PL_V_exact)

PL_V_avgs  = sp.alpha_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                            obs_angle_avg, sp.POWER_LAW,
                            sp.STOKES_V, theta_e, power_law_p,
                            gamma_min, gamma_max, gamma_cutoff,
                            kappa, kappa_width)

print 'PL_V_abs   ', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (PL_V_exact_avg 
                                                  - PL_V_avgs)/PL_V_exact_avg

#-------------------------------kappa_I-------------------------------------------#
kappa_I_exact_avg = 0
kappa_I_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        kappa_I_exact[j][i] = sp.alpha_nu_fit_py(nuratio * nu_c[j][i],
                                             B_mag[j][i],
                                             n_e[j][i], 
                                             obs_angle[j][i],
                                             sp.KAPPA_DIST, 
                                             sp.STOKES_I, theta_e, power_law_p, 
                                             gamma_min, gamma_max, gamma_cutoff, 
                                             kappa, kappa_width
                                            )


#np.savetxt("kappa_I_using_symphony_fits.txt", kappa_I_exact)
#pl.contourf(kappa_I_exact, 100)
#pl.colorbar()
#pl.show()

kappa_I_exact_avg = np.mean(kappa_I_exact)

kappa_I_avgs  = sp.alpha_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                               obs_angle_avg, sp.KAPPA_DIST,
                               sp.STOKES_I, theta_e, power_law_p,
                               gamma_min, gamma_max, gamma_cutoff,
                               kappa, kappa_width)

print 'kappa_I_abs', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (kappa_I_exact_avg 
                                                     - kappa_I_avgs)/kappa_I_exact_avg

#-------------------------------kappa_Q-------------------------------------------#
kappa_Q_exact_avg = 0
kappa_Q_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        kappa_Q_exact[j][i] = sp.alpha_nu_fit_py(nuratio * nu_c[j][i],
                                             B_mag[j][i],
                                             n_e[j][i], 
                                             obs_angle[j][i],
                                             sp.KAPPA_DIST, 
                                             sp.STOKES_Q, theta_e, power_law_p, 
                                             gamma_min, gamma_max, gamma_cutoff, 
                                             kappa, kappa_width
                                            )


#np.savetxt("kappa_Q_using_symphony_fits.txt", kappa_Q_exact)
#pl.contourf(kappa_Q_exact, 100)
#pl.colorbar()
#pl.show()

kappa_Q_exact_avg = np.mean(kappa_Q_exact)

kappa_Q_avgs  = sp.alpha_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                               obs_angle_avg, sp.KAPPA_DIST,
                               sp.STOKES_Q, theta_e, power_law_p,
                               gamma_min, gamma_max, gamma_cutoff,
                               kappa, kappa_width)

print 'kappa_Q_abs', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (kappa_Q_exact_avg 
                                                     - kappa_Q_avgs)/kappa_Q_exact_avg

#-------------------------------kappa_V-------------------------------------------#
kappa_V_exact_avg = 0
kappa_V_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

for j in range(jIndexStart, jIndexEnd):
#    print "j = ", j
    for i in range(0, N1):
        kappa_V_exact[j][i] = sp.alpha_nu_fit_py(nuratio * nu_c[j][i],
                                             B_mag[j][i],
                                             n_e[j][i], 
                                             obs_angle[j][i],
                                             sp.KAPPA_DIST, 
                                             sp.STOKES_V, theta_e, power_law_p, 
                                             gamma_min, gamma_max, gamma_cutoff, 
                                             kappa, kappa_width
                                            )


#np.savetxt("kappa_V_using_symphony_fits.txt", kappa_V_exact)
#pl.contourf(kappa_V_exact, 100)
#pl.colorbar()
#pl.show()

kappa_V_exact_avg = np.mean(kappa_V_exact)

kappa_V_avgs  = sp.alpha_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                               obs_angle_avg, sp.KAPPA_DIST,
                               sp.STOKES_V, theta_e, power_law_p,
                               gamma_min, gamma_max, gamma_cutoff,
                               kappa, kappa_width)

print 'kappa_V_abs', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (kappa_V_exact_avg 
                                                     - kappa_V_avgs)/kappa_V_exact_avg


