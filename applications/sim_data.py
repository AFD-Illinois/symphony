#file to calculate <j_nu()> vs. j_nu(<avgs>)
import sys

#symphony_build_path = '/home/mani/work/symphony/build'
#symphony_build_path = '/home/alex/Documents/Spring_2016/1symphony/symphony/build'
symphony_build_path = '../build'
sys.path.append(symphony_build_path)

import symphonyPy as sp
import numpy as np
import pylab as pl
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()
name = MPI.Get_processor_name()

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


num_skip = 2

N2 = N2 / num_skip
N1 = N1 / num_skip

#rank = 0
#size = 1

comm.Barrier()
if (rank==0):
#    datafiles_path = '/home/mani/work/kunz_data/'
#    datafiles_path = '/home/alex/Documents/Spring_2016/'
    datafiles_path = ''
    B_x = np.loadtxt(datafiles_path + 'mirror_bx.out')[::num_skip, ::num_skip] * B_scale
    B_y = np.loadtxt(datafiles_path + 'mirror_by.out')[::num_skip, ::num_skip] * B_scale
    B_z = np.loadtxt(datafiles_path + 'mirror_bz.out')[::num_skip, ::num_skip] * B_scale
    B_mag = np.sqrt(B_x**2. + B_y**2. + B_z**2.)
    n_e = np.loadtxt(datafiles_path + 'mirror_d.out')[::num_skip, ::num_skip]

    B_x   = np.ascontiguousarray(B_x)
    B_y   = np.ascontiguousarray(B_y)
    B_z   = np.ascontiguousarray(B_z)
    B_mag = np.ascontiguousarray(B_mag)
    n_e   = np.ascontiguousarray(n_e)

else:
    B_x = np.zeros([N2, N1])
    B_y = np.zeros([N2, N1])
    B_z = np.zeros([N2, N1])
    B_mag = np.zeros([N2, N1])
    n_e = np.zeros([N2, N1])

comm.Bcast([B_x, MPI.DOUBLE])
comm.Bcast([B_y, MPI.DOUBLE])
comm.Bcast([B_z, MPI.DOUBLE])
comm.Bcast([B_mag, MPI.DOUBLE])
comm.Bcast([n_e, MPI.DOUBLE])


#still thinking about how to do observer_angle.  Is there a better way than this?
obs_angle = np.arccos(1./(1.) * (0.*B_x + -1.*B_y + 0.*B_z)/B_mag)

#Generate all averages
B_x_avg        = np.mean(B_x)
B_y_avg        = np.mean(B_y)
B_z_avg        = np.mean(B_z)
n_e_avg        = np.mean(n_e)
obs_angle_avg  = np.mean(obs_angle)

print "rank = ", rank, " of ", size, "procs"
#-------------------------------MJ_I-------------------------------------------#
MJ_I_exact_avg = 0
MJ_I_exact = np.zeros([N2, N1])

jIndexStart =  rank    * N2 / size
jIndexEnd   = (rank+1) * N2 / size

nu_c = electron_charge * (B_mag * 1.) / (2. * np.pi * m * c)

for j in range(jIndexStart, jIndexEnd):
    print "j = ", j, ", rank = ", rank
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

MJ_I_exact_gathered = np.zeros([N1, N2])
comm.barrier()
comm.Allgather([MJ_I_exact[jIndexStart:jIndexEnd, :], MPI.DOUBLE],
               [MJ_I_exact_gathered, MPI.DOUBLE]
              )
if (rank==0):
    np.savetxt("MJ_I_using_symphony_fits.txt", MJ_I_exact_gathered)
    pl.contourf(MJ_I_exact_gathered, 100)
    pl.colorbar()
    pl.show()

#MJ_I_avgs  = sp.j_nu_fit_py(nuratio, B_mag_avg, n_e_avg,
#                            obs_angle_avg, sp.MAXWELL_JUETTNER,
#                            sp.STOKES_I, theta_e, power_law_p,
#                            gamma_min, gamma_max, gamma_cutoff,
#                            kappa, kappa_width)
#
#print 'MJ_I', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (MJ_I_exact_avg 
#                                                  - MJ_I_avgs)/MJ_I_exact_avg


##-------------------------------MJ_Q-------------------------------------------#
#x = 0
#y = 0
#MJ_Q_exact_avg = 0
#MJ_Q_exact = [[0 for i in range(1152)] for j in range(1152)]
#
#for x in range(0, 1151):
#  for y in range(0, 1151):
#    MJ_Q_exact[x][y] = sp.j_nu_fit_py(nuratio, B_mag[x][y], n_e[x][y],
#                                      obs_angle[x][y], sp.MAXWELL_JUETTNER,
#                                      sp.STOKES_Q, theta_e, power_law_p,
#                                      gamma_min, gamma_max, gamma_cutoff,
#                                      kappa, kappa_width)
#    MJ_Q_exact_avg = MJ_Q_exact_avg + MJ_Q_exact[x][y]
#
#MJ_Q_exact_avg = MJ_Q_exact_avg/(1152*1152)
#
#MJ_Q_avgs  = sp.j_nu_fit_py(nuratio, B_mag_avg, n_e_avg,
#                            obs_angle_avg, sp.MAXWELL_JUETTNER,
#                            sp.STOKES_Q, theta_e, power_law_p,
#                            gamma_min, gamma_max, gamma_cutoff,
#                            kappa, kappa_width)
#
#print 'MJ_Q', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (MJ_Q_exact_avg
#                                                  - MJ_Q_avgs)/MJ_Q_exact_avg
#
#
##-------------------------------MJ_V-------------------------------------------#
#x = 0
#y = 0
#MJ_V_exact_avg = 0
#MJ_V_exact = [[0 for i in range(1152)] for j in range(1152)]
#
#for x in range(0, 1151):
#  for y in range(0, 1151):
#    MJ_V_exact[x][y] = sp.j_nu_fit_py(nuratio, B_mag[x][y], n_e[x][y],
#                                      obs_angle[x][y], sp.MAXWELL_JUETTNER,
#                                      sp.STOKES_V, theta_e, power_law_p,
#                                      gamma_min, gamma_max, gamma_cutoff,
#                                      kappa, kappa_width)
#    MJ_V_exact_avg = MJ_V_exact_avg + MJ_V_exact[x][y]
#
#MJ_V_exact_avg = MJ_V_exact_avg/(1152*1152)
#
#MJ_V_avgs  = sp.j_nu_fit_py(nuratio, B_mag_avg, n_e_avg,
#                            obs_angle_avg, sp.MAXWELL_JUETTNER,
#                            sp.STOKES_V, theta_e, power_law_p,
#                            gamma_min, gamma_max, gamma_cutoff,
#                            kappa, kappa_width)
#
#print 'MJ_V', '(<j_nu()>-j_nu(<>))/<j_nu()> = ', (MJ_V_exact_avg
#                                                  - MJ_V_avgs)/MJ_V_exact_avg
