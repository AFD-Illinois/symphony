import sys

symphony_build_path = '/home/mani/work/symphony/build'
#symphony_build_path = '/home/alex/Documents/Spring_2016/1symphony/symphony/build'
sys.path.append(symphony_build_path)

import symphonyPy as sp
import numpy as np
import pylab as pl

# Set plot parameters to make beautiful plots
pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 15
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex']     = True
pl.rcParams['axes.linewidth']  = 1.5
pl.rcParams['axes.titlesize']  = 'medium'
pl.rcParams['axes.labelsize']  = 'medium'

pl.rcParams['xtick.major.size'] = 8     
pl.rcParams['xtick.minor.size'] = 4     
pl.rcParams['xtick.major.pad']  = 8     
pl.rcParams['xtick.minor.pad']  = 8     
pl.rcParams['xtick.color']      = 'k'     
pl.rcParams['xtick.labelsize']  = 'medium'
pl.rcParams['xtick.direction']  = 'in'    

pl.rcParams['ytick.major.size'] = 8     
pl.rcParams['ytick.minor.size'] = 4     
pl.rcParams['ytick.major.pad']  = 8     
pl.rcParams['ytick.minor.pad']  = 8     
pl.rcParams['ytick.color']      = 'k'     
pl.rcParams['ytick.labelsize']  = 'medium'
pl.rcParams['ytick.direction']  = 'in'    

#--------------------set constant parameters for the calculation--------------#
m = 9.1093826e-28
c = 2.99792458e10
electron_charge = 4.80320680e-10
theta_e = 10. #TODO: Do we need to do anything about electron temp?
e = 4.80320680e-10
h = 6.6260693e-27
gamma_min = 1.
gamma_max = 1000.
gamma_cutoff = 1e10
power_law_p = 3.5
kappa = 3.5
kappa_width = 10.
B_scale = 30. 

#----------------------set important parameters-------------------------------#

num_skip              = 64                      #sample every nth point
max_nuratio           = 1.e8                   #max nu/nu_c
number_of_points      = 64                      #size of grid
distribution_function = sp.KAPPA_DIST     #distribution function

#---------------------------import data from Dr. Kunz's simulation------------#
rank = 0
size = 1

#datafiles_path = '/home/alex/Documents/Spring_2016/'
datafiles_path = '/home/mani/work/kunz_data/'
B_x = np.loadtxt(datafiles_path + 'mirror_bx.out')[::num_skip, ::num_skip] * B_scale
B_y = np.loadtxt(datafiles_path + 'mirror_by.out')[::num_skip, ::num_skip] * B_scale
B_z = np.loadtxt(datafiles_path + 'mirror_bz.out')[::num_skip, ::num_skip] * B_scale
B_mag = np.sqrt(B_x**2. + B_y**2. + B_z**2.)
n_e = np.loadtxt(datafiles_path + 'mirror_d.out')[::num_skip, ::num_skip]

N1 = B_x.shape[0]
N2 = B_x.shape[1]


nuratio = 1.e0

B_x_avg        = np.mean(B_x)
B_y_avg        = np.mean(B_y)
B_z_avg        = np.mean(B_z)
n_e_avg        = np.mean(n_e)
B_mag_avg      = np.mean(B_mag)

B_vector     = [B_x, B_y, B_z]
B_avg_vector = [B_x_avg, B_y_avg, B_z_avg]
B_avg_vector = B_avg_vector/np.linalg.norm(B_avg_vector)

#-------------------------------define rotation matrix------------------------#
def rotation_matrix(axis, theta):
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


#-------------------------set up nu-theta space scan--------------------------#

relative_difference_I = np.empty([number_of_points, number_of_points])
exact_avg_only_I      = np.empty([number_of_points, number_of_points])
avgs_only_I           = np.empty([number_of_points, number_of_points])
relative_difference_Q = np.empty([number_of_points, number_of_points])
exact_avg_only_Q      = np.empty([number_of_points, number_of_points])
avgs_only_Q           = np.empty([number_of_points, number_of_points])
relative_difference_V = np.empty([number_of_points, number_of_points])
exact_avg_only_V      = np.empty([number_of_points, number_of_points])
avgs_only_V           = np.empty([number_of_points, number_of_points])


nuratio_used        = np.empty([number_of_points])
obs_angle_used      = np.empty([number_of_points])

for x in range(0, number_of_points):
	nuratio = 1. * 10.**(np.log10(max_nuratio) * x/(number_of_points-1.))
        print 100.0*x/number_of_points, '% complete'
	for y in range(0, number_of_points):
		axis_of_rot = [-B_avg_vector[1], B_avg_vector[0], 0]
		theta       = (1.0*y/number_of_points * np.pi/2.)

		obs_vector  = np.dot(rotation_matrix(axis_of_rot, theta), 
 				     B_avg_vector)

		obs_angle = np.arccos((B_x * obs_vector[0] 
                                       + B_y * obs_vector[1]
                                       + B_z * obs_vector[2])
                                      / (np.linalg.norm(obs_vector)*B_mag))

		#generate remaining averages
		obs_angle_avg  = np.mean(obs_angle)
		nu_c           = electron_charge * (B_mag * B_scale) / (2. * np.pi * m * c)
		nu_c_avg       = np.mean(nu_c)
		nu_avg         = nuratio * nu_c_avg


#----------------------scan over simulation------------------------------------#
		exact_avg = 0
		exact_I = np.zeros([N2, N1])
		exact_Q = np.zeros([N2, N1])
		exact_V = np.zeros([N2, N1])

		jIndexStart =  rank    * N2 / size
		jIndexEnd   = (rank+1) * N2 / size

		for j in range(jIndexStart, jIndexEnd):
		    for i in range(0, N1):
			exact_I[j][i] = sp.j_nu_fit_py(nuratio * nu_c[j][i],
	                                             B_mag[j][i],
	                                             n_e[j][i], 
	                                             obs_angle[j][i],
	                                             distribution_function, 
	                                             sp.STOKES_I, 
						     theta_e, 
						     power_law_p, 
	                                             gamma_min, 
						     gamma_max, 
						     gamma_cutoff, 
	                                             kappa, 
						     kappa_width
	                                            )
			exact_Q[j][i] = sp.j_nu_fit_py(nuratio * nu_c[j][i],
                                                     B_mag[j][i],
                                                     n_e[j][i],
                                                     obs_angle[j][i],
                                                     distribution_function,
                                                     sp.STOKES_Q,
                                                     theta_e,
                                                     power_law_p,
                                                     gamma_min,
                                                     gamma_max,
                                                     gamma_cutoff,
                                                     kappa,
                                                     kappa_width
                                                    )
			exact_V[j][i] = sp.j_nu_fit_py(nuratio * nu_c[j][i],
                                                     B_mag[j][i],
                                                     n_e[j][i],
                                                     obs_angle[j][i],
                                                     distribution_function,
                                                     sp.STOKES_V,
                                                     theta_e,
                                                     power_law_p,
                                                     gamma_min,
                                                     gamma_max,
                                                     gamma_cutoff,
                                                     kappa,
                                                     kappa_width
                                                    )


		exact_avg_I = np.mean(exact_I)
		exact_avg_Q = np.mean(exact_Q)
		exact_avg_V = np.mean(exact_V)

		avgs_I      = sp.j_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
		                       obs_angle_avg, distribution_function,
 		                       sp.STOKES_I, theta_e, power_law_p,
 		                       gamma_min, gamma_max, gamma_cutoff,
   	   	                       kappa, kappa_width)
		avgs_Q      = sp.j_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                                       obs_angle_avg, distribution_function,
                                       sp.STOKES_Q, theta_e, power_law_p,
                                       gamma_min, gamma_max, gamma_cutoff,
                                       kappa, kappa_width)
		avgs_V      = sp.j_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
                                       obs_angle_avg, distribution_function,
                                       sp.STOKES_V, theta_e, power_law_p,
                                       gamma_min, gamma_max, gamma_cutoff,
                                       kappa, kappa_width)

		relative_difference_I[y][x] = np.fabs((exact_avg_I - avgs_I)/exact_avg_I)
		exact_avg_only_I[y][x]      = exact_avg_I
		avgs_only_I[y][x]           = avgs_I

		relative_difference_Q[y][x] = np.fabs((exact_avg_Q - avgs_Q)/exact_avg_Q)
                exact_avg_only_Q[y][x]      = exact_avg_Q
                avgs_only_Q[y][x]           = avgs_Q

		relative_difference_V[y][x] = np.fabs((exact_avg_V - avgs_V)/exact_avg_V)
                exact_avg_only_V[y][x]      = exact_avg_V
                avgs_only_V[y][x]           = avgs_V

		if(x == 0):
			obs_angle_used[y] = obs_angle_avg * 180. / np.pi
	nuratio_used[x] = nuratio


#------------------------make contour plot-------------------------------------#
X, Y = np.meshgrid(nuratio_used, obs_angle_used)

figure, ax = pl.subplots(3, 3, figsize=(15, 15))
figure.suptitle("Kappa distribution")
plot1 = ax[0,0].contourf(np.log10(X), Y, exact_avg_only_I, 200)
figure.colorbar(plot1, ax=ax[0,0])
ax[0,0].set_title('$<j_\\nu(n, \mathbf{B})>$')

plot2 = ax[0,1].contourf(np.log10(X), Y, avgs_only_I, 200)
figure.colorbar(plot2, ax=ax[0,1])
ax[0,1].set_title('$j_\\nu(<n>, <\mathbf{B}>)$')

plot3 = ax[0,2].contourf(np.log10(X), Y, relative_difference_I, 200)
figure.colorbar(plot3, ax=ax[0,2])
ax[0,2].set_title('$|\mathrm{Error}|$')

plot4 = ax[1,0].contourf(np.log10(X), Y, exact_avg_only_Q, 200)
figure.colorbar(plot4, ax=ax[1,0])

plot5 = ax[1,1].contourf(np.log10(X), Y, avgs_only_Q, 200)
figure.colorbar(plot5, ax=ax[1,1])

plot6 = ax[1,2].contourf(np.log10(X), Y, relative_difference_Q, 200)
figure.colorbar(plot6, ax=ax[1,2])

plot7 = ax[2,0].contourf(np.log10(X), Y, exact_avg_only_V, 200)
figure.colorbar(plot7, ax=ax[2,0])

plot8 = ax[2,1].contourf(np.log10(X), Y, avgs_only_V, 200)
figure.colorbar(plot8, ax=ax[2,1])

plot9 = ax[2,2].contourf(np.log10(X), Y, relative_difference_V, 200)
figure.colorbar(plot9, ax=ax[2,2])


figure.add_subplot(111, frameon=False)
pl.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
pl.xlabel('$log_{10}(\\nu/\\nu_c)$')
pl.ylabel('$\\theta$ (deg)')
pl.tight_layout()

pl.show()
