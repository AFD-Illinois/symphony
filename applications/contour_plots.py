import sys

#symphony_build_path = '/home/mani/work/symphony/build'
symphony_build_path = '/home/alex/Documents/Spring_2016/1symphony/symphony/build'
sys.path.append(symphony_build_path)

import symphonyPy as sp
import numpy as np
import pylab as pl

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

nuratio_used    = []
obs_angle_used  = []

#----------------------set important parameters-------------------------------#

num_skip              = 64                      #sample every nth point
max_nuratio           = 1.e15                   #max nu/nu_c
number_of_points      = 128                     #size of grid
distribution_function = sp.KAPPA_DIST            #distribution function
polarization          = sp.STOKES_I             #Stokes parameter

#---------------------------import data from Dr. Kunz's simulation------------#
rank = 0
size = 1

datafiles_path = '/home/alex/Documents/Spring_2016/'
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

relative_difference = np.empty([number_of_points, number_of_points])
exact_avg_only      = np.empty([number_of_points, number_of_points])
avgs_only           = np.empty([number_of_points, number_of_points])

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
		exact = np.zeros([N2, N1])

		jIndexStart =  rank    * N2 / size
		jIndexEnd   = (rank+1) * N2 / size

		for j in range(jIndexStart, jIndexEnd):
		    for i in range(0, N1):
			exact[j][i] = sp.j_nu_fit_py(nuratio * nu_c[j][i],
	                                             B_mag[j][i],
	                                             n_e[j][i], 
	                                             obs_angle[j][i],
	                                             distribution_function, 
	                                             polarization, 
						     theta_e, 
						     power_law_p, 
	                                             gamma_min, 
						     gamma_max, 
						     gamma_cutoff, 
	                                             kappa, 
						     kappa_width
	                                            )


		exact_avg = np.mean(exact)

		avgs  = sp.j_nu_fit_py(nu_avg, B_mag_avg, n_e_avg,
		                       obs_angle_avg, distribution_function,
 		                       polarization, theta_e, power_law_p,
 		                       gamma_min, gamma_max, gamma_cutoff,
   	   	                       kappa, kappa_width)

		relative_difference[x][y] = ((exact_avg - avgs)/exact_avg)
		exact_avg_only[x][y]      = exact_avg
		avgs_only[x][y]           = avgs

		if(x == 0):
			obs_angle_used.append(obs_angle_avg * 180. / np.pi)
	nuratio_used.append(nuratio)

obs_angle_used.reverse() #obs_angle_used is reversed
nuratio_used.reverse()   #nuratio_used is also reversed

#------------------------make contour plot-------------------------------------#

X, Y = np.meshgrid(nuratio_used, obs_angle_used)
Z = relative_difference

pl.contourf(np.log10(X), Y, exact_avg_only, 200)
pl.title('<j_nu(exact)>')
pl.xlabel('$log_{10}(\\nu/\\nu_c)$')
pl.ylabel('$\\theta$ (deg)')
pl.colorbar()
pl.show()

pl.contourf(np.log10(X), Y, avgs_only, 200)
pl.title('j_nu(<avgs>)')
pl.xlabel('$log_{10}(\\nu/\\nu_c)$')
pl.ylabel('$\\theta$ (deg)')
pl.colorbar()
pl.show()

pl.contourf(np.log10(X), Y, Z, 200)
pl.xlabel('$log_{10}(\\nu/\\nu_c)$')
pl.ylabel('$\\theta$ (deg)')
pl.colorbar()
pl.show()
