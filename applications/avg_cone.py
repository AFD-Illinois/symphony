import sys
symphony_build_path = '../build'
sys.path.append(symphony_build_path)
import symphonyPy as sp
import numpy as np
import pylab as pl
import numpy.ma
from mpi4py import MPI


#----------------------set important parameters--------------------------------#

num_skip              = 128                     #sample every nth point
distribution_function = sp.MAXWELL_JUETTNER	#distribution function
EMISS                 = True                    #True = j_nu, False = alpha_nu
IN_PLANE              = False		        #True = obs_angle in plane
mask_tolerance        = 1.			#error > tolerance is white
number_of_points      = 100

#--------------------set constant parameters for the calculation--------------#
m = 9.1093826e-28
c = 2.99792458e10
electron_charge = 4.80320680e-10
theta_e = 10.             #TODO: Do we need to do anything about electron temp?
h = 6.6260693e-27
gamma_min = 1.
gamma_max = 1000.
gamma_cutoff = 1e10
power_law_p = 3.
kappa = 2.5
kappa_width = 10.
B_scale = 30.


nu = 230e9
#---------------------------import data from Dr. Kunz's simulation------------#
comm = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()
name = MPI.Get_processor_name()

N1, N2 = 1152/num_skip, 1152/num_skip #TODO: figure out way to avoid hardcoding 1152?

comm.Barrier()
if (rank == 0):
#   datafiles_path = ''
   datafiles_path = '/home/alex/Documents/Spring_2016/'
   B_x = np.loadtxt(datafiles_path + 'mirror_bx.out')[::num_skip, ::num_skip] * B_scale
   B_y = np.loadtxt(datafiles_path + 'mirror_by.out')[::num_skip, ::num_skip] * B_scale
   B_z = np.loadtxt(datafiles_path + 'mirror_bz.out')[::num_skip, ::num_skip] * B_scale


#   B_x = np.ones(np.shape(B_x))
#   B_y = np.zeros(np.shape(B_x))
#   B_z = np.zeros(np.shape(B_x))
#   n_e = np.ones(np.shape(B_x))

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

print "rank = ", rank, " of ", size, "procs"


B_x_avg        = np.mean(B_x)
B_y_avg        = np.mean(B_y)
B_z_avg        = np.mean(B_z)
n_e_avg        = np.mean(n_e)
B_mag_avg      = np.mean(B_mag)

B_vector     = [B_x, B_y, B_z]
B_avg_vector = [B_x_avg, B_y_avg, B_z_avg]
B_avg_vector = B_avg_vector/np.linalg.norm(B_avg_vector)
b_hat        = B_vector / np.linalg.norm(B_vector)


nu_c_avg = (electron_charge * B_mag_avg) / (2. * np.pi * m * c)

cone_angles      = np.zeros([number_of_points])
cone_angles_used = np.zeros([number_of_points])

#-------------------------------define rotation matrix------------------------#
#Euler-Rodrigues rotation matrix, as implemented by user "unutbu" from
# http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
#It matches the literature Euler-Rodrigues formula, and has been shown to
#produce the correct results. This version provides a clockwise rotation, 
#whereas its transpose (as on Wikipedia) produces a counterclockwise rotation.
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

#-----------------choose j_nu or alpha_nu--------------------------------------#
def j_nu_or_alpha_nu(nu, B, n_e, obs_angle, distribution_function,
	             polarization, theta_e, power_law_p, gamma_min,
		     gamma_max, gamma_cutoff, kappa, kappa_width):

#	if(obs_angle < 1e-5):
#                   return 0.

	if(EMISS == True):

		return sp.j_nu_fit_py(nu, B, n_e, obs_angle, distribution_function,
                     polarization, theta_e, power_law_p, gamma_min,
                     gamma_max, gamma_cutoff, kappa, kappa_width)

	else:
		return sp.alpha_nu_fit_py(nu, B, n_e, obs_angle, distribution_function,
                     polarization, theta_e, power_law_p, gamma_min,
                     gamma_max, gamma_cutoff, kappa, kappa_width)  

#-------------------------set up nu-theta space scan--------------------------#

#relative_difference_I = np.empty([number_of_points, number_of_points])
#exact_avg_only_I      = np.empty([number_of_points, number_of_points])
#avgs_only_I           = np.empty([number_of_points, number_of_points])
#relative_difference_Q = np.empty([number_of_points, number_of_points])
#exact_avg_only_Q      = np.empty([number_of_points, number_of_points])
#avgs_only_Q           = np.empty([number_of_points, number_of_points])
#relative_difference_U = np.empty([number_of_points, number_of_points])
#exact_avg_only_U      = np.empty([number_of_points, number_of_points])
#avgs_only_U           = np.empty([number_of_points, number_of_points])
#relative_difference_V = np.empty([number_of_points, number_of_points])
#exact_avg_only_V      = np.empty([number_of_points, number_of_points])
#avgs_only_V           = np.empty([number_of_points, number_of_points])


#nu_used             = np.empty([number_of_points])
#obs_angle_used      = np.empty([number_of_points])

#for x in range(0, number_of_points):
#	nu = nu_c_avg * 10.**(1.*x / (number_of_points - 1) * np.log10(max_nuratio))
#        if (rank == 0):
#	   print 100.0*x/(number_of_points - 1), '% complete'
#	for y in range(0, number_of_points):

if(IN_PLANE == True):
	#rotates observer vector in plane
	axis_of_rot = [0, B_avg_vector[2], -B_avg_vector[1]]
else:
	#rotates observer vector out of plane
	axis_of_rot = [-B_avg_vector[1], B_avg_vector[0], 0]

#desired angle between observer vector and mean B
theta = np.pi/4.
#theta       = (1.0*y/number_of_points * np.pi/2.)

#		print theta * 180./np.pi

#obs_vector  = np.dot(rotation_matrix(axis_of_rot, theta), 
#		     B_avg_vector)

obs_vector = np.array([-0.5, 0.5634654, 1])

obs_angle = np.arccos((B_x * obs_vector[0] 
                       + B_y * obs_vector[1]
                       + B_z * obs_vector[2])
                       / (np.linalg.norm(obs_vector)*B_mag))


#		obs_angle_avg  = np.mean(obs_angle) #not physically correct
angle_to_mean_field = np.arccos(np.dot(B_avg_vector, obs_vector)) 
obs_angle_avg = angle_to_mean_field


#-----------------------avg over cone angle------------------------------------#
for x in range(number_of_points):


	cone_rot = (1.0 * x / number_of_points * 2.*np.pi)

	cone_angles_used[x] = cone_rot

	k_parallel = (np.dot(obs_vector, B_avg_vector) 
	              / np.dot(B_avg_vector, B_avg_vector)) * B_avg_vector

	k_perp = obs_vector - k_parallel

	k_perp_rotated = np.dot(rotation_matrix(k_parallel, cone_rot), 
        	                k_perp)

	obs_vector = k_perp_rotated + k_parallel

#	print obs_vector

#	print k_perp_rotated, np.arccos(np.dot(k_perp_rotated, k_perp)/(np.dot(k_perp, k_perp)*np.dot(k_perp_rotated, k_perp_rotated)))
#	print np.dot((k_parallel + k_perp), B_avg_vector)#, np.dot(obs_vector, B_avg_vector)
#	continue

#-------------------------BEGIN MODIFICATIONS----------------------------------#

	#create NxN array of symphony frame e_alpha axis #TODO: check this
	obs_vector_array = [np.full(np.shape(B_x), obs_vector[0]),
			    np.full(np.shape(B_x), obs_vector[1]),
		            np.full(np.shape(B_x), obs_vector[2])]		
	b_hat_dot_k_hat = b_hat[0] * obs_vector_array[0] + b_hat[1] * obs_vector_array[1] + b_hat[2] * obs_vector_array[2]
	
	e_alpha_unnormed = (b_hat - obs_vector_array * b_hat_dot_k_hat)
	e_alpha_norm = np.sqrt(e_alpha_unnormed[0]**2. + e_alpha_unnormed[1]**2. + e_alpha_unnormed[2]**2.)
	e_alpha_normed = [e_alpha_unnormed[0]/e_alpha_norm, e_alpha_unnormed[1]/e_alpha_norm, e_alpha_unnormed[2]/e_alpha_norm]
	
	#create e_alpha_prime #TODO: check this
	e_alpha_p_not_perp = B_avg_vector
	e_alpha_p_perp = (e_alpha_p_not_perp - obs_vector * np.dot(obs_vector, e_alpha_p_not_perp))
	e_alpha_p_normed = e_alpha_p_perp / np.linalg.norm(e_alpha_p_perp)
	
	

	xi = np.arccos((e_alpha_normed[0] * e_alpha_p_normed[0] +
                       e_alpha_normed[1] * e_alpha_p_normed[1] +
                       e_alpha_normed[2] * e_alpha_p_normed[2]))
	
	xi = np.nan_to_num(xi)


#---------------------------END MODIFICATIONS---------------------------------#


#----------------------scan over simulation------------------------------------#
	exact_avg = 0
	exact_I = np.zeros([N2, N1])
	exact_Q = np.zeros([N2, N1])
	exact_U = np.zeros([N2, N1])
	exact_V = np.zeros([N2, N1])

	jIndexStart =  rank    * N2 / size
	jIndexEnd   = (rank+1) * N2 / size

	for j in range(jIndexStart, jIndexEnd):
	    for i in range(0, N1):
	
		exact_I[j][i] = j_nu_or_alpha_nu(nu, 
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
	
		exact_Q[j][i] = j_nu_or_alpha_nu(nu,
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
	                                        ) * np.cos(2.*xi[j][i]) # * np.cos(alpha_obs[j][i]) PUT IN XI
	
		exact_U[j][i] = j_nu_or_alpha_nu(nu,
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
	                                        ) * -np.sin(2.*xi[j][i]) # * np.cos(alpha_obs[j][i]) PUT IN XI
	
		exact_V[j][i] = j_nu_or_alpha_nu(nu,
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
	
	exact_I_gathered = np.zeros([N2, N1])
	exact_Q_gathered = np.zeros([N2, N1])
	exact_U_gathered = np.zeros([N2, N1])
	exact_V_gathered = np.zeros([N2, N1])
	
	comm.barrier()
	comm.Allgather([exact_I[jIndexStart:jIndexEnd, :], MPI.DOUBLE],
	               [exact_I_gathered, MPI.DOUBLE])
	comm.Allgather([exact_Q[jIndexStart:jIndexEnd, :], MPI.DOUBLE],
	               [exact_Q_gathered, MPI.DOUBLE])
	comm.Allgather([exact_U[jIndexStart:jIndexEnd, :], MPI.DOUBLE],
	               [exact_U_gathered, MPI.DOUBLE])
	comm.Allgather([exact_V[jIndexStart:jIndexEnd, :], MPI.DOUBLE],
	               [exact_V_gathered, MPI.DOUBLE])
	
	exact_I = exact_I_gathered
	exact_Q = exact_Q_gathered
	exact_U = exact_U_gathered
	exact_V = exact_V_gathered
	
	exact_avg_I = np.mean(exact_I)
	exact_avg_Q = np.mean(exact_Q)
	exact_avg_U = np.mean(exact_U)
	exact_avg_V = np.mean(exact_V)
	
			#angles between the assumed B field and the mean B field
	#	  	beta_obs_avg  = np.arccos(1. * B_y_avg / B_mag_avg) + observer_rotation
	#                alpha_obs_avg = np.arccos(1. * B_x_avg / B_mag_avg)
	
	
#	print 'I:', exact_avg_I
#	print 'Q:', exact_avg_Q
#	print 'U:', exact_avg_U
#	print 'V:', exact_avg_V
	cone_angles[x] = exact_avg_U

	print x

pl.plot(cone_angles_used * 180. / np.pi, cone_angles)
mean_cone = np.mean(cone_angles)
print mean_cone
pl.axhline(mean_cone, c='r', lw=2)
pl.xlim([0, 360])
pl.show()
#	figure, ax = pl.subplots(1, 4, figsize=(10, 10))
#	figure, ax = pl.subplots(1, 2)
#	plot1 = ax[0].contourf(exact_I, 100)
#	plot2 = ax[0].contourf(exact_Q, 100)
#	plot3 = ax[1].contourf(exact_U, 100)
#	plot4 = ax[3].contourf(exact_V, 100)
#	figure.colorbar(plot1, ax=ax[0])
#	figure.colorbar(plot2, ax=ax[0])
#	figure.colorbar(plot3, ax=ax[1])
#	figure.colorbar(plot4, ax=ax[3])
#	pl.tight_layout()
#	pl.show()
#	print x
#	figure.savefig(str(x)+'.png')
#	pl.close(figure)

#--------------MORE MODIFICATIONS---------------------------------------------#
#if(np.abs(np.dot(B_avg_vector, obs_vector) - 1) < 1e-7):
#	xi_avg = np.pi/4. #choose it arbitrarily; emission/absorption will be zero anyway		 

#else:
#	e_alpha_avg = (B_avg_vector - obs_vector * np.dot(B_avg_vector, obs_vector))
#	e_alpha_avg = e_alpha_avg / np.linalg.norm(e_alpha_avg)
#	print np.dot(e_alpha_avg, e_alpha_p_normed)
#  	          print 'HERE IT IS', np.dot(e_alpha_avg, e_alpha_avg), np.dot(e_alpha_avg, obs_vector)
#	if(np.dot(e_alpha_avg, e_alpha_p_normed) - 1. > 0 and np.dot(e_alpha_avg, e_alpha_p_normed) - 1. < 1e-6):
#                print 'GOT HERE', np.dot(e_alpha_avg, e_alpha_p_normed)
#		xi_avg = np.arccos(1.)
#	else:
#		xi_avg = np.arccos(np.dot(e_alpha_avg, e_alpha_p_normed))

#---------------END MODIFICATIONS---------------------------------------------#

#	avgs_I      = j_nu_or_alpha_nu(nu, B_mag_avg, n_e_avg, 
#	                       obs_angle_avg, distribution_function,
#	                       sp.STOKES_I, theta_e, power_law_p,
#	                       gamma_min, gamma_max, gamma_cutoff,
#	                       kappa, kappa_width)
#	
#	avgs_Q      = j_nu_or_alpha_nu(nu, B_mag_avg, n_e_avg, 
#	                       obs_angle_avg, distribution_function,
#	                       sp.STOKES_Q, theta_e, power_law_p,
#	                       gamma_min, gamma_max, gamma_cutoff,
#	                       kappa, kappa_width) * np.cos(2.*xi_avg) #* np.cos(alpha_obs_avg)  
#	
#	avgs_U      = j_nu_or_alpha_nu(nu, B_mag_avg, n_e_avg, 
#	                       obs_angle_avg, distribution_function,
#	                       sp.STOKES_Q, theta_e, power_law_p,
#	                       gamma_min, gamma_max, gamma_cutoff,
#	                       kappa, kappa_width) * -np.sin(2.*xi_avg) #* np.cos(alpha_obs_avg) 
#	
#	avgs_V      = j_nu_or_alpha_nu(nu, B_mag_avg, n_e_avg, 
#	                       obs_angle_avg, distribution_function,
#	                       sp.STOKES_V, theta_e, power_law_p,
#	                       gamma_min, gamma_max, gamma_cutoff,
#	                       kappa, kappa_width)

#	relative_difference_I[y][x] = np.fabs((exact_avg_I - avgs_I)/exact_avg_I)
#	exact_avg_only_I[y][x]      = exact_avg_I
#	avgs_only_I[y][x]           = avgs_I

#		relative_difference_Q[y][x] = np.fabs((exact_avg_Q - avgs_Q)/exact_avg_Q)
#                exact_avg_only_Q[y][x]      = exact_avg_Q
#                avgs_only_Q[y][x]           = avgs_Q
#
#		relative_difference_U[y][x] = np.fabs((exact_avg_U - avgs_U)/exact_avg_U)
#                exact_avg_only_U[y][x]      = exact_avg_U
#                avgs_only_U[y][x]           = avgs_U
#
#		relative_difference_V[y][x] = np.fabs((exact_avg_V - avgs_V)/exact_avg_V)
#                exact_avg_only_V[y][x]      = exact_avg_V
#                avgs_only_V[y][x]           = avgs_V


#		if(x == 0):
#			obs_angle_used[y] = obs_angle_avg * 180. / np.pi
#	nu_used[x] = nu / nu_c_avg

#------------------------make contour plot-------------------------------------#

#if (rank == 0):
#   # Set plot parameters to make beautiful plots
#   pl.rcParams['figure.figsize']  = 12, 7.5
#   pl.rcParams['lines.linewidth'] = 1.5
#   pl.rcParams['font.family']     = 'serif'
##   pl.rcParams['font.weight']     = 'bold'
#   pl.rcParams['font.size']       = 15
#   pl.rcParams['font.sans-serif'] = 'serif'
#   pl.rcParams['text.usetex']     = False
#   pl.rcParams['axes.linewidth']  = 1.5
#   pl.rcParams['axes.titlesize']  = 'medium'
#   pl.rcParams['axes.labelsize']  = 'medium'
#   
#   pl.rcParams['xtick.major.size'] = 8
#   pl.rcParams['xtick.minor.size'] = 4
#   pl.rcParams['xtick.major.pad']  = 8
#   pl.rcParams['xtick.minor.pad']  = 8
#   pl.rcParams['xtick.color']      = 'k'
#   pl.rcParams['xtick.labelsize']  = 'medium'
#   pl.rcParams['xtick.direction']  = 'in'
#   
#   pl.rcParams['ytick.major.size'] = 8
#   pl.rcParams['ytick.minor.size'] = 4
#   pl.rcParams['ytick.major.pad']  = 8
#   pl.rcParams['ytick.minor.pad']  = 8
#   pl.rcParams['ytick.color']      = 'k'
#   pl.rcParams['ytick.labelsize']  = 'medium'
#   pl.rcParams['ytick.direction']  = 'in'
#   
#   X, Y = np.meshgrid(nu_used, obs_angle_used)
#
#   figure, ax = pl.subplots(4, 3, figsize=(10, 10))
#   figure.suptitle(figure_title)
#   
#   plot1 = ax[0,0].contourf(np.log10(X), Y, exact_avg_only_I, 200)
#   figure.colorbar(plot1, ax=ax[0,0])
#   ax[0,0].set_title('$<j_\\nu(n, \mathbf{B})>$')
#   
#   
#   plot2 = ax[0,1].contourf(np.log10(X), Y, avgs_only_I, 200)
#   figure.colorbar(plot2, ax=ax[0,1])
#   ax[0,1].set_title('$j_\\nu(<n>, <\mathbf{B}>)$')
#   
#   if(EMISS == False):
#           ax[0,0].set_title('$<\\alpha_\\nu(n, \mathbf{B})>$')
#           ax[0,1].set_title('$\\alpha_\\nu(<n>, <\mathbf{B}>)$')
#   
#   relative_difference_I = np.ma.array(relative_difference_I, mask=relative_difference_I > mask_tolerance)
#   plot3 = ax[0,2].contourf(np.log10(X), Y, relative_difference_I, 200)
#   figure.colorbar(plot3, ax=ax[0,2])
#   ax[0,2].set_title('$|\mathrm{Error}|$')
#   
#   plot4 = ax[1,0].contourf(np.log10(X), Y, exact_avg_only_Q, 200)
#   figure.colorbar(plot4, ax=ax[1,0])
#   
#   plot5 = ax[1,1].contourf(np.log10(X), Y, avgs_only_Q, 200)
#   figure.colorbar(plot5, ax=ax[1,1])
#   
#   
#   relative_difference_Q = np.ma.array(relative_difference_Q, mask=relative_difference_Q > mask_tolerance)
#   plot6 = ax[1,2].contourf(np.log10(X), Y, relative_difference_Q, 200)
#   figure.colorbar(plot6, ax=ax[1,2])
#   
#   plot7 = ax[2,0].contourf(np.log10(X), Y, exact_avg_only_U, 200)
#   figure.colorbar(plot7, ax=ax[2,0])
#   
#   plot8 = ax[2,1].contourf(np.log10(X), Y, avgs_only_U, 200)
#   figure.colorbar(plot8, ax=ax[2,1])
#   
#   relative_difference_U = np.ma.array(relative_difference_U, mask=relative_difference_U > mask_tolerance)
#   plot9 = ax[2,2].contourf(np.log10(X), Y, relative_difference_U, 200)
#   figure.colorbar(plot9, ax=ax[2,2])
#   
#   plot10 = ax[3,0].contourf(np.log10(X), Y, exact_avg_only_V, 200)
#   figure.colorbar(plot10, ax=ax[3,0])
#   
#   plot11 = ax[3,1].contourf(np.log10(X), Y, avgs_only_V, 200)
#   figure.colorbar(plot11, ax=ax[3,1])
#   
#   relative_difference_V = np.ma.array(relative_difference_V, mask=relative_difference_V > mask_tolerance)
#   plot12 = ax[3,2].contourf(np.log10(X), Y, relative_difference_V, 200)
#   figure.colorbar(plot12, ax=ax[3,2])
#   
#   figure.add_subplot(111, frameon=False)
#   pl.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
#   pl.xlabel('$\\log_{10}(\\nu/\\overline{\\nu}_c)$', fontsize='large')
#   pl.ylabel('$\\theta$ (deg)', fontsize='large')
#   pl.tight_layout()
#   pl.show()
#   
#   #print error at position where j_nu or alpha_nu is maximal
#   max_value_I = np.amax(exact_avg_only_I) #TODO: note that I removed abs from this
#   location_I = np.where(exact_avg_only_I == max_value_I) #TODO: removed abs from here too
#   x_I, y_I = location_I[0][0], location_I[1][0]
#   
#   max_value_Q = np.amax(np.abs(exact_avg_only_Q))
#   location_Q = np.where(np.abs(exact_avg_only_Q) == max_value_Q)
#   x_Q, y_Q = location_Q[0][0], location_Q[1][0]
#   
#   max_value_U = np.amax(np.abs(exact_avg_only_U))
#   location_U = np.where(np.abs(exact_avg_only_U) == max_value_U)
#   x_U, y_U = location_U[0][0], location_U[1][0]
#   
#   max_value_V = np.amax(np.abs(exact_avg_only_V))
#   location_V = np.where(np.abs(exact_avg_only_V) == max_value_V)
#   x_V, y_V = location_V[0][0], location_V[1][0]
#
#   print '\nError at location of max j_nu or alpha_nu: '
#   print 'STOKES I: ', relative_difference_I[x_I][y_I]
#   print 'STOKES Q: ', relative_difference_Q[x_Q][y_Q]
#   print 'STOKES U: ', relative_difference_U[x_U][y_U]
#   print 'STOKES V: ', relative_difference_V[x_V][y_V]
