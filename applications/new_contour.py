import sys
symphony_build_path = '../build'
sys.path.append(symphony_build_path)
import symphonyPy as sp
import numpy as np

#PREVENT SSH LOGOFF FROM STOPPING PLOT
#import matplotlib
#matplotlib.use('Agg')

import pylab as pl
import numpy.ma
from mpi4py import MPI


#----------------------set important parameters--------------------------------#

num_skip              = 128                      #sample every nth point
distribution_function = sp.POWER_LAW	#distribution function
EMISS                 = True                    #True = j_nu, False = alpha_nu
IN_PLANE              = True		        #True = obs_angle in plane
mask_tolerance        = 1.			#error > tolerance is white
number_of_points      = 10
cone_resolution       = 5
max_nuratio           = 1e6
figure_title          = ''

#--------------------set constant parameters for the calculation--------------#

m = 9.1093826e-28
c = 2.99792458e10
electron_charge = 4.80320680e-10
theta_e = 10.
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

#generate NxN array where each element is the local B vector
B_vector     = [B_x, B_y, B_z]

#makes <mean B_x, mean B_y, mean B_z>, then normalizes it
B_avg_vector = [B_x_avg, B_y_avg, B_z_avg]
B_avg_vector = B_avg_vector/np.linalg.norm(B_avg_vector)

#normalzed version of B_vector
b_hat        = B_vector / np.sqrt(B_x**2. + B_y**2. + B_z**2.)

#cyclotron frequency nu_c, calculated with average |B|
nu_c_avg = (electron_charge * B_mag_avg) / (2. * np.pi * m * c)

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

	#symphony's integrator breaks down for very low theta; TODO: fix this
	if(obs_angle < 1e-3):
		return 0.
	#symphony's integrator breaks down at theta = pi/2; TODO: fix
	if(np.abs(obs_angle - np.pi/2.) < 1e-3):
		obs_angle = 89.9 * np.pi / 180.
	#for some reason, it also breaks down at exactly theta = pi/3; TODO: fix
	if(np.abs(obs_angle - np.pi/6.) < 1e-3):
		obs_angle = 29.9 * np.pi/180.


	if(EMISS == True):

		#dimensional prefactor from j_nu, to be divided out
		dim_prefactor = (n_e_avg * electron_charge**2. * 
				 electron_charge * B_mag_avg 
				 / (2. * np.pi * m * c) / c)

		return sp.j_nu_fit_py(nu, B, n_e, obs_angle, distribution_function,
                     polarization, theta_e, power_law_p, gamma_min,
                     gamma_max, gamma_cutoff, kappa, kappa_width) / dim_prefactor


	else:
		#dimensional prefactor from alpha_nu, to be divided out
		dim_prefactor = n_e_avg * electron_charge**2. / (nu * 
				m * c)

		return sp.alpha_nu_fit_py(nu, B, n_e, obs_angle, distribution_function,
                     polarization, theta_e, power_law_p, gamma_min,
                     gamma_max, gamma_cutoff, kappa, kappa_width) / dim_prefactor 

#-------------------------set up nu-theta space scan--------------------------#

relative_difference_I = np.zeros([number_of_points, number_of_points])
exact_avg_only_I      = np.zeros([number_of_points, number_of_points])
avgs_only_I           = np.zeros([number_of_points, number_of_points])

relative_difference_lin = np.zeros([number_of_points, number_of_points])
exact_avg_only_lin      = np.zeros([number_of_points, number_of_points])
avgs_only_lin           = np.zeros([number_of_points, number_of_points])

relative_difference_circ = np.zeros([number_of_points, number_of_points])
exact_avg_only_circ      = np.zeros([number_of_points, number_of_points])
avgs_only_circ           = np.zeros([number_of_points, number_of_points])

nu_used             = np.empty([number_of_points])
obs_angle_used      = np.empty([number_of_points])

avgs_I = 0.
avgs_Q = 0.
avgs_U = 0.
avgs_V = 0.

avgs_lin  = 0.
avgs_circ = 0.


for x in range(0, number_of_points):
	nu = nu_c_avg * 10.**(1.*x / (number_of_points - 1) * np.log10(max_nuratio))
	if (rank == 0):
		print 100.0*x/(number_of_points - 1), '% complete'
	for y in range(0, number_of_points):

		if(IN_PLANE == True):
			#rotates observer vector in plane
			axis_of_rot = [0, B_avg_vector[2], -B_avg_vector[1]]
			axis_of_rot = axis_of_rot / np.linalg.norm(axis_of_rot)
		else:
			#rotates observer vector out of plane
			axis_of_rot = [-B_avg_vector[1], B_avg_vector[0], 0]
			axis_of_rot = axis_of_rot / np.linalg.norm(axis_of_rot)
		
		#desired angle between observer vector and mean B
		theta       = (1.0*y/number_of_points * np.pi/2.)

		#vector from sim. point to observer; wavevector		
		obs_vector  = np.dot(rotation_matrix(axis_of_rot, theta), 
				     B_avg_vector)
		
		#The quantity physically available to observers will be the angle
		#between k and the mean B field. Earlier versions naively generated
		#an array of the observer angles to each point in the simulation,
		#and then averaged that; the lines below are more correct.
		angle_to_mean_field = np.arccos(np.dot(B_avg_vector, obs_vector)) 
		obs_angle_avg = angle_to_mean_field
		
		#-----------------------avg over cone angle------------------------------------#
		#Since the simulation represents a small piece within a larger turbulent flow
		#(where only the mean field of the piece is known) there are a set of equivalent
		#wavevectors connecting the observer and simulation, each with the same observer
		#angle to the mean B field.  We average over these equivalent wavevectors.
		for cone_point in range(cone_resolution):
#			if(rank == 0):
#				print 'cone point: ', cone_point			

			#cone rotation angle, varies from 0 to 2*pi	
			cone_rot = (1.0 * cone_point / cone_resolution * 2.*np.pi)
	
			#component of wavevector parallel to mean B field
			k_parallel = (np.dot(obs_vector, B_avg_vector) 
			              / np.dot(B_avg_vector, B_avg_vector)) * B_avg_vector
			
			#component of wavevector perpendicular to mean B field
			k_perp = obs_vector - k_parallel
			
			#produce a rotated version of k_perp about mean B
			k_perp_rotated = np.dot(rotation_matrix(k_parallel, cone_rot), 
		        	                k_perp)
			
			#the new wavevector is now rotated about mean B, but has 
			#the same observer angle to mean B as the previous one
			k_vector = k_perp_rotated + k_parallel

			#generate array of observer angles for each point in the 
			#simulation
			obs_angle = np.arccos((B_x * k_vector[0]
                                       + B_y * k_vector[1]
                                       + B_z * k_vector[2])
                                       / (np.linalg.norm(k_vector)*B_mag))

			#Now we need to perform a rotation to align symphony's assumed
			#frame (local B always along x-hat) to the observer's fixed 
			#frame (primed coordinates)
			
			#create NxN array of symphony frame e_alpha axis
			obs_vector_array = [np.full(np.shape(B_x), k_vector[0]),
			                   np.full(np.shape(B_x), k_vector[1]),
			                   np.full(np.shape(B_x), k_vector[2])]  
		
			#we want the component of the local B field that is
			#in the plane of the sky; we first need b_hat dot k_hat	
			b_hat_dot_k_hat = b_hat[0] * obs_vector_array[0] + b_hat[1] * obs_vector_array[1] + b_hat[2] * obs_vector_array[2]
			
			#now we produce the coordinate e_alpha in the plane of the sky
			#and normalize it
			e_alpha_unnormed = (b_hat - obs_vector_array * b_hat_dot_k_hat)
			e_alpha_norm = np.sqrt(e_alpha_unnormed[0]**2. + e_alpha_unnormed[1]**2. + e_alpha_unnormed[2]**2.)
			e_alpha_normed = [e_alpha_unnormed[0]/e_alpha_norm, e_alpha_unnormed[1]/e_alpha_norm, e_alpha_unnormed[2]/e_alpha_norm]
			
					
			#Now we produce the observer frame (primed) coordinates; we
			#choose to align e_alpha along the component of the mean B
			#field in the plane of the sky
			e_alpha_p_not_perp = B_avg_vector
			e_alpha_p_perp = (e_alpha_p_not_perp - k_vector * np.dot(k_vector, e_alpha_p_not_perp))
			
			#If B_avg_vector has no component in the plane of the sky,
			#this if statement prevents divide-by-zero errors.
			if(np.linalg.norm(e_alpha_p_perp) == 0.):
				e_alpha_p_normed = e_alpha_p_not_perp
			else:
				e_alpha_p_normed = e_alpha_p_perp / np.linalg.norm(e_alpha_p_perp)
				
			#the angle between the symphony and observer frames is denoted
			#xi; the below statement generates an NxN array of xi, one
			#element for each point in the simulation
			xi = np.arccos((e_alpha_normed[0] * e_alpha_p_normed[0] +
		                       e_alpha_normed[1] * e_alpha_p_normed[1] +
		                       e_alpha_normed[2] * e_alpha_p_normed[2]))
			
			#numerical errors can sometimes make e_alpha dot e_alpha_p
			#slightly greater than 1, producing NaNs; the below code
			#changes those NaNs to arccos(1) = 0
			xi = np.nan_to_num(xi)
			
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
				#Note: rotation factor applied to produce correct
				#Stokes Q and U emission/absorption	
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
			                                        ) * np.cos(2.*xi[j][i])
				#Note: again, rotation factor applied below
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
			                                        ) * -np.sin(2.*xi[j][i]) 
			
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

			#prevent divide-by-zero errors
			if(exact_avg_I == 0):
				exact_avg_lin = 0.
				exact_avg_circ = 0.
			else:
				exact_avg_lin = np.sqrt(exact_avg_Q**2. + exact_avg_U**2.)/exact_avg_I
				exact_avg_circ = np.abs(exact_avg_V) / exact_avg_I

			#If theta = 0, emission and absorption are 0; the below
			#code chooses an arbitrary value of xi_avg to prevent
			#divide-by-zero errors (j_nu and alpha_nu are zero no
			#matter what value of xi_avg is chosen in this case)
			if(np.abs(np.dot(B_avg_vector, obs_vector) - 1) < 1e-7):
				xi_avg = np.pi/4.		 
		
			#analogous to the definitions of e_alpha and e_alpha'
			#above, the below code defines e_alpha_avg and
			#e_alpha_prime_avg
			else:
				e_alpha_avg = (B_avg_vector - k_vector * np.dot(B_avg_vector, k_vector))
	
				#prevent divide-by-zero errors when e_alpha gets normalized 
				if(np.linalg.norm(e_alpha_avg) == 0.):
					e_alpha_avg = B_avg_vector
				e_alpha_avg = e_alpha_avg / np.linalg.norm(e_alpha_avg)

				#again, sometimes numerical errors make the dot product of 
				#normalized vectors slightly > 1; prevent NaNs from arccos	
				if(np.dot(e_alpha_avg, e_alpha_p_normed) - 1. > 0 and np.dot(e_alpha_avg, e_alpha_p_normed) - 1. < 1e-6):
					xi_avg = np.arccos(1.)
				else:
					xi_avg = np.arccos(np.dot(e_alpha_avg, e_alpha_p_normed))
		
		#-------------calculate j_nu or alpha_nu of averaged sim. quantities--------#

			avgs_I      = j_nu_or_alpha_nu(nu, B_mag_avg, n_e_avg, 
			                       obs_angle_avg, distribution_function,
			                       sp.STOKES_I, theta_e, power_law_p,
			                       gamma_min, gamma_max, gamma_cutoff,
			                       kappa, kappa_width)
		
			avgs_Q      = j_nu_or_alpha_nu(nu, B_mag_avg, n_e_avg, 
			                       obs_angle_avg, distribution_function,
			                       sp.STOKES_Q, theta_e, power_law_p,
			                       gamma_min, gamma_max, gamma_cutoff,
			                       kappa, kappa_width) * np.cos(2.*xi_avg)   
			
			avgs_U      = j_nu_or_alpha_nu(nu, B_mag_avg, n_e_avg, 
			                       obs_angle_avg, distribution_function,
			                       sp.STOKES_Q, theta_e, power_law_p,
			                       gamma_min, gamma_max, gamma_cutoff,
			                       kappa, kappa_width) * -np.sin(2.*xi_avg)  
			
			avgs_V      = j_nu_or_alpha_nu(nu, B_mag_avg, n_e_avg, 
			                       obs_angle_avg, distribution_function,
			                       sp.STOKES_V, theta_e, power_law_p,
			                       gamma_min, gamma_max, gamma_cutoff,
			                       kappa, kappa_width)

			#prevent divide-by-zero errors
			if(avgs_I == 0):
				avgs_lin  = 0.
				avgs_circ = 0.
			else:
				avgs_lin = np.sqrt(avgs_Q**2. + avgs_U**2.) / avgs_I
				avgs_circ = np.abs(avgs_V) / avgs_I
		
			#calculate 'exact' (calculated point-by-point using local quantities),
			#'avg' (calculated using simulation-averaged quantities), and relative
			#difference abs(exact - avg)/exact
			if(rank == 0):
				relative_difference_I[y][x] += np.fabs((exact_avg_I - avgs_I)/exact_avg_I) / cone_resolution
	                	exact_avg_only_I[y][x]      += exact_avg_I / cone_resolution
	                	avgs_only_I[y][x]           += avgs_I / cone_resolution
	
				relative_difference_lin[y][x] += np.fabs((exact_avg_lin - avgs_lin)/exact_avg_lin) / cone_resolution
                                exact_avg_only_lin[y][x]      += exact_avg_lin / cone_resolution
                                avgs_only_lin[y][x]           += avgs_lin / cone_resolution

				relative_difference_circ[y][x] += np.fabs((exact_avg_circ - avgs_circ)/exact_avg_circ) / cone_resolution
                                exact_avg_only_circ[y][x]      += exact_avg_circ / cone_resolution
                                avgs_only_circ[y][x]           += avgs_circ / cone_resolution

			#keep track of nu and theta used for x and y axes of contour plots
			if(x == 0):
				obs_angle_used[y] = obs_angle_avg * 180. / np.pi
			nu_used[x] = nu / nu_c_avg
	
#------------------------make contour plots------------------------------------#

if (rank == 0):
   # Set plot parameters to make beautiful plots
   pl.rcParams['figure.figsize']  = 12, 7.5
   pl.rcParams['lines.linewidth'] = 1.5
   pl.rcParams['font.family']     = 'serif'
#   pl.rcParams['font.weight']     = 'bold'
   pl.rcParams['font.size']       = 15
   pl.rcParams['font.sans-serif'] = 'serif'
   pl.rcParams['text.usetex']     = False
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
   
   X, Y = np.meshgrid(nu_used, obs_angle_used)

   figure, ax = pl.subplots(3, 3, figsize=(10, 10))
   figure.suptitle(figure_title)
   
   plot1 = ax[0,0].contourf(np.log10(X), Y, exact_avg_only_I, 200)
   figure.colorbar(plot1, ax=ax[0,0])
   ax[0,0].set_title('$<j_\\nu(n, \mathbf{B})>$')
   
   
   plot2 = ax[0,1].contourf(np.log10(X), Y, avgs_only_I, 200)
   figure.colorbar(plot2, ax=ax[0,1])
   ax[0,1].set_title('$j_\\nu(<n>, <\mathbf{B}>)$')
   
   if(EMISS == False):
           ax[0,0].set_title('$<\\alpha_\\nu(n, \mathbf{B})>$')
           ax[0,1].set_title('$\\alpha_\\nu(<n>, <\mathbf{B}>)$')
   
   relative_difference_I = np.ma.array(relative_difference_I, mask=relative_difference_I > mask_tolerance)
   plot3 = ax[0,2].contourf(np.log10(X), Y, relative_difference_I, 200)
   figure.colorbar(plot3, ax=ax[0,2])
   ax[0,2].set_title('$|\mathrm{Error}|$')
   
   v = np.linspace(0., 1., 100)
   plot4 = ax[1,0].contourf(np.log10(X), Y, exact_avg_only_lin, v)
   cbar4 = figure.colorbar(plot4, ax=ax[1,0])

   plot5 = ax[1,1].contourf(np.log10(X), Y, avgs_only_lin, v)
   cbar5 = figure.colorbar(plot5, ax=ax[1,1])

   relative_difference_lin = np.ma.array(relative_difference_lin, mask=relative_difference_lin > mask_tolerance)
   plot6 = ax[1,2].contourf(np.log10(X), Y, relative_difference_lin, 200)
   figure.colorbar(plot6, ax=ax[1,2])

   #errors in the Stokes V formulae allow for Stokes V > Stokes I as theta->0; 
   # in this case Stokes V = Stokes I, so circ. pol. frac. = 1.
   exact_avg_only_circ[exact_avg_only_circ > 1.] = 1.
   avgs_only_circ[avgs_only_circ > 1.] = 1.

   plot10 = ax[2,0].contourf(np.log10(X), Y, exact_avg_only_circ, v)
   cbar10 = figure.colorbar(plot10, ax=ax[2,0])
   
   plot11 = ax[2,1].contourf(np.log10(X), Y, avgs_only_circ, v)
   cbar11 = figure.colorbar(plot11, ax=ax[2,1])

   relative_difference_circ = np.ma.array(relative_difference_circ, mask=relative_difference_circ > mask_tolerance)
   plot12 = ax[2,2].contourf(np.log10(X), Y, relative_difference_circ, 200)
   figure.colorbar(plot12, ax=ax[2,2])
   
   figure.add_subplot(111, frameon=False)
   pl.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
   pl.xlabel('$\\log_{10}(\\nu/\\overline{\\nu}_c)$', fontsize='large')
   pl.ylabel('$\\theta$ (deg)', fontsize='large')
   pl.tight_layout()
   pl.show()
#   figure.savefig('PL_10_10_5.png')
#   pl.close(figure)
