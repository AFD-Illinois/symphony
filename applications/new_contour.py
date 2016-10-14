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
distribution_function = sp.KAPPA_DIST	#distribution function
EMISS                 = True                    #True = j_nu, False = alpha_nu
IN_PLANE              = True		        #True = obs_angle in plane
mask_tolerance        = 1.			#error > tolerance is white
number_of_points      = 32
cone_resolution       = 50
max_nuratio           = 1e6
figure_title          = ''

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
   datafiles_path = ''
#   datafiles_path = '/home/alex/Documents/Spring_2016/'
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

B_vector     = [B_x, B_y, B_z]
B_avg_vector = [B_x_avg, B_y_avg, B_z_avg]
B_avg_vector = B_avg_vector/np.linalg.norm(B_avg_vector)
b_hat        = B_vector / np.linalg.norm(B_vector)


nu_c_avg = (electron_charge * B_mag_avg) / (2. * np.pi * m * c)

cone_angles1     = np.zeros([number_of_points])
cone_angles2     = np.zeros([number_of_points])
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

	if(obs_angle == 0.):
		return 0.


	if(EMISS == True):

		dim_prefactor = (n_e_avg * electron_charge**2. * 
				 electron_charge * B_mag_avg 
				 / (2. * np.pi * m * c) / c)

		return sp.j_nu_fit_py(nu, B, n_e, obs_angle, distribution_function,
                     polarization, theta_e, power_law_p, gamma_min,
                     gamma_max, gamma_cutoff, kappa, kappa_width) / dim_prefactor

	else:

		dim_prefactor = n_e_avg * electron_charge**2. / (nu * 
				m * c)

		return sp.alpha_nu_fit_py(nu, B, n_e, obs_angle, distribution_function,
                     polarization, theta_e, power_law_p, gamma_min,
                     gamma_max, gamma_cutoff, kappa, kappa_width)  

#-------------------------set up nu-theta space scan--------------------------#

#relative_difference_I = np.zeros([number_of_points, number_of_points])
#exact_avg_only_I      = np.zeros([number_of_points, number_of_points])
#avgs_only_I           = np.zeros([number_of_points, number_of_points])
#relative_difference_Q = np.zeros([number_of_points, number_of_points])
#exact_avg_only_Q      = np.zeros([number_of_points, number_of_points])
#avgs_only_Q           = np.zeros([number_of_points, number_of_points])
#relative_difference_U = np.zeros([number_of_points, number_of_points])
#exact_avg_only_U      = np.zeros([number_of_points, number_of_points])
#avgs_only_U           = np.zeros([number_of_points, number_of_points])
#relative_difference_V = np.zeros([number_of_points, number_of_points])
#exact_avg_only_V      = np.zeros([number_of_points, number_of_points])
#avgs_only_V           = np.zeros([number_of_points, number_of_points])


relative_difference_I_TOT = np.zeros([number_of_points, number_of_points])
exact_avg_only_I_TOT      = np.zeros([number_of_points, number_of_points])
avgs_only_I_TOT           = np.zeros([number_of_points, number_of_points])

relative_difference_Q_TOT = np.zeros([number_of_points, number_of_points])
exact_avg_only_Q_TOT      = np.zeros([number_of_points, number_of_points])
avgs_only_Q_TOT           = np.zeros([number_of_points, number_of_points])

relative_difference_U_TOT = np.zeros([number_of_points, number_of_points])
exact_avg_only_U_TOT      = np.zeros([number_of_points, number_of_points])
avgs_only_U_TOT           = np.zeros([number_of_points, number_of_points])

relative_difference_V_TOT = np.zeros([number_of_points, number_of_points])
exact_avg_only_V_TOT      = np.zeros([number_of_points, number_of_points])
avgs_only_V_TOT           = np.zeros([number_of_points, number_of_points])



nu_used             = np.empty([number_of_points])
obs_angle_used      = np.empty([number_of_points])

avgs_I = 0.
avgs_Q = 0.
avgs_U = 0.
avgs_V = 0.


for cone_point in range(cone_resolution):
	if(rank == 0):
		print 1.0*cone_point / cone_resolution * 100., '%'
	for x in range(0, number_of_points):
		nu = nu_c_avg * 10.**(1.*x / (number_of_points - 1) * np.log10(max_nuratio))
#	        if (rank == 0):
#		   print 100.0*x/(number_of_points - 1), '% complete'
		for y in range(0, number_of_points):
	
			if(IN_PLANE == True):
				#rotates observer vector in plane
				axis_of_rot = [0, B_avg_vector[2], -B_avg_vector[1]]
			else:
				#rotates observer vector out of plane
				axis_of_rot = [-B_avg_vector[1], B_avg_vector[0], 0]
			
			#desired angle between observer vector and mean B
			theta       = (1.0*y/number_of_points * np.pi/2.)
			
			#		print theta * 180./np.pi
			
			obs_vector  = np.dot(rotation_matrix(axis_of_rot, theta), 
					     B_avg_vector)
			
			obs_angle = np.arccos((B_x * obs_vector[0] 
			                       + B_y * obs_vector[1]
			                       + B_z * obs_vector[2])
			                       / (np.linalg.norm(obs_vector)*B_mag))
			
			
			#		obs_angle_avg  = np.mean(obs_angle) #not physically correct
			angle_to_mean_field = np.arccos(np.dot(B_avg_vector, obs_vector)) 
			obs_angle_avg = angle_to_mean_field
			
			
			#-----------------------avg over cone angle------------------------------------#
			#for cone_point in range(cone_resolution):
			
			
			cone_rot = (1.0 * cone_point / cone_resolution * 2.*np.pi)
	
			#cone_rot = np.pi/5.		
	
			#cone_angles_used[cone_point] = cone_rot
			
			k_parallel = (np.dot(obs_vector, B_avg_vector) 
			              / np.dot(B_avg_vector, B_avg_vector)) * B_avg_vector
			
			k_perp = obs_vector - k_parallel
			
			k_perp_rotated = np.dot(rotation_matrix(k_parallel, cone_rot), 
		        	                k_perp)
			
			k_vector = k_perp_rotated + k_parallel
	
			#-------------------------BEGIN MODIFICATIONS----------------------------------#
			
				#create NxN array of symphony frame e_alpha axis #TODO: check this
			obs_vector_array = [np.full(np.shape(B_x), k_vector[0]),
			                   np.full(np.shape(B_x), k_vector[1]),
			                   np.full(np.shape(B_x), k_vector[2])]  
			
			b_hat_dot_k_hat = b_hat[0] * obs_vector_array[0] + b_hat[1] * obs_vector_array[1] + b_hat[2] * obs_vector_array[2]
			
			
			e_alpha_unnormed = (b_hat - obs_vector_array * b_hat_dot_k_hat)
			e_alpha_norm = np.sqrt(e_alpha_unnormed[0]**2. + e_alpha_unnormed[1]**2. + e_alpha_unnormed[2]**2.)

			e_alpha_normed = [e_alpha_unnormed[0]/e_alpha_norm, e_alpha_unnormed[1]/e_alpha_norm, e_alpha_unnormed[2]/e_alpha_norm]
			
			
			
			#create e_alpha_prime #TODO: check this
			e_alpha_p_not_perp = B_avg_vector
			
			
			e_alpha_p_perp = (e_alpha_p_not_perp - k_vector * np.dot(k_vector, e_alpha_p_not_perp))
			
			#TODO: CHECK THIS
			if(np.linalg.norm(e_alpha_p_perp) == 0.):
				e_alpha_p_normed = e_alpha_p_not_perp
			else:
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
			                                        ) * np.cos(2.*xi[j][i])
			
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
			
		#--------------MORE MODIFICATIONS---------------------------------------------#
			if(np.abs(np.dot(B_avg_vector, obs_vector) - 1) < 1e-7):
				xi_avg = np.pi/4. #choose it arbitrarily; emission/absorption will be zero anyway		 
		
			else:
				e_alpha_avg = (B_avg_vector - k_vector * np.dot(B_avg_vector, k_vector))

			#TODO: testing this
				if(np.linalg.norm(e_alpha_avg) == 0.):
					e_alpha_avg = B_avg_vector

				e_alpha_avg = e_alpha_avg / np.linalg.norm(e_alpha_avg)
				if(np.dot(e_alpha_avg, e_alpha_p_normed) - 1. > 0 and np.dot(e_alpha_avg, e_alpha_p_normed) - 1. < 1e-6):
					xi_avg = np.arccos(1.)
				else:
					xi_avg = np.arccos(np.dot(e_alpha_avg, e_alpha_p_normed))
		
		#---------------END MODIFICATIONS---------------------------------------------#
		
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
		
#			relative_difference_I[y][x] = np.fabs((exact_avg_I - avgs_I)/exact_avg_I)
#			exact_avg_only_I[y][x]      = exact_avg_I
#			avgs_only_I[y][x]           = avgs_I
#		
#			relative_difference_Q[y][x] = np.fabs((exact_avg_Q - avgs_Q)/exact_avg_Q)
#		        exact_avg_only_Q[y][x]      = exact_avg_Q
#		        avgs_only_Q[y][x]           = avgs_Q
#		
#			relative_difference_U[y][x] = np.fabs((exact_avg_U - avgs_U)/exact_avg_U)
#		        exact_avg_only_U[y][x]      = exact_avg_U
#		        avgs_only_U[y][x]           = avgs_U
#		
#			relative_difference_V[y][x] = np.fabs((exact_avg_V - avgs_V)/exact_avg_V)
#		        exact_avg_only_V[y][x]      = exact_avg_V
#		        avgs_only_V[y][x]           = avgs_V
		

			if(rank == 0):
				relative_difference_I_TOT[y][x] += np.fabs((exact_avg_I - avgs_I)/exact_avg_I) / cone_resolution
                        	exact_avg_only_I_TOT[y][x]      += exact_avg_I / cone_resolution
                        	avgs_only_I_TOT[y][x]           += avgs_I / cone_resolution
	
				relative_difference_Q_TOT[y][x] += np.fabs((exact_avg_Q - avgs_Q)/exact_avg_Q) / cone_resolution
                                exact_avg_only_Q_TOT[y][x]      += exact_avg_Q / cone_resolution
                                avgs_only_Q_TOT[y][x]           += avgs_Q / cone_resolution

				relative_difference_U_TOT[y][x] += np.fabs((exact_avg_U - avgs_U)/exact_avg_U) / cone_resolution
                                exact_avg_only_U_TOT[y][x]      += exact_avg_U / cone_resolution
                                avgs_only_U_TOT[y][x]           += avgs_U / cone_resolution

				relative_difference_V_TOT[y][x] += np.fabs((exact_avg_V - avgs_V)/exact_avg_V) / cone_resolution
                                exact_avg_only_V_TOT[y][x]      += exact_avg_V / cone_resolution
                                avgs_only_V_TOT[y][x]           += avgs_V / cone_resolution

	
			if(x == 0):
				obs_angle_used[y] = obs_angle_avg * 180. / np.pi
			nu_used[x] = nu / nu_c_avg
	
#------------------------make contour plot-------------------------------------#

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

   X, Y = np.meshgrid(nu_used, obs_angle_used)

   figure, ax = pl.subplots(3, 3, figsize=(10, 10))
   figure.suptitle(figure_title)
   
   plot1 = ax[0,0].contourf(np.log10(X), Y, exact_avg_only_I_TOT, 200)
   figure.colorbar(plot1, ax=ax[0,0])
   ax[0,0].set_title('$<j_\\nu(n, \mathbf{B})>$')
   
   
   plot2 = ax[0,1].contourf(np.log10(X), Y, avgs_only_I_TOT, 200)
   figure.colorbar(plot2, ax=ax[0,1])
   ax[0,1].set_title('$j_\\nu(<n>, <\mathbf{B}>)$')
   
   if(EMISS == False):
           ax[0,0].set_title('$<\\alpha_\\nu(n, \mathbf{B})>$')
           ax[0,1].set_title('$\\alpha_\\nu(<n>, <\mathbf{B}>)$')
   
   relative_difference_I_TOT = np.ma.array(relative_difference_I_TOT, mask=relative_difference_I_TOT > mask_tolerance)
   plot3 = ax[0,2].contourf(np.log10(X), Y, relative_difference_I_TOT, 200)
   figure.colorbar(plot3, ax=ax[0,2])
   ax[0,2].set_title('$|\mathrm{Error}|$')
   
   plot4 = ax[1,0].contourf(np.log10(X), Y, exact_avg_only_Q_TOT, 200)
   figure.colorbar(plot4, ax=ax[1,0])
   
   plot5 = ax[1,1].contourf(np.log10(X), Y, avgs_only_Q_TOT, 200)
   figure.colorbar(plot5, ax=ax[1,1])
   
   
   relative_difference_Q_TOT = np.ma.array(relative_difference_Q_TOT, mask=relative_difference_Q_TOT > mask_tolerance)
   plot6 = ax[1,2].contourf(np.log10(X), Y, relative_difference_Q_TOT, 200)
   figure.colorbar(plot6, ax=ax[1,2])
   
#   plot7 = ax[2,0].contourf(np.log10(X), Y, exact_avg_only_U_TOT, 200)
#   figure.colorbar(plot7, ax=ax[2,0])
#   
#   plot8 = ax[2,1].contourf(np.log10(X), Y, avgs_only_U_TOT, 200)
#   figure.colorbar(plot8, ax=ax[2,1])
#   
#   relative_difference_U_TOT = np.ma.array(relative_difference_U_TOT, mask=relative_difference_U_TOT > mask_tolerance)
#   plot9 = ax[2,2].contourf(np.log10(X), Y, relative_difference_U_TOT, 200)
#   figure.colorbar(plot9, ax=ax[2,2])
   
   plot10 = ax[2,0].contourf(np.log10(X), Y, exact_avg_only_V_TOT, 200)
   figure.colorbar(plot10, ax=ax[2,0])
   
   plot11 = ax[2,1].contourf(np.log10(X), Y, avgs_only_V_TOT, 200)
   figure.colorbar(plot11, ax=ax[2,1])
   
   relative_difference_V_TOT = np.ma.array(relative_difference_V_TOT, mask=relative_difference_V_TOT > mask_tolerance)
   plot12 = ax[2,2].contourf(np.log10(X), Y, relative_difference_V_TOT, 200)
   figure.colorbar(plot12, ax=ax[2,2])
   
   figure.add_subplot(111, frameon=False)
   pl.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
   pl.xlabel('$\\log_{10}(\\nu/\\overline{\\nu}_c)$', fontsize='large')
   pl.ylabel('$\\theta$ (deg)', fontsize='large')
   pl.tight_layout()
   pl.show()


   #print error at position where j_nu or alpha_nu is maximal
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

#   print '\nError at location of max j_nu or alpha_nu: '
#   print 'STOKES I: ', relative_difference_I[x_I][y_I]
#   print 'STOKES Q: ', relative_difference_Q[x_Q][y_Q]
#   print 'STOKES U: ', relative_difference_U[x_U][y_U]
#   print 'STOKES V: ', relative_difference_V[x_V][y_V]
