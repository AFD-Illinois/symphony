import numpy as np
import pylab as pl

from scipy.interpolate import interp2d
import scipy.special as special
from scipy.integrate import quad, fixed_quad

#set up spline interpolation datafiles

main_directory = '/home/aapandy2/chi_ij_files_theta/'

gamma = np.loadtxt(main_directory + 'step_array.txt')

gam = gamma
omrat = gamma

gamma_min = gamma[0]
gamma_max = gamma[np.size(gamma)-1]

#load datafiles used to produce spline fits

comp_folder = 'chi_11_real/'

chi_11_real_01 = np.loadtxt(main_directory + comp_folder + 'chi_11_real_mod_step01.txt')
chi_11_real_10 = np.loadtxt(main_directory + comp_folder + 'chi_11_real_mod_step10.txt')
chi_11_real_20 = np.loadtxt(main_directory + comp_folder + 'chi_11_real_mod_step20.txt')
chi_11_real_30 = np.loadtxt(main_directory + comp_folder + 'chi_11_real_mod_step30.txt')
chi_11_real_40 = np.loadtxt(main_directory + comp_folder + 'chi_11_real_mod_step40.txt')
chi_11_real_50 = np.loadtxt(main_directory + comp_folder + 'chi_11_real_mod_step50.txt')
chi_11_real_60 = np.loadtxt(main_directory + comp_folder + 'chi_11_real_mod_step60.txt')
chi_11_real_70 = np.loadtxt(main_directory + comp_folder + 'chi_11_real_mod_step70.txt')
chi_11_real_80 = np.loadtxt(main_directory + comp_folder + 'chi_11_real_mod_step80.txt')
chi_11_real_89 = np.loadtxt(main_directory + comp_folder + 'chi_11_real_mod_step89.txt')

comp_folder = 'chi_11_imag/'

chi_11_imag_01 = np.loadtxt(main_directory + comp_folder + 'chi_11_imag_mod_step01.txt')
chi_11_imag_10 = np.loadtxt(main_directory + comp_folder + 'chi_11_imag_mod_step10.txt')
chi_11_imag_20 = np.loadtxt(main_directory + comp_folder + 'chi_11_imag_mod_step20.txt')
chi_11_imag_30 = np.loadtxt(main_directory + comp_folder + 'chi_11_imag_mod_step30.txt')
chi_11_imag_40 = np.loadtxt(main_directory + comp_folder + 'chi_11_imag_mod_step40.txt')
chi_11_imag_50 = np.loadtxt(main_directory + comp_folder + 'chi_11_imag_mod_step50.txt')
chi_11_imag_60 = np.loadtxt(main_directory + comp_folder + 'chi_11_imag_mod_step60.txt')
chi_11_imag_70 = np.loadtxt(main_directory + comp_folder + 'chi_11_imag_mod_step70.txt')
chi_11_imag_80 = np.loadtxt(main_directory + comp_folder + 'chi_11_imag_mod_step80.txt')
chi_11_imag_89 = np.loadtxt(main_directory + comp_folder + 'chi_11_imag_mod_step89.txt')


comp_folder = 'chi_12_real/'

chi_12_real_01 = np.loadtxt(main_directory + comp_folder + 'chi_12_real_mod_step01.txt')
chi_12_real_10 = np.loadtxt(main_directory + comp_folder + 'chi_12_real_mod_step10.txt')
chi_12_real_20 = np.loadtxt(main_directory + comp_folder + 'chi_12_real_mod_step20.txt')
chi_12_real_30 = np.loadtxt(main_directory + comp_folder + 'chi_12_real_mod_step30.txt')
chi_12_real_40 = np.loadtxt(main_directory + comp_folder + 'chi_12_real_mod_step40.txt')
chi_12_real_50 = np.loadtxt(main_directory + comp_folder + 'chi_12_real_mod_step50.txt')
chi_12_real_60 = np.loadtxt(main_directory + comp_folder + 'chi_12_real_mod_step60.txt')
chi_12_real_70 = np.loadtxt(main_directory + comp_folder + 'chi_12_real_mod_step70.txt')
chi_12_real_80 = np.loadtxt(main_directory + comp_folder + 'chi_12_real_mod_step80.txt')
chi_12_real_89 = np.loadtxt(main_directory + comp_folder + 'chi_12_real_mod_step89.txt')

comp_folder = 'chi_12_imag/'

chi_12_imag_01 = np.loadtxt(main_directory + comp_folder + 'chi_12_imag_mod_step01.txt')
chi_12_imag_10 = np.loadtxt(main_directory + comp_folder + 'chi_12_imag_mod_step10.txt')
chi_12_imag_20 = np.loadtxt(main_directory + comp_folder + 'chi_12_imag_mod_step20.txt')
chi_12_imag_30 = np.loadtxt(main_directory + comp_folder + 'chi_12_imag_mod_step30.txt')
chi_12_imag_40 = np.loadtxt(main_directory + comp_folder + 'chi_12_imag_mod_step40.txt')
chi_12_imag_50 = np.loadtxt(main_directory + comp_folder + 'chi_12_imag_mod_step50.txt')
chi_12_imag_60 = np.loadtxt(main_directory + comp_folder + 'chi_12_imag_mod_step60.txt')
chi_12_imag_70 = np.loadtxt(main_directory + comp_folder + 'chi_12_imag_mod_step70.txt')
chi_12_imag_80 = np.loadtxt(main_directory + comp_folder + 'chi_12_imag_mod_step80.txt')
chi_12_imag_89 = np.loadtxt(main_directory + comp_folder + 'chi_12_imag_mod_step89.txt')


comp_folder = 'chi_13_real/'

chi_13_real_01 = np.loadtxt(main_directory + comp_folder + 'chi_13_real_mod_step01.txt')
chi_13_real_10 = np.loadtxt(main_directory + comp_folder + 'chi_13_real_mod_step10.txt')
chi_13_real_20 = np.loadtxt(main_directory + comp_folder + 'chi_13_real_mod_step20.txt')
chi_13_real_30 = np.loadtxt(main_directory + comp_folder + 'chi_13_real_mod_step30.txt')
chi_13_real_40 = np.loadtxt(main_directory + comp_folder + 'chi_13_real_mod_step40.txt')
chi_13_real_50 = np.loadtxt(main_directory + comp_folder + 'chi_13_real_mod_step50.txt')
chi_13_real_60 = np.loadtxt(main_directory + comp_folder + 'chi_13_real_mod_step60.txt')
chi_13_real_70 = np.loadtxt(main_directory + comp_folder + 'chi_13_real_mod_step70.txt')
chi_13_real_80 = np.loadtxt(main_directory + comp_folder + 'chi_13_real_mod_step80.txt')
chi_13_real_89 = np.loadtxt(main_directory + comp_folder + 'chi_13_real_mod_step89.txt')

comp_folder = 'chi_13_imag/'

chi_13_imag_01 = np.loadtxt(main_directory + comp_folder + 'chi_13_imag_mod_step01.txt')
chi_13_imag_10 = np.loadtxt(main_directory + comp_folder + 'chi_13_imag_mod_step10.txt')
chi_13_imag_20 = np.loadtxt(main_directory + comp_folder + 'chi_13_imag_mod_step20.txt')
chi_13_imag_30 = np.loadtxt(main_directory + comp_folder + 'chi_13_imag_mod_step30.txt')
chi_13_imag_40 = np.loadtxt(main_directory + comp_folder + 'chi_13_imag_mod_step40.txt')
chi_13_imag_50 = np.loadtxt(main_directory + comp_folder + 'chi_13_imag_mod_step50.txt')
chi_13_imag_60 = np.loadtxt(main_directory + comp_folder + 'chi_13_imag_mod_step60.txt')
chi_13_imag_70 = np.loadtxt(main_directory + comp_folder + 'chi_13_imag_mod_step70.txt')
chi_13_imag_80 = np.loadtxt(main_directory + comp_folder + 'chi_13_imag_mod_step80.txt')
chi_13_imag_89 = np.loadtxt(main_directory + comp_folder + 'chi_13_imag_mod_step89.txt')

comp_folder = 'chi_22_real/'

chi_22_real_01 = np.loadtxt(main_directory + comp_folder + 'chi_22_real_mod_step01.txt')
chi_22_real_10 = np.loadtxt(main_directory + comp_folder + 'chi_22_real_mod_step10.txt')
chi_22_real_20 = np.loadtxt(main_directory + comp_folder + 'chi_22_real_mod_step20.txt')
chi_22_real_30 = np.loadtxt(main_directory + comp_folder + 'chi_22_real_mod_step30.txt')
chi_22_real_40 = np.loadtxt(main_directory + comp_folder + 'chi_22_real_mod_step40.txt')
chi_22_real_50 = np.loadtxt(main_directory + comp_folder + 'chi_22_real_mod_step50.txt')
chi_22_real_60 = np.loadtxt(main_directory + comp_folder + 'chi_22_real_mod_step60.txt')
chi_22_real_70 = np.loadtxt(main_directory + comp_folder + 'chi_22_real_mod_step70.txt')
chi_22_real_80 = np.loadtxt(main_directory + comp_folder + 'chi_22_real_mod_step80.txt')
chi_22_real_89 = np.loadtxt(main_directory + comp_folder + 'chi_22_real_mod_step89.txt')

comp_folder = 'chi_22_imag/'

chi_22_imag_01 = np.loadtxt(main_directory + comp_folder + 'chi_22_imag_mod_step01.txt')
chi_22_imag_10 = np.loadtxt(main_directory + comp_folder + 'chi_22_imag_mod_step10.txt')
chi_22_imag_20 = np.loadtxt(main_directory + comp_folder + 'chi_22_imag_mod_step20.txt')
chi_22_imag_30 = np.loadtxt(main_directory + comp_folder + 'chi_22_imag_mod_step30.txt')
chi_22_imag_40 = np.loadtxt(main_directory + comp_folder + 'chi_22_imag_mod_step40.txt')
chi_22_imag_50 = np.loadtxt(main_directory + comp_folder + 'chi_22_imag_mod_step50.txt')
chi_22_imag_60 = np.loadtxt(main_directory + comp_folder + 'chi_22_imag_mod_step60.txt')
chi_22_imag_70 = np.loadtxt(main_directory + comp_folder + 'chi_22_imag_mod_step70.txt')
chi_22_imag_80 = np.loadtxt(main_directory + comp_folder + 'chi_22_imag_mod_step80.txt')
chi_22_imag_89 = np.loadtxt(main_directory + comp_folder + 'chi_22_imag_mod_step89.txt')

comp_folder = 'chi_32_real/'

chi_32_real_01 = np.loadtxt(main_directory + comp_folder + 'chi_32_real_mod_step01.txt')
chi_32_real_10 = np.loadtxt(main_directory + comp_folder + 'chi_32_real_mod_step10.txt')
chi_32_real_20 = np.loadtxt(main_directory + comp_folder + 'chi_32_real_mod_step20.txt')
chi_32_real_30 = np.loadtxt(main_directory + comp_folder + 'chi_32_real_mod_step30.txt')
chi_32_real_40 = np.loadtxt(main_directory + comp_folder + 'chi_32_real_mod_step40.txt')
chi_32_real_50 = np.loadtxt(main_directory + comp_folder + 'chi_32_real_mod_step50.txt')
chi_32_real_60 = np.loadtxt(main_directory + comp_folder + 'chi_32_real_mod_step60.txt')
chi_32_real_70 = np.loadtxt(main_directory + comp_folder + 'chi_32_real_mod_step70.txt')
chi_32_real_80 = np.loadtxt(main_directory + comp_folder + 'chi_32_real_mod_step80.txt')
chi_32_real_89 = np.loadtxt(main_directory + comp_folder + 'chi_32_real_mod_step89.txt')

comp_folder = 'chi_32_imag/'

chi_32_imag_01 = np.loadtxt(main_directory + comp_folder + 'chi_32_imag_mod_step01.txt')
chi_32_imag_10 = np.loadtxt(main_directory + comp_folder + 'chi_32_imag_mod_step10.txt')
chi_32_imag_20 = np.loadtxt(main_directory + comp_folder + 'chi_32_imag_mod_step20.txt')
chi_32_imag_30 = np.loadtxt(main_directory + comp_folder + 'chi_32_imag_mod_step30.txt')
chi_32_imag_40 = np.loadtxt(main_directory + comp_folder + 'chi_32_imag_mod_step40.txt')
chi_32_imag_50 = np.loadtxt(main_directory + comp_folder + 'chi_32_imag_mod_step50.txt')
chi_32_imag_60 = np.loadtxt(main_directory + comp_folder + 'chi_32_imag_mod_step60.txt')
chi_32_imag_70 = np.loadtxt(main_directory + comp_folder + 'chi_32_imag_mod_step70.txt')
chi_32_imag_80 = np.loadtxt(main_directory + comp_folder + 'chi_32_imag_mod_step80.txt')
chi_32_imag_89 = np.loadtxt(main_directory + comp_folder + 'chi_32_imag_mod_step89.txt')

comp_folder = 'chi_33_real/'

chi_33_real_01 = np.loadtxt(main_directory + comp_folder + 'chi_33_real_mod_step01.txt')
chi_33_real_10 = np.loadtxt(main_directory + comp_folder + 'chi_33_real_mod_step10.txt')
chi_33_real_20 = np.loadtxt(main_directory + comp_folder + 'chi_33_real_mod_step20.txt')
chi_33_real_30 = np.loadtxt(main_directory + comp_folder + 'chi_33_real_mod_step30.txt')
chi_33_real_40 = np.loadtxt(main_directory + comp_folder + 'chi_33_real_mod_step40.txt')
chi_33_real_50 = np.loadtxt(main_directory + comp_folder + 'chi_33_real_mod_step50.txt')
chi_33_real_60 = np.loadtxt(main_directory + comp_folder + 'chi_33_real_mod_step60.txt')
chi_33_real_70 = np.loadtxt(main_directory + comp_folder + 'chi_33_real_mod_step70.txt')
chi_33_real_80 = np.loadtxt(main_directory + comp_folder + 'chi_33_real_mod_step80.txt')
chi_33_real_89 = np.loadtxt(main_directory + comp_folder + 'chi_33_real_mod_step89.txt')

comp_folder = 'chi_33_imag/'

chi_33_imag_01 = np.loadtxt(main_directory + comp_folder + 'chi_33_imag_mod_step01.txt')
chi_33_imag_10 = np.loadtxt(main_directory + comp_folder + 'chi_33_imag_mod_step10.txt')
chi_33_imag_20 = np.loadtxt(main_directory + comp_folder + 'chi_33_imag_mod_step20.txt')
chi_33_imag_30 = np.loadtxt(main_directory + comp_folder + 'chi_33_imag_mod_step30.txt')
chi_33_imag_40 = np.loadtxt(main_directory + comp_folder + 'chi_33_imag_mod_step40.txt')
chi_33_imag_50 = np.loadtxt(main_directory + comp_folder + 'chi_33_imag_mod_step50.txt')
chi_33_imag_60 = np.loadtxt(main_directory + comp_folder + 'chi_33_imag_mod_step60.txt')
chi_33_imag_70 = np.loadtxt(main_directory + comp_folder + 'chi_33_imag_mod_step70.txt')
chi_33_imag_80 = np.loadtxt(main_directory + comp_folder + 'chi_33_imag_mod_step80.txt')
chi_33_imag_89 = np.loadtxt(main_directory + comp_folder + 'chi_33_imag_mod_step89.txt')

#Generate splines#

chi_11_real_spline_01 = interp2d(gam, omrat, chi_11_real_01)
chi_11_real_spline_10 = interp2d(gam, omrat, chi_11_real_10)
chi_11_real_spline_20 = interp2d(gam, omrat, chi_11_real_20)
chi_11_real_spline_30 = interp2d(gam, omrat, chi_11_real_30)
chi_11_real_spline_40 = interp2d(gam, omrat, chi_11_real_40)
chi_11_real_spline_50 = interp2d(gam, omrat, chi_11_real_50)
chi_11_real_spline_60 = interp2d(gam, omrat, chi_11_real_60)
chi_11_real_spline_70 = interp2d(gam, omrat, chi_11_real_70)
chi_11_real_spline_80 = interp2d(gam, omrat, chi_11_real_80)
chi_11_real_spline_89 = interp2d(gam, omrat, chi_11_real_89)

chi_11_imag_spline_01 = interp2d(gam, omrat, chi_11_imag_01)
chi_11_imag_spline_10 = interp2d(gam, omrat, chi_11_imag_10)
chi_11_imag_spline_20 = interp2d(gam, omrat, chi_11_imag_20)
chi_11_imag_spline_30 = interp2d(gam, omrat, chi_11_imag_30)
chi_11_imag_spline_40 = interp2d(gam, omrat, chi_11_imag_40)
chi_11_imag_spline_50 = interp2d(gam, omrat, chi_11_imag_50)
chi_11_imag_spline_60 = interp2d(gam, omrat, chi_11_imag_60)
chi_11_imag_spline_70 = interp2d(gam, omrat, chi_11_imag_70)
chi_11_imag_spline_80 = interp2d(gam, omrat, chi_11_imag_80)
chi_11_imag_spline_89 = interp2d(gam, omrat, chi_11_imag_89)

chi_12_real_spline_01 = interp2d(gam, omrat, chi_12_real_01)
chi_12_real_spline_10 = interp2d(gam, omrat, chi_12_real_10)
chi_12_real_spline_20 = interp2d(gam, omrat, chi_12_real_20)
chi_12_real_spline_30 = interp2d(gam, omrat, chi_12_real_30)
chi_12_real_spline_40 = interp2d(gam, omrat, chi_12_real_40)
chi_12_real_spline_50 = interp2d(gam, omrat, chi_12_real_50)
chi_12_real_spline_60 = interp2d(gam, omrat, chi_12_real_60)
chi_12_real_spline_70 = interp2d(gam, omrat, chi_12_real_70)
chi_12_real_spline_80 = interp2d(gam, omrat, chi_12_real_80)
chi_12_real_spline_89 = interp2d(gam, omrat, chi_12_real_89)

chi_12_imag_spline_01 = interp2d(gam, omrat, chi_12_imag_01)
chi_12_imag_spline_10 = interp2d(gam, omrat, chi_12_imag_10)
chi_12_imag_spline_20 = interp2d(gam, omrat, chi_12_imag_20)
chi_12_imag_spline_30 = interp2d(gam, omrat, chi_12_imag_30)
chi_12_imag_spline_40 = interp2d(gam, omrat, chi_12_imag_40)
chi_12_imag_spline_50 = interp2d(gam, omrat, chi_12_imag_50)
chi_12_imag_spline_60 = interp2d(gam, omrat, chi_12_imag_60)
chi_12_imag_spline_70 = interp2d(gam, omrat, chi_12_imag_70)
chi_12_imag_spline_80 = interp2d(gam, omrat, chi_12_imag_80)
chi_12_imag_spline_89 = interp2d(gam, omrat, chi_12_imag_89)

chi_13_real_spline_01 = interp2d(gam, omrat, chi_13_real_01)
chi_13_real_spline_10 = interp2d(gam, omrat, chi_13_real_10)
chi_13_real_spline_20 = interp2d(gam, omrat, chi_13_real_20)
chi_13_real_spline_30 = interp2d(gam, omrat, chi_13_real_30)
chi_13_real_spline_40 = interp2d(gam, omrat, chi_13_real_40)
chi_13_real_spline_50 = interp2d(gam, omrat, chi_13_real_50)
chi_13_real_spline_60 = interp2d(gam, omrat, chi_13_real_60)
chi_13_real_spline_70 = interp2d(gam, omrat, chi_13_real_70)
chi_13_real_spline_80 = interp2d(gam, omrat, chi_13_real_80)
chi_13_real_spline_89 = interp2d(gam, omrat, chi_13_real_89)

chi_13_imag_spline_01 = interp2d(gam, omrat, chi_13_imag_01)
chi_13_imag_spline_10 = interp2d(gam, omrat, chi_13_imag_10)
chi_13_imag_spline_20 = interp2d(gam, omrat, chi_13_imag_20)
chi_13_imag_spline_30 = interp2d(gam, omrat, chi_13_imag_30)
chi_13_imag_spline_40 = interp2d(gam, omrat, chi_13_imag_40)
chi_13_imag_spline_50 = interp2d(gam, omrat, chi_13_imag_50)
chi_13_imag_spline_60 = interp2d(gam, omrat, chi_13_imag_60)
chi_13_imag_spline_70 = interp2d(gam, omrat, chi_13_imag_70)
chi_13_imag_spline_80 = interp2d(gam, omrat, chi_13_imag_80)
chi_13_imag_spline_89 = interp2d(gam, omrat, chi_13_imag_89)

chi_22_real_spline_01 = interp2d(gam, omrat, chi_22_real_01)
chi_22_real_spline_10 = interp2d(gam, omrat, chi_22_real_10)
chi_22_real_spline_20 = interp2d(gam, omrat, chi_22_real_20)
chi_22_real_spline_30 = interp2d(gam, omrat, chi_22_real_30)
chi_22_real_spline_40 = interp2d(gam, omrat, chi_22_real_40)
chi_22_real_spline_50 = interp2d(gam, omrat, chi_22_real_50)
chi_22_real_spline_60 = interp2d(gam, omrat, chi_22_real_60)
chi_22_real_spline_70 = interp2d(gam, omrat, chi_22_real_70)
chi_22_real_spline_80 = interp2d(gam, omrat, chi_22_real_80)
chi_22_real_spline_89 = interp2d(gam, omrat, chi_22_real_89)

chi_22_imag_spline_01 = interp2d(gam, omrat, chi_22_imag_01)
chi_22_imag_spline_10 = interp2d(gam, omrat, chi_22_imag_10)
chi_22_imag_spline_20 = interp2d(gam, omrat, chi_22_imag_20)
chi_22_imag_spline_30 = interp2d(gam, omrat, chi_22_imag_30)
chi_22_imag_spline_40 = interp2d(gam, omrat, chi_22_imag_40)
chi_22_imag_spline_50 = interp2d(gam, omrat, chi_22_imag_50)
chi_22_imag_spline_60 = interp2d(gam, omrat, chi_22_imag_60)
chi_22_imag_spline_70 = interp2d(gam, omrat, chi_22_imag_70)
chi_22_imag_spline_80 = interp2d(gam, omrat, chi_22_imag_80)
chi_22_imag_spline_89 = interp2d(gam, omrat, chi_22_imag_89)

chi_32_real_spline_01 = interp2d(gam, omrat, chi_32_real_01)
chi_32_real_spline_10 = interp2d(gam, omrat, chi_32_real_10)
chi_32_real_spline_20 = interp2d(gam, omrat, chi_32_real_20)
chi_32_real_spline_30 = interp2d(gam, omrat, chi_32_real_30)
chi_32_real_spline_40 = interp2d(gam, omrat, chi_32_real_40)
chi_32_real_spline_50 = interp2d(gam, omrat, chi_32_real_50)
chi_32_real_spline_60 = interp2d(gam, omrat, chi_32_real_60)
chi_32_real_spline_70 = interp2d(gam, omrat, chi_32_real_70)
chi_32_real_spline_80 = interp2d(gam, omrat, chi_32_real_80)
chi_32_real_spline_89 = interp2d(gam, omrat, chi_32_real_89)

chi_32_imag_spline_01 = interp2d(gam, omrat, chi_32_imag_01)
chi_32_imag_spline_10 = interp2d(gam, omrat, chi_32_imag_10)
chi_32_imag_spline_20 = interp2d(gam, omrat, chi_32_imag_20)
chi_32_imag_spline_30 = interp2d(gam, omrat, chi_32_imag_30)
chi_32_imag_spline_40 = interp2d(gam, omrat, chi_32_imag_40)
chi_32_imag_spline_50 = interp2d(gam, omrat, chi_32_imag_50)
chi_32_imag_spline_60 = interp2d(gam, omrat, chi_32_imag_60)
chi_32_imag_spline_70 = interp2d(gam, omrat, chi_32_imag_70)
chi_32_imag_spline_80 = interp2d(gam, omrat, chi_32_imag_80)
chi_32_imag_spline_89 = interp2d(gam, omrat, chi_32_imag_89)

chi_33_real_spline_01 = interp2d(gam, omrat, chi_33_real_01)
chi_33_real_spline_10 = interp2d(gam, omrat, chi_33_real_10)
chi_33_real_spline_20 = interp2d(gam, omrat, chi_33_real_20)
chi_33_real_spline_30 = interp2d(gam, omrat, chi_33_real_30)
chi_33_real_spline_40 = interp2d(gam, omrat, chi_33_real_40)
chi_33_real_spline_50 = interp2d(gam, omrat, chi_33_real_50)
chi_33_real_spline_60 = interp2d(gam, omrat, chi_33_real_60)
chi_33_real_spline_70 = interp2d(gam, omrat, chi_33_real_70)
chi_33_real_spline_80 = interp2d(gam, omrat, chi_33_real_80)
chi_33_real_spline_89 = interp2d(gam, omrat, chi_33_real_89)

chi_33_imag_spline_01 = interp2d(gam, omrat, chi_33_imag_01)
chi_33_imag_spline_10 = interp2d(gam, omrat, chi_33_imag_10)
chi_33_imag_spline_20 = interp2d(gam, omrat, chi_33_imag_20)
chi_33_imag_spline_30 = interp2d(gam, omrat, chi_33_imag_30)
chi_33_imag_spline_40 = interp2d(gam, omrat, chi_33_imag_40)
chi_33_imag_spline_50 = interp2d(gam, omrat, chi_33_imag_50)
chi_33_imag_spline_60 = interp2d(gam, omrat, chi_33_imag_60)
chi_33_imag_spline_70 = interp2d(gam, omrat, chi_33_imag_70)
chi_33_imag_spline_80 = interp2d(gam, omrat, chi_33_imag_80)
chi_33_imag_spline_89 = interp2d(gam, omrat, chi_33_imag_89)

#define constants#

epsilon0  = 1./(4. * np.pi)
e         = 4.80320680e-10
m         = 9.1093826e-28
c         = 2.99792458e10
epsilon   = -1.

#  Each spline above is at fixed value of observer angle theta.  The below
#functions do linear interpolation to be valid for all theta in [0, pi/2].

def chi_11_real_spline(gamma, omratio, angle):
    angle_deg = angle / np.pi * 180.
    if(angle_deg <= 10):
        y0 = chi_11_real_spline_01(gamma, omratio)
        y1 = chi_11_real_spline_10(gamma, omratio) 
        x0 = 01.
        x1 = 10.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(10. < angle_deg <= 20.):
        y0 = chi_11_real_spline_10(gamma, omratio)
        y1 = chi_11_real_spline_20(gamma, omratio) 
        x0 = 10.
        x1 = 20.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(20. < angle_deg <= 30.):
        y0 = chi_11_real_spline_20(gamma, omratio)
        y1 = chi_11_real_spline_30(gamma, omratio) 
        x0 = 20.
        x1 = 30.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(30. < angle_deg <= 40.):
        y0 = chi_11_real_spline_30(gamma, omratio)
        y1 = chi_11_real_spline_40(gamma, omratio) 
        x0 = 30.
        x1 = 40.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(40. < angle_deg <= 50.):
        y0 = chi_11_real_spline_40(gamma, omratio)
        y1 = chi_11_real_spline_50(gamma, omratio) 
        x0 = 40.
        x1 = 50.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(50. < angle_deg <= 60.):
        y0 = chi_11_real_spline_50(gamma, omratio)
        y1 = chi_11_real_spline_60(gamma, omratio) 
        x0 = 50.
        x1 = 60.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(60. < angle_deg <= 70.):
        y0 = chi_11_real_spline_60(gamma, omratio)
        y1 = chi_11_real_spline_70(gamma, omratio) 
        x0 = 60.
        x1 = 70.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(70. < angle_deg <= 80.):
        y0 = chi_11_real_spline_70(gamma, omratio)
        y1 = chi_11_real_spline_80(gamma, omratio) 
        x0 = 70.
        x1 = 80.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(80. < angle_deg <= 89.):
        y0 = chi_11_real_spline_80(gamma, omratio)
        y1 = chi_11_real_spline_89(gamma, omratio) 
        x0 = 80.
        x1 = 89.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    else:
        print 'ERROR ENCOUNTERED IN chi_11_real_spline'
        return 0.
    
def chi_11_imag_spline(gamma, omratio, angle):
    angle_deg = angle / np.pi * 180.
    if(angle_deg <= 10):
        y0 = chi_11_imag_spline_01(gamma, omratio)
        y1 = chi_11_imag_spline_10(gamma, omratio) 
        x0 = 01.
        x1 = 10.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(10. < angle_deg <= 20.):
        y0 = chi_11_imag_spline_10(gamma, omratio)
        y1 = chi_11_imag_spline_20(gamma, omratio) 
        x0 = 10.
        x1 = 20.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(20. < angle_deg <= 30.):
        y0 = chi_11_imag_spline_20(gamma, omratio)
        y1 = chi_11_imag_spline_30(gamma, omratio) 
        x0 = 20.
        x1 = 30.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(30. < angle_deg <= 40.):
        y0 = chi_11_imag_spline_30(gamma, omratio)
        y1 = chi_11_imag_spline_40(gamma, omratio) 
        x0 = 30.
        x1 = 40.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(40. < angle_deg <= 50.):
        y0 = chi_11_imag_spline_40(gamma, omratio)
        y1 = chi_11_imag_spline_50(gamma, omratio) 
        x0 = 40.
        x1 = 50.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(50. < angle_deg <= 60.):
        y0 = chi_11_imag_spline_50(gamma, omratio)
        y1 = chi_11_imag_spline_60(gamma, omratio) 
        x0 = 50.
        x1 = 60.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(60. < angle_deg <= 70.):
        y0 = chi_11_imag_spline_60(gamma, omratio)
        y1 = chi_11_imag_spline_70(gamma, omratio) 
        x0 = 60.
        x1 = 70.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(70. < angle_deg <= 80.):
        y0 = chi_11_imag_spline_70(gamma, omratio)
        y1 = chi_11_imag_spline_80(gamma, omratio) 
        x0 = 70.
        x1 = 80.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(80. < angle_deg <= 89.):
        y0 = chi_11_imag_spline_80(gamma, omratio)
        y1 = chi_11_imag_spline_89(gamma, omratio) 
        x0 = 80.
        x1 = 89.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    else:
        print 'ERROR ENCOUNTERED IN chi_11_imag_spline'
        return 0.

def chi_12_real_spline(gamma, omratio, angle):
    angle_deg = angle / np.pi * 180.
    if(angle_deg <= 10):
        y0 = chi_12_real_spline_01(gamma, omratio)
        y1 = chi_12_real_spline_10(gamma, omratio) 
        x0 = 01.
        x1 = 10.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(10. < angle_deg <= 20.):
        y0 = chi_12_real_spline_10(gamma, omratio)
        y1 = chi_12_real_spline_20(gamma, omratio) 
        x0 = 10.
        x1 = 20.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(20. < angle_deg <= 30.):
        y0 = chi_12_real_spline_20(gamma, omratio)
        y1 = chi_12_real_spline_30(gamma, omratio) 
        x0 = 20.
        x1 = 30.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(30. < angle_deg <= 40.):
        y0 = chi_12_real_spline_30(gamma, omratio)
        y1 = chi_12_real_spline_40(gamma, omratio) 
        x0 = 30.
        x1 = 40.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(40. < angle_deg <= 50.):
        y0 = chi_12_real_spline_40(gamma, omratio)
        y1 = chi_12_real_spline_50(gamma, omratio) 
        x0 = 40.
        x1 = 50.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(50. < angle_deg <= 60.):
        y0 = chi_12_real_spline_50(gamma, omratio)
        y1 = chi_12_real_spline_60(gamma, omratio) 
        x0 = 50.
        x1 = 60.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(60. < angle_deg <= 70.):
        y0 = chi_12_real_spline_60(gamma, omratio)
        y1 = chi_12_real_spline_70(gamma, omratio) 
        x0 = 60.
        x1 = 70.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(70. < angle_deg <= 80.):
        y0 = chi_12_real_spline_70(gamma, omratio)
        y1 = chi_12_real_spline_80(gamma, omratio) 
        x0 = 70.
        x1 = 80.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(80. < angle_deg <= 89.):
        y0 = chi_12_real_spline_80(gamma, omratio)
        y1 = chi_12_real_spline_89(gamma, omratio) 
        x0 = 80.
        x1 = 89.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    else:
        print 'ERROR ENCOUNTERED IN chi_12_real_spline'
        return 0.
    
def chi_12_imag_spline(gamma, omratio, angle):
    angle_deg = angle / np.pi * 180.
    if(angle_deg <= 10):
        y0 = chi_12_imag_spline_01(gamma, omratio)
        y1 = chi_12_imag_spline_10(gamma, omratio) 
        x0 = 01.
        x1 = 10.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(10. < angle_deg <= 20.):
        y0 = chi_12_imag_spline_10(gamma, omratio)
        y1 = chi_12_imag_spline_20(gamma, omratio) 
        x0 = 10.
        x1 = 20.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(20. < angle_deg <= 30.):
        y0 = chi_12_imag_spline_20(gamma, omratio)
        y1 = chi_12_imag_spline_30(gamma, omratio) 
        x0 = 20.
        x1 = 30.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(30. < angle_deg <= 40.):
        y0 = chi_12_imag_spline_30(gamma, omratio)
        y1 = chi_12_imag_spline_40(gamma, omratio) 
        x0 = 30.
        x1 = 40.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(40. < angle_deg <= 50.):
        y0 = chi_12_imag_spline_40(gamma, omratio)
        y1 = chi_12_imag_spline_50(gamma, omratio) 
        x0 = 40.
        x1 = 50.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(50. < angle_deg <= 60.):
        y0 = chi_12_imag_spline_50(gamma, omratio)
        y1 = chi_12_imag_spline_60(gamma, omratio) 
        x0 = 50.
        x1 = 60.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(60. < angle_deg <= 70.):
        y0 = chi_12_imag_spline_60(gamma, omratio)
        y1 = chi_12_imag_spline_70(gamma, omratio) 
        x0 = 60.
        x1 = 70.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(70. < angle_deg <= 80.):
        y0 = chi_12_imag_spline_70(gamma, omratio)
        y1 = chi_12_imag_spline_80(gamma, omratio) 
        x0 = 70.
        x1 = 80.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(80. < angle_deg <= 89.):
        y0 = chi_12_imag_spline_80(gamma, omratio)
        y1 = chi_12_imag_spline_89(gamma, omratio) 
        x0 = 80.
        x1 = 89.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    else:
        print 'ERROR ENCOUNTERED IN chi_12_imag_spline'
        return 0.

def chi_13_real_spline(gamma, omratio, angle):
    angle_deg = angle / np.pi * 180.
    if(angle_deg <= 10):
        y0 = chi_13_real_spline_01(gamma, omratio)
        y1 = chi_13_real_spline_10(gamma, omratio) 
        x0 = 01.
        x1 = 10.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(10. < angle_deg <= 20.):
        y0 = chi_13_real_spline_10(gamma, omratio)
        y1 = chi_13_real_spline_20(gamma, omratio) 
        x0 = 10.
        x1 = 20.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(20. < angle_deg <= 30.):
        y0 = chi_13_real_spline_20(gamma, omratio)
        y1 = chi_13_real_spline_30(gamma, omratio) 
        x0 = 20.
        x1 = 30.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(30. < angle_deg <= 40.):
        y0 = chi_13_real_spline_30(gamma, omratio)
        y1 = chi_13_real_spline_40(gamma, omratio) 
        x0 = 30.
        x1 = 40.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(40. < angle_deg <= 50.):
        y0 = chi_13_real_spline_40(gamma, omratio)
        y1 = chi_13_real_spline_50(gamma, omratio) 
        x0 = 40.
        x1 = 50.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(50. < angle_deg <= 60.):
        y0 = chi_13_real_spline_50(gamma, omratio)
        y1 = chi_13_real_spline_60(gamma, omratio) 
        x0 = 50.
        x1 = 60.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(60. < angle_deg <= 70.):
        y0 = chi_13_real_spline_60(gamma, omratio)
        y1 = chi_13_real_spline_70(gamma, omratio) 
        x0 = 60.
        x1 = 70.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(70. < angle_deg <= 80.):
        y0 = chi_13_real_spline_70(gamma, omratio)
        y1 = chi_13_real_spline_80(gamma, omratio) 
        x0 = 70.
        x1 = 80.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(80. < angle_deg <= 89.):
        y0 = chi_13_real_spline_80(gamma, omratio)
        y1 = chi_13_real_spline_89(gamma, omratio) 
        x0 = 80.
        x1 = 89.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    else:
        print 'ERROR ENCOUNTERED IN chi_13_real_spline'
        return 0.
    
def chi_13_imag_spline(gamma, omratio, angle):
    angle_deg = angle / np.pi * 180.
    if(angle_deg <= 10):
        y0 = chi_13_imag_spline_01(gamma, omratio)
        y1 = chi_13_imag_spline_10(gamma, omratio) 
        x0 = 01.
        x1 = 10.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(10. < angle_deg <= 20.):
        y0 = chi_13_imag_spline_10(gamma, omratio)
        y1 = chi_13_imag_spline_20(gamma, omratio) 
        x0 = 10.
        x1 = 20.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(20. < angle_deg <= 30.):
        y0 = chi_13_imag_spline_20(gamma, omratio)
        y1 = chi_13_imag_spline_30(gamma, omratio) 
        x0 = 20.
        x1 = 30.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(30. < angle_deg <= 40.):
        y0 = chi_13_imag_spline_30(gamma, omratio)
        y1 = chi_13_imag_spline_40(gamma, omratio) 
        x0 = 30.
        x1 = 40.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(40. < angle_deg <= 50.):
        y0 = chi_13_imag_spline_40(gamma, omratio)
        y1 = chi_13_imag_spline_50(gamma, omratio) 
        x0 = 40.
        x1 = 50.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(50. < angle_deg <= 60.):
        y0 = chi_13_imag_spline_50(gamma, omratio)
        y1 = chi_13_imag_spline_60(gamma, omratio) 
        x0 = 50.
        x1 = 60.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(60. < angle_deg <= 70.):
        y0 = chi_13_imag_spline_60(gamma, omratio)
        y1 = chi_13_imag_spline_70(gamma, omratio) 
        x0 = 60.
        x1 = 70.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(70. < angle_deg <= 80.):
        y0 = chi_13_imag_spline_70(gamma, omratio)
        y1 = chi_13_imag_spline_80(gamma, omratio) 
        x0 = 70.
        x1 = 80.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(80. < angle_deg <= 89.):
        y0 = chi_13_imag_spline_80(gamma, omratio)
        y1 = chi_13_imag_spline_89(gamma, omratio) 
        x0 = 80.
        x1 = 89.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    else:
        print 'ERROR ENCOUNTERED IN chi_13_imag_spline'
        return 0.

def chi_22_real_spline(gamma, omratio, angle):
    angle_deg = angle / np.pi * 180.
    if(angle_deg <= 10):
        y0 = chi_22_real_spline_01(gamma, omratio)
        y1 = chi_22_real_spline_10(gamma, omratio) 
        x0 = 01.
        x1 = 10.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(10. < angle_deg <= 20.):
        y0 = chi_22_real_spline_10(gamma, omratio)
        y1 = chi_22_real_spline_20(gamma, omratio) 
        x0 = 10.
        x1 = 20.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(20. < angle_deg <= 30.):
        y0 = chi_22_real_spline_20(gamma, omratio)
        y1 = chi_22_real_spline_30(gamma, omratio) 
        x0 = 20.
        x1 = 30.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(30. < angle_deg <= 40.):
        y0 = chi_22_real_spline_30(gamma, omratio)
        y1 = chi_22_real_spline_40(gamma, omratio) 
        x0 = 30.
        x1 = 40.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(40. < angle_deg <= 50.):
        y0 = chi_22_real_spline_40(gamma, omratio)
        y1 = chi_22_real_spline_50(gamma, omratio) 
        x0 = 40.
        x1 = 50.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(50. < angle_deg <= 60.):
        y0 = chi_22_real_spline_50(gamma, omratio)
        y1 = chi_22_real_spline_60(gamma, omratio) 
        x0 = 50.
        x1 = 60.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(60. < angle_deg <= 70.):
        y0 = chi_22_real_spline_60(gamma, omratio)
        y1 = chi_22_real_spline_70(gamma, omratio) 
        x0 = 60.
        x1 = 70.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(70. < angle_deg <= 80.):
        y0 = chi_22_real_spline_70(gamma, omratio)
        y1 = chi_22_real_spline_80(gamma, omratio) 
        x0 = 70.
        x1 = 80.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(80. < angle_deg <= 89.):
        y0 = chi_22_real_spline_80(gamma, omratio)
        y1 = chi_22_real_spline_89(gamma, omratio) 
        x0 = 80.
        x1 = 89.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    else:
        print 'ERROR ENCOUNTERED IN chi_22_real_spline'
        return 0.
    
def chi_22_imag_spline(gamma, omratio, angle):
    angle_deg = angle / np.pi * 180.
    if(angle_deg <= 10):
        y0 = chi_22_imag_spline_01(gamma, omratio)
        y1 = chi_22_imag_spline_10(gamma, omratio) 
        x0 = 01.
        x1 = 10.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(10. < angle_deg <= 20.):
        y0 = chi_22_imag_spline_10(gamma, omratio)
        y1 = chi_22_imag_spline_20(gamma, omratio) 
        x0 = 10.
        x1 = 20.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(20. < angle_deg <= 30.):
        y0 = chi_22_imag_spline_20(gamma, omratio)
        y1 = chi_22_imag_spline_30(gamma, omratio) 
        x0 = 20.
        x1 = 30.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(30. < angle_deg <= 40.):
        y0 = chi_22_imag_spline_30(gamma, omratio)
        y1 = chi_22_imag_spline_40(gamma, omratio) 
        x0 = 30.
        x1 = 40.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(40. < angle_deg <= 50.):
        y0 = chi_22_imag_spline_40(gamma, omratio)
        y1 = chi_22_imag_spline_50(gamma, omratio) 
        x0 = 40.
        x1 = 50.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(50. < angle_deg <= 60.):
        y0 = chi_22_imag_spline_50(gamma, omratio)
        y1 = chi_22_imag_spline_60(gamma, omratio) 
        x0 = 50.
        x1 = 60.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(60. < angle_deg <= 70.):
        y0 = chi_22_imag_spline_60(gamma, omratio)
        y1 = chi_22_imag_spline_70(gamma, omratio) 
        x0 = 60.
        x1 = 70.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(70. < angle_deg <= 80.):
        y0 = chi_22_imag_spline_70(gamma, omratio)
        y1 = chi_22_imag_spline_80(gamma, omratio) 
        x0 = 70.
        x1 = 80.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(80. < angle_deg <= 89.):
        y0 = chi_22_imag_spline_80(gamma, omratio)
        y1 = chi_22_imag_spline_89(gamma, omratio) 
        x0 = 80.
        x1 = 89.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    else:
        print 'ERROR ENCOUNTERED IN chi_22_imag_spline'
        return 0.

def chi_32_real_spline(gamma, omratio, angle):
    angle_deg = angle / np.pi * 180.
    if(angle_deg <= 10):
        y0 = chi_32_real_spline_01(gamma, omratio)
        y1 = chi_32_real_spline_10(gamma, omratio) 
        x0 = 01.
        x1 = 10.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(10. < angle_deg <= 20.):
        y0 = chi_32_real_spline_10(gamma, omratio)
        y1 = chi_32_real_spline_20(gamma, omratio) 
        x0 = 10.
        x1 = 20.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(20. < angle_deg <= 30.):
        y0 = chi_32_real_spline_20(gamma, omratio)
        y1 = chi_32_real_spline_30(gamma, omratio) 
        x0 = 20.
        x1 = 30.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(30. < angle_deg <= 40.):
        y0 = chi_32_real_spline_30(gamma, omratio)
        y1 = chi_32_real_spline_40(gamma, omratio) 
        x0 = 30.
        x1 = 40.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(40. < angle_deg <= 50.):
        y0 = chi_32_real_spline_40(gamma, omratio)
        y1 = chi_32_real_spline_50(gamma, omratio) 
        x0 = 40.
        x1 = 50.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(50. < angle_deg <= 60.):
        y0 = chi_32_real_spline_50(gamma, omratio)
        y1 = chi_32_real_spline_60(gamma, omratio) 
        x0 = 50.
        x1 = 60.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(60. < angle_deg <= 70.):
        y0 = chi_32_real_spline_60(gamma, omratio)
        y1 = chi_32_real_spline_70(gamma, omratio) 
        x0 = 60.
        x1 = 70.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(70. < angle_deg <= 80.):
        y0 = chi_32_real_spline_70(gamma, omratio)
        y1 = chi_32_real_spline_80(gamma, omratio) 
        x0 = 70.
        x1 = 80.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(80. < angle_deg <= 89.):
        y0 = chi_32_real_spline_80(gamma, omratio)
        y1 = chi_32_real_spline_89(gamma, omratio) 
        x0 = 80.
        x1 = 89.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    else:
        print 'ERROR ENCOUNTERED IN chi_32_real_spline'
        return 0.
    
def chi_32_imag_spline(gamma, omratio, angle):
    angle_deg = angle / np.pi * 180.
    if(angle_deg <= 10):
        y0 = chi_32_imag_spline_01(gamma, omratio)
        y1 = chi_32_imag_spline_10(gamma, omratio) 
        x0 = 01.
        x1 = 10.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(10. < angle_deg <= 20.):
        y0 = chi_32_imag_spline_10(gamma, omratio)
        y1 = chi_32_imag_spline_20(gamma, omratio) 
        x0 = 10.
        x1 = 20.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(20. < angle_deg <= 30.):
        y0 = chi_32_imag_spline_20(gamma, omratio)
        y1 = chi_32_imag_spline_30(gamma, omratio) 
        x0 = 20.
        x1 = 30.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(30. < angle_deg <= 40.):
        y0 = chi_32_imag_spline_30(gamma, omratio)
        y1 = chi_32_imag_spline_40(gamma, omratio) 
        x0 = 30.
        x1 = 40.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(40. < angle_deg <= 50.):
        y0 = chi_32_imag_spline_40(gamma, omratio)
        y1 = chi_32_imag_spline_50(gamma, omratio) 
        x0 = 40.
        x1 = 50.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(50. < angle_deg <= 60.):
        y0 = chi_32_imag_spline_50(gamma, omratio)
        y1 = chi_32_imag_spline_60(gamma, omratio) 
        x0 = 50.
        x1 = 60.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(60. < angle_deg <= 70.):
        y0 = chi_32_imag_spline_60(gamma, omratio)
        y1 = chi_32_imag_spline_70(gamma, omratio) 
        x0 = 60.
        x1 = 70.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(70. < angle_deg <= 80.):
        y0 = chi_32_imag_spline_70(gamma, omratio)
        y1 = chi_32_imag_spline_80(gamma, omratio) 
        x0 = 70.
        x1 = 80.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(80. < angle_deg <= 89.):
        y0 = chi_32_imag_spline_80(gamma, omratio)
        y1 = chi_32_imag_spline_89(gamma, omratio) 
        x0 = 80.
        x1 = 89.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    else:
        print 'ERROR ENCOUNTERED IN chi_32_imag_spline'
        return 0.


def chi_33_real_spline(gamma, omratio, angle):
    angle_deg = angle / np.pi * 180.
    if(angle_deg <= 10):
        y0 = chi_33_real_spline_01(gamma, omratio)
        y1 = chi_33_real_spline_10(gamma, omratio) 
        x0 = 01.
        x1 = 10.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(10. < angle_deg <= 20.):
        y0 = chi_33_real_spline_10(gamma, omratio)
        y1 = chi_33_real_spline_20(gamma, omratio) 
        x0 = 10.
        x1 = 20.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(20. < angle_deg <= 30.):
        y0 = chi_33_real_spline_20(gamma, omratio)
        y1 = chi_33_real_spline_30(gamma, omratio) 
        x0 = 20.
        x1 = 30.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(30. < angle_deg <= 40.):
        y0 = chi_33_real_spline_30(gamma, omratio)
        y1 = chi_33_real_spline_40(gamma, omratio) 
        x0 = 30.
        x1 = 40.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(40. < angle_deg <= 50.):
        y0 = chi_33_real_spline_40(gamma, omratio)
        y1 = chi_33_real_spline_50(gamma, omratio) 
        x0 = 40.
        x1 = 50.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(50. < angle_deg <= 60.):
        y0 = chi_33_real_spline_50(gamma, omratio)
        y1 = chi_33_real_spline_60(gamma, omratio) 
        x0 = 50.
        x1 = 60.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(60. < angle_deg <= 70.):
        y0 = chi_33_real_spline_60(gamma, omratio)
        y1 = chi_33_real_spline_70(gamma, omratio) 
        x0 = 60.
        x1 = 70.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(70. < angle_deg <= 80.):
        y0 = chi_33_real_spline_70(gamma, omratio)
        y1 = chi_33_real_spline_80(gamma, omratio) 
        x0 = 70.
        x1 = 80.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(80. < angle_deg <= 89.):
        y0 = chi_33_real_spline_80(gamma, omratio)
        y1 = chi_33_real_spline_89(gamma, omratio) 
        x0 = 80.
        x1 = 89.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    else:
        print 'ERROR ENCOUNTERED IN chi_33_real_spline'
        return 0.
    
def chi_33_imag_spline(gamma, omratio, angle):
    angle_deg = angle / np.pi * 180.
    if(angle_deg <= 10):
        y0 = chi_33_imag_spline_01(gamma, omratio)
        y1 = chi_33_imag_spline_10(gamma, omratio) 
        x0 = 01.
        x1 = 10.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(10. < angle_deg <= 20.):
        y0 = chi_33_imag_spline_10(gamma, omratio)
        y1 = chi_33_imag_spline_20(gamma, omratio) 
        x0 = 10.
        x1 = 20.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(20. < angle_deg <= 30.):
        y0 = chi_33_imag_spline_20(gamma, omratio)
        y1 = chi_33_imag_spline_30(gamma, omratio) 
        x0 = 20.
        x1 = 30.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(30. < angle_deg <= 40.):
        y0 = chi_33_imag_spline_30(gamma, omratio)
        y1 = chi_33_imag_spline_40(gamma, omratio) 
        x0 = 30.
        x1 = 40.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(40. < angle_deg <= 50.):
        y0 = chi_33_imag_spline_40(gamma, omratio)
        y1 = chi_33_imag_spline_50(gamma, omratio) 
        x0 = 40.
        x1 = 50.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(50. < angle_deg <= 60.):
        y0 = chi_33_imag_spline_50(gamma, omratio)
        y1 = chi_33_imag_spline_60(gamma, omratio) 
        x0 = 50.
        x1 = 60.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(60. < angle_deg <= 70.):
        y0 = chi_33_imag_spline_60(gamma, omratio)
        y1 = chi_33_imag_spline_70(gamma, omratio) 
        x0 = 60.
        x1 = 70.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(70. < angle_deg <= 80.):
        y0 = chi_33_imag_spline_70(gamma, omratio)
        y1 = chi_33_imag_spline_80(gamma, omratio) 
        x0 = 70.
        x1 = 80.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    elif(80. < angle_deg <= 89.):
        y0 = chi_33_imag_spline_80(gamma, omratio)
        y1 = chi_33_imag_spline_89(gamma, omratio) 
        x0 = 80.
        x1 = 89.
        x = angle_deg
        y   = y0 * (1 - (x - x0) / (x1 - x0) ) + y1 * ( (x - x0) / (x1 - x0) )
        return y
    else:
        print 'ERROR ENCOUNTERED IN chi_33_imag_spline'
        return 0.

#function to choose which component's spline to use.

def spline_selector(real, component, gamma, omratio, angle):    
    if(real == 1):
        if(component == 11):
            spline_term = chi_11_real_spline(gamma, omratio, angle)
        elif(component == 12):
            spline_term = chi_12_real_spline(gamma, omratio, angle)
        elif(component == 13):
            spline_term = chi_13_real_spline(gamma, omratio, angle)
        elif(component == 22):
            spline_term = chi_22_real_spline(gamma, omratio, angle)
        elif(component == 32):
            spline_term = chi_32_real_spline(gamma, omratio, angle)
        elif(component == 33):
            spline_term = chi_33_real_spline(gamma, omratio, angle)
    else:
        if(component == 11):
            spline_term = chi_11_imag_spline(gamma, omratio, angle)
        elif(component == 12):
            spline_term = chi_12_imag_spline(gamma, omratio, angle)
        elif(component == 13):
            spline_term = chi_13_imag_spline(gamma, omratio, angle)
        elif(component == 22):
            spline_term = chi_22_imag_spline(gamma, omratio, angle)
        elif(component == 32):
            spline_term = chi_32_imag_spline(gamma, omratio, angle)
        elif(component == 33):
            spline_term = chi_33_imag_spline(gamma, omratio, angle)
            
    return spline_term

# gamma integrand for each component.  The term that changes between components
#is chosen via the function spline_selector, above.

def kappa_unnormalized(gamma, kappa, kappa_width):
    beta = np.sqrt(1. - 1./gamma**2.)
    body = (1. + (gamma - 1.)/(kappa * kappa_width))**(-1. - kappa)
    d3p  = 4. * np.pi * gamma*gamma * beta
    ans  = body * d3p
    return ans

def kappa_normalized(gamma, kappa, kappa_width):
    norm_constant = 1./quad(lambda gamma: kappa_unnormalized(gamma, kappa, kappa_width), 1., np.inf)[0]
    body = (1. + (gamma - 1.)/(kappa * kappa_width))**(-2. - kappa)*(-1. - kappa) / (kappa_width * kappa)
    ans = norm_constant * body
    return ans

def chi_ij_integrand(gamma, omratio, angle, dist_func, real, component, theta_e, power_law_p, gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width):
    beta = np.sqrt(1. - 1./gamma**2.)

    if(dist_func == 0):
        dist = -np.exp(-gamma/theta_e) / (4. * np.pi * theta_e**2. * special.kn(2, 1./theta_e)) * gamma**3. * beta**3.
    elif(dist_func == 1):
        dist = -( (power_law_p - 1.) * (-1 + 2. * gamma**2. + power_law_p * (gamma**2. - 1.))
                / (4. * np.pi * (gamma_min**(-1. - power_law_p)
                - gamma_max**(-1. - power_law_p)) * 1. * 1.)
                * gamma**(-3. - power_law_p) ) * gamma
    elif(dist_func == 2):
       dist = kappa_normalized(gamma, kappa, kappa_width) * gamma**3. * beta**3.

    gam_term = dist

    #NOTE: due to a bug, all datafiles contain a *-1 error; we correct it below.
    #the bug is corrected throughout the code.
    spline_term = -spline_selector(real, component, gamma, omratio, angle)

    ans = gam_term * spline_term
    return ans

def chi_ij(nu,
           magnetic_field,
           electron_density,
           observer_angle,
           distribution,
           real_part,
           theta_e,
           power_law_p,
           gamma_min,
           gamma_max,
           gamma_cutoff,
           kappa,
           kappa_width,
           component):

    """returns chi_ij(nu, magnetic_field, electron_density, observer_angle,
                      distribution, real_part, theta_e, power_law_p, gamma_min
                      gamma_max, gamma_cutoff, kappa, kappa_width, component) """
    
    omega = 2. * np.pi * nu
    omega_p = np.sqrt(electron_density * e*e / (m * epsilon0))
    omega_c = e * magnetic_field / (m * c)

    prefactor = 2. * np.pi * omega_p*omega_p / (omega * omega)
    ans = quad(lambda gamma: np.vectorize(chi_ij_integrand)(gamma, omega / omega_c, 
                                                            observer_angle, distribution, real_part,
							    component, theta_e, power_law_p, gamma_min,
							    gamma_max, gamma_cutoff, kappa, kappa_width), 
               1., 1e4, epsabs=0, epsrel=1e-7, limit=5000)[0] * prefactor #NOTE: limits are bounds of interpolated data in gamma
    
    return ans
