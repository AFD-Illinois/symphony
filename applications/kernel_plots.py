import numpy as np
import pylab as pl
from matplotlib import colors, ticker, cm
import matplotlib.patches as patches

import sys
sys.path.insert(0, '../src/susceptibility_tensor')
sys.path.insert(0, '../build')

from susceptibility_interpolator import chi_ij, spline_selector
import susceptibility_tensor.susceptibilityPy as susp

# Set plot parameters to make beautiful plots
pl.rcParams['figure.figsize']  = 10, 10
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 20  
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex']     = True
pl.rcParams['axes.linewidth']  = 1.5
pl.rcParams['axes.titlesize']  = 'medium'
pl.rcParams['axes.labelsize']  = 'large'

pl.rcParams['xtick.major.size'] = 8     
pl.rcParams['xtick.minor.size'] = 4     
pl.rcParams['xtick.major.pad']  = 8     
pl.rcParams['xtick.minor.pad']  = 8     
pl.rcParams['xtick.color']      = 'k'     
pl.rcParams['xtick.labelsize']  = 'large'
pl.rcParams['xtick.direction']  = 'in'    

pl.rcParams['ytick.major.size'] = 8     
pl.rcParams['ytick.minor.size'] = 4     
pl.rcParams['ytick.major.pad']  = 8     
pl.rcParams['ytick.minor.pad']  = 8     
pl.rcParams['ytick.color']      = 'k'     
pl.rcParams['ytick.labelsize']  = 'large'
pl.rcParams['ytick.direction']  = 'in'

pl.rcParams['image.cmap'] = 'jet'

#define constants
epsilon0  = 1./(4. * np.pi)
e         = 4.80320680e-10
m         = 9.1093826e-28
c         = 2.99792458e10
epsilon   = -1.

def nu_c(magnetic_field):
    ans = e * magnetic_field / (2. * np.pi * m * c)
    return ans

real = 1
angle = np.pi/3.
#print spline_selector(real, component, gamma, omratio, angle)

Gamma   = np.logspace(0., 3., 50)
Nuratio = np.logspace(0., 3., 50)

G, N = np.meshgrid(Gamma, Nuratio)

figure, ax = pl.subplots(nrows=3, ncols=2, sharex=True, sharey=True)

plot1 = ax[0,0].contourf(G, N, (-np.vectorize(spline_selector)(real, 11, G, N, angle)), levels=np.logspace(-4., 0., 100), locator=ticker.LogLocator())
cbar1 = figure.colorbar(plot1, ax=ax[0,0], ticks=[1e-4, 1e-3, 1e-2, 1e-1, 1])
ax[0,0].annotate('$-\\chi_{11}$', xy=(0.7, 0.1), xycoords="axes fraction")

plot2 = ax[1,0].contourf(G, N, (np.vectorize(spline_selector)(real, 12, G, N, angle)), levels=np.logspace(-8., 1., 100), locator=ticker.LogLocator())
cbar2 = figure.colorbar(plot2, ax=ax[1,0], ticks=[1e-8, 1e-6, 1e-4, 1e-2, 1])
ax[1,0].annotate('$\\chi_{12}$', xy=(0.8, 0.1), xycoords="axes fraction")

plot3 = ax[2,0].contourf(G, N, (-np.vectorize(spline_selector)(real, 13, G, N, angle)), levels=np.logspace(-4., 0., 100), locator=ticker.LogLocator())
cbar3 = figure.colorbar(plot3, ax=ax[2,0], ticks=[1e-4, 1e-3, 1e-2, 1e-1, 1])
ax[2,0].annotate('$-\\chi_{13}$', xy=(0.7, 0.1), xycoords="axes fraction")

plot4 = ax[0,1].contourf(G, N, (-np.vectorize(spline_selector)(real, 22, G, N, angle)), levels=np.logspace(-4., 0., 100), locator=ticker.LogLocator())
cbar4 = figure.colorbar(plot4, ax=ax[0,1], ticks=[1e-4, 1e-3, 1e-2, 1e-1, 1])
ax[0,1].annotate('$-\\chi_{22}$', xy=(0.7, 0.1), xycoords="axes fraction")

arr32 = (np.vectorize(spline_selector)(real, 32, G, N, angle))
arr32[arr32 < 1e-8] = 1e-8
plot5 = ax[1,1].contourf(G, N, arr32, levels=np.logspace(-8., 1., 100), locator=ticker.LogLocator())
cbar5 = figure.colorbar(plot5, ax=ax[1,1], ticks=[1e-8, 1e-6, 1e-4, 1e-2, 1])
ax[1,1].annotate('$\\chi_{32}$', xy=(0.8, 0.1), xycoords="axes fraction")

plot6 = ax[2,1].contourf(G, N, (-np.vectorize(spline_selector)(real, 33, G, N, angle)), levels=np.logspace(-4., 0., 100), locator=ticker.LogLocator())
cbar6 = figure.colorbar(plot6, ax=ax[2,1], ticks=[1e-4, 1e-3, 1e-2, 1e-1, 1])
ax[2,1].annotate('$-\\chi_{33}$', xy=(0.7, 0.1), xycoords="axes fraction")

figure.add_subplot(111, frameon=False)
pl.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
pl.xlabel('$\\gamma$', fontsize='large')
#pl.ylabel('$\\omega/\\omega_c$', fontsize='large')

pl.setp(ax, xticks=[300, 600, 900], yticks=[300, 600, 900])

pl.text(-0.17, 0.5, '$\\omega/\\omega_c$', fontsize=25)

pl.tight_layout()
pl.show()

#----------make imaginary part plot-------------------------------------------#
real = 0

figure, ax = pl.subplots(nrows=3, ncols=2, sharex=True, sharey=True)

arr11 = (np.vectorize(spline_selector)(real, 11, G, N, angle))
arr11[arr11 < 1e-4] = 1e-4
plot1 = ax[0,0].contourf(G, N, arr11, levels=np.logspace(-4., 0., 100), locator=ticker.LogLocator())
cbar1 = figure.colorbar(plot1, ax=ax[0,0], ticks=[1e-4, 1e-3, 1e-2, 1e-1, 1])
ax[0,0].annotate('$-\\chi_{11}$', xy=(0.7, 0.1), xycoords="axes fraction")

arr12 = (np.vectorize(spline_selector)(real, 12, G, N, angle))
arr12[arr12 < 1e-8] = 1e-8
plot2 = ax[1,0].contourf(G, N, arr12, levels=np.logspace(-8., 0., 100), locator=ticker.LogLocator())
cbar2 = figure.colorbar(plot2, ax=ax[1,0], ticks=[1e-8, 1e-6, 1e-4, 1e-2, 1])
ax[1,0].annotate('$\\chi_{12}$', xy=(0.8, 0.1), xycoords="axes fraction")

arr13 = (np.vectorize(spline_selector)(real, 13, G, N, angle))
arr13[arr13 < 1e-4] = 1e-4
plot3 = ax[2,0].contourf(G, N, arr13, levels=np.logspace(-4., 0., 100), locator=ticker.LogLocator())
cbar3 = figure.colorbar(plot3, ax=ax[2,0], ticks=[1e-4, 1e-3, 1e-2, 1e-1, 1])
ax[2,0].annotate('$-\\chi_{13}$', xy=(0.7, 0.1), xycoords="axes fraction")

arr22 = (np.vectorize(spline_selector)(real, 22, G, N, angle))
arr22[arr22 < 1e-8] = 1e-8
plot4 = ax[0,1].contourf(G, N, arr22, levels=np.logspace(-8., 0., 100), locator=ticker.LogLocator())
cbar4 = figure.colorbar(plot4, ax=ax[0,1], ticks=[1e-8, 1e-6, 1e-4, 1e-2, 1])
ax[0,1].annotate('$-\\chi_{22}$', xy=(0.7, 0.1), xycoords="axes fraction")

arr32 = (np.vectorize(spline_selector)(real, 32, G, N, angle))
arr32[arr32 < 1e-8] = 1e-8
plot5 = ax[1,1].contourf(G, N, arr32, levels=np.logspace(-8., 0., 100), locator=ticker.LogLocator())
cbar5 = figure.colorbar(plot5, ax=ax[1,1], ticks=[1e-8, 1e-6, 1e-4, 1e-2, 1])
ax[1,1].annotate('$\\chi_{32}$', xy=(0.8, 0.1), xycoords="axes fraction")

arr33 = (np.vectorize(spline_selector)(real, 33, G, N, angle))
arr33[arr33 < 1e-4] = 1e-4
plot6 = ax[2,1].contourf(G, N, arr33, levels=np.logspace(-4., 0., 100), locator=ticker.LogLocator())
cbar6 = figure.colorbar(plot6, ax=ax[2,1], ticks=[1e-4, 1e-3, 1e-2, 1e-1, 1])
ax[2,1].annotate('$-\\chi_{33}$', xy=(0.7, 0.1), xycoords="axes fraction")

figure.add_subplot(111, frameon=False)
pl.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
pl.xlabel('$\\gamma$', fontsize='large')
#pl.ylabel('$\\omega/\\omega_c$', fontsize='large')

pl.setp(ax, xticks=[300, 600, 900], yticks=[300, 600, 900])

pl.text(-0.17, 0.5, '$\\omega/\\omega_c$', fontsize=25)

pl.tight_layout()
pl.show()
