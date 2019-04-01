import numpy as np
import matplotlib.pyplot as plt
from mpdaf.obj import Cube
from functions_mgii import *
from mpdaf.obj import Spectrum
from scipy.optimize import bisect, leastsq
import scipy.stats

cube=Cube('/Users/antonia/taller-de-investigacion/SGASJ1226/muse/r4_cube_mgii_mock.fits')

wave0=2796.352
z=0.77113
wratio=1.0025672375
wavep = wave0*(1+z)
wavep_2 = 2803.531*(1+z)

arr=cube.data

ncoords = arr.count()
coord = np.zeros((ncoords,2),dtype=np.int)

aux = 0

# j es x, i es y
for i in range(43, 44):
    for j in range(45,46):
        coord[aux] = i,j
        aux = aux + 1
ncoords=aux

contador = 0

for k in range(ncoords):
    y = coord[k,1] #q[i]
    x = coord[k,0] #p[j]
    spe = cube[:,x,y]
    var = spe.var
    sigma = np.sqrt(var)
    spe=normalize(spe, wavep, z, wratio=wratio, show=False)
    #plt.close()
    #plt.figure()
    #plt.plot(spe.wave.coord(),spe.data)
    #plt.show()
    A_1_fit, err_A_1, A_2_fit, err_A_2, mu_fit, err_mu, flux_1_fit, flux_2_fit = fit_doublet(spe, z)
    flux_1, err_flux_1 = flux(A_1_fit, err_A_1)
    flux_2, err_flux_2 = flux(A_2_fit, err_A_2)
    EW_1, err_EW_1 = eq_width(flux_1,err_flux_1,z)
    EW_2, err_EW_2 = eq_width(flux_2,err_flux_2,z)
    if EW_1/err_EW_1 >= 2.5:
    	print "succesful fit at: "+str(x)+", "+str(y)
        contador = contador +1
    else:
        print ":("

best_fit = double_model((A_1_fit, A_2_fit, mu_fit), np.linspace(spe.wave.coord()[0], spe.wave.coord()[-1], num=100000))

plt.figure()
plt.plot(np.linspace(spe.wave.coord()[0], spe.wave.coord()[-1], num=100000), best_fit)
plt.plot(spe.wave.coord(), spe.data)
plt.show()
#print contador
