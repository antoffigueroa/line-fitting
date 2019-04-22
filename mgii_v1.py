import numpy as np
import matplotlib.pyplot as plt
from mpdaf.obj import Cube
from functions_mgii import *
from mpdaf.obj import Spectrum
from scipy.optimize import bisect, leastsq
import scipy.stats
import aplpy
from astropy.io import fits

fits_file = '/Users/antonia/taller-de-investigacion/SGASJ1226/muse/r4_cube_mgii_mock.fits'

cube = Cube(fits_file)

wave0 = 2796.352
z = 0.77113
wratio = 1.0025672375
wavep = wave0*(1+z)
wavep_2 = 2803.531*(1+z)

wavep_oii = 3727.092*(1+z)
wratio_oii = 3729.875/wavep_oii

arr = cube.data

ncoords = arr.count()
coord = np.zeros((ncoords, 2), dtype=np.int)

aux = 0

# j es x, i es y
# very bright spaxel, used for experimentation
# for i in range(43, 44):
#    for j in range(45,46):
#        coord[aux] = i,j
#        aux = aux + 1
# ncoords=aux

# select a rectangular region to iterate

#for i in range(32, 52):
#    for j in range(32, 48):
#        coord[aux] = i, j
#        aux = aux + 1
#
#ncoords = aux

# apply mask

north = np.load('/Users/antonia/taller-de-investigacion/SGASJ1226/other-scripts/mask_N.npy')
galaxy = np.load('/Users/antonia/taller-de-investigacion/SGASJ1226/other-scripts/mask_G1.npy')
south = np.load('/Users/antonia/taller-de-investigacion/SGASJ1226/other-scripts/mask_S.npy')

mask_1 = np.logical_or(north, galaxy)
mask_final = np.logical_or(mask_1, south)

mask_shape = mask_final.shape

for i in range(mask_shape[0]):
    for j in range(mask_shape[1]):
        if mask_final[i,j]:
            coord[aux] = i, j
            aux = aux + 1
ncoords = aux

contador = 0

EW = np.array([])
x_succes = np.array([])
y_succes = np.array([])
vel = np.array([])

EW_oii = np.array([])
x_succes_oii = np.array([])
y_succes_oii = np.array([])
vel_oii = np.array([])
SN_array = np.array([])

plt.close()
fig = plt.figure(figsize = (6.4,10))

for k in range(ncoords):
    y = coord[k, 1]  # q[i]
    x = coord[k, 0]  # p[j]
    spe = cube[:, x, y]
    var = spe.var
    sigma = np.sqrt(var)
    spe = normalize(spe, wavep, z, wratio=wratio, show=False)
    # plt.close()
    # plt.figure()
    # plt.plot(spe.wave.coord(),spe.data)
    # plt.show()
    A_1_fit, err_A_1, A_2_fit, err_A_2, mu_fit, err_mu, flux_1_fit, flux_2_fit = fit_doublet(spe, z)
    flux_1, err_flux_1 = flux(A_1_fit, err_A_1)
    flux_2, err_flux_2 = flux(A_2_fit, err_A_2)
    EW_1, err_EW_1 = eq_width(flux_1, err_flux_1, z)
    EW_2, err_EW_2 = eq_width(flux_2, err_flux_2, z)
    SN = EW_1/err_EW_1
    if SN >= 2. and A_1_fit > 0 and A_2_fit > 0 and EW_1 < 5:
        print "succesful mgii fit at: "+str(x)+", "+str(y)+' with EW: '+str(EW_1)
        contador = contador + 1
        velocidad = velocity(mu_fit, z)
        EW = np.append(EW, EW_1)
        x_succes = np.append(x_succes, x)
        y_succes = np.append(y_succes, y)
        vel = np.append(vel, velocidad)
        SN_array = np.append(SN_array, SN)
        best_fit = double_model((A_1_fit, A_2_fit, mu_fit), np.linspace(spe.wave.coord()[0], spe.wave.coord()[-1], num=100000))
        ax = fig.add_subplot(28,1, contador)
        ax.plot(np.linspace(spe.wave.coord()[0], spe.wave.coord()[-1], num=100000), best_fit)
        ax.plot(spe.wave.coord(), spe.data)
        ax.set_xlim(4920,4980)
    else:
        print ":("
    """
    spe = cube[:, x, y]
    spe = normalize(spe, wavep_oii, z, wratio=wratio_oii, show=True)
    A_1_fit_oii, err_A_1_oii, A_2_fit_oii, err_A_2_oii, mu_fit_oii, err_mu_oii, flux_1_fit_oii, flux_2_fit_oii = fit_doublet(spe, z, how='emi')
    flux_1_oii, err_flux_1_oii = flux(A_1_fit_oii, err_A_1_oii)
    flux_2_oii, err_flux_2_oii = flux(A_2_fit_oii, err_A_2_oii)
    EW_1_oii, err_EW_1_oii = eq_width(flux_1_oii, err_flux_1_oii, z)
    EW_2_oii, err_EW_2_oii = eq_width(flux_2_oii, err_flux_2_oii, z)
    if EW_1_oii/err_EW_1_oii >= 2.:
        print "succesful oii fit at: "+str(x)+", "+str(y)+' with EW: '+str(EW_1_oii)
        velocidad = velocity(mu_fit_oii, z, how = 'emi')
        EW_oii = np.append(EW_oii, EW_1_oii)
        x_succes_oii = np.append(x_succes_oii, x)
        y_succes_oii = np.append(y_succes_oii, y)
        vel_oii = np.append(vel_oii, velocidad)
        # best_fit = double_model((A_1_fit, A_2_fit, mu_fit), np.linspace(spe.wave.coord()[0], spe.wave.coord()[-1], num=100000))
        # plt.figure()
        # plt.plot(np.linspace(spe.wave.coord()[0], spe.wave.coord()[-1], num=100000), best_fit)
        # plt.plot(spe.wave.coord(), spe.data)
        # plt.show()
    else:
        print ":("
    """


best_fit = double_model((A_1_fit, A_2_fit, mu_fit), np.linspace(spe.wave.coord()[0], spe.wave.coord()[-1], num=100000))

EW_matrix=np.zeros((cube.data.shape[1],cube.data.shape[2]))
EW_matrix[EW_matrix==0] = np.nan

v_matrix=np.zeros((cube.data.shape[1],cube.data.shape[2]))
v_matrix[v_matrix==0] = np.nan

SN_matrix = np.zeros((cube.data.shape[1],cube.data.shape[2]))
SN_matrix[SN_matrix==0] = np.nan

for i in range(len(x_succes)):
    EW_matrix[int(x_succes[i]), int(y_succes[i])]=EW[i]
    v_matrix[int(x_succes[i]), int(y_succes[i])] = vel[i]
    SN_matrix[int(x_succes[i]), int(y_succes[i])] = SN_array[i]

"""
EW_matrix_oii = np.zeros((cube.data.shape[1],cube.data.shape[2]))
EW_matrix_oii[EW_matrix_oii==0] = np.nan

v_matrix_oii = np.zeros((cube.data.shape[1],cube.data.shape[2]))
v_matrix_oii[v_matrix_oii==0] = np.nan

for i in range(len(x_succes)):
	EW_matrix[int(x_succes[i]), int(y_succes[i])]=EW[i]
	v_matrix[int(x_succes[i]), int(y_succes[i])] = vel[i]

for i in range(len(x_succes_oii)):
	EW_matrix_oii[int(x_succes_oii[i]), int(y_succes_oii[i])]=EW_oii[i]
	v_matrix_oii[int(x_succes_oii[i]), int(y_succes_oii[i])] = vel_oii[i]
"""

image_data = fits.getdata(fits_file, ext=0)

plt.figure()
plt.contour(np.log(np.mean(image_data,axis=0)), levels=1)
plt.imshow(EW_matrix, cmap = 'jet')
plt.xlim(0,image_data.shape[2])
plt.ylim(image_data.shape[1],0)
plt.gca().invert_yaxis()
cbar= plt.colorbar()
cbar.set_label('Equivalent Width [$\AA$]')
plt.title('Mg II')
plt.legend().draggable()
plt.savefig('Figure 1.pdf')
#plt.show()

plt.figure()
plt.contour(np.log(np.mean(image_data,axis=0)), levels=1)
plt.imshow(v_matrix, cmap = 'jet')
plt.xlim(0,image_data.shape[2])
plt.ylim(image_data.shape[1],0)
plt.gca().invert_yaxis()
cbar= plt.colorbar()
cbar.set_label('Velocity [km/s]')
plt.title('Mg II')
plt.legend().draggable()
plt.savefig('Figure 2.pdf')
#plt.show()

plt.figure()
plt.contour(np.log(np.mean(image_data,axis=0)), levels=1)
plt.imshow(SN_matrix, cmap = 'jet')
plt.xlim(0,image_data.shape[2])
plt.ylim(image_data.shape[1],0)
plt.gca().invert_yaxis()
cbar= plt.colorbar()
cbar.set_label('SN')
plt.title('Mg II')
plt.legend().draggable()
plt.savefig('Figure 3.pdf')
#plt.show()

"""
plt.figure()
plt.contour(np.log(np.mean(image_data,axis=0)), levels=1)
plt.imshow(EW_matrix_oii)
plt.xlim(0,image_data.shape[2])
plt.ylim(image_data.shape[1],0)
plt.gca().invert_yaxis()
cbar= plt.colorbar()
cbar.set_label('Equivalent Width [$\AA$]')
plt.title('O II')
plt.legend().draggable()
plt.show()

plt.figure()
plt.contour(np.log(np.mean(image_data,axis=0)), levels=1)
plt.imshow(v_matrix_oii)
plt.xlim(0,image_data.shape[2])
plt.ylim(image_data.shape[1],0)
plt.gca().invert_yaxis()
cbar= plt.colorbar()
cbar.set_label('Velocity [km/s]')
plt.title('O II')
plt.legend().draggable()
plt.show()
"""

#plt.show()

# plt.figure()
# plt.plot(np.linspace(spe.wave.coord()[0], spe.wave.coord()[-1], num=100000), best_fit)
# plt.plot(spe.wave.coord(), spe.data)
# plt.show()
print contador
