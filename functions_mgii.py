import numpy as np
import matplotlib.pyplot as plt
from mpdaf.obj import Spectrum
from scipy.optimize import bisect, leastsq
import scipy.stats


def create_mask():
    """
    """


def cut_spectra(spe, w_ini, w_fin):
    spe_new = spe
    flux = spe_new.data
    wavelength = spe_new.wave.coord()
    p1 = np.int(spe_new.wave.pixel(w_ini))
    p2 = np.int(spe_new.wave.pixel(w_fin))
    cut_flux = flux[p1:p2]
    cut_wavelength = wavelength[p1:p2]
    return cut_flux, cut_wavelength


def normalize(spe, w_0, z, wratio=None, show=False):
    """
    spe: Spectrum object from mpdaf
    w_0: wavelength of the line you are working with
    wratio: only used if you know there is another line close to the line you
    are working with
    """
    if wratio is None:
        wratio_new = 1
    else:
        wratio_new = wratio*(1+z)
    r1 = np.int(spe.wave.pixel(w_0+10*wratio_new))
    r2 = np.int(spe.wave.pixel(w_0+10*wratio_new+20))
    l1 = np.int(spe.wave.pixel(w_0-10*wratio_new-20))
    l2 = np.int(spe.wave.pixel(w_0-10*wratio_new))
    flux = spe.data
    right_part = flux[r1:r2]
    left_part = flux[l1:l2]
    mean_right = np.mean(right_part)
    mean_left = np.mean(left_part)
    mean = (mean_right+mean_left)/2.0
    spe.data = spe.data/mean
    if show:
        plt.figure()
        plt.plot(spe.wave.coord(), spe.data)
        plt.plot(w_0+10*wratio, 1, '>')
        plt.plot(w_0+10*wratio+20, 1, '<')
        plt.plot(w_0-10*wratio-20, 1, '>')
        plt.plot(w_0-10*wratio, 1, '<')
        plt.show()
    return spe


def normalize_poly(spe, z):
    """
    """
    flux = spe.data
    wavelength = spe.wave.coord()
    flux_cut, wavelength_cut = cut_spectra(
        spe, w1_mgii*(1+z)-20, w1_mgii*wratio_mgii*(1+z)+20)
    poly = np.polyfit(wavelength_cut, flux_cut, 1)
    line_poly = line((poly[1], poly[0]), wavelength)
    spe_new = spe
    spe_new.data = flux/line_poly
    return spe_new


sigma_fixed = 2.7/(2*np.sqrt(2*np.log(2)))
err_sigma_fixed = 0


def gaussian(parameters, x):
    """
    """
    if len(parameters) == 2:
        mu, A = parameters
        sigma = sigma_fixed
    else:
        mu, A, sigma = parameters
    g = A * sigma * np.sqrt(2*np.pi) * scipy.stats.norm(loc=mu,
                                                        scale=sigma).pdf(x)
    return g


def line(parameters, x):
    m, n = parameters
    y = m * x + n
    return y


def simple_model(parameters, x):
    """
    """
    mu, A, m, n = parameters
    return 1-gaussian((mu, A), x)


w1_mgii = 2796.352
wratio_mgii = 1.0025672375


def double_model(parameters, x):
    """
    """
    A1, A2, mu = parameters
    return 1-gaussian((mu, A1), x)-gaussian((mu*wratio_mgii, A2), x)


w1_oii = 3727.092
wratio_oii = 3729.875/w1_oii


def double_model_em(parameters, x):
    """
    """
    A1, A2, mu = parameters
    return 1+gaussian((mu, A1), x) + gaussian((mu*wratio_emission), x)


def chi_cuadrado_abs(parameters, x, y):
    """
    """
    return y-double_model(parameters, x)


def chi_cuadrado_em(parameters, x, y):
    """
    """
    return y-double_model_em(parameters, x)


def fit_doublet(spe, z, how='abs', fwhm=None):
    """
    how: 'abs' for absorption
         'em' for emission
    """
    if how == 'abs':
        w1 = w1_mgii
        wratio = wratio_mgii
        function = chi_cuadrado_abs
    elif how == 'em':
        w1 = w1_oii
        wratio = wratio_oii
        function = chi_cuadrado_em
    # cut the spectra
    flux, wavelength = cut_spectra(spe, w1*(1+z)-20, w1*wratio*(1+z)+20)
    # print flux
    # print wavelength
    # define priors
    A_1 = abs(max(flux) - min(flux))
    A_2 = abs(max(flux) - min(flux))
    mu = w1*(1+z)
    # print line_prior[0]
    parameters = A_1, A_2, mu
    # make fits
    v, covar, info, mesg, success = leastsq(function, parameters,
                                            args=(wavelength, flux),
                                            full_output=1,  maxfev=100000)
    A_1_fit, A_2_fit, mu_fit = v[0], v[1], v[2]
    # calculate fluxes
    flux_1_fit = A_1_fit*sigma
    flux_2_fit = A_2_fit*sigma
    # calculate errors
    chisq = sum(info["fvec"] * info["fvec"])
    dof = len(info["fvec"]) - len(v)
    if covar is not None:
        err = np.array([np.sqrt(np.abs(covar[i, i])) *
                        np.sqrt(np.abs(chisq / dof)) for i in range(len(v))])
    else:
        err = None
    if err is not None:
        err_A_1 = err[0]
        err_A_2 = err[1]
        err_mu = err[2]
    else:
        err_A_1 = np.NAN
        err_A_2 = np.NAN
        err_mu = np.NAN
    return A_1_fit, err_A_1, A_2_fit, err_A_2, mu_fit, err_mu, flux_1_fit, flux_2_fit


def fit_line(spe, w1):
    """
    """
    # cut the spectra
    flux, wavelength = cut_spectra(spe, w1-20, w1+20)
    # define priors
    line = np.polyfit(wavelength, flux, 1)
    A = abs(max(flux) - min(flux))
    mu = np.mean(wavelength)
    parameters = mu, A, line[1], line[0]
    # make fits
    fit_gaussiano = leastsq(chi_cuadrado, parametros, args=(wavelength, flux))
    mu_fit, A_fit, n_fit, m = fit_gaussiano[0][0], fit_gaussiano[0][1], fit_gaussiano[0][2], fit_gaussiano[0][3]
    # calculate flux
    flux_fit = A_fit*sigma
    return mu_fit, A_fit, flux_fit, n_fit, m_fit


def flux(A, err_A, sigma=sigma_fixed):
    """
    """
    flux = A*sigma*np.sqrt(2*np.pi)
    err_flux = A*sigma*((err_A/A)**2+(err_sigma/sigma)
                        ** 2)**0.5*np.sqrt(2*np.pi)
    return flux, err_flux


def eq_width(flux, err_flux, z):
    """
    takes flux and turns it into equivalent width
    """
    ew = flux/(1+z)
    err_ew = err_flux/(1+z)
    return ew, err_ew


def velocity(w_obs, z, how='abs'):
    """
    """
    if how == 'abs':
        w1 = w1_mgii
    elif how == 'emi':
        w1 = w1_oii
    c = 299792.458
    z_obs = w_obs/w1 - 1
    vel = (z_obs-z)/(1+z)*c
    return vel


def def_zero(vel_matrix):
    """
    Recieves a velocity matrix and returns the position closer to 0
    """
    result = np.where(vel_matrix == np.nanmin(vel_matrix))
    return (result[0][0], result[1][0])


def calculate_distance(x1, y1, x2, y2, scale):
    """
    Recieves two point coordinates and calculates their distance
    """
    distance = np.sqrt((x1-x2)**2+(y1-y2)**2)*scale
    return distance


def scale_pixel(scale, distance):
    """
    turns pixel distances into kpc
    """
    kpc = scale * distance
    return kpc
