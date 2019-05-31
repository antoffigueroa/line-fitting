import numpy as np
import matplotlib.pyplot as plt
from mpdaf.obj import Spectrum
from scipy.optimize import bisect, leastsq
import scipy.stats
from astropy.coordinates import SkyCoord


def cut_spectra(spe, w_ini, w_fin):
    spe_new = spe
    flux = spe_new.data
    wavelength = spe_new.wave.coord()
    p1 = np.int(spe_new.wave.pixel(w_ini))
    p2 = np.int(spe_new.wave.pixel(w_fin))
    cut_flux = flux[p1:p2]
    cut_wavelength = wavelength[p1:p2]
    return cut_flux, cut_wavelength

def normalize(spe, w_0, z, how='constant', wratio=None, show=False):
    if how == 'constant':
        return normalize_constant(spe, w_0, z, wratio, show)
    if how == 'line':
        return normalize_line(spe, w_0, z, wratio, show)
    if how == 'parabole':
        return normalize_parabole(spe, w_0, z, wratio, show)


def normalize_constant(spe, w_0, z, wratio=None, show=False):
    """
    spe: Spectrum object from mpdaf
    w_0: wavelength of the line you are working with
    wratio: only used if you know there is another line close to the line you
    are working with
    """
    if wratio is None:
        wratio_new = 1
    else:
        wratio_new = wratio
    r1 = np.int(spe.wave.pixel(w_0+20*wratio_new))
    r2 = np.int(spe.wave.pixel(w_0+20*wratio_new+20))
    l1 = np.int(spe.wave.pixel(w_0-20*wratio_new-20))
    l2 = np.int(spe.wave.pixel(w_0-20*wratio_new))
    flux = spe.data
    wavelength = spe.wave.coord()
    right_part = flux[r1:r2]
    left_part = flux[l1:l2]
    mean_right = np.mean(right_part)
    mean_left = np.mean(left_part)
    mean = (mean_right+mean_left)/2.0
    spe_new = spe.copy()
    spe_new.data = spe_new.data/mean
    if show:
        plt.figure()
        plt.plot(velocity(wavelength, z), spe.data, label='Spectrum')
        plt.plot(velocity(wavelength, z), line((0, mean), wavelength), label='Fitted continuum')
        plt.axvspan(velocity(w_0+20*wratio_new, z), velocity(w_0+20*wratio_new+20, z),
                    facecolor='#2ca02c', alpha=0.5, label='Window used to fit continuum')
        plt.axvspan(velocity(w_0-20*wratio_new-20, z), velocity(w_0 -
                                                                20*wratio_new, z), facecolor='#2ca02c', alpha=0.5)
        plt.xlim(-3000, 3000)
        plt.legend().draggable()
    return spe_new


def normalize_line(spe, w_0, z, wratio=None, show=False):
    """
    """
    if wratio is None:
        wratio_new = 1
    else:
        wratio_new = wratio
    flux = spe.data
    wavelength = spe.wave.coord()
    flux_cut_1, wavelength_cut_1 = cut_spectra(
        spe, w_0+20*wratio_new, w_0+20*wratio_new+20)
    flux_cut_2, wavelength_cut_2 = cut_spectra(
        spe, w_0-20*wratio_new-20, w_0-20*wratio_new)
    flux_cut = np.concatenate((flux_cut_1, flux_cut_2))
    wavelength_cut = np.concatenate((wavelength_cut_1, wavelength_cut_2))
    poly = np.polyfit(wavelength_cut, flux_cut, 1)
    line_poly = poly[0]*wavelength + poly[1]
    spe_new = spe.copy()
    spe_new.data = flux/line_poly
    if show:
        plt.figure()
        plt.plot(velocity(wavelength, z), spe.data, label='Spectrum')
        plt.plot(velocity(wavelength, z), line_poly, label='Fitted continuum')
        plt.axvspan(velocity(w_0+20*wratio_new, z), velocity(w_0+20*wratio_new+20, z),
                    facecolor='#2ca02c', alpha=0.5, label='Window used to fit continuum')
        plt.axvspan(velocity(w_0-20*wratio_new-20, z), velocity(w_0 -
                                                                20*wratio_new, z), facecolor='#2ca02c', alpha=0.5)
        plt.xlim(-3000, 3000)
        plt.legend().draggable()
    return spe_new

def normalize_parabole(spe, w_0, z, wratio=None, show=False):
    if wratio is None:
        wratio_new = 1
    else:
        wratio_new = wratio
    flux = spe.data
    wavelength = spe.wave.coord()
    flux_cut_1, wavelength_cut_1 = cut_spectra(
        spe, w_0+20*wratio_new, w_0+20*wratio_new+20)
    flux_cut_2, wavelength_cut_2 = cut_spectra(
        spe, w_0-20*wratio_new-20, w_0-20*wratio_new)
    flux_cut = np.concatenate((flux_cut_1, flux_cut_2))
    wavelength_cut = np.concatenate((wavelength_cut_1, wavelength_cut_2))
    poly = np.polyfit(wavelength_cut, flux_cut, 2)
    parabole = poly[0]*wavelength**2 + poly[1]*wavelength + poly[2]
    spe_new = spe.copy()
    spe_new.data = flux/parabole
    if show:
        plt.figure()
        plt.plot(velocity(wavelength, z), spe.data, label='Spectrum')
        plt.plot(velocity(wavelength, z), parabole, label='Fitted continuum')
        plt.axvspan(velocity(w_0+20*wratio_new, z), velocity(w_0+20*wratio_new+20, z),
                    facecolor='#2ca02c', alpha=0.5, label='Window used to fit continuum')
        plt.axvspan(velocity(w_0-20*wratio_new-20, z), velocity(w_0 -
                                                                20*wratio_new, z), facecolor='#2ca02c', alpha=0.5)
        plt.xlim(-3000, 3000)
        plt.legend().draggable()
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
    mu, A, sigma, m, n = parameters
    return gaussian((mu, A, sigma), x) + m*x + n


def chi_cuadrado_g(parameters, x, y):
    return y-simple_model(parameters, x)


w1_mgii = 2796.352
wratio_mgii = 1.0025672375


def double_model(parameters, x):
    """
    """
    if len(parameters) == 3:
        A1, A2, mu = parameters
        return 1-gaussian((mu, A1), x)-gaussian((mu*wratio_mgii, A2), x)
    else:
        A1, A2, mu, sigma_1, sigma_2 = parameters
        return 1-gaussian((mu, A1, sigma_1), x)-gaussian((mu*wratio_mgii, A2, sigma_2), x)


w1_oii = 3727.092
wratio_oii = 3729.875/w1_oii


def double_model_em(parameters, x):
    """
    """
    if len(parameters) == 3:
        A1, A2, mu = parameters
        return 1+gaussian((mu, A1), x)+gaussian((mu*wratio_mgii, A2), x)
    else:
        A1, A2, mu, sigma_1, sigma_2 = parameters
        return 1+gaussian((mu, A1, sigma_1), x)+gaussian((mu*wratio_mgii, A2, sigma_2), x)


def chi_cuadrado_abs(parameters, x, y):
    """
    """
    return y-double_model(parameters, x)


def chi_cuadrado_em(parameters, x, y):
    """
    """
    return y-double_model_em(parameters, x)


def fit_doublet(spe, z, how='abs', fwhm=2.7):
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
    sigma_1 = sigma_fixed
    sigma_2 = sigma_fixed
    if fwhm != 2.7:
        parameters = A_1, A_2, mu, sigma_1, sigma_2
    else:
        parameters = A_1, A_2, mu
    # make fits
    v, covar, info, mesg, success = leastsq(function, parameters,
                                            args=(wavelength, flux),
                                            full_output=1,  maxfev=100000)
    A_1_fit, A_2_fit, mu_fit = v[0], v[1], v[2]
    if fwhm != 2.7:
        sigma_1_fit, sigma_2_fit = v[3], v[4]
        fitted_parameters = [A_1_fit, A_2_fit,
                             mu_fit, sigma_1_fit, sigma_2_fit]
    else:
        fitted_parameters = [A_1_fit, A_2_fit, mu_fit]
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
        if fwhm != 2.7:
            err_sigma_1 = err[3]
            err_sigma_2 = err[4]
    else:
        err_A_1 = np.NAN
        err_A_2 = np.NAN
        err_mu = np.NAN
        if fwhm != 2.7:
            err_sigma_1 = np.NAN
            err_sigma_2 = np.NAN
    if fwhm != 2.7:
        error_parameters = [err_A_1, err_A_2, err_mu, err_sigma_1, err_sigma_2]
    else:
        error_parameters = [err_A_1, err_A_2, err_mu]
    return fitted_parameters, error_parameters


def fit_gaussian(wavelength, flux):
    """
    """
    # define priors
    A = abs(max(flux) - min(flux))
    mu = np.mean(wavelength)
    sigma = np.std(flux)
    m = 0
    n = 1
    parameters = mu, A, sigma, m, n
    # make fits
    fit_gaussiano = leastsq(chi_cuadrado_g, parameters, args=(
        wavelength, flux), full_output=1,  maxfev=100000)
    mu_fit, A_fit, sigma_fit, m_fit, n_fit = fit_gaussiano[0][0], fit_gaussiano[0][1], fit_gaussiano[0][2], fit_gaussiano[0][3], fit_gaussiano[0][4]
    return mu_fit, A_fit, sigma_fit, m_fit, n_fit


def flux(A, err_A, sigma=sigma_fixed, err_sigma=err_sigma_fixed):
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


def calculate_distance(c1, c2):
    """
    Recieves two point coordinates and calculates their distance
    """
    distance = c1.separation(c2)
    return distance


def scale_pixel(scale, distance):
    """
    turns pixel distances into kpc
    """
    kpc = scale * distance
    return kpc

def copy_header(headered_file, unheadered_file):
    """
    Copies the header of one cube and the data of another cube and creates a new
    cube with the header and the data
    """
    fits_header = fits.open(headered_file)
    w = wcs.WCS(fits_header[1].header)
    header = w.to_header()
    fits_header.close()
    fits_datos = fits.open(unheadered_file)
    data = fits_datos[0].data
    fits_datos.close()
    hdu = fits.PrimaryHDU(data, header=header)
    return hdu

def copy_header_npy(headered_file, numpy_array):
    """
    Copies the header of one cube and the data of another cube and creates a new
    cube with the header and the data
    """
    fits_header = fits.open(headered_file)
    w = wcs.WCS(fits_header[1].header)
    header = w.to_header()
    fits_header.close()
    data = numpy_array
    hdu = fits.PrimaryHDU(data, header=header)
    return hdu
