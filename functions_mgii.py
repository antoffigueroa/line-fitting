import numpy as np
import matplotlib.pyplot as plt
from mpdaf.obj import Spectrum
from scipy.optimize import bisect, leastsq
import scipy.stats


def create_mask():
    """
    """


def cut_spectra(spe, w_ini, w_fin):
    flux = spe.data
    wavelength = spe.wave.coord()
    spe_nuevo = spe
    p1 = np.int(spe.wave.pixel(w_ini))
    p2 = np.int(spe.wave.pixel(w_fin))
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
    wratio = wratio*(1+z)
    if wratio is None:
        wratio = 1
    r1 = np.int(spe.wave.pixel(w_0+10*wratio))
    r2 = np.int(spe.wave.pixel(w_0+10*wratio+20))
    l1 = np.int(spe.wave.pixel(w_0-10*wratio-20))
    l2 = np.int(spe.wave.pixel(w_0-10*wratio))
    print r1, r2, l1, l2
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


def gaussian(parameters, x):
    """
    """
    mu, sigma, A = parameters
    g = A * scipy.stats.norm(loc=mu, scale=sigma).pdf(x)
    return g


def line(parameters, x):
    """
    """
    m, n = parameters
    y = m * x + n
    return y


def simple_model(parameters, x):
    """
    """
    mu, sigma, A, m, n = parameters
    return line((m, n), x)-gaussian((mu, sigma, A), x)


def double_model(parameters, x):
    """
    """
    mu1, sigma1, A1, mu2, sigma2, A2, m, n = parameters
    return line((m, n), x)-gaussian((mu1, sigma1, A1), x)-gaussian((mu2, sigma2, A2), x)


def chi_cuadrado(parameters, x, y):
    """
    """
    return y-double_model(parameters, x)


def fit_doublet(spe, w1, w2):
    """
    """
    # cut the spectra
    flux, wavelength = cut_spectra(spe, w1-20, w2+20)
    # print flux
    # print wavelength
    # define priors
    line_prior = np.polyfit(wavelength, flux, 1)
    print (3*w1-w2)/2.0, (w1+w2)/2.0
    flux_cut1, wl_cut1 = cut_spectra(spe, (3*w1-w2)/2.0, (w1+w2)/2.0)
    A_1 = abs(max(flux_cut1) - min(flux_cut1))
    mu_1 = np.mean(wl_cut1)
    sigma_1 = np.std(wl_cut1)
    flux_cut2, wl_cut2 = cut_spectra(spe, (w1+w2)/2.0, (3*w2-w1)/2.0)
    A_2 = abs(max(flux_cut2) - min(flux_cut2))
    mu_2 = np.mean(wl_cut2)
    sigma_2 = np.std(wl_cut2)
    # print line_prior[0]
    parameters = mu_1, sigma_1, A_1, mu_2, sigma_2, A_2, line_prior[1], line_prior[0]
    # make fits
    v, covar, info, mesg, success = leastsq(chi_cuadrado, parameters, args=(wavelength, flux), full_output=1,  maxfev=100000)
    mu_1_fit, sigma_1_fit, A_1_fit, mu_2_fit, sigma_2_fit, A_2_fit, n_fit, m_fit = v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]
    # calculate fluxes
    flux_1_fit = A_1_fit*sigma_1_fit
    flux_2_fit = A_2_fit*sigma_2_fit
    # calculate errors
    chisq = sum(info["fvec"] * info["fvec"])
    dof = len(info["fvec"]) - len(v)
    if covar is not None:
        err = np.array([np.sqrt(np.abs(covar[i, i])) * np.sqrt(np.abs(chisq / dof)) for i in range(len(v))])
    else:
        err = None
    if err is not None:
        err_mu_1 = err[0]
        err_sigma_1 = err[1]
        err_A_1 = err[2]
        err_mu_2 = err[3]
        err_sigma_2 = err[4]
        err_A_2 = err[5]
        err_n = err[6]
        err_m = err[7]
    else:
        err_mu_1 = np.NAN
        err_sigma_1 = np.NAN
        err_A_1 = np.NAN
        err_mu_2 = np.NAN
        err_sigma_2 = np.NAN
        err_A_2 = np.NAN
        err_n = np.NAN
        err_m = np.NAN
    return mu_1_fit, err_mu_1, mu_2_fit, err_mu_2, sigma_1_fit, err_sigma_1, sigma_2_fit, err_sigma_2, A_1_fit, err_A_1, A_2_fit, err_A_2, flux_1_fit, flux_2_fit, n_fit, m_fit


def fit_line(spe, w1):
    """
    """
    # cut the spectra
    flux, wavelength = cut_spectra(spe, w1-20, w1+20)
    # define priors
    line = np.polyfit(wavelength, flux, 1)
    A = abs(max(flux) - min(flux))
    mu = np.mean(wavelength)
    sigma = np.std(wavelength)
    parameters = mu, sigma, A, line[1], line[0]
    # make fits
    fit_gaussiano = leastsq(chi_cuadrado, parametros, args=(wavelength, flux))
    mu_fit, sigma_fit, A_fit, n_fit, m = fit_gaussiano[0][0], fit_gaussiano[0][1], fit_gaussiano[0][2], fit_gaussiano[0][3], fit_gaussiano[0][4]
    # calculate flux
    flux_fit = A_fit*sigma_fit
    return mu_fit, sigma_fit, A_fit, flux_fit, n_fit, m_fit


def flux(A, sigma, err_A, err_sigma):
    """
    """
    flux = A*sigma
    err_flux = A*sigma*((err_A/A)**2+(err_sigma/sigma)**2)**0.5
    return flux, err_flux


def eq_width(flux, err_flux, z):
    """
    takes flux and turns it into equivalent width
    """
    ew = -1.0*flux/(1+z)
    err_ew = err_flux/(1+z)
    return ew, err_ew
