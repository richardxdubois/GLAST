import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import quad

def SED_function(E, A, alpha, E_cut):
    E = np.float64(E)
    A = np.float64(A)
    alpha = np.float64(alpha)
    E_cut = np.float64(E_cut)
    #E_0 = np.float64(E_0)

    return E*E*A * (E / 1800.)**(-alpha) * np.exp(-E / E_cut)

def flux_function(E, A, alpha, E_cut):
    E = np.float64(E)
    A = np.float64(A)
    alpha = np.float64(alpha)
    E_cut = np.float64(E_cut)
    #E_0 = np.float64(E_0)
    return A * (E/1800.)**(-alpha) * np.exp(-E/E_cut)

def fit_SED(E, flux, errors, initial_guesses):

    # Fit the model
    params, covariance = curve_fit(
        SED_function,
        E,
        flux,
        p0=initial_guesses,
        sigma=errors,
        absolute_sigma=True, method="trf", bounds=(0, [1e-4, 3., 15000.])
        )

    return params, covariance


def fit_SED_errors(E_min, E_max, params, covariance):

    i_params = []
    integral_diffs = []

    for p in range(len(params)):

        cov = covariance[p][p]
        err = []
        for i in range(-1,2):
            i_params = list(params)
            i_params[p] += i*np.sqrt(cov)
            integrated_fits, int_error = quad(flux_function, E_min, E_max,
                                              args=(i_params[0], i_params[1], i_params[2]))
            err.append(integrated_fits)
        integral_diffs.append(abs(err[1]-err[0]))
        integral_diffs.append(abs(err[2]-err[1]))

    return max(integral_diffs)
