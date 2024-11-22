import numpy as np
from scipy.optimize import curve_fit


def SED_function(E, A, alpha, E_cut, beta, E_0):
    return E*E*A * (E / E_0)**(-alpha) * np.exp(-E / E_cut) + beta


def fit_SED(E, flux, errors, initial_guesses):

    # Fit the model
    params, covariance = curve_fit(
        SED_function,
        E,
        flux,
        p0=initial_guesses,
        sigma=errors,
        absolute_sigma=True
    )


    # Extract fitted parameters
    A, alpha, E_cut, beta, E_0 = params

    return params, covariance
