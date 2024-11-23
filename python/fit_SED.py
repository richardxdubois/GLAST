import numpy as np
from scipy.optimize import curve_fit


def SED_function(E, A, alpha, E_cut, E_0):
    E = np.float64(E)
    A = np.float64(A)
    alpha = np.float64(alpha)
    E_cut = np.float64(E_cut)
    E_0 = np.float64(E_0)

    return E*E*A * (E / E_0)**(-alpha) * np.exp(-E / E_cut)


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
    A, alpha, E_cut, E_0 = params

    return params, covariance
