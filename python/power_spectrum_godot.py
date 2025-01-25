import pylab as pl
import numpy as np
from scipy.stats import chi2,norm
import os

from godot import core

in_dir = "/sdf/home/r/richard/fermi-user/LSI61303/periods/periodicity/godot/diffrsp/fits/"

# Get the list of files in the directory
ft1_u = [f for f in os.listdir(in_dir) if os.path.isfile(os.path.join(in_dir, f))]
ft1 = sorted(ft1_u)

ft2 = ["/sdf/home/r/richard/fermi-user/LSI61303/fssc_data/L24082417075904476C3F57_SC00.fits"]

ra = 40.143
dec = 61.229
spectrum = lambda E: (E/1000)**-2.1

data = core.Data(ft1, ft2, ra, dec, weight_col="4FGL J0240.5+6113", base_spectrum=spectrum, zenith_cut=90)

print("done so far")

"""
ts = data.get_cells(tcell=300,time_series_only=True,
            trim_zero_exposure=False,use_barycenter=True)
    f,window = core.power_spectrum_fft(ts,exp_only=True)
    scale = 50./40000
    f,dlogl_nobg,dlogl,dlogl_null = core.power_spectrum_fft(ts)
    fday = f*86400

    forb = 2.963145573933919e-06
    fprec = 2.1777777777777778e-07
    freqs = np.asarray([fprec,forb,2*forb])
    corr,pows = core.get_orbital_modulation(ts,freqs)
    f2,dlogl_nobg2,dlogl2,dlogl_null2 = core.power_spectrum_fft(ts,
            exposure_correction=corr)

    add_power = np.zeros_like(dlogl_nobg2)
    for freq,p in zip(freqs,pows):
        idx = np.argmin(np.abs(f[1:]-freq))
        add_power[idx] = p

"""
