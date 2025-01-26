import pylab as pl
import numpy as np
from scipy.stats import chi2,norm
import os

from godot import core

from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column, gridplot
from bokeh.models import (Label, Span, LinearAxis, Range1d, Whisker, ColumnDataSource, BasicTicker, Tabs,
                          TabPanel, RangeSlider, CustomJS, Button, NumeralTickFormatter, BasicTickFormatter)
from bokeh.models.widgets import Div
from bokeh.palettes import Plasma256 as palette
from bokeh.transform import linear_cmap, log_cmap

in_dir = "/sdf/home/r/richard/fermi-user/LSI61303/periods/periodicity/godot/diffrsp/fits/"

# Get the list of files in the directory
ft1_u = [os.path.join(in_dir, f) for f in os.listdir(in_dir) if os.path.isfile(os.path.join(in_dir, f))]
ft1 = sorted(ft1_u)
print("Input files:", ft1)

ft2 = ["/sdf/home/r/richard/fermi-user/LSI61303/fssc_data/L24082417075904476C3F57_SC00.fits"]

ra = 40.143
dec = 61.229
spectrum = lambda E: (E/1000)**-2.1

data = core.Data(ft1, ft2, ra, dec, weight_col="4FGL J0240.5+6113", base_spectrum=spectrum, zenith_cut=90)

print("data prep done so far")

ts = data.get_cells(tcell=300,time_series_only=True, trim_zero_exposure=False, use_barycenter=True)
f,window = core.power_spectrum_fft(ts,exp_only=True)
scale = 50./40000
f,dlogl_nobg,dlogl,dlogl_null = core.power_spectrum_fft(ts)
fday = f*86400

forb = 2.963145573933919e-06
#forb = 4.3676e-7
fprec = 2.1777777777777778e-07
freqs = np.asarray([fprec,forb,2*forb])
corr,pows = core.get_orbital_modulation(ts,freqs)
f2,dlogl_nobg2,dlogl2,dlogl_null2 = core.power_spectrum_fft(ts, exposure_correction=corr)

add_power = np.zeros_like(dlogl_nobg2)
for freq,p in zip(freqs,pows):
    idx = np.argmin(np.abs(f[1:]-freq))
    add_power[idx] = p

print("power spectrum done")

fmask = fday < 0.1

pday = np.array([1./f for f in fday if f != 0])
print("pday: ", len(pday), max(pday), min(pday))

pmask = [i for i in range(len(pday)) if pday[i] > 23 and pday[i] < 29]
smask = [i for i in range(len(pday)) if pday[i] > 1000 and pday[i] < 3000]

fig1 = figure(title="frequency: dlogl", y_axis_type="log", width=800, height=640)
fig1.line(fday[fmask], dlogl_nobg[fmask])

fig2 = figure(title="period: dlogl", y_axis_type="log", width=800, height=640)
fig2.line(pday[pmask], dlogl[pmask])

fig3 = figure(title="period: dlogl", width=800, height=640)
fig3.line(pday[pmask], dlogl[pmask])

fig3a = figure(title="period: dlogl_nobg", width=800, height=640)
fig3a.line(pday[pmask], dlogl_nobg[pmask])

fig3b = figure(title="period: dlogl_nobg + add_power", width=800, height=640)
fig3b.line(pday[pmask], (dlogl+add_power)[pmask])

fig4 = figure(title="period: dlogl", width=800, height=640)
fig4.line(pday[smask], dlogl[smask])

fig4a = figure(title="period: dlogl_nobg", width=800, height=640)
fig4a.line(pday[smask], dlogl_nobg[smask])

fig12 = figure(title="frequency: dlogl2", y_axis_type="log", width=800, height=640)
fig12.line(fday[fmask], dlogl2[fmask])

fig22 = figure(title="period: dlogl", y_axis_type="log", width=800, height=640)
fig22.line(pday[pmask], dlogl2[pmask])

fig32 = figure(title="period: dlogl2", width=800, height=640)
fig32.line(pday[pmask], dlogl2[pmask])

fig32a = figure(title="period: dlogl_nobg2", width=800, height=640)
fig32a.line(pday[pmask], dlogl_nobg2[pmask])

fig32b = figure(title="period: dlogl_nobg2 + add_power", width=800, height=640)
fig32b.line(pday[pmask], (dlogl_nobg2+add_power)[pmask])

fig42 = figure(title="period: dlogl2", width=800, height=640)
fig42.line(pday[smask], dlogl2[smask])

fig42a = figure(title="period: dlogl_nobg2", width=800, height=640)
fig42a.line(pday[smask], dlogl_nobg2[smask])

l = layout(fig1, row(fig2, fig3), row(fig3a, fig3b), row(fig4, fig4a),
           fig12, row(fig22, fig32), row(fig32a, fig32b), row(fig42, fig42a))

output_file("power_spectrum_godot.html")
save(l)
