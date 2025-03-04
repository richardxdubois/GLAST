import pylab as pl
import numpy as np
from scipy.stats import chi2,norm
from scipy.signal import find_peaks
import yaml
import argparse
from datetime import datetime, date
from astropy.io import fits

import os
from pathlib import Path

from godot import core

from bokeh.models.widgets import Div
from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column, gridplot
from bokeh.models import Label, Span

parser = argparse.ArgumentParser()

parser.add_argument('--app_config',
                    default="process_exposure_config.yaml",
                    help="overall app config file")
args = parser.parse_args()

with open(args.app_config, "r") as f:
    data = yaml.safe_load(f)

source = data["source"]
source_FGL = data["source_FGL"]
try:
    in_file = data["in_file"]
except KeyError:
    in_file = None

in_dir = data["in_dir"]
ft2 = data["ft2"]

ra = data["ra"]
dec = data["dec"]
emin = data["emin"]
emax = data["emax"]

try:
    tmin = data["tmin"]
    tmax = data["tmax"]
except KeyError:
    tmax = None
    tmin = None

porb = data["porb"]
psuper = data["psuper"]
forb = data["forb"]
fprec = data["fprec"]

orb_low = data["orb_low"]
orb_high = data["orb_high"]
super_low = data["super_low"]
super_high = data["super_high"]

html_file = data["html_file"]
html_title = data["html_title"]

if in_file is None:
    # Get the list of files in the directory
    ft1_u = [os.path.join(in_dir, f) for f in os.listdir(in_dir) if os.path.isfile(os.path.join(in_dir, f))]
    ft1 = sorted(ft1_u)
else:
    ft1 = [in_file]

print("Input files:", ft1)

for f in ft1:
    hdu = fits.open(f)
    hp = hdu[0].header
    hg = hdu['gti']
    fpath = Path(f)
    print(fpath.name, hp.get('TSTART'), hp.get('TSTOP'), hg.data['start'][0], hg.data['stop'][-1])
    hdu.close()

spectrum = lambda E: (E/1000)**-2.1

data = core.Data(ft1, ft2, ra, dec, weight_col=source_FGL, base_spectrum=spectrum, zenith_cut=90, emin=emin,
                 emax=emax, tstart=tmin, tstop=tmax)

print("data prep done so far")

ts = data.get_cells(tcell=300,time_series_only=True, trim_zero_exposure=False, use_barycenter=True)
f,window = core.power_spectrum_fft(ts,exp_only=True)
scale = 50./40000
f,dlogl_nobg,dlogl,dlogl_null = core.power_spectrum_fft(ts)
fday = f*86400

peaks_ls, props_ls = find_peaks(dlogl_nobg, height=0.5 * max(dlogl_nobg))
pk_days = (1. / f[peaks_ls] / 86400.)
print("found power peaks (days)", pk_days)

freqs = np.asarray([fprec,forb,2*forb])
corr,pows = core.get_orbital_modulation(ts,freqs)
f2,dlogl_nobg2,dlogl2,dlogl_null2 = core.power_spectrum_fft(ts, exposure_correction=corr)

add_power = np.zeros_like(dlogl_nobg2)
for freq,p in zip(freqs,pows):
    idx = np.argmin(np.abs(f[1:]-freq))
    add_power[idx] = p

print("power spectrum done")

fmask = [i for i in range(len(fday)) if fday[i] > 1./orb_high and fday[i] < 1./orb_low]

pday = np.array([1./f for f in fday if f != 0])
print("pday: ", len(pday), max(pday), min(pday))

pmask = [i for i in range(len(pday)) if pday[i] > orb_low and pday[i] < orb_high]
smask = [i for i in range(len(pday)) if pday[i] > super_low and pday[i] < super_high]

vline_p1 = Span(location=porb, dimension='height', line_color='red', line_width=2,
                line_dash='dashed')
vline_p2 = Span(location=psuper, dimension='height', line_color='blue', line_width=2,
                line_dash='dashed')

res_label = Label(x=23.5, y=400., text_font_size="8pt",
                  text="Peak : " + str('{0:.3f}'.format(pk_days[0])) + " days")


fig1 = figure(title="frequency: dlogl_nobg", y_axis_type="log", width=800, height=640)
fig1.line(fday[fmask], dlogl_nobg[fmask])

fig1a = figure(title="frequency: dlogl_nobg", width=800, height=640)
fig1a.line(fday[fmask], dlogl_nobg[fmask])

fig2 = figure(title="period: dlogl_nobg", y_axis_type="log", width=800, height=640)
fig2.line(pday[pmask], dlogl_nobg[pmask])
fig2.add_layout(vline_p1)

fig3a = figure(title="period: dlogl_nobg", width=800, height=640)
fig3a.line(pday[pmask], dlogl_nobg[pmask])
fig3a.add_layout(vline_p1)
fig3a.add_layout(res_label)

fig4a = figure(title="period: dlogl_nobg", width=800, height=640)
fig4a.line(pday[smask], dlogl_nobg[smask])
fig4a.add_layout(vline_p2)

d_text = source + " Run on: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S")
if in_file is None:
    d_text += " " + in_dir
else:
    d_text += " " + in_file

c_text = "<br>Emin " + str(emin) + " Emax " + str(emax) + " Tmin " + str(tmin) + " Tmax " + str(tmax)
div_text = d_text + c_text
del_div = Div(text=div_text)

l = layout(del_div, row(fig1, fig1a), row(fig2, fig3a), fig4a)

output_file(html_file)
save(l, title=html_title)
