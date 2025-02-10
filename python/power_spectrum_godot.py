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
    print(fpath.name, hp.get('TSTART'), hp.get('TSTOP'), hg.data['start'], hg.data['stop'])

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

fmask = fday < 0.1

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


fig1 = figure(title="frequency: dlogl", y_axis_type="log", width=800, height=640)
fig1.line(fday[fmask], dlogl_nobg[fmask])

fig2 = figure(title="period: dlogl", y_axis_type="log", width=800, height=640)
fig2.line(pday[pmask], dlogl[pmask])
fig2.add_layout(vline_p1)

fig3 = figure(title="period: dlogl", width=800, height=640)
fig3.line(pday[pmask], dlogl[pmask])
fig3.add_layout(vline_p1)

fig3a = figure(title="period: dlogl_nobg", width=800, height=640)
fig3a.line(pday[pmask], dlogl_nobg[pmask])
fig3a.add_layout(vline_p1)
fig3a.add_layout(res_label)

fig3b = figure(title="period: dlogl_nobg + add_power", width=800, height=640)
fig3b.line(pday[pmask], (dlogl+add_power)[pmask])
fig3b.add_layout(vline_p1)

fig4 = figure(title="period: dlogl", width=800, height=640)
fig4.line(pday[smask], dlogl[smask])
fig4.add_layout(vline_p2)

fig4a = figure(title="period: dlogl_nobg", width=800, height=640)
fig4a.line(pday[smask], dlogl_nobg[smask])
fig4a.add_layout(vline_p2)

fig12 = figure(title="frequency: dlogl2", y_axis_type="log", width=800, height=640)
fig12.line(fday[fmask], dlogl2[fmask])

fig22 = figure(title="period: dlogl", y_axis_type="log", width=800, height=640)
fig22.line(pday[pmask], dlogl2[pmask])
fig22.add_layout(vline_p1)

fig32 = figure(title="period: dlogl2", width=800, height=640)
fig32.line(pday[pmask], dlogl2[pmask])
fig32.add_layout(vline_p1)

fig32a = figure(title="period: dlogl_nobg2", width=800, height=640)
fig32a.line(pday[pmask], dlogl_nobg2[pmask])
fig32a.add_layout(vline_p1)

fig32b = figure(title="period: dlogl_nobg2 + add_power", width=800, height=640)
fig32b.line(pday[pmask], (dlogl_nobg2+add_power)[pmask])
fig32b.add_layout(vline_p1)

fig42 = figure(title="period: dlogl2", width=800, height=640)
fig42.line(pday[smask], dlogl2[smask])
fig42.add_layout(vline_p2)

fig42a = figure(title="period: dlogl_nobg2", width=800, height=640)
fig42a.line(pday[smask], dlogl_nobg2[smask])
fig42a.add_layout(vline_p2)

d_text = source + " Run on: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S")
if in_file is None:
    d_text += " " + in_dir
else:
    d_text += " " + in_file

del_div = Div(text=d_text)

l = layout(del_div, fig1, row(fig2, fig3), row(fig3a, fig3b), row(fig4, fig4a),
           fig12, row(fig22, fig32), row(fig32a, fig32b), row(fig42, fig42a))

output_file(html_file)
save(l, title=html_title)
