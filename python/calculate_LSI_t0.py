from datetime import datetime
from astropy.io import fits
import math
from gti import Gti  # gift from Matthew Kerr: https://github.com/kerrm/godot/blob/master/gti.py

from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column

import numpy as np

t0_MJD = 43366.275
daysecs = 86400.
orb = 26.4960 * daysecs

#infile = '/Users/richarddubois/Code/GLAST/tmp/dbg/L24082417075904476C3F57_PH00.fits'
infile = '/Users/richarddubois/Code/GLAST/tmp/dbg/ft1_0-0.fits'
h = fits.open(infile)
h.info()
gti_hdr = h[2]
print(gti_hdr.columns)

# Access the primary HDU (Header/Data Unit)
primary_hdu = h[0]
header = primary_hdu.header

t_max = np.float128(header["TSTOP"])
t_min = np.float128(header["TSTART"])

gti_starts = gti_hdr.data["START"]
gti_stops = gti_hdr.data["STOP"]
file_gti = Gti(gti_starts, gti_stops)

num_periods = (t_max - t_min) / orb

date_t0_MJD = datetime.strptime("1977-08-11 06:36:00.000", "%Y-%m-%d %H:%M:%S.%f")

start_MET = datetime.strptime("2001-01-01 00:00:00.000", "%Y-%m-%d %H:%M:%S.%f")
aug_2008 = datetime.strptime("2008-08-15 00:00:00.000", "%Y-%m-%d %H:%M:%S.%f")

delta = np.float128((start_MET - date_t0_MJD).total_seconds())
delta_08 = np.float128((aug_2008 - start_MET).total_seconds())
delta_sum = delta_08 + delta
mod_t0 = np.mod(delta_sum, orb)

t0_post_MET = np.float128(delta_08 - mod_t0)
check = delta_sum - mod_t0
check_mod = np.mod(check, orb)

print(num_periods, delta, t0_post_MET, check_mod)

gti = {}
new_gti = {}

for b in np.arange(10.):

    gti[b] = []
    p = b/10.
    t = np.float128(t0_post_MET - math.ceil((t_min - t0_post_MET)/orb) * orb) # start on orbit boundary before first time measure

    while t < t_max + orb:
        gti[b].append(np.array([t+p*orb, t+(p+0.1)*orb]))
        t += orb

    t_start = np.array(list(zip(*gti[b]))[0])
    t_stop = np.array(list(zip(*gti[b]))[1])

    new_gti[b] = Gti(t_start, t_stop)

    rc = new_gti[b].intersection(file_gti)
    merged_starts = new_gti[b].get_edges(starts=True)
    merged_stops = new_gti[b].get_edges(starts=False)

    gti_hdr.header["NAXIS2"] = len(merged_starts)
    t0 = fits.Column(name='start', array=merged_starts, format='D', unit='s')
    t1 = fits.Column(name='stop', array=merged_stops, format='D', unit='s')
    cols = fits.ColDefs([t0, t1])
    hdu = fits.BinTableHDU.from_columns(cols)
    gti_hdr = hdu
    #gti_hdr.data["START"] = merged_starts
    #gti_hdr.data["STOP"] = merged_stops
    h.writeto(infile, overwrite=True)

    int_diffs = np.diff(gti_starts)/daysecs
    int_lens = np.array((gti_stops - gti_starts))/daysecs
    print(len(int_diffs), len(int_lens))

    q_hist = figure(title="stop-start",
                    x_axis_label='Time (days MET)', y_axis_label='counts',
                    width=750)
    counts_hist, counts_edges = np.histogram(int_lens, bins=200)
    q_hist.vbar(top=counts_hist, x=counts_edges[1:], width=counts_edges[1]-counts_edges[0], fill_color='red',
                fill_alpha=0.2, bottom=0)

    r_hist = figure(title="start - start",
                    x_axis_label='Time (days MET)', y_axis_label='counts',
                    width=750)
    diff_hist, diff_edges = np.histogram(int_diffs, bins=200, range=(0., 0.5))
    r_hist.vbar(top=diff_hist, x=diff_edges[1:], width=diff_edges[1]-diff_edges[0], fill_color='red',
                fill_alpha=0.2, bottom=0)

    print(len(gti[b]))
    break

print("done")

output_file("../../../../Home/HomeStuff/python/test_phases.html")
l = layout(q_hist, r_hist)
save(l)
