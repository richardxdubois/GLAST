from astropy.io import fits
import numpy as np
from datetime import datetime, date, timedelta
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

from bokeh.models.widgets import Tabs, Panel, DataTable, TableColumn, Div, NumberFormatter
from bokeh.models.formatters import DatetimeTickFormatter
from bokeh.models import ColumnDataSource, Span, LinearAxis, Range1d, LinearColorMapper
from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column

source = "LS5039"
html_name = source + ".html"

# Open the FITS file
if source == "LSI61303":
    #hdul = fits.open('/Users/richarddubois/Code/GLAST/tmp/LSI61303_1_deg_500s.fits')
    hdul = fits.open('/Users/richarddubois/Code/GLAST/tmp/LSI61303_3_deg_10800s.fits')
    f_start = 1./40./86400.
    f_stop = 1./5./86400.
    nom_freq = 1/26.5/86400.
elif source == "LS5039":
    #hdul = fits.open('/Users/richarddubois/Code/GLAST/tmp/LS5039_1_deg_500s.fits')
    hdul = fits.open('/Users/richarddubois/Code/GLAST/tmp/LS5039_3_deg_10800s.fits')
    f_start = 1./10./86400.
    f_stop = 1./2./86400.
    nom_freq = 1/3.9/86400.

# Print information about the FITS file
hdul.info()

# Access the primary HDU (Header/Data Unit)
primary_hdu = hdul[0]

# Get the data and header
print(hdul[1].columns)
header = primary_hdu.header
data = hdul[1].data

indices = np.where(data.COUNTS > 0.)

suppress_zero = True   # suppress entries with zero counts

if suppress_zero:
    d_non_zero_times = data.TIME[indices]
    exp = data.EXPOSURE[indices]
    d_weighted = data.COUNTS[indices]
else:
    d_non_zero_times = data.TIME
    exp = data.EXPOSURE
    d_weighted = data.COUNTS

r_weighted = []


mean_counts = np.average(d_weighted)
mean_rate = np.sum(d_weighted) / np.sum(exp)
mean_sub = d_weighted - mean_counts

for i, e in enumerate(exp):
    if e == 0.:
        r_weighted.append(0.)
    else:
        ave = mean_rate/np.sqrt(e)  # "predicted" count error
        r_weighted.append(mean_sub[i] / ave)

d_weighted = np.array(r_weighted)
counts_hist, counts_edges = np.histogram(d_weighted, bins=100)
q_hist = figure(title="Counts",
                x_axis_label='Counts', y_axis_label='counts',
                width=750)

q_hist.vbar(top=counts_hist, x=counts_edges[1:], width=counts_edges[1]-counts_edges[0], fill_color='red',
            fill_alpha=0.2, bottom=0)

exp_hist, exp_edges = np.histogram(exp, bins=100)
r_hist = figure(title="Exposure",
                x_axis_label='cm^2 s', y_axis_label='counts',
                width=750)

r_hist.vbar(top=exp_hist, x=exp_edges[1:], width=exp_edges[1]-exp_edges[0], fill_color='red', fill_alpha=0.2, bottom=0)

t_start = data.TIME[0]
t_end = data.TIME[-1]

print("t_start", t_start, "t_end", t_end, "num orig bins", len(data.TIME), "num non-zero bins", len(d_weighted))

frequency = np.linspace(f_start, f_stop, 100000)  # for orbital
freq_days = 1./frequency/86400.

power = LombScargle(d_non_zero_times, d_weighted).power(frequency)
#power = LombScargle(data.TIME, data.COUNTS).power(frequency)
#power = LombScargle(data.TIME, data.EXPOSURE).power(frequency)
print(max(frequency), max(power))
print(min(frequency), min(power))

peaks_ls, props_ls = find_peaks(power, height=0.1 * max(power))
print(peaks_ls, props_ls)
print(frequency[peaks_ls])
print(1./frequency[peaks_ls]/86400.)

f1 = figure(title="full time span: power vs frequency",
                        x_axis_label='period (days)', y_axis_label='power',
                        width=750)
f1.line(freq_days, power, line_width=2)

timespan = data.TIME[-1] - data.TIME[0]
yr_bins = 16
delta = timespan / yr_bins    # ~ 1 yr
yr_figs = []
yr_orb_power = []

for yr in range(yr_bins):
    tmin = data.TIME[0] + yr * delta
    tmax = tmin + delta

    yr_ind = [i for i, t in enumerate(d_non_zero_times) if (t >= tmin and t < tmax)]
    yr_times = d_non_zero_times[yr_ind]
    yr_weights = d_weighted[yr_ind]
    sum_weights = np.sum(yr_weights)

    yr_power = LombScargle(yr_times, yr_weights).power(frequency)
    yr_figs.append(figure(title="time span: " + str(yr) + " " + str(tmin) + ", " + str(tmax) + " power vs frequency",
                        x_axis_label='period (days)', y_axis_label='power',
                        width=750))
    yr_figs[yr].line(freq_days, yr_power, line_width=2)

    yr_peaks_ls, yr_props_ls = find_peaks(yr_power, height=0.1 * max(yr_power))
    print("yr", yr, "tmin", tmin, "tmax", tmax, "sum_weights", sum_weights, "num bins", len(yr_weights))
    print(yr_peaks_ls, yr_props_ls)
    print(frequency[yr_peaks_ls])
    print(1. / frequency[yr_peaks_ls] / 86400.)

    orb_power = 0.
    for p in np.arange(len(yr_peaks_ls)):
        i = yr_peaks_ls[p]
        if abs(frequency[i] - nom_freq) < 0.05 * nom_freq:
            orb_power = yr_props_ls["peak_heights"][p]
            print(i, frequency[i], orb_power)
            break
    yr_orb_power.append(orb_power)

f2 = figure(title="orb period power vs time bin",
                        x_axis_label='time bin', y_axis_label='power',
                        width=750)
f2.line(np.arange(yr_bins), yr_orb_power, line_width=2)
f2.circle(np.arange(yr_bins), yr_orb_power, size=6, fill_color="white")


del_div = Div(text=source + " Run on: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

#output_file("LS5039.html")
output_file(html_name)
l = layout(column(del_div, f1, f2, q_hist, r_hist, column(yr_figs)))
save(l, title=source + " Power vs Frequency")

# Close the FITS file
hdul.close()

# Print the data and header
print(header)
