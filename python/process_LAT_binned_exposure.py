from astropy.io import fits
import numpy as np
from datetime import datetime, date, timedelta
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

from bokeh.models.widgets import DataTable, TableColumn, Div, NumberFormatter
from bokeh.models.formatters import DatetimeTickFormatter
from bokeh.models import ColumnDataSource, Span, LinearAxis, Range1d, LinearColorMapper
from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column

source = "LSI61303"
html_name = source + ".html"
fn = ""

# Open the FITS file
if source == "LSI61303":
    fn = '/Users/richarddubois/Code/GLAST/tmp/LSI61303_1_deg_mkt_500s.fits'
    #fn ='/Users/richarddubois/Code/GLAST/tmp/LSI61303_1_deg_mkt_225000s.fits'
    #fn = '/Users/richarddubois/Code/GLAST/tmp/LSI61303_3_deg_mkt_10800s.fits'

    hdul = fits.open(fn)

    f_start = 1./27.5/86400.  # 40.
    f_stop = 1./23.5/86400.   # 5.
    nom_period = 26.496
    nom_freq = 1/nom_period/86400.

elif source == "LS5039":
    #hdul = fits.open('/Users/richarddubois/Code/GLAST/tmp/LS5039_1_deg_500s.fits')
    hdul = fits.open('/Users/richarddubois/Code/GLAST/tmp/LS5039_3_deg_10800s.fits')
    f_start = 1./10./86400.
    f_stop = 1./2./86400.
    nom_period = 3.9
    nom_freq = 1/nom_period/86400.

# Print information about the FITS file
hdul.info()

# Access the primary HDU (Header/Data Unit)
primary_hdu = hdul[0]

# Get the data and header
print(hdul[1].columns)
header = primary_hdu.header
data = hdul[1].data

indices = np.where(data.COUNTS > 0.)
timedel = data.TIMEDEL

suppress_zero = True   # suppress entries with zero counts
do_weighting = True

if suppress_zero:
    d_non_zero_times = data.TIME[indices] + timedel[0]/2.
    exp = data.EXPOSURE[indices]
    d_weighted = data.COUNTS[indices]
else:
    d_non_zero_times = data.TIME + timedel[0]/2.
    exp = data.EXPOSURE
    d_weighted = data.COUNTS

r_weighted = []
weights = []

mean_counts = np.average(d_weighted)
mean_rate = np.sum(d_weighted) / np.sum(exp)
mean_sub = d_weighted - mean_counts

# best weighting so far:
#  ave - mean_rate*np.sqrt(e)
#  r_weighted.append(d_weighted[i]/ave)

for i, e in enumerate(exp):
    if e == 0.:
        r_weighted.append(0.)
        print("skipped over zero exposure at ", i, d_weighted[i])
    else:
        ave = mean_rate*np.sqrt(e)  # "predicted" count error - optimal
        if not do_weighting:
            ave = 1.
        #r_weighted.append(mean_sub[i] / e)
        r_weighted.append(d_weighted[i]/ave)
        weights.append(ave)

r_weighted = np.array(r_weighted)
weights = np.array(weights)
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

f3 = figure(title="Exposure vs Counts",  x_axis_label='Exposure', y_axis_label='Counts', width=750)
f3.scatter(x=exp, y=d_weighted, size=3)

t_start = data.TIME[0]
t_end = data.TIME[-1]

print("t_start", t_start, "t_end", t_end, "num orig bins", len(data.TIME), "num non-zero bins", len(d_weighted))

frequency = np.linspace(f_start, f_stop, 10000)  # for orbital 100000
freq_days = 1./frequency/86400.

power = LombScargle(t=d_non_zero_times, y=r_weighted).power(frequency)
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
vline_p1 = Span(location=nom_period, dimension='height', line_color='red', line_width=2, line_dash='dashed')
f1.add_layout(vline_p1)

timespan = data.TIME[-1] - data.TIME[0]
yr_bins = 4
delta = timespan / yr_bins    # ~ 1 yr
yr_figs = []
yr_exp = []
yr_c = []
yr_orb_power = []

for yr in range(yr_bins):
    tmin = data.TIME[0] + yr * delta
    tmax = tmin + delta

    yr_ind = [i for i, t in enumerate(d_non_zero_times) if (t >= tmin and t < tmax)]
    yr_times = d_non_zero_times[yr_ind]
    yr_counts_d = d_weighted[yr_ind]
    yr_counts = r_weighted[yr_ind]
    yr_weights = weights[yr_ind]
    yr_exposure = exp[yr_ind]
    sum_weights = np.sum(yr_weights)

    yr_power = LombScargle(t=yr_times, y=yr_counts).power(frequency)
    yr_figs.append(figure(title="time span: " + str(yr) + " " + str(tmin) + ", " + str(tmax) + " power vs frequency",
                        x_axis_label='period (days)', y_axis_label='power',
                        width=750))
    yr_figs[yr].line(freq_days, yr_power, line_width=2)
    yr_figs[yr].add_layout(vline_p1)

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

    counts_hist, counts_edges = np.histogram(yr_counts_d, bins=100)

    yr_c.append(figure(title="Counts",
                    x_axis_label='Counts', y_axis_label='counts',
                    width=750))

    yr_c[yr].vbar(top=counts_hist, x=counts_edges[1:], width=counts_edges[1] - counts_edges[0], fill_color='red',
                fill_alpha=0.2, bottom=0)

    exp_hist, exp_edges = np.histogram(yr_exposure, bins=100)
    yr_exp.append(figure(title="Exposure",
                    x_axis_label='cm^2 s', y_axis_label='counts',
                    width=750))

    yr_exp[yr].vbar(top=exp_hist, x=exp_edges[1:], width=exp_edges[1] - exp_edges[0], fill_color='red', fill_alpha=0.2,
                bottom=0)


f2 = figure(title="orb period power vs time bin",
                        x_axis_label='time bin', y_axis_label='power',
                        width=750)
f2.line(np.arange(yr_bins), yr_orb_power, line_width=2)
f2.scatter(np.arange(yr_bins), yr_orb_power, size=6, fill_color="white")


del_div = Div(text=source + " Run on: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + " for " + fn)

#output_file("LS5039.html")
output_file(html_name)
l = layout(column(del_div, f1, f2, q_hist, r_hist, f3, column(yr_figs), column(yr_c), column(yr_exp)))
save(l, title=source + " Power vs Frequency")

# Close the FITS file
hdul.close()

# Print the data and header
print(header)
