from astropy.io import fits
from collections import OrderedDict
import numpy as np
from datetime import datetime, date, timedelta
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

from bokeh.models.widgets import Tabs, Panel, DataTable, TableColumn, Div, NumberFormatter
from bokeh.models.formatters import DatetimeTickFormatter
from bokeh.models import ColumnDataSource, Span, LinearAxis, Range1d, LinearColorMapper
from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column

source = "LSI61303"
html_name = source + "_phase.html"
infile = OrderedDict()

# Open the FITS file
if source == "LSI61303":
    infile['/Users/richarddubois/Code/GLAST/tmp/dbg/L24082417075904476C3F57_PH00.fits'] = []
    #infile['/Users/richarddubois/Code/GLAST/tmp/LSI61303_1_deg_500s.fits'] = []
    infile['/Users/richarddubois/Code/GLAST/tmp/dbg/LSI61303_file_3_phase_1.fits'] = []
    infile['/Users/richarddubois/Code/GLAST/tmp/dbg/LSI61303_file_0_phase_0.fits'] = []
    infile['/Users/richarddubois/Code/GLAST/tmp/dbg/ft1_00-00.fits'] = []
    infile['/Users/richarddubois/Code/GLAST/tmp/dbg/ft1_00-01.fits'] = []
    infile['/Users/richarddubois/Code/GLAST/tmp/dbg/ft1_00-02.fits'] = []
    infile['/Users/richarddubois/Code/GLAST/tmp/dbg/ft1_00-09.fits'] = []
elif source == "LS5039":
    #infile = fits.open('/Users/richarddubois/Code/GLAST/tmp/LS5039_1_deg_500s.fits')
    infile.append('/Users/richarddubois/Code/GLAST/tmp/LS5039_3_deg_10800s.fits')


for f in infile:

    del_div = Div(text="Run on: " + datetime.now().strftime("%Y-%m-%d") + " for " + f, width=400)
    infile[f].append(del_div)

    # Print information about the FITS file
    h = fits.open(f)
    h.info()

    # Access the primary HDU (Header/Data Unit)
    primary_hdu = h[0]
    gti = h[2]

    # Get the data and header
    print(h[1].columns)
    print(gti.columns)
    header = primary_hdu.header
    data = h[1].data

    d_non_zero_times = data.TIME / 86400.
    energy = data.ENERGY
    ra = data.RA
    dec = data.DEC
    zenith_angle = data.ZENITH_ANGLE

    print(len(d_non_zero_times))

    p_hist = figure(title="Time",
                    x_axis_label='Time (days MET)', y_axis_label='counts',
                    width=750)
    infile[f].append(p_hist)

    t_min = min(d_non_zero_times) - 30.
    counts_hist, counts_edges = np.histogram(d_non_zero_times, range=(t_min, t_min+300.), bins=100)
    p_hist.vbar(top=counts_hist, x=counts_edges[1:], width=counts_edges[1]-counts_edges[0], fill_color='red',
                fill_alpha=0.2, bottom=0)

    q_hist = figure(title="Time",
                    x_axis_label='Time (days MET)', y_axis_label='counts',
                    width=750)
    infile[f].append(q_hist)

    #counts_hist, counts_edges = np.histogram(d_non_zero_times, range=(2780, 3580), bins=1000)
    counts_hist, counts_edges = np.histogram(d_non_zero_times, bins=200)
    q_hist.vbar(top=counts_hist, x=counts_edges[1:], width=counts_edges[1]-counts_edges[0], fill_color='red',
                fill_alpha=0.2, bottom=0)

    r_hist = figure(title="Energy",
                    x_axis_label='Energy (MeV)', y_axis_label='counts', x_axis_type='log',
                    width=750)
    infile[f].append(r_hist)

    #counts_hist, counts_edges = np.histogram(d_non_zero_times, range=(2780, 3580), bins=1000)
    e_hist, e_edges = np.histogram(energy, bins=100)
    r_hist.vbar(top=e_hist, x=e_edges[1:], width=e_edges[1]-e_edges[0], fill_color='red',
                fill_alpha=0.2, bottom=0)

    s_hist = figure(title="RA",
                    x_axis_label='RA (deg)', y_axis_label='counts',
                    width=750)
    infile[f].append(s_hist)

    #counts_hist, counts_edges = np.histogram(d_non_zero_times, range=(2780, 3580), bins=1000)
    ra_hist, ra_edges = np.histogram(ra, bins=100)
    s_hist.vbar(top=ra_hist, x=ra_edges[1:], width=ra_edges[1]-ra_edges[0], fill_color='red',
                fill_alpha=0.2, bottom=0)

    t_hist = figure(title="DEC",
                    x_axis_label='DEC (deg)', y_axis_label='counts',
                    width=750)
    infile[f].append(t_hist)

    #counts_hist, counts_edges = np.histogram(d_non_zero_times, range=(2780, 3580), bins=1000)
    dec_hist, dec_edges = np.histogram(dec, bins=100)
    t_hist.vbar(top=dec_hist, x=dec_edges[1:], width=dec_edges[1]-dec_edges[0], fill_color='red',
                fill_alpha=0.2, bottom=0)

    u_hist = figure(title="Zenith",
                    x_axis_label='Zenith (deg)', y_axis_label='counts',
                    width=750)
    infile[f].append(u_hist)


    #counts_hist, counts_edges = np.histogram(d_non_zero_times, range=(2780, 3580), bins=1000)
    z_hist, z_edges = np.histogram(zenith_angle, bins=100)
    u_hist.vbar(top=z_hist, x=z_edges[1:], width=z_edges[1]-z_edges[0], fill_color='red',
                fill_alpha=0.2, bottom=0)

    h.close()

output_file(html_name)
l = layout(children=infile.values())
save(l, title=source + " Phase selected times")
