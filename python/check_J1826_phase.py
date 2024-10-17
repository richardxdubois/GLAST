from astropy.io import fits
import numpy as np
import argparse


from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column
from bokeh.models.widgets import Div

# Command line arguments
parser = argparse.ArgumentParser(description='run phased analysis')

parser.add_argument('--infile', default='/Users/richarddubois/Code/GLAST/tmp/LS5039_J1826_5deg.fits',
                    help="events file")

args = parser.parse_args()


infile = args.infile

h = fits.open(infile)
h.info()

primary_hdu = h[0]
# Get the data and header
print(h[1].columns)

data = h[1].data
pulse_phase = data.PULSE_PHASE

print(len(pulse_phase))

q_hist = figure(title="Pulse Phase",
                x_axis_label='Phase', y_axis_label='counts',
                width=750)

# counts_hist, counts_edges = np.histogram(d_non_zero_times, range=(2780, 3580), bins=1000)
counts_hist, counts_edges = np.histogram(pulse_phase, bins=200)
q_hist.vbar(top=counts_hist, x=counts_edges[1:], width=counts_edges[1] - counts_edges[0], fill_color='red',
            fill_alpha=0.2, bottom=0)

output_file("LS5039_J1826_phase.html")
l = layout(q_hist)
save(l, title=" LS5039-J1826 Phase selected times")

