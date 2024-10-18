from astropy.io import fits
import numpy as np
import argparse
import yaml
from datetime import datetime

from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column
from bokeh.models import Label, Span, LinearAxis, Range1d, Whisker, ColumnDataSource
from bokeh.models.widgets import Div

# Command line arguments
parser = argparse.ArgumentParser(description='plot SEDs')

parser.add_argument('--app_config',
                    default="/Users/richarddubois/Code/GLAST/tmp/dbg/plot_eflux_config.yaml",
                    help="overall app config file")

args = parser.parse_args()

with open(args.app_config, "r") as f:
    data = yaml.safe_load(f)

source_name = data["source"]
fgl_source = data["4FGL_source"]

num_pickles = data["num_pickles"]
p_bins = np.arange(num_pickles)
html = data["html"]
page_title = data["page_title"]

print(num_pickles, p_bins)

base_fn = data["base_fn"]

seds = []

for phase_bin in np.arange(num_pickles):

    infile = base_fn + str(phase_bin) + "/pickle_" + str(phase_bin) + "_" + fgl_source + "_sed.npy"

    p = np.load(infile, allow_pickle=True).flat[0]

    eflux = p["e2dnde"]
    eflux_err_lo = p["e2dnde_err_lo"]
    eflux_err_hi = p["e2dnde_err_hi"]
    loge_ctr = np.power(10, p["loge_ctr"])

    b_upper = [x+e for x, e in zip(eflux, eflux_err_hi)]
    b_lower = [x-e for x, e in zip(eflux, eflux_err_lo)]

    source = ColumnDataSource(data=dict(x=loge_ctr, y=eflux, upper=b_upper, lower=b_lower))

    # Log-log plot
    p_fig = figure(x_axis_type="log", y_axis_type="log",
                   title="Phase bin " + str(phase_bin) + ": E^2 dN/dE vs Energy", x_axis_label='E (MeV)',
                   y_axis_label='E^2 dN/dE', tooltips=[('Eflux', '@y'), ('Energy', '@x')], x_range=(100, 100000),
                   y_range=(5e-7, 5.e-4))

    p_fig.line(source=source, x="x", y="y", legend_label="E^2 dN/dE", line_width=2)
    p_fig.scatter(loge_ctr, eflux, size=8, fill_color="white")
    p_fig.add_layout(Whisker(source=source, base="x", upper="upper", lower="lower", level="overlay"))

    seds.append(p_fig)

del_div = Div(text=source_name + " Run on: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

output_file(html)
l = layout(column(seds))
save(l, title=page_title)
