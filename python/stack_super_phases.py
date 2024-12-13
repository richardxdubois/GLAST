import numpy as np
from datetime import datetime
import pickle
import argparse
import yaml

from bokeh.models.widgets import Div, NumberFormatter
from bokeh.models.formatters import DatetimeTickFormatter
from bokeh.models import ColumnDataSource, Span, LinearAxis, Range1d, LinearColorMapper, Whisker
from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column

parser = argparse.ArgumentParser()

parser.add_argument('--app_config', default="", help="overall app config file")
args = parser.parse_args()


with open(args.app_config, "r") as f:
    data = yaml.safe_load(f)


source_name = data["source_name"]
source_4fgl = data["source_4fgl"]

num_pickles = data["num_pickles"]

p_bins = np.arange(num_pickles)
print(num_pickles, p_bins)

base_fn = data["base_fn"]
base_super = data["base_super"]

plot_kind = data["plot_kind"]

try:
    pickle_base = data["pickle_base"]
except KeyError:
    pickle_base = "pickle_"

html_file = data["html_file"]

p_bins = np.arange(num_pickles)
q_bins = p_bins
q_bins = np.append(q_bins, p_bins)
r_bins = np.arange(len(q_bins))
dict_ticker = {}
for i in r_bins:
    dict_ticker[r_bins[i]] = str(q_bins[i]/float(num_pickles))


d_hist = []

all_x = []
all_y = []
all_flux = []

for phase_bin in np.arange(num_pickles):

    norms = []
    norms_errors = []
    alphas = []
    alphas_errors = []
    betas = []
    betas_errors = []
    npreds = []
    tss = []
    fluxs = []
    fluxs_errors = []
    efluxs = []
    efluxs_errors = []

    for super_bin in np.arange(num_pickles):

            p_name = (base_fn + str(phase_bin) + "/" +
                      base_super + str(super_bin) + "/" + pickle_base + str(super_bin) + ".npy")

            p = np.load(p_name, allow_pickle=True).flat[0]

            #LSI = p["sources"]['4FGL J0240.5+6113']
            LSI = p["sources"][source_4fgl]
            p_values = LSI["param_values"]
            p_errors = LSI["param_errors"]

            norm = p_values[0]
            norm_error = p_errors[0]
            alpha = p_values[1]
            alpha_error = p_errors[1]
            beta = p_values[2]
            beta_error = p_errors[2]
            flux = LSI["flux"]
            flux_error =LSI["flux_err"]
            eflux = LSI["eflux"]
            eflux_error =LSI["eflux_err"]
            npred = LSI["npred"]
            ts = LSI["ts"]

            norms.append(norm)
            norms_errors.append(norm_error)
            alphas.append(alpha)
            alphas_errors.append(alpha_error)
            betas.append(beta)
            betas_errors.append(beta_error)
            fluxs.append(flux)
            fluxs_errors.append(flux_error)
            efluxs.append(eflux)
            efluxs_errors.append(eflux_error)
            npreds.append(npred)
            tss.append(ts)

            print(norm, norm_error)

    fluxs.extend(fluxs)
    fluxs_errors.extend(fluxs_errors)

    title = "flux - " + plot_kind + " phase bin:" + str(phase_bin)
    u_hist = figure(title=title, x_axis_label='Phase bin', width=750)
    #u_hist.line(p_bins, fluxs, line_width=2)
    color = "red"
    if phase_bin < 5:
        color = "green"
    u_hist.vbar(top=fluxs, x=r_bins, width=1., fill_color=color, fill_alpha=0.05, bottom=0)
    u_hist.circle(r_bins, fluxs, size=6, fill_color="white")
    f_upper = [x+e for x,e in zip(fluxs, fluxs_errors)]
    f_lower = [x-e for x,e in zip(fluxs, fluxs_errors)]
    f_source = ColumnDataSource(data=dict(groups=r_bins, counts=fluxs, upper=f_upper, lower=f_lower))
    u_hist.add_layout(Whisker(source=f_source, base="groups", upper="upper", lower="lower", level="overlay"))
    u_hist.y_range.start = 1.e-8
    u_hist.y_range = Range1d(0., 1.2e-6)
    
    if phase_bin != 9:
        u_hist.xaxis.visible = False  # Hide x-axis on the first plot
        u_hist.xgrid.visible = False  # Optionally hide x-grid lines
        u_hist.min_border_top = 0  # Remove top border padding
        u_hist.min_border_bottom = 0  # Remove bottom border padding
    else:
        u_hist.xaxis.ticker = r_bins
        u_hist.xaxis.major_label_overrides = dict_ticker

    d_hist.append(u_hist)
    all_flux.extend(fluxs)
    all_x.extend(r_bins)
    
for r in range(num_pickles):
    all_y.extend(np.full(2*num_pickles, r))

# Open a file in binary write mode
with open('heatmap.pkl', 'wb') as file:
    # Dump the data into the file
    pickle.dump([all_x, all_y, all_flux], file)
    
output_file(html_file)
del_div = Div(text=source_name + " Run on: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

l = layout(del_div, column(d_hist, sizing_mode="stretch_width", spacing=0))
save(l, title="LS super sliced fit params")
print("done")
