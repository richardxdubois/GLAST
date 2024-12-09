import numpy as np
from datetime import datetime
import yaml
import argparse

from bokeh.models.widgets import Div, NumberFormatter
from bokeh.models.formatters import DatetimeTickFormatter
from bokeh.models import ColumnDataSource, Span, LinearAxis, Range1d, LinearColorMapper, Whisker
from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column


def shift_list(lst, n):
    n = n % len(lst)  # Ensure n is within the bounds of the list length
    return lst[n:] + lst[:n]


parser = argparse.ArgumentParser()

parser.add_argument('--app_config',
                    default="/Users/richarddubois/Code/GLAST/tmp/dbg/fetch_fpy_config.yaml",
                    help="overall app config file")
args = parser.parse_args()


with open(args.app_config, "r") as f:
    data = yaml.safe_load(f)

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

source_name = data["source"]

num_pickles = data["num_pickles"]

try:
    phase_offset = int(data["phase_offset"])
except KeyError:
    phase_offset = 0

p_bins = np.arange(num_pickles)
print(num_pickles, p_bins)

#super = True
super = data["super"]

base_fn = data["base_fn"]

try:
    pickle_base = data["pickle_base"]
except KeyError:
    pickle_base = "pickle_"

for phase_bin in np.arange(num_pickles):

    p_name = base_fn + str(phase_bin) + "/" + pickle_base + str(phase_bin) + ".npy"

    p = np.load(p_name, allow_pickle=True).flat[0]

    LSI = p["sources"][data["4FGL_source"]]

    p_values = LSI["param_values"]
    p_errors = LSI["param_errors"]
    norm = p_values[0]
    norm_error = p_errors[0]
    alpha = p_values[1]
    alpha_error = p_errors[1]

    try:
        beta = p_values[2]
        beta_error = p_errors[2]
    except KeyError:  # PL fit
        beta = -999.
        beta_error = -999.

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

if phase_offset != 0:
    fluxs = shift_list(fluxs, phase_offset)
    fluxs_errors = shift_list(fluxs_errors, phase_offset)
    alphas = shift_list(alphas, phase_offset)
    betas = shift_list(betas, phase_offset)
    efluxs = shift_list(efluxs, phase_offset)
    efluxs_errors = shift_list(fluxs_errors, phase_offset)
    npreds = shift_list(npreds, phase_offset)
    tss = shift_list(tss, phase_offset)

fluxs.extend(fluxs)
fluxs_errors.extend(fluxs_errors)
q_bins = p_bins
q_bins = np.append(q_bins, p_bins)
r_bins = np.arange(len(q_bins))
dict_ticker = {}
for i in r_bins:
    dict_ticker[r_bins[i]] = str(q_bins[i]/num_pickles)

q_hist = figure(title="N0", x_axis_label='Phase bin', y_axis_label='ph/MeV/cm^2/s', width=750)
q_hist.line(p_bins, norms, line_width=2)
q_hist.circle(p_bins, norms, size=6, fill_color="white")
upper = [x+e for x,e in zip(norms, norms_errors)]
lower = [x-e for x,e in zip(norms, norms_errors)]
source = ColumnDataSource(data=dict(groups=p_bins, counts=norms, upper=upper, lower=lower))
q_hist.add_layout(Whisker(source=source, base="groups", upper="upper", lower="lower", level="overlay"))
q_hist.y_range.start = 1.e-12

r_hist = figure(title="Alpha", x_axis_label='Phase bin', width=750)
r_hist.line(p_bins, alphas, line_width=2)
r_hist.circle(p_bins, alphas, size=6, fill_color="white")
a_upper = [x+e for x,e in zip(alphas, alphas_errors)]
a_lower = [x-e for x,e in zip(alphas, alphas_errors)]
a_source = ColumnDataSource(data=dict(groups=p_bins, counts=alphas, upper=a_upper, lower=a_lower))
r_hist.add_layout(Whisker(source=a_source, base="groups", upper="upper", lower="lower", level="overlay"))
r_hist.y_range.start = 0

s_hist = figure(title="Beta", x_axis_label='Phase bin', width=750)
s_hist.line(p_bins, betas, line_width=2)
s_hist.circle(p_bins, betas, size=6, fill_color="white")
b_upper = [x+e for x,e in zip(betas, betas_errors)]
b_lower = [x-e for x,e in zip(betas, betas_errors)]
b_source = ColumnDataSource(data=dict(groups=p_bins, counts=betas, upper=b_upper, lower=b_lower))
s_hist.add_layout(Whisker(source=b_source, base="groups", upper="upper", lower="lower", level="overlay"))
s_hist.y_range.start = 0

t_hist = figure(title="Npred (black) TS (red)", x_axis_label='Phase bin', width=750)
t_hist.line(p_bins, npreds, line_width=2)
t_hist.circle(p_bins, npreds, size=6, fill_color="white")
t_hist.line(p_bins, tss, line_width=2, color="red")
t_hist.circle(p_bins, tss, size=6, fill_color="white")
t_hist.y_range.start = 0

u_hist = figure(title="flux", x_axis_label='Phase', width=750)
#u_hist.line(p_bins, fluxs, line_width=2)
u_hist.vbar(top=fluxs, x=r_bins, width=1., fill_color='red', fill_alpha=0.05, bottom=0)
u_hist.circle(r_bins, fluxs, size=6, fill_color="white")
f_upper = [x+e for x,e in zip(fluxs, fluxs_errors)]
f_lower = [x-e for x,e in zip(fluxs, fluxs_errors)]
f_source = ColumnDataSource(data=dict(groups=r_bins, counts=fluxs, upper=f_upper, lower=f_lower))
u_hist.add_layout(Whisker(source=f_source, base="groups", upper="upper", lower="lower", level="overlay"))
u_hist.y_range.start = 1.e-8
u_hist.xaxis.ticker = r_bins
u_hist.xaxis.major_label_overrides = dict_ticker
u_hist.xaxis.major_label_orientation = 0.7
vline_p1 = Span(location=data["periastron"]*num_pickles, dimension='height', line_color='red', line_width=2,
                line_dash='dashed')
vline_a1 = Span(location=data["apastron"]*num_pickles, dimension='height', line_color='blue', line_width=2,
                line_dash='dashed')

vline_p2 = Span(location=(1.+data["periastron"])*num_pickles, dimension='height', line_color='red', line_width=2,
                line_dash='dashed')
vline_a2 = Span(location=(1.+data["apastron"])*num_pickles, dimension='height', line_color='blue', line_width=2,
                line_dash='dashed')

if not super:
    u_hist.add_layout(vline_p1)
    u_hist.add_layout(vline_a1)
    u_hist.add_layout(vline_p2)
    u_hist.add_layout(vline_a2)

v_hist = figure(title="eflux", x_axis_label='Phase bin', width=750)
v_hist.line(p_bins, efluxs, line_width=2)
v_hist.circle(p_bins, efluxs, size=6, fill_color="white")
ef_upper = [x+e for x,e in zip(efluxs, efluxs_errors)]
ef_lower = [x-e for x,e in zip(efluxs, efluxs_errors)]
ef_source = ColumnDataSource(data=dict(groups=p_bins, counts=efluxs, upper=ef_upper, lower=ef_lower))
v_hist.add_layout(Whisker(source=ef_source, base="groups", upper="upper", lower="lower", level="overlay"))
v_hist.y_range.start = 1.e-5

output_file(data["html"])
del_div = Div(text=source_name + " Run on: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

l = layout(del_div, u_hist, v_hist, q_hist, r_hist, s_hist, t_hist)
#l = layout(del_div, u_hist, v_hist)
save(l, title=source_name + " fit params")
print("done")
