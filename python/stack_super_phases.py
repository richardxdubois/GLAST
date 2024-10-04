import numpy as np
from datetime import datetime

from bokeh.models.widgets import Div, NumberFormatter
from bokeh.models.formatters import DatetimeTickFormatter
from bokeh.models import ColumnDataSource, Span, LinearAxis, Range1d, LinearColorMapper, Whisker
from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column

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

source_name = "LSI61303"

num_pickles = 10
p_bins = np.arange(num_pickles)
print(num_pickles, p_bins)

d_hist = []

for phase_bin in np.arange(num_pickles):

    for super_bin in np.arange(num_pickles):

            p_name = ("/Users/richarddubois/Code/GLAST/tmp/dbg/phase" + str(phase_bin) +
                      "super_bin" + str(super_bin) + ".npy")

            p = np.load(p_name, allow_pickle=True).flat[0]

            LSI = p["sources"]['4FGL J0240.5+6113']
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

    norms.extend(norms)
    norms_errors.extend(norms_errors)
    q_bins = p_bins
    q_bins = np.append(q_bins, p_bins)
    r_bins = np.arange(len(q_bins))
    dict_ticker = {}
    for i in r_bins:
        dict_ticker[str(r_bins[i])] = str(q_bins[i])

    u_hist = figure(title="flux", x_axis_label='Phase bin', width=750)
    #u_hist.line(p_bins, fluxs, line_width=2)
    color = "red"
    if phase_bin < 5:
        color = "green"
    u_hist.vbar(top=fluxs, x=p_bins, width=1., fill_color=color, fill_alpha=0.05, bottom=0)
    u_hist.circle(p_bins, fluxs, size=6, fill_color="white")
    f_upper = [x+e for x,e in zip(fluxs, fluxs_errors)]
    f_lower = [x-e for x,e in zip(fluxs, fluxs_errors)]
    f_source = ColumnDataSource(data=dict(groups=p_bins, counts=fluxs, upper=f_upper, lower=f_lower))
    u_hist.add_layout(Whisker(source=f_source, base="groups", upper="upper", lower="lower", level="overlay"))
    u_hist.y_range.start = 1.e-8
    vline_p = Span(location=2.3, dimension='height', line_color='red', line_width=2, line_dash='dashed')
    vline_a = Span(location=7.75, dimension='height', line_color='blue', line_width=2, line_dash='dashed')
    u_hist.add_layout(vline_p)
    u_hist.add_layout(vline_a)
    if phase_bin != 9:
        u_hist.xaxis.visible = False  # Hide x-axis on the first plot
        u_hist.xgrid.visible = False  # Optionally hide x-grid lines
        u_hist.min_border_top = 0  # Remove top border padding
        u_hist.min_border_bottom = 0  # Remove bottom border padding
    d_hist.append(u_hist)

output_file("LS_params_super_sliced.html")
del_div = Div(text=source_name + " Run on: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

l = layout(del_div, column(d_hist, sizing_mode="stretch_width", spacing=0))
save(l, title="LS super sliced fit params")
print("done")
