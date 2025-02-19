import numpy as np
from datetime import datetime
import yaml
import argparse

from bokeh.models.widgets import Div, NumberFormatter
from bokeh.models.formatters import DatetimeTickFormatter
from bokeh.models import ColumnDataSource, Span, LinearAxis, Range1d, LinearColorMapper, Whisker
from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column

super_base = "/sdf/home/r/richard/fermi-user/LSI61303/periods/phased_1/analyses/super_only_p/super_only_1/"

super_bins = "phase"
orb_bins = "orb_bins"
pickle = "pickle_"

super_flux_1 = np.zeros(10)
super_error_1 = np.zeros(10)
super_flux_2 = np.zeros(10)
super_error_2 = np.zeros(10)

for s in range(0,10):
    for o in range(0,10):
        fn = super_base + super_bins + str(s) + "/" + orb_bins + str(o) + "/" + pickle + str(o) + ".npy"

        p = np.load(fn, allow_pickle=True).flat[0]

        LSI = p["sources"]["4FGL J0240.5+6113"]

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
        flux_error = LSI["flux_err"]
        eflux = LSI["eflux"]
        eflux_error = LSI["eflux_err"]
        npred = LSI["npred"]
        ts = LSI["ts"]

        if s < 3:
            super_flux_1[o] += flux
            super_error_1[o] += flux_error**2
        else:
            super_flux_2[o] += flux
            super_error_2[o] += flux_error**2

super_error_1 = np.sqrt(super_error_1)
super_error_2 = np.sqrt(super_error_2)

bins = np.arange(10)

tooltips = [('phase', 'orbital: @groups'), ('flux', '@counts')]

f_upper = [x+e for x,e in zip(super_flux_1, super_error_1)]
f_lower = [x-e for x,e in zip(super_flux_1, super_error_1)]
f_source = ColumnDataSource(data=dict(groups=bins, counts=super_flux_1, upper=f_upper, lower=f_lower))

u_hist = figure(title="super 0-2 flux vs phase", x_axis_label='Orbital Phase', width=750, tooltips=tooltips)

u_hist.vbar(top="counts", x="groups", width=1., fill_color='red', fill_alpha=0.05, bottom=0)
u_hist.scatter(bins, super_flux_1, size=6, fill_color="white")

u_hist.add_layout(Whisker(source=f_source, base="groups", upper="upper", lower="lower", level="overlay"))

v_hist = figure(title="super 3-9 flux vs phase", x_axis_label='Orbital Phase', width=750, tooltips=tooltips)

v_upper = [x+e for x,e in zip(super_flux_2, super_error_2)]
v_lower = [x-e for x,e in zip(super_flux_2, super_error_2)]
v_source = ColumnDataSource(data=dict(groups=bins, counts=super_flux_2, upper=v_upper, lower=v_lower))

v_hist.vbar(top="counts", x="groups", width=1., fill_color='red', fill_alpha=0.05, bottom=0)
v_hist.scatter(bins, super_flux_2, size=6, fill_color="white")

v_hist.add_layout(Whisker(source=v_source, base="groups", upper="upper", lower="lower", level="overlay"))

output_file("LSI61303_super_sums.html")
del_div = Div(text="LSI61303 Run on: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

l = layout(del_div, u_hist, v_hist)
#l = layout(del_div, u_hist, v_hist)
save(l, title="LSI61303 super sums fit params")
