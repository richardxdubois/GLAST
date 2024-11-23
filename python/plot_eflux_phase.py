from astropy.io import fits
import numpy as np
import argparse
import yaml
from datetime import datetime

from fit_SED import SED_function, fit_SED

from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column, gridplot
from bokeh.models import Label, Span, LinearAxis, Range1d, Whisker, ColumnDataSource
from bokeh.models.widgets import Div


class plot_eflux_phase():

    def __init__(self, app_config):


        with open(app_config, "r") as f:
            data = yaml.safe_load(f)

        self.source_name = data["source"]
        self.fgl_source = data["4FGL_source"]

        self.num_pickles = data["num_pickles"]
        self.p_bins = np.arange(self.num_pickles)
        self.html = data["html"]
        self.page_title = data["page_title"]
        self.fig_height = data["fig_height"]
        self.fig_width = data["fig_width"]

        self.base_fn = data["base_fn"]
        self.sed_prefix = data["sed_prefix"]

        self.seds = {}
        self.eflux = None
        self.eflux_err_lo = None
        self.eflux_err_hi = None
        self.loge_ctr = None
        self.b_lower = None
        self.b_upper = None

        self.type_1 = None
        self.type_2 = None

        self.num_pickles_2 = 0
        self.base_fn_2 = None
        self.p_bins_2 = []
        self.failed_fits = 0
        self.initial_guesses = [1e-4, 1., 1e3, 1e3]
        self.backup_guesses = [1e-4, -0.5, 1e3, 1e3]

        try:
            self.num_pickles_2 = data["num_pickles_2"]
            self.base_fn_2 = data["base_fn_2"]
            self.p_bins_2 = np.arange(self.num_pickles_2)
            self.type_1 = data["type_1"]
            self.type_2 = data["type_2"]
        except KeyError:
            pass

    def loop_over_bins(self):

        for phase_bin in self.p_bins:

            self.seds[phase_bin] = []

            for phase_bin_2 in self.p_bins_2:

                if self.num_pickles_2 == 0:
                    infile = self.base_fn + str(phase_bin) + "/" + self.sed_prefix + \
                              str(phase_bin) + "_" + self.fgl_source + "_sed.npy"
                else:
                    infile = (self.base_fn + str(phase_bin) + "/" + self.base_fn_2 +
                              str(phase_bin_2) + "/" + self.sed_prefix + str(phase_bin_2) + "_" +
                              self.fgl_source + "_sed.npy")

                print("working on ",infile)
                p = np.load(infile, allow_pickle=True).flat[0]

                self.eflux = p["e2dnde"]
                self.eflux_err_lo = p["e2dnde_err_lo"]
                self.eflux_err_hi = p["e2dnde_err_hi"]
                self.loge_ctr = np.power(10, p["loge_ctr"])

                self.b_upper = [x+e for x, e in zip(self.eflux, self.eflux_err_hi)]
                self.b_lower = [x-e for x, e in zip(self.eflux, self.eflux_err_lo)]

                rc = self.make_plot(phase_bin1=phase_bin, phase_bin2=phase_bin_2)

    def make_plot(self, phase_bin1, phase_bin2=None):

        source = ColumnDataSource(data=dict(x=self.loge_ctr, y=self.eflux, upper=self.b_upper,
                                            lower=self.b_lower))

        title = self.type_1 + " Phase bin " + str(phase_bin1)
        if phase_bin2 is not None:
            title += " " + self.type_2 + " Phase bin " + str(phase_bin2)
        title += " : E^2 dN/dE vs Energy"

        # Log-log plot
        p_fig = figure(x_axis_type="log", y_axis_type="log", height=self.fig_height, width=self.fig_width,
                       title=title, x_axis_label='E (MeV)',
                       y_axis_label='E^2 dN/dE', tooltips=[('Eflux', '@y'), ('Energy', '@x')], x_range=(100, 100000),
                       y_range=(5e-7, 5.e-4))

        #p_fig.line(source=source, x="x", y="y", legend_label="E^2 dN/dE", line_width=2)
        p_fig.scatter(self.loge_ctr, self.eflux, size=8, fill_color="white")
        p_fig.add_layout(Whisker(source=source, base="x", upper="upper",
                                 lower="lower", level="overlay"))

        flux = []
        errors = []
        E = []

        for i, E_i in enumerate(self.loge_ctr):
            if E_i <= 10000.:
                flux.append(self.eflux[i])
                errors.append(self.eflux_err_hi[i])
                E.append(E_i)

        print("Processing super bin ", phase_bin1, "orb bin ", phase_bin2)

        try:
            params, covariance = fit_SED(E, flux, errors, self.initial_guesses)
            print("Fitted parameters:", params)
            # Generate data for the fit line
            E_fit = np.linspace(1e2, 1e4, 100)  # Energy range for the fit
            flux_fit = SED_function(E_fit, *params)  # Calculate the fitted flux
            p_fig.line(E_fit, flux_fit, color='red', legend_label='Fitted Model')
        except RuntimeError:
            print("fit failed")
            self.failed_fits += 1
            pass

        self.seds[phase_bin1].append(p_fig)

    def output_plot(self):

        del_div = Div(text=self.source_name + " Run on: " +
                           datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        output_file(self.html)

        # Prepare Bokeh figures for grid layout
        plots = []
        for i in self.p_bins:
            plots.append(self.seds[i])
        grid = gridplot(plots)

        l = layout(del_div, grid)
        save(l, title=self.page_title)

if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser(description='plot SEDs')

    parser.add_argument('--app_config',
                        default="/Users/richarddubois/Code/GLAST/tmp/dbg/plot_eflux_config.yaml",
                        help="overall app config file")

    args = parser.parse_args()

    p = plot_eflux_phase(args.app_config)
    rc = p.loop_over_bins()
    print("# failed fits", p.failed_fits)

    rc = p.output_plot()
