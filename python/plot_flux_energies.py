from astropy.io import fits
import numpy as np
import pandas as pd
import argparse
import yaml
import pickle
from datetime import datetime
from scipy.integrate import quad

from fit_SED import SED_function, fit_SED, flux_function, fit_SED_errors

from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column, gridplot
from bokeh.models import (Label, Span, LinearAxis, Range1d, Whisker, ColumnDataSource, BasicTicker, Tabs,
                          TabPanel, RangeSlider, CustomJS, Button, NumeralTickFormatter, BasicTickFormatter)
from bokeh.models.widgets import Div
from bokeh.palettes import Plasma256 as palette
from bokeh.transform import linear_cmap, log_cmap


class plot_flux_energies():
    def __init__(self, app_config=None):

        with open(app_config, "r") as f:
            data = yaml.safe_load(f)
        # Read from the pickle file

        self.source_name = data["source"]
        self.fgl_source = data["4FGL_source"]

        self.num_pickles = data["num_pickles"]
        self.p_bins = np.arange(self.num_pickles)
        self.html = data["html"]
        self.hist_flux_max = data["hist_flux_max"]
        self.params_save_pickle = self.html.split(".")[0] + "_params.pkl"
        self.page_title = data["page_title"]
        self.fig_height = data["fig_height"]
        self.fig_width = data["fig_width"]

        self.base_fn = data["base_fn"]
        self.sed_prefix = data["sed_prefix"]

        self.num_pickles_2 = data["num_pickles_2"]
        self.base_fn_2 = data["base_fn_2"]
        self.p_bins_2 = np.arange(self.num_pickles_2)
        if self.num_pickles_2 == 0:
            self.p_bins_2 = np.arange(1)

        self.type_1 = data["type_1"]
        self.type_2 = data["type_2"]

        self.dict_ticker = {}
        self.phase_h = np.arange(2*self.num_pickles)
        for i, ph in enumerate(self.phase_h):
            self.dict_ticker[i] = str(self.phase_h[i] / self.num_pickles)

        self.orbital = {}
        self.orbital_100_300 = {}
        self.orbital_300_1000 = {}
        self.orbital_1000_10000 = {}

        self.orbital_per_super = {}
        self.orbital_per_super_100_300 = {}
        self.orbital_per_super_300_1000 = {}
        self.orbital_per_super_1000_10000 = {}

        self.super = {}
        self.super_100_300 = {}
        self.super_300_1000 = {}
        self.super_1000_10000 = {}

        self.super_per_orbital = {}
        self.super_per_orbital_100_300 = {}
        self.super_per_orbital_300_1000 = {}
        self.super_per_orbital_1000_10000 = {}

        self.no_energy_overlay = False
        self.energy_index = 0
        self.energy_index_flux = ["fluxs", "fluxs_100", "fluxs_300", "fluxs_1000"]

    def fill_maps(self):

        E_end_bin = [300, 1000, 10000]

        for s in self.p_bins:
            for o in self.p_bins:

                infile_b = (self.base_fn + str(s) + "/" + self.base_fn_2 +
                          str(o) + "/" + self.sed_prefix + str(o))
                infile = infile_b + "_" + self.fgl_source + "_sed.npy"

                p = np.load(infile, allow_pickle=True).flat[0]

                flux_E_bin = np.zeros(4)  # index 3 is total flux
                flux_errors_E_bin = np.zeros(4)
                index_e_bin = 0

                for i, E_i in enumerate(p["flux"]):
                    E_end = p["e_max"][i]
                    if E_end < E_end_bin[2] and not pd.isna(p["e2dnde_err_lo"][i]):
                        if E_end > E_end_bin[index_e_bin]:
                            index_e_bin += 1
                        flux_E_bin[index_e_bin] += p["flux"][i]
                        flux_errors_E_bin[index_e_bin] += p["flux_err"][i] ** 2

                        flux_E_bin[3] += p["flux"][i]
                        flux_errors_E_bin[3] += p["flux_err"][i] ** 2

                rc = self.orbital.setdefault(o, 0.)
                self.orbital[o] += flux_E_bin[3]/self.num_pickles

                rc = self.orbital_100_300.setdefault(o, 0.)
                self.orbital_100_300[o] += flux_E_bin[0]/self.num_pickles

                rc = self.orbital_300_1000.setdefault(o, 0.)
                self.orbital_300_1000[o] += flux_E_bin[1]/self.num_pickles

                rc = self.orbital_1000_10000.setdefault(o, 0.)
                self.orbital_1000_10000[o] += flux_E_bin[2]/self.num_pickles

                rc = self.super.setdefault(s, 0.)
                self.super[s] += flux_E_bin[3]/self.num_pickles

                rc = self.super_100_300.setdefault(s, 0.)
                self.super_100_300[s] += flux_E_bin[0]/self.num_pickles

                rc = self.super_300_1000.setdefault(s, 0.)
                self.super_300_1000[s] += flux_E_bin[1]/self.num_pickles

                rc = self.super_1000_10000.setdefault(s, 0.)
                self.super_1000_10000[s] += flux_E_bin[2]/self.num_pickles

                rc = self.orbital_per_super.setdefault(s, {})
                self.orbital_per_super[s][o] = flux_E_bin[3]

                rc = self.orbital_per_super_100_300.setdefault(s, {})
                self.orbital_per_super_100_300[s][o] = flux_E_bin[0]

                rc = self.orbital_per_super_300_1000.setdefault(s, {})
                self.orbital_per_super_300_1000[s][o] = flux_E_bin[1]

                rc = self.orbital_per_super_1000_10000.setdefault(s, {})
                self.orbital_per_super_1000_10000[s][o] = flux_E_bin[2]

                rc = self.super_per_orbital.setdefault(o, {})
                self.super_per_orbital[o][s] = flux_E_bin[3]

                rc = self.super_per_orbital_100_300.setdefault(o, {})
                self.super_per_orbital_100_300[o][s] = flux_E_bin[0]

                rc = self.super_per_orbital_300_1000.setdefault(o, {})
                self.super_per_orbital_300_1000[o][s] = flux_E_bin[1]

                rc = self.super_per_orbital_1000_10000.setdefault(o, {})
                self.super_per_orbital_1000_10000[o][s] = flux_E_bin[2]

    def top_level_hists(self):

        # orbital phase only

        super_markers = False

        phase = list(self.orbital.keys())
        phase.extend(phase)

        fluxs = list(self.orbital.values())
        fluxs.extend(fluxs)
        fluxs_100 = list(self.orbital_100_300.values())
        fluxs_100.extend(fluxs_100)
        fluxs_300 = list(self.orbital_300_1000.values())
        fluxs_300.extend(fluxs_300)
        fluxs_1000 = list(self.orbital_1000_10000.values())
        fluxs_1000.extend(fluxs_1000)

        u_hist = figure(title="orbital flux", x_axis_label='Phase', width=750)
        #u_hist.line(p_bins, fluxs, line_width=2)
        u_hist.vbar(top=locals().get(self.energy_index_flux[self.energy_index]), x=self.phase_h, width=1.,
                    fill_color='red', fill_alpha=0.05, bottom=0)
        u_hist.scatter(self.phase_h, locals().get(self.energy_index_flux[self.energy_index]), size=6,
                       fill_color="white")
        if not self.no_energy_overlay or self.energy_index == 0:
            u_hist.scatter(self.phase_h, np.array(fluxs_100), size=6, fill_color="black", legend_label="100-300 MeV")
            u_hist.scatter(self.phase_h, np.array(fluxs_300), size=6, fill_color="blue", legend_label="300-1000 MeV")
            u_hist.scatter(self.phase_h, np.array(fluxs_1000), size=6, fill_color="green", legend_label="1000-10000 MeV")
            u_hist.y_range = Range1d(0., self.hist_flux_max)

        u_hist.xaxis.ticker = self.phase_h
        u_hist.xaxis.major_label_overrides = self.dict_ticker
        u_hist.xaxis.major_label_orientation = 0.7

        vline_p1 = Span(location=0.23*self.num_pickles, dimension='height', line_color='red', line_width=2,
                line_dash='dashed')
        vline_a1 = Span(location=0.78*self.num_pickles, dimension='height', line_color='blue', line_width=2,
                line_dash='dashed')

        vline_p2 = Span(location=(1.+0.23)*self.num_pickles, dimension='height', line_color='red', line_width=2,
                line_dash='dashed')
        vline_a2 = Span(location=(1.+0.78)*self.num_pickles, dimension='height', line_color='blue', line_width=2,
                line_dash='dashed')

        if not super_markers:
            u_hist.add_layout(vline_p1)
            u_hist.add_layout(vline_a1)
            u_hist.add_layout(vline_p2)
            u_hist.add_layout(vline_a2)

        fluxs = list(self.super.values())
        fluxs.extend(fluxs)
        fluxs_100 = list(self.super_100_300.values())
        fluxs_100.extend(fluxs_100)
        fluxs_300 = list(self.super_300_1000.values())
        fluxs_300.extend(fluxs_300)
        fluxs_1000 = list(self.super_1000_10000.values())
        fluxs_1000.extend(fluxs_1000)

        v_hist = figure(title="super flux", x_axis_label='Phase', width=750)

        v_hist.vbar(top=locals().get(self.energy_index_flux[self.energy_index]), x=self.phase_h, width=1.,
                    fill_color='red', fill_alpha=0.05, bottom=0)
        v_hist.scatter(self.phase_h, locals().get(self.energy_index_flux[self.energy_index]), size=6,
                       fill_color="white")
        if not self.no_energy_overlay or self.energy_index == 0:
            v_hist.scatter(self.phase_h, fluxs_100, size=6, fill_color="black", legend_label="100-300 MeV")
            v_hist.scatter(self.phase_h, fluxs_300, size=6, fill_color="blue", legend_label="300-1000 MeV")
            v_hist.scatter(self.phase_h, fluxs_1000, size=6, fill_color="green", legend_label="1000-10000 MeV")
            v_hist.y_range = Range1d(0., self.hist_flux_max)

        v_hist.xaxis.ticker = self.phase_h
        v_hist.xaxis.major_label_overrides = self.dict_ticker
        v_hist.xaxis.major_label_orientation = 0.7

        return u_hist, v_hist

    def matrix_hists(self):
        # orbital per super phase

        fig_orb_by_super = []
        for s_phase in self.orbital_per_super:

            fluxs = list(self.orbital_per_super[s_phase].values())
            fluxs.extend(fluxs)
            fluxs_100 = list(self.orbital_per_super_100_300[s_phase].values())
            fluxs_100.extend(fluxs_100)
            fluxs_300 = list(self.orbital_per_super_300_1000[s_phase].values())
            fluxs_300.extend(fluxs_300)
            fluxs_1000 = list(self.orbital_per_super_1000_10000[s_phase].values())
            fluxs_1000.extend(fluxs_1000)

            title = "orbital flux per super bin " + str(s_phase)
            a_hist = figure(title=title, x_axis_label='Phase', width=750)
            a_hist.vbar(top=locals().get(self.energy_index_flux[self.energy_index]), x=self.phase_h, width=1.,
                        fill_color='red', fill_alpha=0.05, bottom=0)
            a_hist.scatter(self.phase_h, locals().get(self.energy_index_flux[self.energy_index]), size=6,
                           fill_color="white")
            a_hist.y_range = Range1d(0., self.hist_flux_max)

            if not self.no_energy_overlay:
                # Create a second y-axis
                a_hist.extra_y_ranges = {"y2": Range1d(start=0, end=1.e-7)}
                #a_hist.add_layout(a_hist.yaxis[0], 'left')  # Attach the first y-axis
                a_hist.scatter(self.phase_h, fluxs_100, size=6, fill_color="black", legend_label="100-300 MeV")
                a_hist.scatter(self.phase_h, fluxs_300, size=6, fill_color="blue", legend_label="300-1000 MeV")

                a_hist.square(self.phase_h, fluxs_1000, size=6, fill_color="green",
                               legend_label="1000-10000 MeV", y_range_name="y2")
                a_hist.add_layout(LinearAxis(y_range_name="y2", axis_label='1000-10000 MeV'), 'right')

            a_hist.xaxis.ticker = self.phase_h
            a_hist.xaxis.major_label_overrides = self.dict_ticker
            a_hist.xaxis.major_label_orientation = 0.7

            vline_p1 = Span(location=0.23*self.num_pickles, dimension='height', line_color='red', line_width=2,
                            line_dash='dashed')
            vline_a1 = Span(location=0.78*self.num_pickles, dimension='height', line_color='blue', line_width=2,
                            line_dash='dashed')

            vline_p2 = Span(location=(1.+0.23)*self.num_pickles, dimension='height', line_color='red', line_width=2,
                            line_dash='dashed')
            vline_a2 = Span(location=(1.+0.78)*self.num_pickles, dimension='height', line_color='blue', line_width=2,
                            line_dash='dashed')

            a_hist.add_layout(vline_p1)
            a_hist.add_layout(vline_a1)
            a_hist.add_layout(vline_p2)
            a_hist.add_layout(vline_a2)

            fig_orb_by_super.append(a_hist)

        # orbital per super phase

        fig_super_by_orb = []
        for o_phase in self.super_per_orbital:

            fluxs = list(self.super_per_orbital[o_phase].values())
            fluxs.extend(fluxs)
            fluxs_100 = list(self.super_per_orbital_100_300[o_phase].values())
            fluxs_100.extend(fluxs_100)
            fluxs_300 = list(self.super_per_orbital_300_1000[o_phase].values())
            fluxs_300.extend(fluxs_300)
            fluxs_1000 = list(self.super_per_orbital_1000_10000[o_phase].values())
            fluxs_1000.extend(fluxs_1000)

            title = "super flux per orbital bin " + str(o_phase)
            a_hist = figure(title=title, x_axis_label='Phase', width=750)

            a_hist.vbar(top=locals().get(self.energy_index_flux[self.energy_index]), x=self.phase_h, width=1.,
                        fill_color='red', fill_alpha=0.05, bottom=0)

            a_hist.scatter(self.phase_h, locals().get(self.energy_index_flux[self.energy_index]), size=6,
                           fill_color="white")
            a_hist.y_range = Range1d(0., self.hist_flux_max)

            if not self.no_energy_overlay:
                # Create a second y-axis
                a_hist.extra_y_ranges = {"y2": Range1d(start=0, end=1.e-7)}
                a_hist.add_layout(LinearAxis(y_range_name="y2", axis_label='1000-10000 MeV'), 'right')

                a_hist.scatter(self.phase_h, fluxs_100, size=6, fill_color="black", legend_label="100-300 MeV")
                a_hist.scatter(self.phase_h, fluxs_300, size=6, fill_color="blue", legend_label="300-1000 MeV")
                a_hist.square(self.phase_h, fluxs_1000, size=6, fill_color="green",
                               legend_label="1000-10000 MeV", y_range_name="y2")

            a_hist.xaxis.ticker = self.phase_h
            a_hist.xaxis.major_label_overrides = self.dict_ticker
            a_hist.xaxis.major_label_orientation = 0.7

            fig_super_by_orb.append(a_hist)

        return fig_orb_by_super, fig_super_by_orb

    def make_plots(self):

        rc = self.fill_maps()
        u_hist, v_hist = self.top_level_hists()
        fig_orb_by_super, fig_super_by_orb = self.matrix_hists()
        e_orb_super = {}
        e_super_orb = {}
        self.no_energy_overlay = True
        e_panel = []
        for e in range(4):
            self.energy_index = e
            u_h, v_h = self.top_level_hists()
            f_o, f_s = self.matrix_hists()
            e_orb_super[e] = f_o
            e_super_orb[e] = f_s
            l_e = row(column(u_h, column(f_o)), column(v_h, column(f_s)))
            e_panel.append(TabPanel(child=l_e, title=self.energy_index_flux[e]))

        l_hists = row(column(u_hist, column(fig_orb_by_super)), column(v_hist, column(fig_super_by_orb)))
        del_div = Div(text=self.source_name + " " + self.sed_prefix + " Run on: " +
                           datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        output_file(self.html)
        save(column(del_div, l_hists), title=self.page_title)


if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser(description='plot SEDs')

    parser.add_argument('--config', default="", help="config yaml file")

    args = parser.parse_args()

    p = plot_flux_energies(app_config=args.config)
    p.make_plots()