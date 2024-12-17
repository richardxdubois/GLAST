from astropy.io import fits
import numpy as np
import pandas as pd
import argparse
import yaml
import pickle
from datetime import datetime
from scipy.integrate import quad

from fit_SED import SED_function, fit_SED, flux_function

from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column, gridplot
from bokeh.models import (Label, Span, LinearAxis, Range1d, Whisker, ColumnDataSource, BasicTicker, Tabs,
                          TabPanel, RangeSlider, CustomJS, Button, NumeralTickFormatter, BasicTickFormatter)
from bokeh.models.widgets import Div
from bokeh.palettes import Plasma256 as palette
from bokeh.transform import linear_cmap, log_cmap


class plot_eflux_phase():

    def __init__(self, app_config):


        with open(app_config, "r") as f:
            data = yaml.safe_load(f)

        self.source_name = data["source"]
        self.fgl_source = data["4FGL_source"]

        self.num_pickles = data["num_pickles"]
        self.p_bins = np.arange(self.num_pickles)
        self.html = data["html"]
        self.params_save_pickle = self.html.split(".")[0] + "_params.pkl"
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
        self.initial_guesses = [1e-4, 1., 1e3]
        self.backup_guesses = [1e-4, -0.5, 1e3, 1e3]

        self.all_x = []
        self.all_y = []
        self.all_A = []
        self.all_alpha = []
        self.all_E_0 = []
        self.all_E_cut = []
        self.integrated_fits = []
        self.covariance = []

        self.fermipy_fit = []
        self.fermipy_flux = []
        self.fermipy_alpha = []
        self.fermipy_TS = []

        try:
            self.num_pickles_2 = data["num_pickles_2"]
            self.base_fn_2 = data["base_fn_2"]
            self.p_bins_2 = np.arange(self.num_pickles_2)
            if self.num_pickles_2 == 0:
                self.p_bins_2 = np.arange(1)

            self.type_1 = data["type_1"]
            self.type_2 = data["type_2"]
        except KeyError:
            pass

        try:
            self.phase_offset = int(data["phase_offset"])
        except KeyError:
            self.phase_offset = 0

        try:
            self.no_comp = data["no_comp"]
        except KeyError:
            self.no_comp = False

    def shift_list(self, lst, n):
        n = n % len(lst)  # Ensure n is within the bounds of the list length
        return lst[n:] + lst[:n]

    def shift_map(self, map, n):

        keys = list(map.keys())  # Get the keys as a list
        n = n % len(keys)  # Ensure n is within the bounds of the keys length
        shifted_keys = keys[n:] + keys[:n]  # Shift the keys
        # Create a new dictionary with shifted keys
        shifted_dict = {key: map[key] for key in shifted_keys}

        return shifted_dict


    def loop_over_bins(self):

        for phase_bin in self.p_bins:

            self.seds[phase_bin] = []

            for phase_bin_2 in self.p_bins_2:

                self.all_y.append(phase_bin)
                self.all_x.append(phase_bin_2)

                if self.num_pickles_2 == 0:
                    infile_b = self.base_fn + str(phase_bin) + "/" + self.sed_prefix + \
                              str(phase_bin)
                    infile = infile_b + "_" + self.fgl_source + "_sed.npy"
                else:
                    infile_b = (self.base_fn + str(phase_bin) + "/" + self.base_fn_2 +
                              str(phase_bin_2) + "/" + self.sed_prefix + str(phase_bin_2))
                    infile = infile_b + "_" + self.fgl_source + "_sed.npy"
                infile_f = infile_b + ".npy"

                print("working on ", infile)
                p = np.load(infile, allow_pickle=True).flat[0]

                self.eflux = p["e2dnde"]
                self.eflux_err_lo = p["e2dnde_err_lo"]
                self.eflux_err_hi = p["e2dnde_err_hi"]
                self.loge_ctr = np.power(10, p["loge_ctr"])

                self.b_upper = [x+e for x, e in zip(self.eflux, self.eflux_err_hi)]
                self.b_lower = [x-e for x, e in zip(self.eflux, self.eflux_err_lo)]

                p_fermipy = np.load(infile_f, allow_pickle=True).flat[0]

                f_name = self.fgl_source.replace("_", " ").upper()
                self.fermipy_fit = p_fermipy["sources"][f_name]

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

        flux_i = []
        errors_i = []
        E_ii = []

        for i, E_i in enumerate(self.loge_ctr):
            if E_i <= 10000.:
                if pd.isna(self.eflux_err_lo[i]):
                    continue
                flux_i.append(self.eflux[i])
                errors_i.append(self.eflux_err_hi[i])
                E_ii.append(E_i)

        flux = np.array(flux_i, dtype=np.float64)
        errors = np.array(errors_i, dtype=np.float64)
        E = np.array(E_ii, dtype=np.float64)

        print("Processing super bin ", phase_bin1, "orb bin ", phase_bin2)
        print(flux, errors, E)

        self.fermipy_flux.append(self.fermipy_fit["flux"])
        self.fermipy_TS.append(self.fermipy_fit["ts"])
        self.fermipy_alpha.append(-1.*self.fermipy_fit["param_values"][1])
        try:
            params, covariance = fit_SED(E, flux, errors, self.initial_guesses)
            print("Fitted parameters:", params)
            self.all_A.append(params[0])
            self.all_alpha.append(params[1])
            self.all_E_cut.append(params[2])
            integrated_fits, int_error = quad(flux_function, 100., 10000., args=tuple(params))
            self.integrated_fits.append(integrated_fits)
            self.covariance.append(covariance)

            # Generate data for the fit line
            E_fit = np.linspace(1e2, 1e4, 100)  # Energy range for the fit
            flux_fit = SED_function(E_fit, *params)  # Calculate the fitted flux
            p_fig.line(E_fit, flux_fit, color='red', legend_label='Fitted Model')
            p_fig.legend.location = "bottom_left"
        except:
            print("fit failed")
            self.failed_fits += 1
            self.all_A.append(-999.)
            self.all_alpha.append(-999.)
            self.all_E_cut.append(-999.)
            self.integrated_fits.append(-999.)
            self.covariance.append(-999.)
            pass

        self.seds[phase_bin1].append(p_fig)

    def output_plot(self):

        if self.phase_offset != 0:
            self.all_x = self.shift_list(self.all_x, self.phase_offset)
            self.all_y = self.shift_list(self.all_y, self.phase_offset)
            self.all_A = self.shift_list(self.all_A, self.phase_offset)
            self.all_alpha = self.shift_list(self.all_alpha, self.phase_offset)
            self.all_E_cut = self.shift_list(self.all_E_cut, self.phase_offset)
            self.integrated_fits = self.shift_list(self.integrated_fits, self.phase_offset)
            self.covariance = self.shift_list(self.covariance, self.phase_offset)
            self.seds = self.shift_map(self.seds, self.phase_offset)
            self.fermipy_flux = self.shift_map(self.fermipy_flux, self.phase_offset)
            self.fermipy_alpha = self.shift_map(self.fermipy_alpha, self.phase_offset)
            self.fermipy_TS = self.shift_map(self.fermipy_TS, self.phase_offset)

        print("shifted phase bins by", self.phase_offset)

        min_flux = min(self.integrated_fits)
        min_flux_index = self.integrated_fits.index(min_flux)
        min_flux_params = [self.all_A[min_flux_index], self.all_alpha[min_flux_index],
                                self.all_E_cut[min_flux_index]]
        max_flux = max(self.integrated_fits)
        max_flux_index = self.integrated_fits.index(max_flux)
        max_flux_params = [self.all_A[max_flux_index], self.all_alpha[max_flux_index],
                                self.all_E_cut[max_flux_index]]

        flux_diffs = []   #  comparing lowest to highest flux bins
        for f in self.integrated_fits:
            f_diff_ratio = (f - min_flux) / (max_flux - min_flux)
            flux_diffs.append(f_diff_ratio)

        if self.source_name == "LSI61303" and not self.no_comp:
            E_fit = np.linspace(1e2, 1e4, 100)  # Energy range for the fit
            for i_s, s in enumerate(self.all_y):
                o = self.all_x[i_s]
                flux_fit = ((1.-flux_diffs[i_s])*SED_function(E_fit, *min_flux_params) + flux_diffs[i_s] *
                            SED_function(E_fit, self.all_A[i_s], self.all_alpha[i_s], self.all_E_cut[i_s]))  # Calculate the fitted flux
                self.seds[s][o].line(E_fit, flux_fit, color='blue', legend_label='pulsarness')

        all_lists_params = [self.all_x, self.all_y, self.all_A, self.all_alpha, self.all_E_cut,
                            self.integrated_fits, self.covariance]

        # Write to a pickle file
        with open(self.params_save_pickle, 'wb') as file:
            pickle.dump(all_lists_params, file)

        del_div = Div(text=self.source_name + " Run on: " +
                           datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        output_file(self.html)

        # Prepare Bokeh figures for grid layout
        plots = []
        for i in self.p_bins:
            plots.append(self.seds[i])
        grid = gridplot(plots)

        l = column(del_div, grid)

        # do heatmaps of fit parameters

        source = ColumnDataSource(data=dict(x=self.all_x, y=self.all_y, A=self.all_A, alpha=self.all_alpha,
                                            E_cut=self.all_E_cut, int_f=self.integrated_fits,
                                            fpy_flux=self.fermipy_flux, fpy_alpha=self.fermipy_alpha,
                                            TS=self.fermipy_TS))

        # this is the colormap from the original NYTimes plot
        colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]

        TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

        heatmap_figs = []
        dict_ticker = {}
        for i, x in enumerate(self.all_x):
            dict_ticker[self.all_x[i]] = str(x / 10.)

        tooltips = [[('phases', 'super: @y orbital: @x'), ('A', '@A')],
                    [('phases', 'super: @y orbital: @x'), ('alpha', '@alpha')],
                    [('phases', 'super: @y orbital: @x'), ('E_cut', '@E_cut')],
                    [('phases', 'super: @y orbital: @x'), ('int_f', '@int_f')],
                    [('phases', 'super: @y orbital: @x'), ('fpy_flux', '@fpy_flux')],
                    [('phases', 'super: @y orbital: @x'), ('fpy_alpha', '@fpy_alpha')],
                    [('phases', 'super: @y orbital: @x'), ('TS', '@TS')]
                    ]
        title = ["A", "alpha", "E_cut", "int_f", "fpy_flux", "fpy_alpha", "TS"]
        high = 1.01*np.array([max(self.all_A), max(self.all_alpha), max(self.all_E_cut), max(self.integrated_fits),
                              max(self.fermipy_flux), max(self.fermipy_alpha), max(self.fermipy_TS)])
        low = 0.99*np.array([max(0.,min(self.all_A)), max(0.,min(self.all_alpha)),
                             max(0.,min(self.all_E_cut)), max(0.,min(self.integrated_fits)),
                             max(0.,min(self.fermipy_flux)), max(0.,min(self.fermipy_alpha)),
                             max(0.,min(self.fermipy_TS))])

        for h in range(len(title)):

            p = figure(title=title[h],
                       x_axis_location="above", width=900, height=900,
                       tools=TOOLS, toolbar_location='below', y_axis_label="super phase", x_axis_label="orbital phase",
                       tooltips=tooltips[h])

            p.grid.grid_line_color = None
            p.axis.axis_line_color = None
            p.axis.major_tick_line_color = None
            p.axis.major_label_text_font_size = "12px"
            p.axis.major_label_standoff = 0
            p.xaxis.major_label_orientation = np.pi / 3

            if h == 2:  # E_cut
                fill_color = log_cmap(title[h], palette=palette, low=low[h], high=high[h])
            else:
                fill_color = linear_cmap(title[h], palette=palette, low=low[h], high=high[h])

            r = p.rect(x="x", y="y", width=1, height=1, source=source,
                       fill_color=fill_color,
                       line_color=None, )
            p.xaxis.ticker = self.all_x
            p.xaxis.major_label_overrides = dict_ticker
            p.yaxis.ticker = self.all_y
            p.yaxis.major_label_overrides = dict_ticker

            p.add_layout(r.construct_color_bar(
                major_label_text_font_size="12px",
                ticker=BasicTicker(desired_num_ticks=len(colors)),

                label_standoff=6,
                border_line_color=None,
                padding=5,
            ), 'right')

            heatmap_figs.append(p)

        # histograms

        fpy_alpha_hist, fpy_alpha_edges = np.histogram(self.fermipy_alpha, bins=25)
        alpha_hist = figure(title="fpy Alpha",
                            x_axis_label='Index', y_axis_label='counts',
                            width=750)

        alpha_hist.vbar(top=fpy_alpha_hist, x=fpy_alpha_edges[1:], width=fpy_alpha_edges[1] - fpy_alpha_edges[0],
                        fill_color='red', fill_alpha=0.1, bottom=0)

        fpy_flux_hist, fpy_flux_edges = np.histogram(self.fermipy_flux, bins=25)
        flux_hist = figure(title="fpy Flux",
                            x_axis_label='Flux', y_axis_label='counts',
                            width=750)

        flux_hist.vbar(top=fpy_flux_hist, x=fpy_flux_edges[1:], width=fpy_flux_edges[1] - fpy_flux_edges[0],
                        fill_color='green', fill_alpha=0.1, bottom=0)

        flux_alpha_scat = figure(title="fpy Flux vs Alpha",
                            x_axis_label='Alpha', y_axis_label='Flux',
                            width=750)
        flux_alpha_scat.scatter(x=fpy_alpha_hist, y=fpy_flux_hist, color="blue", size=2)

        # create sliders
        steps = (high - low)/20.

        slider_TS = RangeSlider(start=low[6], end=high[6], value=(low[6], high[6]), step=steps[6], title="TS")
        slider_A = RangeSlider(start=low[0], end=high[0], value=(low[0], high[0]), step=steps[0], title="A")
        slider_A.format = BasicTickFormatter(use_scientific=True)  # Set format to scientific notation
        slider_alpha = RangeSlider(start=low[1], end=high[1], value=(low[1], high[1]), step=steps[1], title="alpha")

        slider_E_cut = RangeSlider(start=low[2], end=high[2], value=(low[2], high[2]), step=steps[2], title="E_cut")


        # CustomJS callback to update rectangle colors based on slider values
        callback = CustomJS(
            args=dict(source=source, slider0=slider_TS, slider1=slider_A, slider2=slider_alpha,
                      slider3=slider_E_cut), code="""
            //if (Object.keys(original_data).length == 0) {
            if (typeof original_data == 'undefined') {
                window.original_data = new Map().constructor;
                original_data.TS_orig = source.data['TS'].slice();
                original_data.A_orig = source.data['A'].slice();
                original_data.alpha_orig = source.data['alpha'].slice();
                original_data.E_cut_orig = source.data['E_cut'].slice();
                original_data.int_f = source.data['int_f'].slice();
                original_data.fpy_flux = source.data['fpy_flux'].slice();
                original_data.fpy_alpha = source.data['fpy_alpha'].slice();
                }
            const data = source.data;
            const y = data['y'];
            const x = data['x'];
            const TS = data['TS'];
            const A = data['A'];
            const alpha = data['alpha'];
            const E_cut = data['E_cut'];
            const int_f = data['int_f'];
            const fpy_flux = data['fpy_flux'];
            const fpy_alpha = data['fpy_alpha'];
            
            const lowerLimit_TS = slider0.value[0];
            const upperLimit_TS = slider0.value[1];
            const lowerLimit_A = slider1.value[0];
            const upperLimit_A = slider1.value[1];
            const lowerLimit_alpha = slider2.value[0];
            const upperLimit_alpha = slider2.value[1];
            const lowerLimit_E_cut = slider3.value[0];
            const upperLimit_E_cut = slider3.value[1];

            for (let i = 0; i < y.length; i++) {
                
                if (original_data.A_orig[i] < lowerLimit_A || original_data.A_orig[i] > upperLimit_A 
                || original_data.alpha_orig[i] < lowerLimit_alpha 
                || original_data.alpha_orig[i] > upperLimit_alpha 
                || original_data.TS_orig[i] < lowerLimit_TS 
                || original_data.TS_orig[i] > upperLimit_TS
                ||original_data.E_cut_orig[i] < lowerLimit_E_cut || original_data.E_cut_orig[i] > upperLimit_E_cut) {
                    TS[i] = -1.;
                    A[i] = -1.;
                    alpha[i] = -1.;
                    E_cut[i] = -1.;
                    int_f[i] = -1.
                    fpy_flux[i] = -1.;
                    fpy_alpha[i] = -1.;
                } else {
                    source.data["TS"][i] = original_data.TS_orig[i];
                    source.data["A"][i] = original_data.A_orig[i];
                    source.data["alpha"][i] = original_data.alpha_orig[i];
                    source.data["E_cut"][i] = original_data.E_cut_orig[i];
                    source.data["int_f"][i] = original_data.int_f[i];
                    source.data["fpy_flux"][i] = original_data.fpy_flux[i];
                    source.data["fpy_alpha"][i] = original_data.fpy_alpha[i];
                }
            }
            source.change.emit();
        """)
        # Create a button to execute the custom JavaScript callback
        button = Button(label="Save Original Data before using sliders", button_type="danger")

        # CustomJS callback to save the original data in JavaScript
        button.js_on_click(CustomJS(args=dict(source=source, button=button), code="""
            // Save the original data in a JavaScript variable
        
            original_data.TS_orig = source.data['TS'].slice();
            original_data.A_orig = source.data['A'].slice();
            original_data.alpha_orig = source.data['alpha'].slice();
            original_data.E_cut_orig = source.data['E_cut'].slice();
            original_data.int_f = source.data['int_f'].slice();
            button.button_type = 'success'
            button.label = 'Data saved'
            
        """))

        # Attach the callback to the sliders
        slider_TS.js_on_change('value', callback)
        slider_A.js_on_change('value', callback)
        slider_alpha.js_on_change('value', callback)
        slider_E_cut.js_on_change('value', callback)

        # Layout the sliders and the plot - remove Button from layout. At some point, remove it from code.
        s = column(slider_TS, slider_A, slider_alpha, slider_E_cut)
        h_layout = row(column(del_div, s, column(heatmap_figs)), column(alpha_hist, flux_hist, flux_alpha_scat))

        if self.source_name == "LSI61303" and not self.no_comp:

            # do heatmaps of fit parameters

            source = ColumnDataSource(data=dict(x=self.all_x, y=self.all_y, flux_d=flux_diffs))

            # this is the colormap from the original NYTimes plot
            colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]

            TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

            dict_ticker = {}
            for i, x_i in enumerate(self.all_x[0:self.num_pickles]):
                dict_ticker[self.all_x[i]] = str(x_i / self.num_pickles)

            tooltips = [('phases', 'super: @y orbital: @x'), ('flux_d', '@flux_d')]

            pc = figure(title="Flux ratio",
                       x_axis_location="above", width=900, height=900,
                       tools=TOOLS, toolbar_location='below', y_axis_label="super phase", x_axis_label="orbital phase",
                       tooltips=tooltips)

            pc.grid.grid_line_color = None
            pc.axis.axis_line_color = None
            pc.axis.major_tick_line_color = None
            pc.axis.major_label_text_font_size = "12px"
            pc.axis.major_label_standoff = 0
            pc.xaxis.major_label_orientation = np.pi / 3

            fill_color = linear_cmap("flux_d", palette=palette, low=0, high=1)

            r = pc.rect(x="x", y="y", width=1, height=1, source=source,
                       fill_color=fill_color,
                       line_color=None, )
            pc.xaxis.ticker = self.all_x[0:self.num_pickles]
            pc.xaxis.major_label_overrides = dict_ticker
            pc.yaxis.ticker = self.all_x[0:self.num_pickles_2]
            pc.yaxis.major_label_overrides = dict_ticker

            pc.add_layout(r.construct_color_bar(
                major_label_text_font_size="12px",
                ticker=BasicTicker(desired_num_ticks=len(colors)),

                label_standoff=6,
                border_line_color=None,
                padding=5,
            ), 'right')


        panel1 = TabPanel(child=h_layout, title="Parameter heatmaps")
        panel2 = TabPanel(child=l, title="SED matrix")

        if self.source_name == "LSI61303" and not self.no_comp:
            panel3 = TabPanel(child=pc, title="LSI61303 pulsarness")
            tabs = Tabs(tabs=[panel1, panel2, panel3])
        else:
            tabs = Tabs(tabs=[panel1, panel2])

        save(tabs, title=self.page_title)

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
