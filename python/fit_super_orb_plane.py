import numpy as np
import pandas as pd
import argparse
import yaml
import pickle
from datetime import datetime

from scipy.ndimage import gaussian_filter, label
from scipy.ndimage import maximum_filter

from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column, gridplot
from bokeh.models import (Label, Span, LinearAxis, Range1d, Whisker, ColumnDataSource, BasicTicker, Tabs,
                          TabPanel, RangeSlider, CustomJS, Button, NumeralTickFormatter, BasicTickFormatter)
from bokeh.models.widgets import Div
from bokeh.palettes import Plasma256 as palette
from bokeh.transform import linear_cmap, log_cmap


class fit_super_orb_plane():
    def __init__(self, in_pickle=None):

        # Read from the pickle file

        self.in_pickle = in_pickle
        self.html_name = self.in_pickle.split('.')[0] + ".html"

        with open(in_pickle, 'rb') as file:
            loaded_lists = pickle.load(file)

        self.data_dict = {
            "x": loaded_lists[0],  # "all_x"
            "y": loaded_lists[1],  # "all_y"
            "A": loaded_lists[2],  # "all_A"
            "alpha": loaded_lists[3],       # "all_alpha"
            "E_cut": loaded_lists[4], # "all_E_cut
            "int_f": loaded_lists[5],  # "all_int_f"
            "fpy_flux": loaded_lists[6],
            "fpy_flux_err": loaded_lists[7],
            "fpy_alpha": loaded_lists[8],
            "fpy_TS": loaded_lists[9],
            "cov": loaded_lists[10]  # "parameter covariance"
        }
        self.num_bins = 10
        self.maxima = None
        self.peaks = None
        self.num_peaks = None
        self.labeled_peaks = None
        self.z_2d = None
        self.z2d_err = None
        self.x_2d = None
        self.y_2d = None
        self.z_2d_2d = None

        print("loaded data_dict")

    def find_peaks(self):

        self.x_2d = np.empty((2*self.num_bins, 2*self.num_bins))
        self.y_2d = np.empty((2*self.num_bins, 2*self.num_bins))
        self.z_2d = np.empty((2*self.num_bins, 2*self.num_bins))

        for iy in range(self.num_bins):
            for ix in range(self.num_bins):
                d_i = self.num_bins * iy + ix
                self.x_2d[iy, ix] = np.array(self.data_dict["x"][d_i])
                self.x_2d[iy, ix + self.num_bins] = np.array(self.data_dict["x"][d_i])
                self.x_2d[iy + self.num_bins, ix] = np.array(self.data_dict["x"][d_i]) + self.num_bins
                self.x_2d[iy + self.num_bins, ix + self.num_bins] = np.array(self.data_dict["x"][d_i]) + self.num_bins

                self.y_2d[iy, ix] = np.array(self.data_dict["y"][d_i])
                self.y_2d[iy, ix + self.num_bins] = np.array(self.data_dict["y"][d_i]) + self.num_bins
                self.y_2d[iy + self.num_bins, ix] = np.array(self.data_dict["y"][d_i])
                self.y_2d[iy + self.num_bins, ix + self.num_bins] = np.array(self.data_dict["y"][d_i]) + self.num_bins

                self.z_2d[iy, ix] = np.array(self.data_dict["fpy_flux"][d_i])
                self.z_2d[iy, ix + self.num_bins] = np.array(self.data_dict["fpy_flux"][d_i])
                self.z_2d[iy + self.num_bins, ix] = np.array(self.data_dict["fpy_flux"][d_i])
                self.z_2d[iy + self.num_bins, ix + self.num_bins] = np.array(self.data_dict["fpy_flux"][d_i])

        # Apply a Gaussian filter to smooth the data
        Z_smooth = gaussian_filter(self.z_2d, sigma=0.5)

        # Find local maxima (peaks) with a maximum filter
        self.maxima = maximum_filter(Z_smooth, size=3, mode='constant')

        # Create a boolean array where Z_smooth equals to maxima (indicates peaks)
        self.peaks = (Z_smooth == self.maxima)

        # Label the peaks
        self.labeled_peaks, self.num_peaks = label(self.peaks)
        print("found ", self.num_peaks, " peaks")

        for s in range(self.num_bins):
            for o in range(self.num_bins):
                if self.peaks[s, o]:
                    print("peak: o=", o, "s=", s, self.maxima[s, o], self.z_2d[s, o])

    # Define the 2D Gaussian function for fitting

    def gaussian_2d(self, xy, amplitude, x0, y0, sigma_x, sigma_y, offset):
        x, y = xy
        return offset + amplitude * np.exp(
            -(((x - x0) ** 2) / (2 * sigma_x ** 2) + ((y - y0) ** 2) / (2 * sigma_y ** 2)))

    def fit_surface(self):

        # Fit a 2D Gaussian for each detected peak
        fitted_params = []

        x_flat = self.x_2d.flatten()
        y_flat = self.y_2d.flatten()
        z_flat = self.z_2d.flatten()

    def make_plots(self):

        peak_positions = np.argwhere(self.peaks)

        tooltips = [('phases', 'super: @y orbital: @x'), ('fpy_flux', '@z')]

        # Create a ColumnDataSource for the contour data
        source = ColumnDataSource(data=dict(x=self.x_2d.ravel(), y=self.y_2d.ravel(), z=self.z_2d.ravel()))

        # Create Bokeh figure
        p = figure(title="2D Peak Detection with Bokeh", width=600, height=600, x_axis_label='Orbital',
                   y_axis_label='Super', tooltips=tooltips)
        p.x_range.start = -1

        # Add a heatmap of the smoothed data using image
        fill_color = linear_cmap("z", palette=palette, low=min(self.z_2d.ravel()),
                                 high=max(self.z_2d.ravel()))
        r = p.rect(x="x", y="y", width=1, height=1, source=source,
                   fill_color=fill_color,
                   line_color=None, )

        # Overlay detected peaks
        peak_x = peak_positions[:, 1]  # X coordinates of peaks
        peak_y = peak_positions[:, 0]  # Y coordinates of peaks

        # Add peaks as circle glyphs
        p.scatter(peak_x, peak_y, size=10, color='red', alpha=0.7, legend_label='Detected Peaks')

        # histogram by orbital phase

        orb = np.zeros(2*self.num_bins)

        for o in range(2*self.num_bins):
            for s in range(self.num_bins):
                orb[o] += self.z_2d[s, o]/self.num_bins

        hist_orb = figure(title="orbital phase", width=600, height=600, x_axis_label='Orbital Phase')
        hist_orb.vbar(x=np.arange(20),top=orb, width=1., fill_color="blue", alpha=0.2)

        # histogram by super phase

        super = np.zeros(2*self.num_bins)

        for s in range(2*self.num_bins):
            for o in range(self.num_bins):
                super[s] += self.z_2d[s, o]/self.num_bins

        hist_super = figure(title="super phase", width=600, height=600, y_axis_label='Super Phase')
        hist_super.hbar(y=np.arange(20), right=super, fill_color="blue", alpha=0.2)

        div_text = " Run on: " + self.in_pickle + " " + datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        del_div = Div(text=div_text)


        output_file(self.html_name)
        save(column(del_div, row(column(p, hist_orb), hist_super)))


if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser(description='plot SEDs')

    parser.add_argument('--input',
                        default="/Users/richarddubois/Code/GLAST/tmp/LSI61303_super_orb_sed_PL_fixedP_roi21_params.pkl",
                        help="overall input pickle file")

    args = parser.parse_args()

    f = fit_super_orb_plane(in_pickle=args.input)
    f.find_peaks()
    f.make_plots()
