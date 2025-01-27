from astropy.io import fits
import numpy as np

import yaml
import argparse

from datetime import datetime, date, timedelta
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

from bokeh.models.widgets import DataTable, TableColumn, Div, NumberFormatter
from bokeh.models.formatters import DatetimeTickFormatter
from bokeh.models import Label, Span, LinearAxis, Range1d, Whisker, ColumnDataSource
from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column


class process_LAT_binned_exposure():
    def __init__(self, input_yaml=None, source=None, suppress_zero=True, do_weights=True):

        with open(input_yaml, "r") as f:
            data = yaml.safe_load(f)

        self.source = data["source"]
        self.html_name = data["html"]
        self.fn = data["input_file"]
        self.suppress_zero = data["suppress_zero"]
        self.do_weights = data["do_weights"]

        self.f_start = 1./data["f_start"]/86400.
        self.f_stop = 1./data["f_stop"]/86400.
        self.nom_period = data["nom_period"]
        self.power_threshold = data["power_threshold"]   # 0.4

        self.hdul = None
        self.time = None
        self.counts = None
        self.timedel = None
        self.exposure = None

        try:
            self.super_period = data["super_period"]
            self.super_start = data["super_start"]
            self.s_start = 1. / self.super_start / 86400.
        except KeyError:
            self.super_period = 0.

        self.nom_freq = 1/self.nom_period/86400.

        try:
            self.original_ft1 = data["original_ft1"]
            self.do_src_prob = data["do_src_prob"]
            self.fgl_src = data["fgl_src"]
        except KeyError:
            self.original_ft1 = None
            self.do_src_prob = False
            self.fgl_src = None

    def get_data(self):
        # Print information about the FITS file
        self.hdul = fits.open(self.fn)
        self.hdul.info()

        # Access the primary HDU (Header/Data Unit)
        primary_hdu = self.hdul[0]

        # Get the data and header
        print(self.hdul[1].columns)
        header = primary_hdu.header
        data = self.hdul[1].data

        indices = np.where(data.COUNTS > 0.)
        self.timedel = data.TIMEDEL

        if self.suppress_zero:
            self.exposure = data.EXPOSURE[indices]
            self.time = data.TIME[indices] #+ self.timedel[0]/2.
            self.exposure = data.EXPOSURE[indices]
            self.counts = data.COUNTS[indices]

            indices_exp = np.where(self.exposure != 0)

            self.time = self.time[indices_exp] #+ self.timedel[0]/2.
            self.exposure = self.exposure[indices_exp]
            self.counts = self.counts[indices_exp]
        else:
            self.time = data.TIME + self.timedel[0]/2.
            self.exposure = data.EXPOSURE
            self.counts = data.COUNTS

    def Robin_exp_weight(self, timeslice=None):

        r_weighted = []
        weights = []

        cnts = None
        exp = None
        if timeslice is None:
            cnts = self.counts
            exp = self.exposure
            t = self.time
        else:
            cnts = list(self.counts[timeslice])
            exp = list(self.exposure[timeslice])
            t = list(self.time[timeslice])

        print(" len cnts, exp, t", len(cnts), len(exp), len(t))
        mean_rate = np.sum(cnts) / np.sum(exp)

        for i, e in enumerate(exp):

            pop_counts = mean_rate * e
            perr = np.sqrt(pop_counts)

            cerr1 = 0.5 + np.sqrt(cnts[i] + 0.25)
            cerr2 = -0.5 + np.sqrt(max(0., cnts[i] - 0.25))
            meanerr = (cerr1 + cerr2)/2.
            rmserr = np.sqrt((cerr1*cerr1 + cerr2*cerr2)/2.)

            if self.do_weights:
                rate = cnts[i] / e
                err_rate = rmserr / e

                r_weighted.append(rate)
                weights.append(err_rate)
            else:
                r_weighted.append(cnts[i])
                weights.append(1.)

        r_weighted = np.array(r_weighted)
        weights = np.array(weights)

        return r_weighted, weights

    def get_src_probs(self):

        t = []
        s_prob = []

        for f in self.original_ft1:
            h = fits.open(f)
            h.info()

            data = h[1].data

            t_f = data.TIME
            s_prob_f = data.field(self.fgl_src)  # source prob column

            t.append(t_f)
            s_prob.append(s_prob_f)

            h.close()

        edges = [self.time[0] + i*self.timedel for i in range(len(self.time)+1)]
        bin_i = np.digitize(t, edges)
        binned_sums = np.zeros(range(self.time))

        for i in range(len(t)):
            if 0 <= bin_i[i] < len(binned_sums):
                binned_sums[bin_i[i]] += 1./s_prob[i]

        return binned_sums

    def calc_peak_error(self, frequency=None, power=None, peak_index=None):

        pk_freq_idx = peak_index
        max_power = power[pk_freq_idx]
        max_freq = frequency[pk_freq_idx]

        # find half point above peak

        idx = pk_freq_idx
        while power[idx] > max_power/2.:
            p = power[idx]
            f = frequency[idx]
            idx += 1

        delta_up = (1./frequency[pk_freq_idx] - 1./frequency[idx])/86400.

        idx = pk_freq_idx
        while power[idx] > max_power / 2.:
            idx -= 1

        delta_down = (1./frequency[idx] - 1./frequency[pk_freq_idx])/86400.

        std_err = (delta_up + delta_down)/2.

        return std_err

    def make_plots(self):

        r_weighted, weights = self.Robin_exp_weight()
        if self.do_src_prob:
            src_weights = self.get_src_probs()
            weights = weights * src_weights

        counts_hist, counts_edges = np.histogram(self.counts, bins=100)
        q_hist = figure(title="Counts",
                        x_axis_label='Counts', y_axis_label='counts',
                        width=750)

        q_hist.vbar(top=counts_hist, x=counts_edges[1:], width=counts_edges[1]-counts_edges[0], fill_color='red',
                    fill_alpha=0.2, bottom=0)

        exp_hist, exp_edges = np.histogram(self.exposure, bins=100)
        r_hist = figure(title="Exposure",
                        x_axis_label='cm^2 s', y_axis_label='counts',
                        width=750)

        r_hist.vbar(top=exp_hist, x=exp_edges[1:], width=exp_edges[1]-exp_edges[0], fill_color='red', fill_alpha=0.2, bottom=0)

        ave_hist, ave_edges = np.histogram(weights, bins=100)
        s_hist = figure(title="Weights",
                        x_axis_label='Weights', y_axis_label='counts',
                        width=750)

        s_hist.vbar(top=ave_hist, x=ave_edges[1:], width=ave_edges[1]-ave_edges[0], fill_color='red', fill_alpha=0.2, bottom=0)

        f3 = figure(title="Exposure vs Counts",  x_axis_label='Exposure', y_axis_label='Counts', width=750)
        f3.scatter(x=self.exposure, y=self.counts, size=3)

        t_start = self.time[0]
        t_end = self.time[-1]

        f_hists = []
        print("t_start", t_start, "t_end", t_end, "num orig bins", len(self.time), "num non-zero bins",
              len(self.counts))

        frequency = np.linspace(self.f_start, self.f_stop, 1000)  # for orbital 100000

        power = LombScargle(t=self.time, y=r_weighted, dy=weights).power(frequency)


        print(max(frequency), max(power))
        print(min(frequency), min(power))
        freq_days = 1./frequency/86400.

        peaks_ls, props_ls = find_peaks(power, height=0.1 * max(power))
        pk_days = (1./frequency[peaks_ls]/86400.)
        close_idx = (np.abs(pk_days - self.nom_period)).argmin()
        peak_idx = peaks_ls[close_idx]

        pk_error = self.calc_peak_error(frequency=frequency, power=power, peak_index=peak_idx)

        print(peaks_ls, props_ls)
        print(frequency[peaks_ls])
        print(pk_days)
        print("peak error", pk_error)

        f1 = figure(title="full time span: power vs frequency",
                                x_axis_label='period (days)', y_axis_label='power',
                                width=750)
        f1.line(freq_days, power, line_width=2)
        vline_p1 = Span(location=self.nom_period, dimension='height', line_color='red', line_width=2, line_dash='dashed')
        f1.add_layout(vline_p1)
        res_label = Label(x=min(freq_days), y=props_ls["peak_heights"][close_idx]/2., text_font_size="8pt",
                          text="Peak : " + str('{0:.3f}'.format(pk_days[close_idx])) + "+/- " +
                               str('{0:.3f}'.format(pk_error) + " days"))
        f1.add_layout(res_label)
        f_hists.append(f1)

        if self.super_period != 0.:  # larger frequency range to include super orbital period

            s_frequency = np.linspace(self.s_start, self.f_stop, 100000)  # for orbital 100000

            power = LombScargle(t=self.time, y=r_weighted, dy=weights).power(s_frequency)

            print(max(frequency), max(power))
            print(min(frequency), min(power))
            freq_days = 1. / s_frequency / 86400.

            peaks_ls, props_ls = find_peaks(power, height=0.1 * max(power))
            pk_days = (1. / s_frequency[peaks_ls] / 86400.)
            close_idx = (np.abs(pk_days - self.nom_period)).argmin()
            peak_idx = peaks_ls[close_idx]

            pk_error = self.calc_peak_error(frequency=s_frequency, power=power, peak_index=peak_idx)

            print(peaks_ls, props_ls)
            print(s_frequency[peaks_ls])
            print(pk_days)
            print("peak error", pk_error)

            fs = figure(title="full super time span: power vs frequency", x_axis_type="log",
                        x_axis_label='period (days)', y_axis_label='power',
                        width=750)
            fs.line(freq_days, power, line_width=2)
            vline_p1 = Span(location=self.nom_period, dimension='height', line_color='red', line_width=2,
                            line_dash='dashed')
            fs.add_layout(vline_p1)
            vline_p2 = Span(location=self.super_period, dimension='height', line_color='blue', line_width=2,
                            line_dash='dashed')
            fs.add_layout(vline_p1)
            fs.add_layout(vline_p2)
            res_label = Label(x=min(freq_days), y=props_ls["peak_heights"][close_idx] / 2., text_font_size="8pt",
                              text="Peak : " + str('{0:.3f}'.format(pk_days[close_idx])) + "+/- " +
                                   str('{0:.3f}'.format(pk_error) + " days"))
            fs.add_layout(res_label)
            f_hists.append(fs)

        timespan = self.time[-1] - self.time[0]
        yr_bins = 4
        delta = timespan / yr_bins    # ~ 1 yr
        yr_figs = []
        yr_exp = []
        yr_c = []
        yr_orb_power = []
        yr_peak = []
        yr_peak_error = []

        for yr in range(yr_bins):
            tmin = self.time[0] + yr * delta
            tmax = tmin + delta

            yr_ind = [i for i, t in enumerate(self.time) if (t >= tmin and t < tmax)]

            r_weighted, weights = self.Robin_exp_weight(timeslice=yr_ind)

            yr_times = self.time[yr_ind]
            yr_counts_d = self.counts[yr_ind]
            yr_counts = r_weighted
            yr_weights = weights
            yr_exposure = self.exposure[yr_ind]
            sum_weights = np.sum(yr_weights)

            yr_power = LombScargle(t=yr_times, y=yr_counts, dy=yr_weights).power(frequency)

            freq_days = 1. / frequency / 86400.

            yr_figs.append(figure(title="time span: " + str(yr) + " " + str(tmin) + ", " + str(tmax) + " power vs frequency",
                                x_axis_label='period (days)', y_axis_label='power',
                                width=750))
            yr_figs[yr].line(freq_days, yr_power, line_width=2)
            yr_figs[yr].add_layout(vline_p1)

            yr_peaks_ls, yr_props_ls = find_peaks(yr_power, height=self.power_threshold * max(yr_power))

            pk_days = (1. / frequency[yr_peaks_ls] / 86400.)
            close_idx = (np.abs(pk_days - self.nom_period)).argmin()
            peak_idx = yr_peaks_ls[close_idx]
            pk_error = self.calc_peak_error(frequency=frequency, power=yr_power, peak_index=peak_idx)
            yr_peak.append(pk_days[close_idx])
            yr_peak_error.append(pk_error)

            res_label = Label(x=min(freq_days), y=yr_props_ls["peak_heights"][close_idx] / 2., text_font_size="8pt",
                              text="Peak : " + str('{0:.3f}'.format(pk_days[close_idx])) + "+/- " +
                                   str('{0:.3f}'.format(pk_error) + " days"))
            yr_figs[yr].add_layout(res_label)

            print("yr", yr, "tmin", tmin, "tmax", tmax, "t[0]", yr_times[0], "t[1]", yr_times[-1],
                  "sum_weights", sum_weights, "num w bins", len(yr_weights), "num t bins", len(yr_counts))
            print(yr_peaks_ls, yr_props_ls)
            print(frequency[yr_peaks_ls])
            print(pk_days)
            print("peak error", pk_error)

            orb_power = 0.
            for p in np.arange(len(yr_peaks_ls)):
                i = yr_peaks_ls[p]
                if abs(frequency[i] - self.nom_freq) < 0.05 * self.nom_freq:
                    orb_power = yr_props_ls["peak_heights"][p]
                    print(i, frequency[i], orb_power)
                    break
            yr_orb_power.append(orb_power)

            counts_hist, counts_edges = np.histogram(yr_counts_d, bins=100)

            yr_c.append(figure(title="Counts",
                            x_axis_label='Counts', y_axis_label='counts',
                            width=750))

            yr_c[yr].vbar(top=counts_hist, x=counts_edges[1:], width=counts_edges[1] - counts_edges[0], fill_color='red',
                        fill_alpha=0.2, bottom=0)

            exp_hist, exp_edges = np.histogram(yr_exposure, bins=100)
            yr_exp.append(figure(title="Exposure",
                            x_axis_label='cm^2 s', y_axis_label='counts',
                            width=750))

            yr_exp[yr].vbar(top=exp_hist, x=exp_edges[1:], width=exp_edges[1] - exp_edges[0], fill_color='red', fill_alpha=0.2,
                        bottom=0)

        f2 = figure(title="orb results vs time bin",
                                x_axis_label='time bin', y_axis_label='power',
                                width=750)
        f2.line(np.arange(yr_bins), yr_orb_power, line_width=2, legend_label="Power (left)")

        f2.scatter(np.arange(yr_bins), yr_orb_power, size=6, fill_color="white")
        f2.y_range = Range1d(0.8*min(yr_orb_power), 1.2*max(yr_orb_power))
        f2.extra_y_ranges = {"y2": Range1d(start=0.96*min(yr_peak), end=1.04*max(yr_peak))}
        f2.add_layout(LinearAxis(y_range_name="y2"), 'right')
        f2.line(np.arange(yr_bins), yr_peak, line_width=2, color="blue", legend_label="Period (right)",
                y_range_name="y2")
        f2.scatter(np.arange(yr_bins), yr_peak, size=6, fill_color="white", y_range_name="y2")
        f_upper = [x + e for x, e in zip(yr_peak, yr_peak_error)]
        f_lower = [x - e for x, e in zip(yr_peak, yr_peak_error)]
        f_source = ColumnDataSource(data=dict(groups=np.arange(yr_bins), counts=yr_peak, upper=f_upper, lower=f_lower))
        f2.add_layout(Whisker(source=f_source, base="groups", upper="upper", lower="lower", level="overlay",
                      y_range_name="y2"))

        del_div = Div(text=self.source + " Run on: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + " for " + self.fn)

        output_file(self.html_name)
        l = layout(column(del_div, column(f_hists), f2, q_hist, r_hist, s_hist, f3, column(yr_figs), column(yr_c), column(yr_exp)))
        save(l, title=self.source + " Power vs Frequency")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--app_config',
                        default="process_exposure_config.yaml",
                        help="overall app config file")
    args = parser.parse_args()

    p = process_LAT_binned_exposure(input_yaml=args.app_config)

    p.get_data()
    p.make_plots()

# Close the FITS file
    p.hdul.close()

# Print the data and header
    print(p.hdul[0].header)
