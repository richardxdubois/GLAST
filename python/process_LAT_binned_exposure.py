from astropy.io import fits
import numpy as np
from datetime import datetime, date, timedelta
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

from bokeh.models.widgets import DataTable, TableColumn, Div, NumberFormatter
from bokeh.models.formatters import DatetimeTickFormatter
from bokeh.models import ColumnDataSource, Span, LinearAxis, Range1d, LinearColorMapper
from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column


class process_LAT_binned_exposure():
    def __init__(self, source=None, suppress_zero=True, do_weights=True):

        self.source = source
        self.html_name = source + "_periodicity.html"
        self.fn = ""
        self.hdul = None
        self.time = None
        self.counts = None
        self.timedel = None
        self.exposure = None
        self.suppress_zero = suppress_zero
        self.do_weights = do_weights

# Open the FITS file
        if self.source == "LSI61303":
            #self.fn = '/Users/richarddubois/Code/GLAST/tmp/LSI61303_1_deg_mkt_500s.fits'
            #self.fn = '/Users/richarddubois/Code/GLAST/tmp/LSI61303_1_deg_mkt_86400s.fits'
            self.fn ='/Users/richarddubois/Code/GLAST/tmp/LSI61303_1_deg_mkt_225000s.fits'
            #fn = '/Users/richarddubois/Code/GLAST/tmp/LSI61303_3_deg_mkt_10800s.fits'

            self.f_start = 1./28.5/86400.  # 40.
            self.f_stop = 1./24.5/86400.   # 5.
            self.nom_period = 26.496

        elif self.source == "LS5039":
            self.fn = '/Users/richarddubois/Code/GLAST/tmp/LS5039_1_deg_500s.fits'
            #seld.fn = '/Users/richarddubois/Code/GLAST/tmp/LS5039_3_deg_10800s.fits'
            self.f_start = 1./10./86400.
            self.f_stop = 1./2./86400.
            self.nom_period = 3.9

        self.nom_freq = 1/self.nom_period/86400.

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
            self.time = data.TIME[indices] #+ self.timedel[0]/2.
            self.exposure = data.EXPOSURE[indices]
            self.counts = data.COUNTS[indices]
        else:
            self.time = data.TIME + self.timedel[0]/2.
            self.exposure = data.EXPOSURE
            self.counts = data.COUNTS

    def kludge_exp_weight(self):

        r_weighted = []
        weights = []

        mean_counts = np.average(self.counts)
        mean_rate = np.sum(self.counts) / np.sum(self.exposure)
        mean_sub = self.counts - mean_counts

        print("mean_counts", mean_counts, "mean_rate", mean_rate)

        # best weighting so far:
        #  err - np.sqrt(mean_rate*e)
        #  r_weighted.append(self.counts[i]/ave)

        for i, e in enumerate(self.exposure):
            if e == 0.:
                r_weighted.append(0.)
                print("skipped over zero exposure at ", i, self.counts[i])
            else:
                pred = mean_rate*e
                err = np.sqrt(pred)
                #err = mean_rate*np.sqrt(e)  # "predicted" count error - optimal
                if not self.do_weights:
                    ave = 1.
                #r_weighted.append(mean_sub[i] / e)
                r_weighted.append(self.counts[i]/err)
                weights.append(err)

        r_weighted = np.array(r_weighted)
        #weights = np.array(weights)
        weights = None

        return r_weighted, weights

    def Robin_exp_weight(self, timeslice=None):

        r_weighted = []
        weights = []

        cnts = None
        exp = None
        if timeslice is None:
            cnts = self.counts
            exp = self.exposure
        else:
            cnts = self.counts[timeslice]
            exp = self.exposure[timeslice]

        mean_rate = np.sum(cnts) / np.sum(exp)

        for i, e in enumerate(exp):
            pop_counts = mean_rate * e
            perr = np.sqrt(pop_counts)

            cerr1 = 0.5 + np.sqrt(cnts[i] + 0.25)
            cerr2 = -0.5 + np.sqrt(cnts[i] - 0.25)
            meanerr = (cerr1 + cerr2)/2.
            rmserr = np.sqrt((cerr1*cerr1 + cerr2*cerr2)/2.)

            rate = cnts[i]/e
            err_rate = rmserr/e

            r_weighted.append(rate)

            weights.append(err_rate)

        r_weighted = np.array(r_weighted)
        weights = np.array(weights)

        return r_weighted, weights

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

        print("t_start", t_start, "t_end", t_end, "num orig bins", len(self.time), "num non-zero bins",
              len(self.counts))

        frequency = np.linspace(self.f_start, self.f_stop, 100000)  # for orbital 100000

        if weights is None:
            power = LombScargle(t=self.time, y=r_weighted).power(frequency)
        else:
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

        timespan = self.time[-1] - self.time[0]
        yr_bins = 4
        delta = timespan / yr_bins    # ~ 1 yr
        yr_figs = []
        yr_exp = []
        yr_c = []
        yr_orb_power = []

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

            if weights is None:
                yr_power = LombScargle(t=yr_times, y=yr_counts).power(frequency)
            else:
                yr_power = LombScargle(t=yr_times, y=yr_counts, dy=yr_weights).power(frequency)

            freq_days = 1. / frequency / 86400.

            yr_figs.append(figure(title="time span: " + str(yr) + " " + str(tmin) + ", " + str(tmax) + " power vs frequency",
                                x_axis_label='period (days)', y_axis_label='power',
                                width=750))
            yr_figs[yr].line(freq_days, yr_power, line_width=2)
            yr_figs[yr].add_layout(vline_p1)

            yr_peaks_ls, yr_props_ls = find_peaks(yr_power, height=0.1 * max(yr_power))

            pk_days = (1. / frequency[yr_peaks_ls] / 86400.)
            close_idx = (np.abs(pk_days - self.nom_period)).argmin()
            peak_idx = yr_peaks_ls[close_idx]
            pk_error = self.calc_peak_error(frequency=frequency, power=yr_power, peak_index=peak_idx)

            print("yr", yr, "tmin", tmin, "tmax", tmax, "t[0]", yr_times[0], "t[1]", yr_times[-1],
                  "sum_weights", sum_weights, "num w bins", len(yr_weights), "num t bins", len(yr_counts))
            print(yr_peaks_ls, yr_props_ls)
            print(frequency[yr_peaks_ls])
            print(1. / frequency[yr_peaks_ls] / 86400.)
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

        f2 = figure(title="orb period power vs time bin",
                                x_axis_label='time bin', y_axis_label='power',
                                width=750)
        f2.line(np.arange(yr_bins), yr_orb_power, line_width=2)
        f2.scatter(np.arange(yr_bins), yr_orb_power, size=6, fill_color="white")

        del_div = Div(text=self.source + " Run on: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + " for " + self.fn)

        output_file(self.html_name)
        l = layout(column(del_div, f1, f2, q_hist, r_hist, s_hist, f3, column(yr_figs), column(yr_c), column(yr_exp)))
        save(l, title=self.source + " Power vs Frequency")


if __name__ == "__main__":

    p = process_LAT_binned_exposure("LSI61303")

    p.get_data()
    r_weighted, weights = p.Robin_exp_weight()
    p.make_plots()

# Close the FITS file
    p.hdul.close()

# Print the data and header
    print(p.hdul[0].header)
