from datetime import datetime
from astropy.io import fits
import math
from gti import Gti  # gift from Matthew Kerr: https://github.com/kerrm/godot/blob/master/gti.py
import numpy as np


class prepare_gti():

    def __init__(self):

        self.t0_MJD = None
        self.daysecs = 86400.
        self.orb = None

        self.infile = None
        self.h = None
        self.gti_hdr = None
        self.infile_gti = None
        self.new_gti = None
        self.t_min = None
        self.t_max = None
        self.phase_bin = None
        self.t0_MET = None

    def get_infile_stuff(self, infile):

        self.infile = infile

        self.h = fits.open(infile)
        self.h.info()
        self.gti_hdr = self.h[2]
        print(self.gti_hdr.columns)

        # Access the primary HDU (Header/Data Unit)
        primary_hdu = self.h[0]
        header = primary_hdu.header

        self.t_max = np.float128(header["TSTOP"])
        self.t_min = np.float128(header["TSTART"])

        gti_starts = self.gti_hdr.data["START"]
        gti_stops = self.gti_hdr.data["STOP"]
        self.infile_gti = Gti(gti_starts, gti_stops)

    def extrapolate_t0(self, t0_MJD, orb):

        self.t0_MJD = t0_MJD
        self.orb = orb

        if self.t0_MJD is None:
            self.t0_MET = self.t_min
            return self.t_min

        date_t0_MJD = datetime.strptime(self.t0_MJD, "%Y-%m-%d %H:%M:%S.%f")

        start_MET = datetime.strptime("2001-01-01 00:00:00.000", "%Y-%m-%d %H:%M:%S.%f")
        aug_2008 = datetime.strptime("2008-08-15 00:00:00.000", "%Y-%m-%d %H:%M:%S.%f")

        delta = np.float128((start_MET - date_t0_MJD).total_seconds())
        delta_08 = np.float128((aug_2008 - start_MET).total_seconds())
        delta_sum = delta_08 + delta
        mod_t0 = np.mod(delta_sum, orb)

        t0_post_MET = np.float128(delta_08 - mod_t0)
        self.t0_MET = t0_post_MET
        print("t0_post_MET", t0_post_MET)

        return t0_post_MET

    def make_new_gti(self, phase_bin, num_bins):

        gti = []
        p = int(phase_bin)/num_bins
        dp = 1./num_bins

        # start on orbit boundary before first time measure
        t = np.float128(self.t0_MET - math.ceil((self.t_min - self.t0_MET)/self.orb) * self.orb)

        while t < self.t_max + self.orb:
            gti.append(np.array([t+p*self.orb, t+(p+dp)*self.orb]))
            t += self.orb

        t_start = np.array(list(zip(*gti))[0])
        t_stop = np.array(list(zip(*gti))[1])

        self.new_gti = Gti(t_start, t_stop)

        return 0

    def merge_gtis(self):

        rc = self.new_gti.intersection(self.infile_gti)
        merged_starts = self.new_gti.get_edges(starts=True)
        merged_stops = self.new_gti.get_edges(starts=False)

        self.gti_hdr.header["NAXIS2"] = len(merged_starts)
        t0 = fits.Column(name='start', array=merged_starts, format='D', unit='s')
        t1 = fits.Column(name='stop', array=merged_stops, format='D', unit='s')
        cols = fits.ColDefs([t0, t1])
        hdu = fits.BinTableHDU.from_columns(cols)
        self.h[2] = hdu
        self.h[2].name = "GTI"
        self.h.writeto(self.infile, overwrite=True)
