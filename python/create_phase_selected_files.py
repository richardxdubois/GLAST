# read list of input files and loop over them with ftcopy to apply a filter selecting a phase bin
import numpy as np
import os
import yaml
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--app_config',
                    default="/Users/richarddubois/Code/GLAST/tmp/dbg/phase_config.yaml",
                    help="overall app config file")
args = parser.parse_args()

with open(args.app_config, "r") as f:
    data = yaml.safe_load(f)

source = data["source"]

n_phase_bins = data["n_phase_bins"]
phase_bin = data["phase_bin"]
t0 = data["t0"]

period = data["period"]

infile = data["infile"]

delta = 1./n_phase_bins
filter_string = []

t1 = phase_bin/n_phase_bins
t2 = str(t1 + delta)
t1 = str(t1)
o_file = "ft1_phase.fits"
f_string = "ftcopy '" + infile + "[EVENTS]" \
            "[(((((TIME-" + str(t0) + ")/86400)%" + period + ")/" + period + " >=" + t1 + ") .and. " \
            "((((TIME-" + str(t0) + ")/86400.)%" + period + ")/" + period + "<" + t2 + "))]' outfile=" + o_file
print(f_string)
res = os.system(f_string)
