# read list of input files and loop over them with ftcopy to apply a filter selecting a phase bin
import numpy as np
import os

source = "LSI61303"
period = "26.4960"
file_list = "/Users/richarddubois/Code/GLAST/tmp/binned_events.txt"
n_phase_bins = 10
delta = 1./n_phase_bins
filter_string = []

f = open(file_list, "r")
flist = []
for input in f:
    flist.append(input.rstrip())

for p in np.arange(n_phase_bins):
    g_name = source + "_phase_" + str(p) + ".txt"
    g = open(g_name, "w")
    i = 0
    for infile in flist:
        t1 = p/n_phase_bins
        t2 = str(t1 + delta)
        t1 = str(t1)
        o_file = source + "_file_" + str(i) + "_phase_" + str(p) + ".fits"
        f_string = "ftcopy '" + infile + "[EVENTS]" \
                    "[((((TIME/86400)%" + period + ")/" + period + " >=" + t1 + ") .and. " \
                    "(((TIME/86400.)%" + period + ")/" + period + "<" + t2 + "))]' outfile=" + o_file
        filter_string.append(f_string)
        res = os.system(f_string)
        g.write(o_file + "\n")
        i += 1
    g.close()



print(len(filter_string))
