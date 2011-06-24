#!/usr/bin/python

# Usage:
# In this directory, type:
#    ./viewSAXS.py -rxxxx -m 400 (optional) -x 800 (optional)
# For details, type 
#	 python viewSAXS.py --help
# where rxxxx is the run number of hits and nonhits found using the hitfinder executable. 
# By default, this script looks into the h5 files that are in the appropriate rxxxx directory
#

import os
import sys
import string
import re
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", 
					help="run number you wish to view", metavar="rxxxx", default="")
parser.add_option("-m", "--min", action="store", type="float", dest="min_value", 
					help="ignore intensities below this q-value", metavar="min_value", default="0")
parser.add_option("-x", "--max", action="store", type="float", dest="max_value", 
					help="ignore intensities above this q-value", default="100000")

(options, args) = parser.parse_args()

import numpy as N
import h5py as H

import matplotlib
import matplotlib.pyplot as P

########################################################
# Edit this variable accordingly
# Files are read for source_dir/runtag and
# written to write_dir/runtag.
# Be careful of the trailing "/"; 
# ensure you have the necessary read/write permissions.
########################################################
source_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
write_dir = ""

runtag = "r%s"%(options.runNumber)
print source_dir+runtag+"/"+runtag+"-SAXS.h5"
f = H.File(source_dir+runtag+"/"+runtag+"-SAXS.h5","r")
d = N.array(f['/data/data'])
q = d[0]
saxs = d[1]
f.close()

q[N.isnan(saxs)] = 0
q *= (q > options.min_value)
q *= (q < options.max_value)
q = q[q != 0]
saxs = saxs[q != 0]

print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

fig = P.figure()
P.plot(q, saxs)
canvas = fig.add_subplot(111)
canvas.set_title(runtag + "_SAXS")
P.xlabel("Q")
P.ylabel("I(Q)")
P.show()
