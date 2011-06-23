#!/usr/bin/python

# Usage:
# In this directory, type:
#    ./viewEnergyHistograms.py -rxxxx -m 400 (optional) -x 800 (optional)
# For details, type 
#	 python viewEnergyHistograms.py --help
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
print source_dir+runtag+"/"+runtag+"-energy_histograms.h5"
f = H.File(source_dir+runtag+"/"+runtag+"-energy_histograms.h5","r")
d = N.array(f['/data/data'])
energies = d[0]
xe = d[1]
wavelengths = d[2]
xw = d[3]
f.close()

energies[N.isnan(energies)] = 0
wavelengths[N.isnan(wavelengths)] = 0

print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

fig = P.figure()
P.plot(xe, energies)
canvas = fig.add_subplot(111)
canvas.set_title(runtag + "_energies")
P.xlabel("E (eV)")
P.ylabel("events")
P.draw()

fig = P.figure()
P.plot(xw, wavelengths)
canvas = fig.add_subplot(111)
canvas.set_title(runtag + "_wavelengths")
P.xlabel("lambda (A)")
P.ylabel("events")
P.show()
