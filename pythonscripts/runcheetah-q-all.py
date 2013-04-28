#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/python

# Usage:
# In this directory, type:
#    ./runcheetah-q-all.py -d FILENAME (optional)
# For details, type 
#	 python runcheetah-q-all.py --help
# where FILENAME is the name of the filelist where all the darkcals are listed
# By default, this script looks into the hdf5 files that are in the appropriate darkcal directory
#

import numpy as N
import h5py as H
import glob as G
import matplotlib
import matplotlib.pyplot as P
from pylab import *
import scipy
import scipy.interpolate as I
from scipy import *
from scipy import optimize
import sys, os, re, shutil, subprocess, time
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-s", "--start", action="store", type="int", dest="startrun",
		  help="start run of experiment (default: 1)", metavar="RUNNUMBER", default=1)
parser.add_option("-e", "--end", action="store", type="int", dest="endrun",
		  help="end run of experiment (default: 270)", metavar="RUNNUMBER", default=270)
parser.add_option("-d", "--dark", action="store", type="string", dest="darklist",
		  help="input filelist of darkcals (default: filelist.dat)", metavar="FILENAME", default="filelist.dat")
parser.add_option("-o", "--omit", action="store", type="string", dest="skiplist",
		  help="input filelist of runs to omit", metavar="FILENAME", default="")
parser.add_option("-b", "--break", action="store", type="int", dest="submit_break",
		  help="number of runs after which the script takes a break from submitting new jobs (default: never)", metavar="RUNS", default=0)
parser.add_option("-t", "--time", action="store", type="float", dest="submit_break_time",
		  help="number of minutes that the break lasts (default: 20 min)", metavar="MINUTES", default=20.0)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
		  help="prints out additional information", default=False)
(options, args) = parser.parse_args()

########################################################
# Edit this variable accordingly
# Files are read for dark_dir and ini_dir,
# run from run_dir and written to write_dir.
# Be careful of the trailing "/"; 
# ensure you have the necessary read/write permissions.
########################################################
dark_dir = "/reg/d/psdm/cxi/cxi65913/res/config/darkcal/"
ini_dir = "/reg/d/psdm/cxi/cxi65913/res/config/ini/"
run_dir = "/reg/d/psdm/cxi/cxi65913/res/cheetah_scripts/"
write_dir = "/reg/d/psdm/cxi/cxi65913/ftc/cleaned_hdf5/"

runs = range(options.startrun, options.endrun + 1)
if options.verbose:
	print "Will process r%04d - r%04d" % (options.startrun, options.endrun)
darkruns = []
darkrunnumber = []

print "Reading darkcals from %s ..." % (dark_dir+options.darklist)
f = open(dark_dir+options.darklist, 'r')
for line in f:
	if line[0] != '#':
		line = re.sub("\n", '', line)
		darkrun = int(re.sub("r", '', line))
		darkruns.append(line)
		darkrunnumber.append(darkrun)

f.close()

skipruns = []
skiprunnumber = []
if options.skiplist != "":
	print "Reading runs to skip from %s ..." % (options.skiplist)
	f = open(options.skiplist, 'r')
	for line in f:
		if line[0] != '#':
			line = re.sub("\n", '', line)
			skiprun = int(re.sub("r", '', line))
			skipruns.append(line)
			skiprunnumber.append(skiprun)
	f.close()

print "Reading configs from %s ..." % (ini_dir)
#files = G.glob( ini_dir + '*' )
#files = [os.path.basename(f) for f in files]
#start, end = os.path.basename(f).split('.')[0].split('-')
searchstring="cheetah_r+[0-9]+-+[0-9]+.ini"
inipattern = re.compile(searchstring)
inifiles = [inipattern.findall(x) for x in os.listdir(ini_dir)]
inifiles = [items for sublists in inifiles for items in sublists]

sDarkruns = set(darkruns)
sSkipruns = set(skipruns)
submitted_jobs = 0
ini_last_used = None
for run in runs:
	runtag = "r%04d" % run
	if (runtag not in sDarkruns) and (runtag not in sSkipruns):
		ini_to_use = None
		for ini in inifiles:
			start,end = ini.split('.')[0].split('_')[1].split('-')
			startnumber = int(re.sub('r','',start))
			endnumber = int(end)
			if startnumber > endnumber:
				print "r%04d must be smaller or equal to r%04d, ignoring %s." % (startnumber, endnumber, ini)
			else:
				if (run >= startnumber) and (run <= endnumber):
					os.system("cp " + ini_dir + ini + " " + run_dir + "cheetah.ini")
					ini_to_use = ini
					break
		if ini_to_use == None:
			print "No config file present for %s, aborting." % runtag
			sys.exit(1)
		elif ini_to_use != ini_last_used:
			print "Using config file: %s" % ini_to_use
		cmd = "./runcheetah-q %s" % runtag
		print cmd
		os.system(cmd)
		submitted_jobs += 1
		ini_last_used = ini_to_use
		
		output = ""
		while output != ">-------- Start of job --------<":
			time.sleep(2)
			if os.path.exists(write_dir+runtag+"/log.txt"):
				cmd = ["tail", write_dir+runtag+"/log.txt"]
				output = subprocess.Popen( cmd, stdout=subprocess.PIPE ).communicate()[0].split('\n')[1]
				if options.verbose:
					print output
				if (options.submit_break != 0) and (submitted_jobs % options.submit_break == 0):
					print "Submitted %d runs, taking a %d min break..." % (submitted_jobs, int(options.submit_break_time))
					time.sleep(60*options.submit_break_time)
					print "\t... resuming script after break"
	elif options.verbose:
		print "Skipping %s" % runtag


