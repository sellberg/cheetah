#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/python

# Usage:
# In this directory, type:
#    ./viewRun.py -rxxxx
# where rxxxx is the run number of hits and nonhits found using the hitfinder executable. 
# By default, this script looks into the h5 files that are in the appropriate rxxxx directory
#

import numpy as N
import h5py as H
import matplotlib
import matplotlib.pyplot as P
import sys
import os
import re
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", 
					help="run number you wish to view", metavar="rxxxx", default="")

(options, args) = parser.parse_args()


########################################################
# Edit this variable accordingly
# Files are read for source_dir/runtag and
# written to write_dir/runtag.
# Be careful of the trailing "/"; 
# ensure you have the necessary read/write permissions.
########################################################
source_dir = "/reg/d/psdm/cxi/cxi74613/scratch/cleaned_hdf5/"
#source_dir = "/reg/neh/home3/sellberg/NML-2013/analysis/cheetah_scripts/test_runs/"
#source_dir = "/reg/d/psdm/cxi/cxi74613/ftc/cleaned_hdf5/"
#write_dir = "/reg/d/psdm/cxi/cxi74613/ftc/sellberg/figures/"
write_dir = "/reg/neh/home3/sellberg/NML-2013/analysis/cheetah_scripts/figures/"


if(options.runNumber is not ""):
	print "Now examining H5 files in %sr%s/ ..."%(source_dir,options.runNumber)

runtag = "r%s"%(options.runNumber)

########################################################
# Search specified directory for *.h5 files
########################################################
searchstring="[a-zA-Z0-9\_]+"+runtag+"[a-z0-9\_]+.h5"
h5pattern = re.compile(searchstring)
h5files = [h5pattern.findall(x) for x in os.listdir(source_dir+runtag)]
h5files = [items for sublists in h5files for items in sublists]

colmax = 500
colmin = 0

########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
########################################################
class img_class (object):
	def __init__(self, inarr, filename):
		self.inarr = inarr*(inarr>0)
		for i in range(len(inarr)):
			self.inarr[i] = self.inarr[i][::-1]
		self.filename = filename
		global colmax
		global colmin
	
	def on_keypress(self,event):
		global colmax
		global colmin
		if event.key in ['1', '2', '3', '4', '5', '6', '7','8', '9', '0']:
			if not os.path.exists(write_dir + runtag):
				os.mkdir(write_dir + runtag)
			recordtag = write_dir + runtag + "/" + runtag + "_" + event.key + ".txt"
			print "recording filename in " + recordtag
			f = open(recordtag, 'a+')
			f.write(self.filename+"\n")
			f.close()
		if event.key == 'p':
			if not os.path.exists(write_dir + runtag):
				os.mkdir(write_dir + runtag)
			pngtag = write_dir + runtag + "/%s.png" % (self.filename)	
			print "saving image as " + pngtag 
			P.savefig(pngtag)
		if event.key == 'e':
			if not os.path.exists(write_dir + runtag):
				os.mkdir(write_dir + runtag)
			epstag = write_dir + runtag + "/%s.eps" % (self.filename)	
			print "saving image as " + epstag 
			P.savefig(epstag, format='eps')
		if event.key == 'r':
			colmin = self.inarr.min()
			colmax = self.inarr.max()
			P.clim(colmin, colmax)
			P.draw()


	def on_click(self, event):
		global colmax
		global colmin
		if event.inaxes:
			lims = self.axes.get_clim()
			colmin = lims[0]
			colmax = lims[1]
			range = colmax - colmin
			value = colmin + event.ydata * range
			if event.button is 1 :
				if value > colmin and value < colmax :
					colmin = value
			elif event.button is 2 :
				colmin = self.inarr.min()
				colmax = self.inarr.max()
			elif event.button is 3 :
				if value > colmin and value < colmax:
					colmax = value
			P.clim(colmin, colmax)
			P.draw()
				

	def draw_img(self):
		global colmax
		global colmin
		fig = P.figure()
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(111)
		canvas.set_title(self.filename)
		P.rc('image',origin='lower')
		self.axes = P.imshow(self.inarr, vmax = colmax, vmin = colmin)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		P.show() 

print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Press 'p' to save PNG of image (with the current colorscales) in the appropriately named folder."
print "Press any single digit '0-9' to save the H5 filename of current image to the appropriate file (e.g. r0079/r0079_1.txt) ."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

########################################################
# Loop to display all H5 files found. 
########################################################
for fname in h5files:
	f = H.File(source_dir+runtag+"/"+fname, 'r')
	d = N.array(f['/data/data'])
	f.close()
	currImg = img_class(d, fname)
	currImg.draw_img()
