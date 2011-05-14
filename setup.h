/*
 *  setup.h
 *  cheetah
 *
 *  Created by Anton Barty on 7/2/11.
 *  Copyright 2011 CFEL. All rights reserved.
 *
 *	You can modify this software under the terms of the GNU General Public License 
 *	as published by the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 */

#ifndef _setup_h
#define _setup_h

#include <stdio.h>
#include <pthread.h>

/*
 *	Structure for hitfinder parameters
 */

class cHitfinder {

public:
	int			use;			 // specify whether to use the hitfinder	
	int			Algorithm;		 // 1,2,3 see the commented ini file or the code hitfinder.cpp
	int			ADC;			 // ADC threshold used in the hitfinders
	int			NAT;			 // Threshold on number of pixels above hitfinderADC (Algorithm 1)
	int			Npeaks;		 	 // Number of peaks above hitfinderADC (algorithm 3)
	int			NpeaksMax;
	int			MinPixCount;	 	 // Min number of pixels in peaks required for hit (Algorithm 2 and 3)
	int			MaxPixCount;		 // Max number of pixels in peaks required for hit (Algorithm 3 only)
	int			Cluster;		 // minimum number of pixels to define a cluster (Algorithm 2)
	int			UsePeakmask;		 // set if you want to look for peaks in the region defined by the peak mask
	int 		savehits;		 // set if you want to save the hits
	char		peaksearchFile[1024];			 // the name of the file containing the peak mask (Raw format)
	int16_t		*peakmask;		//stores the peakmask from the file peakmaskFile
	FILE		*cleanedfp;		// file name where the hits of this hitfinder are written.
};


/*
 *	Global variables
 */
class cGlobal {
	
public:
	
	/*
	 *	Various switches and processing options
	 */
	// ini file to read
	char		configFile[1024];		// name of the config file (e.g. cspad-cryst.ini)
	
	// Real-space geometry
	char		geometryFile[1024];		// File containing pixelmap (X,Y coordinate of each pixel in raw data stream)
	float		pixelSize;			// in micrometers
	
	// Bad pixel masks
	int			useBadPixelMask;	// to specify pixels that you know are bad
	char		badpixelFile[1024];		// file containing the map of bad pixel locations (raw format)
	
	// Static dark calibration (static offsets on each pixel to be subtracted)
	char		darkcalFile[1024];		// File containing dark calibration
	int			useDarkcalSubtraction;	// Subtract the darkcal
	int			generateDarkcal;		// Flip this on to generate a darkcal (auto-turns-on appropriate other options)
	
	// Common mode and pedastal subtraction
	int			cmModule;				// Subtract common mode from each ASIC
	int			cmSubModule;			// Subtract common mode from subsets of each ASIC (currently 16 sub-portions)
	float		cmFloor;				// Use lowest x% of values as the offset to subtract (typically lowest 2%)

	// Gain correction
	int			useGaincal;			// whether to read to gain map from a file
	int			invertGain;			// whether to invert the gain map????
	char		gaincalFile[1024];			// Name of file containing the gain map
	
	// Running background subtraction
	int			useSubtractPersistentBackground;  // if set a running background will be calculated and subtracted. 
	int			subtractBg;			  // what is this parameter?? It's also not in the ini file
	int			scaleBackground;		  // scale the running background for each shot to account for intensity fluctuations
	float		bgMemory;				  // number of frames to use for determining the running background
	
	// Kill persistently hot pixels
	int			useAutoHotpixel;		 // determine the hot pixels on the fly
	int			hotpixADC;			 // threshold above which to count as hot pixels
	int			hotpixMemory;			 // number of frames to look for hot pixels in
	float		hotpixFreq;				 // hot often a pixel needs to be above the threshold to be regarded as hot.
	
	// Attenuation correction
	int			useAttenuationCorrection;		// Whether to correct each event's intensity with the calculated attenuation
	char		attenuationFile[1024];			// Name of the file containing the attenuation list
	
	int			startFrames;			 // number of frames to use for forming initial background and hot pixel estimate (no frames outputed; digesting)
	
	// Hitfinding
	cHitfinder 	hitfinder;				// instance of the standard hitfinder
	cHitfinder	icefinder;				// hitfinder to finder patterns with ice peaks
	cHitfinder	waterfinder;				// hitfinder to identify shots with water
	cHitfinder	backgroundfinder;			// identify shots that are definitely part of the background (to make a more conservative background estimate)	

	
	// Powder pattern generation
	int			powdersum;			 // set to calculate powder pattern
	int			powderthresh;			 // pixels with an ADC value above this threshold will be added to the powder
	int			saveInterval;			 // powder pattern is repeatedly saved according to this interval
    
    // Correlation analysis
    int         useCorrelation;     // set to nonzero to turn on angular cross correlation module, also controls what correlation algorithm to be used, 1: regular, 2: fast
	int			autoCorrelationOnly;	// set to nonzero to only calculate autocorrelation (q1=q2)
	
	// Saving options
	int			saveRaw;			 // set to save each hit in raw format, in addition to assembled format. Powders are only saved in raw if saveRaw and powdersum are enabled
	int			hdf5dump;			 // set to write every frame to h5 format 
	
	// Verbosity
	int			debugLevel;			 // doesn't seem to be implemented
	
	// Log files
	char		logfile[1024];
	char		framefile[1024];
	char		cleanedfile[1024];
	char		icefile[1024];
	char 		waterfile[1024];
	char		backgroundfile[1024];
	
	
	/*
	 *	Stuff used for managing the program execution
	 */
	// Run information
	unsigned	runNumber;

	// Log file pointers
	FILE		*framefp;
	//FILE		*cleanedfp;
	
	// Thread management
	long			nThreads;
	long			nActiveThreads;
	long			threadCounter;
	pthread_t		*threadID;
	pthread_mutex_t	nActiveThreads_mutex;		// there should be one mutex variable for each global variable which threads write to.
	pthread_mutex_t	hotpixel_mutex;
	pthread_mutex_t	selfdark_mutex;
	pthread_mutex_t	powdersum1_mutex;
	pthread_mutex_t	powdersum2_mutex;
	pthread_mutex_t correlation_mutex;
    pthread_mutex_t	nhits_mutex;
	pthread_mutex_t	framefp_mutex;
	pthread_mutex_t	watersumassembled_mutex;
	pthread_mutex_t	watersumraw_mutex;
	pthread_mutex_t	icesumassembled_mutex;
	pthread_mutex_t	icesumraw_mutex;	
	
	
	// Detector geometry
	long			pix_nx;
	long			pix_ny;
	long			pix_nn;
	float			*pix_x;
	float			*pix_y;
	float			*pix_z;
	float			pix_dx;
	unsigned		module_rows;
	unsigned		module_cols;
	long			image_nx;
	long			image_nn;
	
	
	// Attenuation variables
	unsigned		nFilters;		// Counter for Si filters in XRT
	unsigned		nThicknesses;	// Counter for number of possible thicknesses
	unsigned		*filterThicknesses;		// Pointer to array of filter thicknesses
	unsigned		*possibleThicknesses;	// Pointer to array of all possible combinations of filter thicknesses
	double			*possibleAttenuations;	// Pointer to array of all possible attenuations obtained from possibleThicknesses
	double			*attenuations;	// Pointer to dynamic array of all calculated attenuations during the run
	unsigned		*changedAttenuationEvents;		// Pointer to dynamic array of all events where the attenuation changed during the run
	unsigned		*totalThicknesses;		// Pointer to dynamic array of all total thicknesses from used Si filters during the run
	unsigned		attenuationCapacity;	// Allocated size of dynamic attenuation array
	unsigned		nAttenuations;	// Number of attenuations saved in attenuation array
	int				attenuationOffset;		// Integer to compensate for the offset of nevents w.r.t. the recorded attenuations

	
	// Common variables
	int32_t			*darkcal;		//stores darkcal from the file darkcalFile
	double			*powderRaw;		//stores powder pattern in raw format, WAS int64_t
	double			*powderAssembled;	//stores the assembled powder pattern, WAS int64_t
	double			*waterRaw;		//stores powder pattern in raw format, WAS int64_t
	double			*waterAssembled;	//stores the assembled powder pattern, WAS int64_t
	double			*iceRaw;		//stores powder pattern in raw format, WAS int64_t
	double			*iceAssembled;		//stores the assembled powder pattern, WAS int64_t
	int16_t			*badpixelmask;		//stores the bad pixel mask from the file badpixelmaskFile
	float			*hotpixelmask;		//stores the hot pixel mask calculated by the auto hot pixel finder
	float			*selfdark;		//stores the background calculated by the running (persistant) background subtraction
	float			*gaincal;		//stores the gain map read from the gaincalFile
	float			avgGMD;			// what is this?
	long			npowder;		// number of frames in the powder??
	long			nprocessedframes;	// number of frames that have been processed by the worker program
	long			nhits;			// number of hits that have been found
	double			detectorZ;		// position (mm) of the detector along the beam direction
	
	clock_t			lastclock;		// variables for keeping track of time the program has been running
	timeval			lasttime;	
	float			datarate;
	time_t			tstart, tend;

	
	
public:
	void defaultConfiguration(void);		// sets default parameter values
	void parseConfigFile(char *);			// reads the config file
	void parseCommandLineArguments(int, char**);	// reads command line arguments
	void setup(void);				// function that sets parameter values from default/config file/command line arguments
	void readDetectorGeometry(char *);		// following functions read the h5 files in raw format
	void readDarkcal(char *);
	void readGaincal(char *);
	void readPeakmask(char *);
	void readBadpixelMask(char *);
	void readIcemask(char *);
	void readWatermask(char *);
	void readBackgroundmask(char *);
	void readAttenuations(char *);
	void expandAttenuationCapacity();

	void writeInitialLog(void);			// functions to write the log file
	void updateLogfile(void);			
	void writeFinalLog(void);

	
private:
	void parseConfigTag(char*, char*);		// parses the ini file or command line options. New input commands must be added to this function.
	
};

#endif
