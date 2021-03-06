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
#include <string>
#include <vector>

/*
 *	Structure for hitfinder parameters
 */

class cHitfinder {

public:
	int			use;			 // specify whether to use the hitfinder	
	int			Algorithm;		 // 1,2,3 see the commented ini file or the code hitfinder.cpp
	int			ADC;			 // ADC threshold used in the hitfinders
	int			NAT;			 // Threshold on number of pixels above hitfinderADC (Algorithm 1)
	int			Npeaks;		 	 // Min number of peaks above hitfinderADC (algorithm 3)
	int			NpeaksMax;			 // Max number of peaks above hitfinderADC (algorithm 3)
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
	float		pixelSize;			// in meters
	double		detectorOffset;		// offset in nominal detector position (mm), follows CXI coordinate system; positive value (+Z) moves detector further downstream
	double		detectorZpos;			// nominal detectorZ position if it fails to read from data stream
	char		detectorZpvname[1024];	// PV name of detectorZ motor
	int			useCenterCorrection;	// set to nonzero to shift 0 of pixel array to the center of the square hole
	float		pixelCenterX;		// option to manually override the calculated center of the square hole (if useCenterCorrection is enabled) and set the center to X in pixels if pixelCenterX is nonzero
	float		pixelCenterY;		// option to manually override the calculated center of the square hole (if useCenterCorrection is enabled) and set the center to Y in pixels if pixelCenterY is nonzero
	int			calculateCenterCorrectionPowder;	// set to nonzero to calculate center of pixel array from powder pattern using a Hough transform, and applies it to powder averages after the correction has been calculated
	int			calculateCenterCorrectionHit;	// set to nonzero to calculate center of pixel array from each hit using a Hough transform, and applies it after the correction has been calculated
	int			calculateCenterCorrectionQuad;	// set to nonzero to calculate center of each quad in pixel array from powder pattern using a Hough transform, and applies it after the correction has been calculated, also controls what correction algorithm to be used: 1 = unified center, 2 = unified radius
	double		centerCorrectionThreshold;	// determines the threshold for which below pixel intensities are set to zero
	double		centerCorrectionMaxR;	// maximum radius value in pixels that the hough transform will be performed for, the calibration ring should lie BELOW this value.
	double		centerCorrectionMinR;	// minimum radius value in pixels that the hough transform will be performed for, the calibration ring should lie ABOVE this value.
	double		centerCorrectionDeltaR;	// step length in pixels for the radius
	double		centerCorrectionMaxC;	// maximum value in pixels of the shift in center that the hough transform will be performed for, the true center should lie inside +/- this value.
	double		centerCorrectionDeltaC;	// step length in pixels for the center
	int			useMetrologyRefinement;	// set to nonzero to shift X/Y values of the quads
	float		quad0DX;
	float		quad0DY;
	float		quad1DX;
	float		quad1DY;
	float		quad2DX;
	float		quad2DY;
	float		quad3DX;
	float		quad3DY;
	int			refineMetrology;	// refine metrology by translating/rotating quads w.r.t. each other, also switches between algorithms: 1 = angular average calculated for each quad separately, 2 = angular average calculated for all quads together
	float		refinementMaxC;			// maximum value in pixels of the shift in the quad center that the refinement will be performed for, the true quad center should lie inside +/- this value.
	float		refinementDeltaC;		// step length in pixels for the quad center
		
	// Bad pixel masks
	int			useBadPixelMask;	// to specify pixels that you know are bad
	char		badpixelFile[1024];		// file containing the map of bad pixel locations (raw format)
	
	// Static dark calibration (static offsets on each pixel to be subtracted)
	char		darkcalFile[1024];		// File containing dark calibration
	int			useDarkcalSubtraction;	// Subtract the darkcal
	int			generateDarkcal;		// Flip this on to generate a darkcal (auto-turns-off appropriate other options)
	int			manualDarkcalGenerationControl;			// Flip this on to manually control the other settings when generating a darkcal
	
	// Common mode and pedastal subtraction
	int			cmModule;				// Subtract common mode from each ASIC, can be done by finding the zero-photon peak (algorithm 1) in the intensity histogram or by using the cmFloor-weighted median (algorithm 2)
	int			cmSubModule;			// Subtract common mode from subsets of each ASIC (only implemented using algorithm 2), the value of cmSubModule defines how many times each dimension of each ASIC will be divided by and should be 1 or a multiple of 2
	int			cmStart;				// Algorithm 1: intensity (ADU) from which the peakfinding should start in the histogram
	int			cmStop;					// Algorithm 1: intensity (ADU) at which the peakfinding should stop in the histogram
	float		cmDelta;				// Algorithm 1: noise threshold intensity (ADU) over which the peakfinding should consider as true peaks in the histogram
	float		cmFloor;				// Algorithm 2: use lowest x% of values as the offset to subtract (typically lowest 2%)
	int			cmSaveHistograms;		// Save intensity histograms for each 2x1 Module. Histograms are saved into separate files for each event. The median defined by cmFloor is saved in the first element unless terminal output states otherwise
	
	// Gain correction
	int			useGaincal;			// whether to read to gain map from a file
	int			invertGain;			// inverts the gain map, this is the standard since gain maps from flat fields should be divided with the measured intensity
	int			normalizeGain;		// normalizes the gain map, so that a flat-field average can be used as gain calibration
	char		gaincalFile[1024];			// Name of file containing the gain map
	
	// Running background subtraction
	int			useSubtractPersistentBackground;  // if set a running background will be calculated and subtracted. 
	int			scaleBackground;		  // scale the running background for each shot to account for intensity fluctuations
	float		bgMemory;				  // number of frames to use for determining the running background
	
	// Kill persistently hot pixels
	int			useAutoHotpixel;		 // determine the hot pixels on the fly
	int			hotpixADC;			 // threshold above which to count as hot pixels
	int			hotpixMemory;			 // number of frames to look for hot pixels in
	float		hotpixFreq;				 // hot often a pixel needs to be above the threshold to be regarded as hot.
	
	// Polarization correction
	int			usePolarizationCorrection;	// set to nonzero to calculate and apply polarization correction to each hit (or saved event)
	double		horizontalPolarization;		// determines the ratio of horizontal polarization versus vertical polarization of the FEL beam, valid values are [0,1]
	
	// Solid angle correction
	int			useSolidAngleCorrection;	// set to nonzero to calculate and apply solid angle correction to each hit (or saved event), only the theta-dependent part of the correction is applied to the 2D scattering pattern, also controls what correction algorithm to be used: 1 = rigorous, 2 = azimuthally symmetric
	
	// Attenuation correction
	int			useAttenuationCorrection;		// Whether to correct each event's intensity with the calculated attenuation, this also toggles if the intensity should be corrected for the FEE gas detectors readouts (in mJ)
	float		filterPositionIn;				// Nominal position for silicon filters in beam
	float		filterPositionOut;				// Nominal position for silicon filters out of beam
	char		attenuationFile[1024];			// Name of the file containing the attenuation list
	
	// Energy calibration
	int			useEnergyCalibration;		// Save histogram of energies and wavelengths for energy calibration
	
	// Average intensity statitiscs
	int			useIntensityStatistics;		// Save histogram of avg intensities
	
	// Single-pixel statitiscs
	int			usePixelStatistics;		// Save intensities for randomly chosen single pixels (currently 26 pixels with high, medium, and low variance)
	char		pixelFile[1024];		// Name of the file containing the pixel list
	
	// Save and flush settings
	int			saveInterval;			 // powder pattern and log is repeatedly saved according to this interval
	int			flushInterval;			// HDF5 data and memory is repeatedly flushed according to this interval
	int			startFrames;			 // number of frames to use for forming initial background and hot pixel estimate (no frames outputed; digesting)
	
	// Hitfinding
	cHitfinder 	hitfinder;				// instance of the standard hitfinder
	cHitfinder	icefinder;				// hitfinder to find patterns with ice peaks
	cHitfinder	waterfinder;				// hitfinder to identify shots with water
	cHitfinder	backgroundfinder;			// identify shots that are definitely part of the background (to make a more conservative background estimate)	
	
	// Listfinding
	cHitfinder	listfinder;				// instance of the hitfinder for using a list as input criteria
	char		listfinderFile[1024];	// Name of the file containing the hit list
	
	// Powder pattern generation
	int			powdersum;			 // set to calculate powder pattern
	int			powderthresh;			 // pixels with an ADC value above this threshold will be added to the powder
	    
	// Angular averages
	int			powderAngularAvg;			// set to nonzero to calculate angular averages of the powder patterns
	int			hitAngularAvg;			// set to nonzero to calculate angular averages of each hit
	double		angularAvgStartQ;		// start position for angular average calculated in pixel from the center
	double		angularAvgStopQ;			// end position for angular average calculated in pixel from the center
	double		angularAvgDeltaQ;			// binning of angular averages, DeltaQ currently determines step size in pixels
	
    // Correlation analysis
    int         useCorrelation;     // set to nonzero to turn on angular cross-correlation module, also controls what correlation algorithm to be used, 1: regular, 2: fast
	int			sumCorrelation;		// set to nonzero to sum cross-correlation patterns for different hits
	int			autoCorrelateOnly;		// set to nonzero to only calculate autocorrelation (q1=q2)
	int			correlationNormalization;	// switch between normalization algorithms: 1 = intensity, 2 = variance
	int			correlationQScale;		// customize the scale for the Q-range of the correlation algorithms
    double 		correlationStartQ;		// 1 = startQ/stopQ are in units of radial distance [pixels]
    double 		correlationStopQ;		// 2 = startQ/stopQ are in absolute units of momentum transfer [Å-1], 3 = startQ/stopQ are in absolute units of perpendicualar momentum transfer [Å-1]
    int 		correlationNumQ;		// number of q values between start and stop
    double 		correlationStartPhi;	// start angle in degrees, default: 0
    double 		correlationStopPhi;		// stop angle in degrees, default: 360
    int 		correlationNumPhi;		// number of angular steps, default: 256 (attention: if possible, use powers of 2, that makes FFT especially fast)
    int 		correlationNumDelta;	// number of angular lag steps, default: 0 (it then calculates suitable number from correlationNumPhi with the same step length)
	int		 	*correlationLUT;		// lookup table (LUT) needed for the fast correlation
	int			correlationLUTdim1;		// dim1 of LUT
	int			correlationLUTdim2;		// dim2 of LUT
	int			correlationOutput;		// switch between output formats: 1 = hdf5, 2 = bin, 3 = hdf5+bin, 4 = tiff, 5 = tiff+hdf5, 6 = tiff+bin, 7 = tiff+hdf5+bin
	
	// Saving options
	int			saveRaw;			 // set to save each hit in raw format, in addition to assembled format. Powders are only saved in raw if saveRaw and powdersum are enabled
	int			hdf5dump;			 // set to write every frame to h5 format 
	
	// Verbosity
	int			debugLevel;			 // set to 0 for regular, 1 for debug mode, and 2 for extra verbose
	
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
	pthread_mutex_t	nActiveThreads_mutex;		// there should be one mutex variable for each global variable which threads write to.
	pthread_mutex_t	hotpixel_mutex;
	pthread_mutex_t	selfdark_mutex;
	pthread_mutex_t	intensities_mutex;
	pthread_mutex_t	powdersumraw_mutex;
	pthread_mutex_t	powdersumassembled_mutex;
	pthread_mutex_t	powdersumcorrelation_mutex;
	pthread_mutex_t	powdersumvariance_mutex;
	pthread_mutex_t	watersumassembled_mutex;
	pthread_mutex_t	watersumraw_mutex;
	pthread_mutex_t	watersumcorrelation_mutex;
	pthread_mutex_t	icesumassembled_mutex;
	pthread_mutex_t	icesumraw_mutex;
	pthread_mutex_t	icesumcorrelation_mutex;
	pthread_mutex_t correlation_mutex;
	pthread_mutex_t correlationFFT_mutex;
	pthread_mutex_t pixelcenter_mutex;
	pthread_mutex_t image_mutex;
    pthread_mutex_t	nhits_mutex;
	pthread_mutex_t	framefp_mutex;
	
	
	// Detector geometry
	long			pix_nx;
	long			pix_ny;
	long			pix_nn;
	float			*pix_x;
	float			*pix_y;
	float			*pix_z;
	double			*pix_r;
	float			pix_xmax;
	float			pix_xmin;
	float			pix_ymax;
	float			pix_ymin;
	float			pix_rmax;
	unsigned		module_rows;
	unsigned		module_cols;
	long			image_nx;
	long			image_nn;
	float			*quad_dx;
	float			*quad_dy;
	
	
	// Angular average variables
	unsigned	angularAvg_nn;	// length of angular average arrays
	int			*angularAvg_i;	// stores index for each pixel in angular average array
	double		*angularAvgQ;	// stores q-values (in pixels) for angular averages	
	double		*angularAvgQcal;	// stores absolute q-values calibrated from the mean photon energy of the run
	
	
	// Polarization correction variables
	double		*phi;
	
	
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
	
	
	// Energy calibration variables
	double			*energies;	// Dynamic array of photon energies (eV) for all events
	double			*wavelengths;	// Dynamic array of wavelengths (Å) for all events
	unsigned		energyCapacity;	// Starting capacity of dynamic arrays
	unsigned		nEnergies;	// Number of energies saved in dynamic array
	double			Emin;	// Lowest photon energy
	double			Emax;	// Highest photon energy
	double			Emean;	// Mean photon energy
	double			Lmin;	// Lowest wavelength
	double			Lmax;	// Highest wavelength
	double			Lmean;	// Mean wavelength
	unsigned		*Ehist;	// Histogram of energies
	unsigned		*Lhist;	// Histogram of wavelengths
	
	
	// Intensity statistics variables
	double			*intensities;		// Dynamic array of avg intensities for all events
	bool			*hits;		// Dynamic array of hits for all events
	unsigned		intensityCapacity;	// Allocated size of dynamic intensity array
	unsigned		nIntensities;		// Number of avg intensities saved in dynamic intensity array
	unsigned		*Ihist;	// Histogram of intensities
	unsigned		*IHhist;// Histogram of intensities of hits
	double			Imin;	// Lowest avg intensity
	double			Imax;	// Highest avg intensity
	double			Imean;	// Mean avg intensity
	
	
	// Pixel statistics variables
	unsigned		*pixels;		// Pointer to dynamic array of all pixel indices in the raw data format
	double			*pixelXYList;	// Pointer to dynamic array of all pixel x/y values imported from python scripts in (1480,1552) format
	unsigned		pixelCapacity;	// Allocated size of dynamic pixel array
	unsigned		nPixels;		// Number of pixels read in dynamic pixel array
	
	
	// Listfinding variables
	bool			eventIsHit;				// keeps track of whether the event is a hit contained in the list
	std::vector<std::string>	hitlist;	// list of all hits as output string names
	
	// Common variables
	float			*darkcal;		// stores darkcal from the file darkcalFile in float format
	double			*powderRaw;		// stores powder pattern in raw format
	double			*powderAssembled;	// stores the assembled powder pattern
	double			*powderAverage;		// stores angular average of powder pattern
	double			*powderCorrelation;	// stores correlation sum of regular hits
	double			*powderVariance;	// stores the variance of the powder pattern
	double			*waterRaw;		// stores powder pattern of water hits in raw format
	double			*waterAssembled;	// stores the assembled powder pattern of water hits
	double			*waterAverage;		// stores angular average of powder pattern of water hits
	double			*waterCorrelation;	// stores correlation sum of water hits
	double			*iceRaw;		// stores powder pattern of ice hits in raw format
	double			*iceAssembled;		// stores the assembled powder pattern of ice hits
	double			*iceAverage;		// stores angular average of powder pattern	of ice hits
	double			*iceCorrelation;	// stores correlation sum of ice hits
	int16_t			*badpixelmask;		// stores the bad pixel mask from the file badpixelmaskFile
	float			*hotpixelmask;		// stores the hot pixel mask calculated by the auto hot pixel finder
	float			*selfdark;		// stores the background calculated by the running (persistant) background subtraction
	float			*gaincal;		// stores the gain map read from the gaincalFile
	float			avgGMD;			// what is this?
	long			npowder;		// number of frames in the powder
	long			nwater;		// number of frames in the water powder
	long			nice;		// number of frames in the ice powder
	long			nprocessedframes;	// number of frames that have been processed by the worker program
	long			nhits;			// number of hits that have been found
	long			correlation_nn;	// length of global cross-correlation arrays
	double			detectorZ;		// position (mm) of the detector along the beam direction
	float			detposold;		// the detector position of the second last event, this makes sure the detector doesn't artificially 'jump' between events due to bug in PV readout
	
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
	float pixelCenter(float pixel_array[]); // help function for readDetectorGeometry to calculate center of pixel array
	void shiftQuads(float xarray[], float dx[], float yarray[], float dy[]); // help function for readDetectorGeometry to shift quads w.r.t. each other
	void readDarkcal(char *);
	void readGaincal(char *);
	void readPeakmask(char *);
	void readBadpixelMask(char *);
	void readIcemask(char *);
	void readWatermask(char *);
	void readBackgroundmask(char *);
	void readHits(char *);
	void readAttenuations(char *);		// read in list of attenuations
	void expandAttenuationCapacity();
	void expandEnergyCapacity();
	void expandIntensityCapacity();
	void readPixels(char *);			// read in list of pixels to be analyzed on a single-pixel basis
	void expandPixelCapacity();
	void createLookupTable();			// create lookup table (LUT) needed for the fast correlation algorithm	

	void writeInitialLog(void);			// functions to write the log file
	void updateLogfile(void);			
	void writeFinalLog(void);

	
private:
	void parseConfigTag(char*, char*);		// parses the ini file or command line options. New input commands must be added to this function.
	
};

#endif
