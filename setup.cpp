/*
 *  setup.cpp
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
 

#include "myana/myana.hh"
#include "myana/main.hh"
#include "myana/XtcRun.hh"
#include "release/pdsdata/cspad/ConfigV1.hh"
#include "release/pdsdata/cspad/ConfigV2.hh"
#include "release/pdsdata/cspad/ConfigV3.hh"
#include "release/pdsdata/cspad/ElementHeader.hh"
#include "release/pdsdata/cspad/ElementIterator.hh"
#include "cspad-gjw/CspadTemp.hh"
#include "cspad-gjw/CspadCorrector.hh"
#include "cspad-gjw/CspadGeometry.hh"


#include <string.h>
#include <sys/time.h>
#include <cmath>
#include <hdf5.h>
#include <fenv.h>
#include <stdlib.h>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::string;

#include <fstream>

#include "setup.h"
#include "worker.h"
#include "data2d.h"
#include "attenuation.h"
#include "arrayclasses.h"


/*
 *	Default settings/configuration
 */
void cGlobal::defaultConfiguration(void) {

	// ini file to use
	strcpy(configFile, "cheetah.ini");

	// Geometry
	strcpy(geometryFile, "geometry/cspad_pixelmap.h5");
	pixelSize = 110e-6;
	useCenterCorrection = 0;
	pixelCenterX = 0;
	pixelCenterY = 0;
	calculateCenterCorrectionPowder = 0;
	calculateCenterCorrectionHit = 0;
	centerCorrectionThreshold = 150;
	centerCorrectionMaxR = 600;
	centerCorrectionMinR = 400;
	centerCorrectionDeltaR = 1;
	centerCorrectionMaxC = 50;
	centerCorrectionDeltaC = 1;
	refineMetrology = 0;
	refinementMaxC = 20;
	refinementDeltaC = 1;
	
	// Bad pixel mask
	strcpy(badpixelFile, "badpixels.h5");
	useBadPixelMask = 0;

	// Static dark calibration (electronic offsets)
	strcpy(darkcalFile, "darkcal.h5");
	useDarkcalSubtraction = 1;
	generateDarkcal = 0;
	
	// Common mode subtraction from each ASIC
	cmModule = 0;
	cmSubModule = 0;
	cmFloor = 0.1;

	// Gain calibration correction
	strcpy(gaincalFile, "gaincal.h5");
	useGaincal = 0;
	invertGain = 0;
	
	// Subtraction of running background (persistent photon background) 
	useSubtractPersistentBackground = 0;
	subtractBg = 0;
	bgMemory = 50;
	startFrames = 0;
	scaleBackground = 0;
	
	// Kill persistently hot pixels
	useAutoHotpixel = 1;
	hotpixFreq = 0.9;
	hotpixADC = 1000;
	hotpixMemory = 50;
	
	// Attenuation correction
	useAttenuationCorrection = 0;
	strcpy(attenuationFile, "attenuations.dat");
	
	// Energy calibration
	useEnergyCalibration = 0;
	
	// Hitfinding
	
	hitfinder.use = 0;
	hitfinder.ADC = 100;
	hitfinder.NAT = 100;
	hitfinder.Npeaks = 50;
	hitfinder.NpeaksMax = 100000;
	hitfinder.Algorithm = 3;
	hitfinder.MinPixCount = 3;
	hitfinder.MaxPixCount = 20;
	hitfinder.UsePeakmask = 0;
	hitfinder.savehits = 0;
	strcpy(hitfinder.peaksearchFile, "peakmask.h5");

	waterfinder.use = 0;
	waterfinder.ADC = 100;
	waterfinder.NAT = 100;
	waterfinder.Npeaks = 50;
	waterfinder.NpeaksMax = 100000;
	waterfinder.Algorithm = 3;
	waterfinder.MinPixCount = 3;
	waterfinder.MaxPixCount = 20;
	waterfinder.UsePeakmask = 0;
	waterfinder.savehits = 0;
	strcpy(waterfinder.peaksearchFile, "watermask.h5");

	icefinder.use = 0;
	icefinder.ADC = 100;
	icefinder.NAT = 100;
	icefinder.Npeaks = 50;
	icefinder.NpeaksMax = 100000;
	icefinder.Algorithm = 3;
	icefinder.MinPixCount = 3;
	icefinder.MaxPixCount = 20;
	icefinder.UsePeakmask = 0;
	icefinder.savehits = 0;
	strcpy(icefinder.peaksearchFile, "icemask.h5");

	backgroundfinder.use = 0;
	backgroundfinder.ADC = 100;
	backgroundfinder.NAT = 100;
	backgroundfinder.Npeaks = 50;
	backgroundfinder.NpeaksMax = 100000;
	backgroundfinder.Algorithm = 3;
	backgroundfinder.MinPixCount = 3;
	backgroundfinder.MaxPixCount = 20;
	backgroundfinder.UsePeakmask = 0;
	backgroundfinder.savehits = 0;
	strcpy(backgroundfinder.peaksearchFile, "backgroundmask.h5");

	// Powder pattern generation
	powdersum = 1;
	powderthresh = 0;
	
	// Angular averages of powder patterns
	powderSAXS = 0;
	powderStartQ = 0;
	powderStopQ = 0;
	deltaqSAXS = 1;
    
	// Correlation analysis
	useCorrelation = 0;
	sumCorrelation = 0;
	autoCorrelationOnly = 1;
    correlationStartQ = 100;
    correlationStopQ = 600;
    correlationNumQ = 51;
    correlationStartPhi = 0;
    correlationStopPhi = 360;
    correlationNumPhi = 256;
	correlationNumDelta = 0;
	correlationLUTdim1 = 100;
	correlationLUTdim2 = 100;
	
	// Saving options
	saveRaw = 0;
	hdf5dump = 0;
	saveInterval = 500;
	
	// Verbosity
	debugLevel = 0;
	
	// Default to only a few threads
	nThreads = 16;
	
	// Log files
	strcpy(logfile, "log.txt");
	strcpy(framefile, "frames.txt");
	strcpy(cleanedfile, "cleaned.txt");
	strcpy(icefile, "icehits.txt");
	strcpy(waterfile, "waterhits.txt");
	strcpy(backgroundfile, "backgroundhits.txt");
	
}



/*
 *	Setup stuff to do with thread management, settings, etc.
 */
void cGlobal::setup() {
	/*
	 *	Set up arrays for remembering powder data, background, etc.
	 */
	
	//initially, set everything to NULL, allocate only those that are actually needed below
	powderRaw = NULL;
	powderAssembled = NULL;
	powderAverage = NULL;
	powderQ = NULL;
	iceRaw = NULL;
	iceAssembled = NULL;
	iceAverage = NULL;
	waterRaw = NULL;
	waterAssembled = NULL;		
	waterAverage = NULL;
	correlation_nn = 0;
	powderCorrelation = NULL;
	iceCorrelation = NULL;
	waterCorrelation = NULL;
	quad_dx = NULL;
	quad_dy = NULL;
	
	selfdark = (float*) calloc(pix_nn, sizeof(float));
	
	if (useAutoHotpixel) 
		hotpixelmask = (float*) calloc(pix_nn, sizeof(float));
	
	if (powdersum) {
		if (hitfinder.use) 
			powderAssembled = (double*) calloc(image_nn, sizeof(double));
		if (icefinder.use) 
			iceAssembled = (double*) calloc(image_nn, sizeof(double));
		if (waterfinder.use) 
			waterAssembled = (double*) calloc(image_nn, sizeof(double));
	
		if (saveRaw) {
			if (hitfinder.use) 
				powderRaw = (double*) calloc(pix_nn, sizeof(double));
			if (icefinder.use) 
				iceRaw = (double*) calloc(pix_nn, sizeof(double));
			if (waterfinder.use) 
				waterRaw = (double*) calloc(pix_nn, sizeof(double));
		}
		if (powderSAXS) {
			powderQ = (double*) calloc(powder_nn, sizeof(double));
			if (hitfinder.use) 
				powderAverage = (double*) calloc(powder_nn, sizeof(double));
			if (icefinder.use) 
				iceAverage = (double*) calloc(powder_nn, sizeof(double));
			if (waterfinder.use) 
				waterAverage = (double*) calloc(powder_nn, sizeof(double));
		}
	}//powdersum end
		
	if (useCorrelation) {
		if (!correlationNumDelta) 
			correlationNumDelta = (int) ceil(correlationNumPhi/2.0+1);
		if (autoCorrelationOnly) {
			correlation_nn = correlationNumQ*correlationNumDelta;
		} else {
			correlation_nn = correlationNumQ*correlationNumQ*correlationNumDelta;
		}
		if (sumCorrelation) {
			if (hitfinder.use) 
				powderCorrelation = (double*) calloc(correlation_nn, sizeof(double));
			if (icefinder.use) 
				iceCorrelation = (double*) calloc(correlation_nn, sizeof(double));
			if (waterfinder.use) 
				waterCorrelation = (double*) calloc(correlation_nn, sizeof(double));
		}
	}// useCorrelation end
	
	if (refineMetrology) {
		quad_dx = (float *) calloc(4, sizeof(float));
		quad_dy = (float *) calloc(4, sizeof(float));
	}
	
	
	/*
	 *	Set up thread management
	 */
	nActiveThreads = 0;
	threadCounter = 0;
	pthread_mutex_init(&nActiveThreads_mutex, NULL);
	pthread_mutex_init(&hotpixel_mutex, NULL);
	pthread_mutex_init(&selfdark_mutex, NULL);
	pthread_mutex_init(&powdersumraw_mutex, NULL);
	pthread_mutex_init(&powdersumassembled_mutex, NULL);
	pthread_mutex_init(&powdersumcorrelation_mutex, NULL);
	pthread_mutex_init(&icesumraw_mutex, NULL);
	pthread_mutex_init(&icesumassembled_mutex, NULL);
	pthread_mutex_init(&icesumcorrelation_mutex, NULL);
	pthread_mutex_init(&watersumraw_mutex, NULL);
	pthread_mutex_init(&watersumassembled_mutex, NULL);
	pthread_mutex_init(&watersumcorrelation_mutex, NULL);
	pthread_mutex_init(&correlation_mutex, NULL);
	pthread_mutex_init(&correlationFFT_mutex, NULL);
	pthread_mutex_init(&pixelcenter_mutex, NULL);
	pthread_mutex_init(&image_mutex, NULL);
    pthread_mutex_init(&nhits_mutex, NULL);
	pthread_mutex_init(&framefp_mutex, NULL);

	threadID = (pthread_t*) calloc(nThreads, sizeof(pthread_t));

	
	/*
	 *	Trap specific configurations and mutually incompatible options
	 */
	if(generateDarkcal) {
		cmModule = 0;
		cmSubModule = 0;
		subtractBg = 0;
		useDarkcalSubtraction = 0;
		useSubtractPersistentBackground = 0;
		hitfinder.use = 0;
		waterfinder.use = 0;
		icefinder.use = 0;
		backgroundfinder.use = 0;
		hitfinder.savehits = 0;
		waterfinder.savehits = 0;
		icefinder.savehits = 0;
		backgroundfinder.savehits = 0;		
		hdf5dump = 0;
		saveRaw = 0;
		useAutoHotpixel = 0;
		startFrames = 0;
		powderthresh = 0;
		powderSAXS = 0;
		calculateCenterCorrectionPowder = 0;
		calculateCenterCorrectionHit = 0;
		useAttenuationCorrection = -1;
		useEnergyCalibration = 0;
		useCorrelation = 0;
	}
	
	/*
	 *	Other stuff
	 */
	npowder = 0;
	nwater = 0;
	nice = 0;
	nprocessedframes = 0;
	nhits = 0;
	lastclock = clock()-10;
	gettimeofday(&lasttime, NULL);
	datarate = 1;
	detectorZ = 0;
	runNumber = getRunNumber();
	time(&tstart);
	avgGMD = 0;
	
	/*
	 *	Setup global attenuation variables
	 */
	if (useAttenuationCorrection >= 0) {
		nFilters = 10; // Counter for Si filters in XRT
		filterThicknesses = new unsigned[nFilters]; // Array of filter thicknesses
		for (int i=0; i<nFilters; i++) {
			filterThicknesses[i] = int(20*pow((float)2.,i));
		}
		nThicknesses = 1; // Counter for number of possible thicknesses (0 um => counter starts at 1)
		for (int i=1; i<=nFilters; i++) {
			nThicknesses += factorial(nFilters)/(factorial(nFilters-i)*factorial(i));
		} // 1023 combinations added
		possibleThicknesses = new unsigned[nThicknesses]; // Array of all possible combinations of filter thicknesses
		for (int i=0; i<=nThicknesses; i++) {
			possibleThicknesses[i] = 20*i;
		}
		possibleAttenuations = new double[nThicknesses]; // Array of all possible attenuations obtained from possibleThicknesses
		attenuationCapacity = 100; // Starting capacity of dynamic arrays
		attenuations = new double[attenuationCapacity]; // Dynamic array of all calculated attenuations during the run
		changedAttenuationEvents = new unsigned[attenuationCapacity]; // Dynamic array of all events where the attenuation changed during the run
		totalThicknesses = new unsigned[attenuationCapacity]; // Dynamic array of all total thicknesses from used Si filters during the run
		nAttenuations = 0; // Number of attenuations saved in attenuation array
		attenuationOffset = 0; // Integer to compensate for the offset of nevents w.r.t. the recorded attenuations
	} else {
		filterThicknesses = NULL;
		possibleThicknesses = NULL;
		possibleAttenuations = NULL;
		attenuations = NULL;
		changedAttenuationEvents = NULL;
		totalThicknesses = NULL;
		nFilters = 0;
		nThicknesses = 0;
		attenuationCapacity = 0;
		nAttenuations = 0;
		attenuationOffset = 0;
	}
	
	/*
	 *	Setup global energy calibration variables
	 */
	if (useEnergyCalibration) {
		energyCapacity = 1000; // Starting capacity of dynamic arrays
		energies = new double[energyCapacity]; // Dynamics array of all photon energies (eV)
		wavelengths = new double[energyCapacity]; // Dynamics array of all wavelengths (Å)
		nEnergies = 0;
		Emin = 100000;	// Lowest photon energy
		Emax = 0;	// Highest photon energy
		Lmin = 100000;	// Lowest wavelength
		Lmax = 0;	// Highest wavelength
		Ehist = NULL;	// Histograms are allocated in makeEnergyHistograms()
		Lhist = NULL;	// Histograms are allocated in makeEnergyHistograms()
	} else {
		energies = NULL;
		wavelengths = NULL;
		Ehist = NULL;
		Lhist = NULL;
		energyCapacity = 0;
		nEnergies = 0;
		Emin = 100000;
		Emax = 0;
		Lmin = 100000;
		Lmax = 0;
	}
	
	/*
	 *	Setup global cross correlation variables
	 */
	correlationLUT = new int[correlationLUTdim1*correlationLUTdim2];
	
	// Make sure to use SLAC timezone!
	setenv("TZ","US/Pacific",1);
	
}




/*
 *	Parse command line arguments 
 */
void cGlobal::parseCommandLineArguments(int argc, char **argv) {
	
	// No arguments specified = ask for help
	if (argc == 1) {
		printf("No configuration file spcified\n");
		printf("\t--> using default settings\n");
		return;
	}
	
	// First argument is always the configuration file
	parseConfigFile(argv[1]);
	
	// Other arguments are optional switches but take same syntax prefixed by an '-'
	if (argc > 2) {
		for (long i=2; i<argc; i++) {
			if (argv[i][0] == '-' && i+1 < argc) {
				parseConfigTag(argv[i]+1, argv[++i]);
			}
		}
	}
}




/*
 *	Read and process configuration file
 */
void cGlobal::parseConfigFile(char* filename) {
	char		cbuf[cbufsize];
	char		tag[cbufsize];
	char		value[cbufsize];
	char		*cp;
	FILE		*fp;
	
	
	/*
	 *	Open configuration file for reading
	 */
	printf("Parsing input configuration file:\t%s\n",filename);
	
	fp = fopen(filename,"r");
	if (fp == NULL) {
		printf("\tCould not open configuration file \"%s\"\n",filename);
		printf("\tUsing default values\n");
		return;
	}
	
	/*
	 *	Loop through configuration file until EOF 
	 *	Ignore lines beginning with a '#' (comments)
	 *	Split each line into tag and value at the '=' sign
	 */
	while (feof(fp) == 0) {
		
		cp = fgets(cbuf, cbufsize, fp);
		if (cp == NULL) 
			break;
		
		if (cbuf[0] == '#')
			continue;
		
		cp = strpbrk(cbuf, "=");
		if (cp == NULL)
			continue;
		
		*(cp) = '\0';
		sscanf(cp+1,"%s",value);
		sscanf(cbuf,"%s",tag);
		
		parseConfigTag(tag, value);
	}
	
	fclose(fp);
	
}

/*
 *	Process tags for both configuration file and command line options
 */
void cGlobal::parseConfigTag(char *tag, char *value) {
	
	/*
	 *	Convert to lowercase
     *  --> config file is not case sensitive!
	 */
	for(int i=0; i<strlen(tag); i++) 
		tag[i] = tolower(tag[i]);
	
	/*
	 *	Parse known tags
	 */
	if (!strcmp(tag, "nthreads")) {
		nThreads = atoi(value);
	}
	else if (!strcmp(tag, "geometry")) {
		strcpy(geometryFile, value);
	}
	else if (!strcmp(tag, "darkcal")) {
		strcpy(darkcalFile, value);
	}
	else if (!strcmp(tag, "gaincal")) {
		strcpy(gaincalFile, value);
	}
	else if (!strcmp(tag, "peakmask")) {
		strcpy(hitfinder.peaksearchFile, value);
	}
	else if (!strcmp(tag, "badpixelmask")) {
		strcpy(badpixelFile, value);
	}
	else if (!strcmp(tag, "attenuation")) {
		strcpy(attenuationFile, value);
	}
	
	// Processing options
	else if (!strcmp(tag, "usecentercorrection")) {
		useCenterCorrection = atoi(value);
	}
	else if (!strcmp(tag, "pixelcenterx")) {
		pixelCenterX = atof(value);
	}
	else if (!strcmp(tag, "pixelcentery")) {
		pixelCenterY = atof(value);
	}
	else if (!strcmp(tag, "calculatecentercorrectionpowder")) {
		calculateCenterCorrectionPowder = atoi(value);
	}
	else if (!strcmp(tag, "calculatecentercorrectionhit")) {
		calculateCenterCorrectionHit = atoi(value);
	}
	else if (!strcmp(tag, "centercorrectionthreshold")) {
		centerCorrectionThreshold = atof(value);
	}
	else if (!strcmp(tag, "centercorrectionmaxr")) {
		centerCorrectionMaxR = atof(value);
	}
	else if (!strcmp(tag, "centercorrectionminr")) {
		centerCorrectionMinR = atof(value);
	}
	else if (!strcmp(tag, "centercorrectiondeltar")) {
		centerCorrectionDeltaR = atof(value);
	}
	else if (!strcmp(tag, "centercorrectionmaxc")) {
		centerCorrectionMaxC = atof(value);
	}
	else if (!strcmp(tag, "centercorrectiondeltac")) {
		centerCorrectionDeltaC = atof(value);
	}
	else if (!strcmp(tag, "refinemetrology")) {
		refineMetrology = atoi(value);
	}
	else if (!strcmp(tag, "refinementmaxc")) {
		refinementMaxC = atof(value);
	}
	else if (!strcmp(tag, "refinementdeltac")) {
		refinementDeltaC = atof(value);
	}
	else if (!strcmp(tag, "subtractcmmodule")) {
		cmModule = atoi(value);
	}
	else if (!strcmp(tag, "cmmodule")) {
		cmModule = atoi(value);
	}
	else if (!strcmp(tag, "usegaincal")) {
		useGaincal = atoi(value);
	}
	else if (!strcmp(tag, "invertgain")) {
		invertGain = atoi(value);
	}
	else if (!strcmp(tag, "subtractcmsubmodule")) {
		cmSubModule = atoi(value);
	}
	else if (!strcmp(tag, "generatedarkcal")) {
		generateDarkcal = atoi(value);
	}
	else if (!strcmp(tag, "subtractbg")) {
		subtractBg = atoi(value);
	}
	else if (!strcmp(tag, "usebadpixelmask")) {
		useBadPixelMask = atoi(value);
	}
	else if (!strcmp(tag, "usedarkcalsubtraction")) {
		useDarkcalSubtraction = atoi(value);
	}
	else if (!strcmp(tag, "hitfinder")) {
		hitfinder.use = atoi(value);
	}
	else if (!strcmp(tag, "savehits")) {
		hitfinder.savehits = atoi(value);
	}
	else if (!strcmp(tag, "powdersum")) {
		powdersum = atoi(value);
	}
	else if (!strcmp(tag, "powdersaxs")) {
		powderSAXS = atoi(value);
	}
	else if (!strcmp(tag, "powderstartq")) {
		powderStartQ = atof(value);
	}
	else if (!strcmp(tag, "powderstopq")) {
		powderStopQ = atof(value);
	}
	else if (!strcmp(tag, "deltaqsaxs")) {
		deltaqSAXS = atof(value);
	}
	else if (!strcmp(tag, "usecorrelation")) {
		useCorrelation = atoi(value);
	}
	else if (!strcmp(tag, "sumcorrelation")) {
		sumCorrelation = atoi(value);
	}
	else if (!strcmp(tag, "autocorrelationonly")) {
		autoCorrelationOnly = atoi(value);
	}
    else if (!strcmp(tag, "correlationstartq")) {
		correlationStartQ = atof(value);
	}
    else if (!strcmp(tag, "correlationstopq")) {
		correlationStopQ = atof(value);
	}
    else if (!strcmp(tag, "correlationnumq")) {
		correlationNumQ = atoi(value);
	}
    else if (!strcmp(tag, "correlationstartphi")) {
		correlationStartPhi = atof(value);
	}
    else if (!strcmp(tag, "correlationstopphi")) {
		correlationStopPhi = atof(value);
	}
    else if (!strcmp(tag, "correlationnumphi")) {
		correlationNumPhi = atoi(value);
	}
    else if (!strcmp(tag, "correlationnumdelta")) {
		correlationNumDelta = atoi(value);
	}
    else if (!strcmp(tag, "correlationlutdim1")) {
		correlationLUTdim1 = atoi(value);
	}
    else if (!strcmp(tag, "correlationlutdim2")) {
		correlationLUTdim2 = atoi(value);
	}
	else if (!strcmp(tag, "saveraw")) {
		saveRaw = atoi(value);
	}
	else if (!strcmp(tag, "hdf5dump")) {
		hdf5dump = atoi(value);
	}
	else if (!strcmp(tag, "saveinterval")) {
		saveInterval = atoi(value);
	}
	else if (!strcmp(tag, "useautohotpixel")) {
		useAutoHotpixel = atoi(value);
	}
	else if (!strcmp(tag, "useselfdarkcal")) {
		useSubtractPersistentBackground = atoi(value);
	}
	else if (!strcmp(tag, "usesubtractpersistentbackground")) {
		useSubtractPersistentBackground = atoi(value);
	}
	else if (!strcmp(tag, "useattenuationcorrection")) {
		useAttenuationCorrection = atoi(value);
	}
	else if (!strcmp(tag, "useenergycalibration")) {
		useEnergyCalibration = atoi(value);
	}
	
	// Power user settings
	else if (!strcmp(tag, "cmfloor")) {
		cmFloor = atof(value);
	}
	else if (!strcmp(tag, "pixelsize")) {
		pixelSize = atof(value);
	}
	else if (!strcmp(tag, "debuglevel")) {
		debugLevel = atoi(value);
	}
	else if (!strcmp(tag, "hotpixfreq")) {
		hotpixFreq = atof(value);
	}
	else if (!strcmp(tag, "hotpixadc")) {
		hotpixADC = atoi(value);
	}
	else if (!strcmp(tag, "hotpixmemory")) {
		hotpixMemory = atoi(value);
	}
	else if (!strcmp(tag, "powderthresh")) {
		powderthresh = atoi(value);
	}
	else if (!strcmp(tag, "hitfinderadc")) {
		hitfinder.ADC = atoi(value);
	}
	else if (!strcmp(tag, "hitfindernat")) {
		hitfinder.NAT = atoi(value);
	}
	else if (!strcmp(tag, "hitfindercluster")) {
		hitfinder.Cluster = atoi(value);
	}
	else if (!strcmp(tag, "hitfindernpeaks")) {
		hitfinder.Npeaks = atoi(value);
	}
	else if (!strcmp(tag, "hitfindernpeaksmax")) {
		hitfinder.NpeaksMax = atoi(value);
	}
	else if (!strcmp(tag, "hitfinderalgorithm")) {
		hitfinder.Algorithm = atoi(value);
	}
	else if (!strcmp(tag, "hitfinderminpixcount")) {
		hitfinder.MinPixCount = atoi(value);
	}
	else if (!strcmp(tag, "hitfindermaxpixcount")) {
		hitfinder.MaxPixCount = atoi(value);
	}
	else if (!strcmp(tag, "hitfinderusepeakmask")) {
		hitfinder.UsePeakmask = atoi(value);
	}
	else if (!strcmp(tag, "selfdarkmemory")) {
		bgMemory = atof(value);
	}
	else if (!strcmp(tag, "bgmemory")) {
		bgMemory = atof(value);
	}
	else if (!strcmp(tag, "scalebackground")) {
		scaleBackground = atoi(value);
	}
	else if (!strcmp(tag, "scaledarkcal")) {
		scaleBackground = atoi(value);
	}
	else if (!strcmp(tag, "startframes")) {
		startFrames = atoi(value);
	}
	
	/* 	
	 *	Tags for water, ice, background finder
	 */
	else if (!strcmp(tag, "icemask")) {
		strcpy(icefinder.peaksearchFile, value);	
	}	
	else if (!strcmp(tag, "icefinder")) {
		icefinder.use = atoi(value);
	}
	else if (!strcmp(tag, "icefinderadc")) {
		icefinder.ADC = atoi(value);
	}
	else if (!strcmp(tag, "icefindernat")) {
		icefinder.NAT = atoi(value);
	}
	else if (!strcmp(tag, "icefindercluster")) {
		icefinder.Cluster = atoi(value);
	}
	else if (!strcmp(tag, "icefindernpeaks")) {
		icefinder.Npeaks = atoi(value);
	}
	else if (!strcmp(tag, "icefindernpeaksmax")) {
		icefinder.NpeaksMax = atoi(value);
	}
	else if (!strcmp(tag, "icefinderalgorithm")) {
		icefinder.Algorithm = atoi(value);
	}
	else if (!strcmp(tag, "icefinderminpixcount")) {
		icefinder.MinPixCount = atoi(value);
	}
	else if (!strcmp(tag, "icefindermaxpixcount")) {
		icefinder.MaxPixCount = atoi(value);
	}
	else if (!strcmp(tag, "icefinderusepeakmask")) {
		icefinder.UsePeakmask = atoi(value);
	}
	else if (!strcmp(tag, "saveicehits")) {
		icefinder.savehits = atoi(value);
	}

	else if (!strcmp(tag, "watermask")) {
		strcpy(waterfinder.peaksearchFile, value);	
	}
	else if (!strcmp(tag, "waterfinder")) {
		waterfinder.use = atoi(value);
	}
	else if (!strcmp(tag, "waterfinderadc")) {
		waterfinder.ADC = atoi(value);
	}
	else if (!strcmp(tag, "waterfindernat")) {
		waterfinder.NAT = atoi(value);
	}
	else if (!strcmp(tag, "waterfindercluster")) {
		waterfinder.Cluster = atoi(value);
	}
	else if (!strcmp(tag, "waterfindernpeaks")) {
		waterfinder.Npeaks = atoi(value);
	}
	else if (!strcmp(tag, "waterfindernpeaksmax")) {
		waterfinder.NpeaksMax = atoi(value);
	}
	else if (!strcmp(tag, "waterfinderalgorithm")) {
		waterfinder.Algorithm = atoi(value);
	}
	else if (!strcmp(tag, "waterfinderminpixcount")) {
		waterfinder.MinPixCount = atoi(value);
	}
	else if (!strcmp(tag, "waterfindermaxpixcount")) {
		waterfinder.MaxPixCount = atoi(value);
	}
	else if (!strcmp(tag, "waterfinderusepeakmask")) {
		waterfinder.UsePeakmask = atoi(value);
	}
	else if (!strcmp(tag, "savewaterhits")) {
		waterfinder.savehits = atoi(value);
	}

	else if (!strcmp(tag, "backgroundmask")) {
		strcpy(backgroundfinder.peaksearchFile, value);	
	}
	else if (!strcmp(tag, "backgroundfinder")) {
		backgroundfinder.use = atoi(value);
	}
	else if (!strcmp(tag, "backgroundfinderadc")) {
		backgroundfinder.ADC = atoi(value);
	}
	else if (!strcmp(tag, "backgroundfindernat")) {
		backgroundfinder.NAT = atoi(value);
	}
	else if (!strcmp(tag, "backgroundfindercluster")) {
		backgroundfinder.Cluster = atoi(value);
	}
	else if (!strcmp(tag, "backgroundfindernpeaks")) {
		backgroundfinder.Npeaks = atoi(value);
	}
	else if (!strcmp(tag, "backgroundfindernpeaksmax")) {
		backgroundfinder.NpeaksMax = atoi(value);
	}
	else if (!strcmp(tag, "backgroundfinderalgorithm")) {
		backgroundfinder.Algorithm = atoi(value);
	}
	else if (!strcmp(tag, "backgroundfinderminpixcount")) {
		backgroundfinder.MinPixCount = atoi(value);
	}
	else if (!strcmp(tag, "backgroundfindermaxpixcount")) {
		backgroundfinder.MaxPixCount = atoi(value);
	}
	else if (!strcmp(tag, "backgroundfinderusepeakmask")) {
		backgroundfinder.UsePeakmask = atoi(value);
	}
	else if (!strcmp(tag, "savebackgroundhits")) {
		backgroundfinder.savehits = atoi(value);
	}
	
	// Unknown tags
	else {
		printf("\tUnknown tag (ignored): %s = %s\n",tag,value);
	}
}



/*
 *	Read in detector configuration
 */
void cGlobal::readDetectorGeometry(char* filename) {
	

	// Pixel size (measurements in geometry file are in m)
	module_rows = ROWS;
	module_cols = COLS;
	pix_dx = pixelSize;

	
	// Set filename here 
	printf("Reading detector configuration:\n");
	printf("\t%s\n",filename);
	
	
	// Check whether pixel map file exists!
	FILE* fp = fopen(filename, "r");
	if (fp) 	// file exists
		fclose(fp);
	else {		// file doesn't exist
		printf("Error: Detector configuration file does not exist: %s\n",filename);
		exit(1);
	}
	
	
	// Read pixel locations from file
	cData2d		detector_x;
	cData2d		detector_y;
	cData2d		detector_z;
	detector_x.readHDF5(filename, (char *) "x");
	detector_y.readHDF5(filename, (char *) "y");
	detector_z.readHDF5(filename, (char *) "z");
	
	// Sanity check that all detector arrays are the same size (!)
	if (detector_x.nn != detector_y.nn || detector_x.nn != detector_z.nn) {
		printf("readDetectorGeometry: array size mismatch\n");
		exit(1);
	}
	

	// Sanity check that size matches what we expect for cspad (!)
	if (detector_x.nx != 8*ROWS || detector_x.ny != 8*COLS) {
		printf("readDetectorGeometry: array size mismatch\n");
		printf("%ix%i != %ix%i\n", 8*ROWS, 8*COLS, (int)detector_x.nx, (int)detector_x.ny);
		exit(1);
	}
	
	
	// Create local arrays for detector pixel locations
	long	nx = 8*ROWS;
	long	ny = 8*COLS;
	long	nn = nx*ny;
	pix_nx = nx;
	pix_ny = ny;
	pix_nn = nn;
	pix_x = (float *) calloc(nn, sizeof(float));
	pix_y = (float *) calloc(nn, sizeof(float));
	pix_z = (float *) calloc(nn, sizeof(float));
	printf("\tPixel map is %li x %li pixel array\n",nx,ny);
	
	
	// Copy values from 2D array
	for(long i=0;i<nn;i++){
		pix_x[i] = (float) detector_x.data[i];
		pix_y[i] = (float) detector_y.data[i];
		pix_z[i] = (float) detector_z.data[i];
	}
	
	
	// Divide array (in m) by pixel size to get pixel location indices (ijk)
	for(long i=0;i<nn;i++){
		pix_x[i] /= pix_dx;
		pix_y[i] /= pix_dx;
		pix_z[i] /= pix_dx;
	}
	
	
	// Center correct the array w.r.t the square hole created by the quads (assume beam is centered)
	if (useCenterCorrection && !calculateCenterCorrectionPowder && !calculateCenterCorrectionHit) {
		if (debugLevel >= 1) cout << "\tX values:" << endl;
		float x0 = pixelCenter(pix_x);
		if (debugLevel >= 1) cout << "\tY values:" << endl;
		float y0 = pixelCenter(pix_y);
		if (pixelCenterX || pixelCenterY) {
			x0 = pixelCenterX;
			y0 = pixelCenterY;
		}
		
		cout << "\tCorrected center (x,y): (" << x0 << ", " << y0 << ")" << endl;
		
		for (int i=0; i<nn; i++) {
			pix_x[i] -= x0;
			pix_y[i] -= y0;
		}
	}
	
	
	// Find bounds of image array
	float	xmax = -1e9;
	float	xmin =  1e9;
	float	ymax = -1e9;
	float	ymin =  1e9;
	float	rmax = 0;
	for(long i=0;i<nn;i++){
		if (pix_x[i] > xmax) xmax = pix_x[i];
		if (pix_x[i] < xmin) xmin = pix_x[i];
		if (pix_y[i] > ymax) ymax = pix_y[i];
		if (pix_y[i] < ymin) ymin = pix_y[i];
		float rtemp = sqrt(pix_x[i]*pix_x[i] + pix_y[i]*pix_y[i]);
		if (rtemp > rmax) rmax = rtemp;
	}
	
	// Initialize global variables
	pix_xmax = xmax;
	pix_xmin = xmin;
	pix_ymax = ymax;
	pix_ymin = ymin;
	pix_rmax = rmax;
	if (powderStopQ == 0 || powderStopQ < powderStartQ) {
		powder_nn = (unsigned) round(rmax/deltaqSAXS)+1;
		powderStartQ = 0;
	} else {
		powder_nn = (unsigned) round((powderStopQ-powderStartQ)/deltaqSAXS)+1;
	}
	
	//xmax = ceil(xmax);
	//xmin = floor(xmin);
	//ymax = ceil(ymax);
	//ymin = floor(ymin);
		
	fesetround(1);
	xmax = lrint(xmax);
	xmin = lrint(xmin);
	ymax = lrint(ymax);
	ymin = lrint(ymin);
	printf("\tImage bounds:\n");
	printf("\tx range %f to %f\n",xmin,xmax);
	printf("\ty range %f to %f\n",ymin,ymax);
	
	
	// How big must the output image be?
	float max = xmax;
	if(ymax > max) max = ymax;
	if(fabs(xmin) > max) max = fabs(xmin);
	if(fabs(ymin) > max) max = fabs(ymin);
	image_nx = 2*(unsigned)max;
	image_nn = image_nx*image_nx;
	printf("\tImage output array will be %i x %i\n",(int)image_nx,(int)image_nx);
}


/*
 *	Help function for readDetectorGeometry to calculate center of pixel array
 */
float cGlobal::pixelCenter( float *pixel_array ) {
	float center = 0;
	int quads = 4;
	float dq1 = 0;
	float dq2 = 0;
	// Loop over quads and pick out closest pixel to center
	for (int i=0; i<quads; i++) {
		if (debugLevel >= 2) cout << "\tQ" << i << ",S1: " << pixel_array[8*ROWS*(2*COLS-1)+i*2*ROWS] << endl;
		if (debugLevel >= 2) cout << "\tQ" << i << ",S2: " << pixel_array[8*ROWS*2*COLS+i*2*ROWS] << endl;
		center += pixel_array[8*ROWS*(2*COLS-1)+i*2*ROWS];
		if (i == 0) dq1 += pixel_array[8*ROWS*(2*COLS-1)+i*2*ROWS];
		else if (i == 1) dq2 += pixel_array[8*ROWS*(2*COLS-1)+i*2*ROWS];
		else if (i == 2) dq1 -= pixel_array[8*ROWS*(2*COLS-1)+i*2*ROWS];
		else dq2 -= pixel_array[8*ROWS*(2*COLS-1)+i*2*ROWS];
	}
	if (debugLevel >= 1) cout << "\tGap(Q0-Q2) = " << dq1 << endl;
	if (debugLevel >= 1) cout << "\tGap(Q1-Q3) = " << dq2 << endl;
	return center/quads;
}


/*
 *	Read in darkcal file
 */
void cGlobal::readDarkcal(char *filename){
	
	printf("Reading darkcal configuration:\n");
	printf("\t%s\n",filename);
	
	
	// Create memory space and pad with zeros
	darkcal = (int32_t*) calloc(pix_nn, sizeof(int32_t));
	memset(darkcal,0, pix_nn*sizeof(int32_t));
	
	
	
	// Check whether pixel map file exists!
	FILE* fp = fopen(filename, "r");
	if (fp) 	// file exists
		fclose(fp);
	else {		// file doesn't exist
		printf("\tDarkcal file does not exist: %s\n",filename);
		printf("\tDefaulting to all-zero darkcal\n");
		return;
	}
	
	
	// Read darkcal data from file
	cData2d		temp2d;
	temp2d.readHDF5(filename);
	
	// Correct geometry?
	if(temp2d.nx != pix_nx || temp2d.ny != pix_ny) {
		printf("\tGeometry mismatch: %ix%x != %ix%i\n",(int)temp2d.nx, (int)temp2d.ny, (int)pix_nx, (int)pix_ny);
		printf("\tDefaulting to all-zero darkcal\n");
		return;
	} 
	
	// Copy into darkcal array
	for(long i=0;i<pix_nn;i++)
		darkcal[i] = (int32_t) temp2d.data[i];
	
}


/*
 *	Read in gaincal file
 */
void cGlobal::readGaincal(char *filename){
	
	printf("Reading detector gain calibration:\n");
	printf("\t%s\n",filename);
	
	
	// Create memory space and set default gain to 1 everywhere
	gaincal = (float*) calloc(pix_nn, sizeof(float));
	for(long i=0;i<pix_nn;i++)
		gaincal[i] = 1;
	
		
	// Check whether gain calibration file exists!
	FILE* fp = fopen(filename, "r");
	if (fp) 	// file exists
		fclose(fp);
	else {		// file doesn't exist
		printf("\tGain calibration file does not exist: %s\n",filename);
		printf("\tDefaulting to uniform gaincal\n");
		return;
	}
	
	
	// Read darkcal data from file
	cData2d		temp2d;
	temp2d.readHDF5(filename);
	

	// Correct geometry?
	if(temp2d.nx != pix_nx || temp2d.ny != pix_ny) {
		printf("\tGeometry mismatch: %ix%x != %ix%i\n",(int)temp2d.nx, (int)temp2d.ny, (int)pix_nx, (int)pix_ny);
		printf("\tDefaulting to uniform gaincal\n");
		return;
	} 
	
	
	// Copy into gaincal array
	for(long i=0;i<pix_nn;i++)
		gaincal[i] = (float) temp2d.data[i];


	// Invert the gain so we have an array that all we need to do is simple multiplication
	// Pixels with zero gain become dead pixels
	if(invertGain) {
		for(long i=0;i<pix_nn;i++) {
			if(gaincal[i] != 0)
				gaincal[i] = 1.0/gaincal[i];
			else 
				gaincal[i] = 0;
		}
	}
	
}


/*
 *	Read in peaksearch mask
 */
void cGlobal::readPeakmask(char *filename){
	
	printf("Reading peak search mask:\n");
	printf("\t%s\n",filename);
	
	
	// Create memory space and default to searching for peaks everywhere
	hitfinder.peakmask = (int16_t*) calloc(pix_nn, sizeof(int16_t));
	for(long i=0;i<pix_nn;i++)
		hitfinder.peakmask[i] = 1;
	
	
	// Check whether file exists!
	FILE* fp = fopen(filename, "r");
	if (fp) 	// file exists
		fclose(fp);
	else {		// file doesn't exist
		printf("\tPeak search mask does not exist: %s\n",filename);
		printf("\tDefaulting to uniform search mask\n");
		return;
	}
	
	
	// Read peakmask data from file
	cData2d		temp2d;
	temp2d.readHDF5(filename);
	
	
	// Correct geometry?
	if(temp2d.nx != pix_nx || temp2d.ny != pix_ny) {
		printf("\tGeometry mismatch: %ix%x != %ix%i\n",(int)temp2d.nx, (int)temp2d.ny, (int)pix_nx, (int)pix_ny);
		printf("\tDefaulting to uniform peak search mask\n");
		return;
	} 
	
	
	// Copy into peakmask array
	for(long i=0;i<pix_nn;i++)
		hitfinder.peakmask[i] = (int16_t) temp2d.data[i];
}


/*
 *	Read in bad pixel mask
 */
void cGlobal::readBadpixelMask(char *filename){
	
	printf("Reading bad pixel mask:\n");
	printf("\t%s\n",filename);
	
	
	// Create memory space and default to searching for peaks everywhere
	badpixelmask = (int16_t*) calloc(pix_nn, sizeof(int16_t));
	for(long i=0;i<pix_nn;i++)
		badpixelmask[i] = 1;
	
	
	// Check whether file exists!
	FILE* fp = fopen(filename, "r");
	if (fp) 	// file exists
		fclose(fp);
	else {		// file doesn't exist
		printf("\tPeak search mask does not exist: %s\n",filename);
		printf("\tDefaulting to uniform search mask\n");
		return;
	}
	
	
	// Read darkcal data from file
	cData2d		temp2d;
	temp2d.readHDF5(filename);
	
	
	// Correct geometry?
	if(temp2d.nx != pix_nx || temp2d.ny != pix_ny) {
		printf("\tGeometry mismatch: %ix%x != %ix%i\n",(int)temp2d.nx, (int)temp2d.ny, (int)pix_nx, (int)pix_ny);
		printf("\tDefaulting to uniform peak search mask\n");
		return;
	} 
	
	
	// Copy back into array
	for(long i=0;i<pix_nn;i++)
		badpixelmask[i] = (int16_t) temp2d.data[i];
}


/*
 *	Read in list of attenuations
 */
void cGlobal::readAttenuations(char *filename) {
	
	printf("Reading list of attenuations:\n");
	printf("\t%s\n",filename);
	
	std::ifstream infile;
	infile.open(filename);
	if (infile.fail()) {
		cout << "\tUnable to open " << filename << endl;
		infile.clear();
		printf("\tReading default list of attenuations\n");
		infile.open("attenuations.dat");
		if (infile.fail()) {
			cout << "\tUnable to open default file" << endl;
			infile.clear();
			printf("\tDisabling atteunation calculation completely\n");
			useAttenuationCorrection = -1;
			return;
		}
	}

	string line;
	int counter = 0;
	while (true) {
		getline(infile, line);
		if (infile.fail()) break;
		if (line[0] != '#') {
			if (!fromString(possibleAttenuations[counter], line)) {
				cout << "\tConversion of string to double failed" << endl;
			}
			counter++;
		}
	}
}


/*
 *	Expand capacity of dynamic attenuation arrays
 */
void cGlobal::expandAttenuationCapacity() {
	attenuationCapacity *= 2;
	double *oldAttenuations = attenuations;
	attenuations = new double[attenuationCapacity];
	unsigned *oldChangedAttenuations = changedAttenuationEvents;
	changedAttenuationEvents = new unsigned[attenuationCapacity];
	unsigned *oldTotalThicknesses = totalThicknesses;
	totalThicknesses = new unsigned[attenuationCapacity];		
	for (int i=0; i<nAttenuations; i++) {
		attenuations[i] = oldAttenuations[i];
		changedAttenuationEvents[i] = oldChangedAttenuations[i];
		totalThicknesses[i] = oldTotalThicknesses[i];
	}
	delete[] oldAttenuations;
	delete[] oldChangedAttenuations;
	delete[] oldTotalThicknesses;
}


/*
 *	Expand capacity of dynamic attenuation arrays
 */
void cGlobal::expandEnergyCapacity() {
	energyCapacity *= 2;
	double *oldEnergies = energies;
	energies = new double[energyCapacity];
	double *oldWavelengths = wavelengths;
	wavelengths = new double[energyCapacity];
	for (int i=0; i<nEnergies; i++) {
		energies[i] = oldEnergies[i];
		wavelengths[i] = oldWavelengths[i];
	}
	delete[] oldEnergies;
	delete[] oldWavelengths;
}


/*
 *	Create lookup table (LUT) needed for the fast correlation algorithm
 */
void cGlobal::createLookupTable(){	
	//write lookup table for fast cross-correlation
	int lutNx = correlationLUTdim1;
	int lutNy = correlationLUTdim2;
	int lutSize = lutNx*lutNy;
	cout << "Creating lookup table (LUT) of size " << lutNx << " x " << lutNy << " (" << lutSize << " entries)" << endl;
	if (correlationLUT) {		// free memory of old LUT first, if necessary
  		delete[] correlationLUT;
	}
	correlationLUT = new int[lutSize];
	
	//initialize to zero!
	for (int i=0; i<lutSize; i++) {
  		correlationLUT[i] = 0;
	}
	
	//find max and min of the q-calibration arrays
	double qx_range = fabs(pix_xmax - pix_xmin);
    double qx_stepsize = qx_range/(double)(lutNx-1);
	double qy_range = fabs(pix_ymax - pix_ymin);
    double qy_stepsize = qy_range/(double)(lutNy-1);
	cout << "\tLUT qx-values: min=" << pix_xmin << ", max=" << pix_xmax << ", range=" << qx_range << ", stepsize=" << qx_stepsize << endl;
	cout << "\tLUT qy-values: min=" << pix_ymin << ", max=" << pix_ymax << ", range=" << qy_range << ", stepsize=" << qy_stepsize << endl;
	
	int lutFailCount = 0;
	for (int i = 0; i < pix_nn; i++){           //go through the whole the q-calibration array
												//get q-values from qx and qy arrays
												//and determine at what index (ix, iy) to put them in the lookup table
		const double ix = (pix_x[i]-pix_xmin) / qx_stepsize;
		const double iy = (pix_y[i]-pix_ymin) / qy_stepsize;
		
		//fill table at the found coordinates with the data index
		//overwriting whatever value it had before
	    //(the add-one-half->floor trick is to achieve reasonable rounded integers)
		int lutindex = (int) ( floor(ix+0.5) + lutNx*floor(iy+0.5) );
		//int lutindex = (int) ( floor(ix) + lutNx*floor(iy) );

		if (lutindex < 0 || lutindex >= lutSize) {
			cerr << endl;
  			cerr << "Error in cGlobal::createLookupTable! lookup table index out of bounds" << endl;
			cerr << "\t (was trying to set correlationLUT[" << lutindex << "] = " << i << ")" << endl;
			cerr << "\tLUT: pix_x, pix_y = (" << pix_x[i] << ", " << pix_y[i] << ") --> ix, iy = (" << ix << ", " << iy << ")" << endl;
			lutFailCount++;
		} else {
			correlationLUT[lutindex] = i;
		}
			
		/////////////////////////////////////////////////////////////////////////////////////////
		//ATTENTION: THIS METHOD WILL LEAD TO A LOSS OF DATA,
		//ESPECIALLY FOR SMALL TABLE SIZES,
		//BUT IT WILL BUY A LOT OF SPEED IN THE LOOKUP PROCESS
		//--> this should be improved to a more precise version, 
		//    maybe even one that allows the lookup(x,y) to interpolate
		//    for that to work, we need to find the four closest points in the data or so
		//    (for instance, instead of one index, the table could contain 
		//    a vector of all applicable indices)
		/////////////////////////////////////////////////////////////////////////////////////////
		
		if(debugLevel>2 && (i<=100 || pix_nn-i<=200) ){	//print the first and the last  entries to check
			cout << "setting correlationLUT[" << lutindex << "] = " << i << "";
			cout << "   LUT: pix_x, pix_y = (" << pix_x[i] << ", " << pix_y[i] << ") --> ix, iy = (" << ix << ", " << iy << ")" << endl;
		}
	}//for
	//after all is done, set the zero element to zero 
	//to make sure a failure of lookup() doesn't result in an acutal value
	correlationLUT[0] = 0;		
	
	cout << "\tLUT created ";
	cout << "(info: LUT assignment failed in " << lutFailCount << " of " << lutSize << " cases)" << endl;
	
	
	//explicit output to test...
	if (debugLevel>1) {
		std::ostringstream osst;
		osst << "-------------------LUT begin---------------------------------" << endl;
			for (int j = 0; j<lutNy; j++){
			osst << " [";
			for (int i = 0; i<lutNx; i++) {
				osst << " " << correlationLUT[i+lutNx*j];
			}
			osst << "]" << endl;
		}
		osst << "------------------LUT end----------------------------------" << endl;		
		cout << osst.str() << endl;
	}//if	
}




/*
 *	Write initial log file
 */
void cGlobal::writeInitialLog(void){
	FILE *fp;
	
	
	// Start time
	char	timestr[1024];
	time_t	rawtime;
	tm		*timeinfo;
	time(&rawtime);
	timeinfo=localtime(&rawtime);
	strftime(timestr,80,"%c",timeinfo);
	
	
	
	// Logfile name
	printf("Writing initial log file: %s\n", logfile);

	fp = fopen (logfile,"w");
	fprintf(fp, "start time: %s\n",timestr);
	fprintf(fp, ">-------- Start of job --------<\n");
	fclose (fp);
	
	
	// Open a new frame file at the same time
	pthread_mutex_lock(&framefp_mutex);
	
	sprintf(framefile,"r%04u-frames.txt",getRunNumber());
	framefp = fopen (framefile,"w");
	fprintf(framefp, "# FrameNumber, UnixTime, EventName, npeaks\n");

	sprintf(cleanedfile,"r%04u-cleaned.txt",getRunNumber());
	hitfinder.cleanedfp = fopen (cleanedfile,"w");
	fprintf(hitfinder.cleanedfp, "# Filename, npeaks\n");
	
	sprintf(icefile,"r%04u-icehits.txt",getRunNumber());
	icefinder.cleanedfp = fopen (icefile,"w");
	fprintf(icefinder.cleanedfp, "# Filename, npeaks\n");

	sprintf(waterfile,"r%04u-waterhits.txt",getRunNumber());
	waterfinder.cleanedfp = fopen (waterfile,"w");
	fprintf(waterfinder.cleanedfp, "# Filename, npeaks\n");

	sprintf(backgroundfile,"r%04u-backgroundhits.txt",getRunNumber());
	backgroundfinder.cleanedfp = fopen (backgroundfile,"w");
	fprintf(backgroundfinder.cleanedfp, "# Filename, npeaks\n");

	pthread_mutex_unlock(&framefp_mutex);
	
}


/*
 *	Update log file
 */
void cGlobal::updateLogfile(void){
	FILE *fp;
	
	// Calculate hit rate
	float hitrate;
	hitrate = 100.*( nhits / (float) nprocessedframes);
	
	// Elapsed processing time
	double	dtime;
	int		hrs, mins, secs; 
	time(&tend);
	dtime = difftime(tend,tstart);
	hrs = (int) floor(dtime / 3600);
	mins = (int) floor((dtime-3600*hrs)/60);
	secs = (int) floor(dtime-3600*hrs-60*mins);
	
	// Average data rate
	float	fps;
	fps = nprocessedframes / dtime;
	
	
	// Update logfile
	printf("Updating log file: %s\n", logfile);
	fp = fopen (logfile,"a");
	fprintf(fp, "nFrames: %i,  nHits: %i (%2.2f%%), wallTime: %ihr %imin %isec (%2.1f fps)\n", 
		(int)nprocessedframes, (int)nhits, hitrate, hrs, mins, secs, fps);
	fclose (fp);
	
	
	// Flush frame file buffer
	pthread_mutex_lock(&framefp_mutex);
	fclose(framefp);
	framefp = fopen (framefile,"a");
	fclose(hitfinder.cleanedfp);
	hitfinder.cleanedfp = fopen (cleanedfile,"a");
	fclose(icefinder.cleanedfp);
	icefinder.cleanedfp = fopen (icefile,"a");
	fclose(waterfinder.cleanedfp);
	waterfinder.cleanedfp = fopen (waterfile,"a");
	fclose(backgroundfinder.cleanedfp);
	backgroundfinder.cleanedfp = fopen (backgroundfile,"a");
	pthread_mutex_unlock(&framefp_mutex);
	
}

/*
 *	Write final log file
 */
void cGlobal::writeFinalLog(void){

	
	FILE *fp;
	
	// Logfile name
	printf("Writing final log file: %s\n", logfile);
	fp = fopen (logfile,"a");

	
	// Calculate hit rate
	float hitrate;
	hitrate = 100.*( nhits / (float) nprocessedframes);
	

	// End time
	char	timestr[1024];
	time_t	rawtime;
	tm		*timeinfo;
	time(&rawtime);
	timeinfo=localtime(&rawtime);
	strftime(timestr,80,"%c",timeinfo);
	
	
	// Elapsed processing time
	double	dtime;
	int		hrs, mins, secs; 
	time(&tend);
	dtime = difftime(tend,tstart);
	hrs = (int) floor(dtime / 3600);
	mins = (int) floor((dtime-3600*hrs)/60);
	secs = (int) floor(dtime-3600*hrs-60*mins);
	

	// Average data rate
	float	fps;
	fps = nprocessedframes / dtime;
				 
				 
	
	// Save log file
	fprintf(fp, ">-------- End of job --------<\n");
	fprintf(fp, "End time: %s\n",timestr);
	fprintf(fp, "Elapsed time: %ihr %imin %isec\n",hrs,mins,secs);
	fprintf(fp, "Frames processed: %i\n",(int)nprocessedframes);
	fprintf(fp, "nFrames in powder pattern: %i\n",(int)npowder);
	fprintf(fp, "nFrames in water powder pattern: %i\n",(int)nwater);
	fprintf(fp, "nFrames in ice powder pattern: %i\n",(int)nice);
	fprintf(fp, "Number of hits: %i\n",(int)nhits);
	fprintf(fp, "Average hit rate: %2.2f %%\n",hitrate);
	fprintf(fp, "Average data rate: %2.2f fps\n",fps);

	fclose (fp);

	
	// Flush frame file buffer
	pthread_mutex_lock(&framefp_mutex);
	fclose(framefp);
	fclose(hitfinder.cleanedfp);
	fclose(icefinder.cleanedfp);
	fclose(waterfinder.cleanedfp);
	fclose(backgroundfinder.cleanedfp);
	pthread_mutex_unlock(&framefp_mutex);
	
	
}

/*
 *	Read in peaksearch mask for ice finder
 */
void cGlobal::readIcemask(char *filename){
	
	printf("Reading peak search mask:\n");
	printf("\t%s\n",filename);
	
	
	// Create memory space and default to searching for peaks everywhere
	icefinder.peakmask = (int16_t*) calloc(pix_nn, sizeof(int16_t));
	for(long i=0;i<pix_nn;i++)
		icefinder.peakmask[i] = 1;
	
	
	// Check whether file exists!
	FILE* fp = fopen(filename, "r");
	if (fp) 	// file exists
		fclose(fp);
	else {		// file doesn't exist
		printf("\tPeak search mask does not exist: %s\n",filename);
		printf("\tDefaulting to uniform search mask\n");
		return;
	}
	
	
	// Read peakmask data from file
	cData2d		temp2d;
	temp2d.readHDF5(filename);
	
	
	// Correct geometry?
	if(temp2d.nx != pix_nx || temp2d.ny != pix_ny) {
		printf("\tGeometry mismatch: %ix%x != %ix%i\n",(int)temp2d.nx, (int)temp2d.ny, (int)pix_nx, (int)pix_ny);
		printf("\tDefaulting to uniform peak search mask\n");
		return;
	} 
	
	
	// Copy into peakmask array
	for(long i=0;i<pix_nn;i++)
		icefinder.peakmask[i] = (int16_t) temp2d.data[i];
}


/*
 *	Read in peaksearch mask for water finder
 */
void cGlobal::readWatermask(char *filename){
	
	printf("Reading peak search mask:\n");
	printf("\t%s\n",filename);
	
	
	// Create memory space and default to searching for peaks everywhere
	waterfinder.peakmask = (int16_t*) calloc(pix_nn, sizeof(int16_t));
	for(long i=0;i<pix_nn;i++)
		waterfinder.peakmask[i] = 1;
	
	
	// Check whether file exists!
	FILE* fp = fopen(filename, "r");
	if (fp) 	// file exists
		fclose(fp);
	else {		// file doesn't exist
		printf("\tPeak search mask does not exist: %s\n",filename);
		printf("\tDefaulting to uniform search mask\n");
		return;
	}
	
	
	// Read peakmask data from file
	cData2d		temp2d;
	temp2d.readHDF5(filename);
	
	
	// Correct geometry?
	if(temp2d.nx != pix_nx || temp2d.ny != pix_ny) {
		printf("\tGeometry mismatch: %ix%x != %ix%i\n",(int)temp2d.nx, (int)temp2d.ny, (int)pix_nx, (int)pix_ny);
		printf("\tDefaulting to uniform peak search mask\n");
		return;
	} 
	
	
	// Copy into peakmask array
	for(long i=0;i<pix_nn;i++)
		waterfinder.peakmask[i] = (int16_t) temp2d.data[i];
}


/*
 *	Read in peaksearch mask for background finder
 */
void cGlobal::readBackgroundmask(char *filename){
	
	printf("Reading peak search mask:\n");
	printf("\t%s\n",filename);
	
	
	// Create memory space and default to searching for peaks everywhere
	backgroundfinder.peakmask = (int16_t*) calloc(pix_nn, sizeof(int16_t));
	for(long i=0;i<pix_nn;i++)
		backgroundfinder.peakmask[i] = 1;
	
	
	// Check whether file exists!
	FILE* fp = fopen(filename, "r");
	if (fp) 	// file exists
		fclose(fp);
	else {		// file doesn't exist
		printf("\tPeak search mask does not exist: %s\n",filename);
		printf("\tDefaulting to uniform search mask\n");
		return;
	}
	
	
	// Read peakmask data from file
	cData2d		temp2d;
	temp2d.readHDF5(filename);
	
	
	// Correct geometry?
	if(temp2d.nx != pix_nx || temp2d.ny != pix_ny) {
		printf("\tGeometry mismatch: %ix%x != %ix%i\n",(int)temp2d.nx, (int)temp2d.ny, (int)pix_nx, (int)pix_ny);
		printf("\tDefaulting to uniform peak search mask\n");
		return;
	} 
	
	
	// Copy into peakmask array
	for(long i=0;i<pix_nn;i++)
		backgroundfinder.peakmask[i] = (int16_t) temp2d.data[i];
}
