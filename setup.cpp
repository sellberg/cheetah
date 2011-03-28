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
#include "release/pdsdata/cspad/ElementHeader.hh"
#include "release/pdsdata/cspad/ElementIterator.hh"
#include "cspad-gjw/CspadTemp.hh"
#include "cspad-gjw/CspadCorrector.hh"
#include "cspad-gjw/CspadGeometry.hh"

#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>
#include <math.h>
#include <hdf5.h>
#include <fenv.h>
#include <stdlib.h>

#include "setup.h"
#include "worker.h"
#include "data2d.h"


/*
 *	Default settings/configuration
 */
void cGlobal::defaultConfiguration(void) {

	// ini file to use
	strcpy(configFile, "cheetah.ini");

	// Geometry
	strcpy(geometryFile, "geometry/cspad_pixelmap.h5");
	pixelSize = 110e-6;
	
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
	strcpy(backgroundfinder.peaksearchFile, "backgroundmask.h5");

	// Powder pattern generation
	powdersum = 1;
	powderthresh = 0;
	
	// Saving options
	savehits = 0;
	saveRaw = 0;
	hdf5dump = 0;
	saveInterval = 500;
	
	// Verbosity
	debugLevel = 2;
	
	// Default to only a few threads
	nThreads = 16;
	
	// Log files
	strcpy(logfile, "log.txt");
	strcpy(framefile, "frames.txt");
	strcpy(cleanedfile, "cleaned.txt");
	
	
}



/*
 *	Setup stuff to do with thread management, settings, etc.
 */
void cGlobal::setup() {
	/*
	 *	Set up arrays for remembering powder data, background, etc.
	 */
	hotpixelmask = (float*) calloc(pix_nn, sizeof(float));
	selfdark = (float*) calloc(pix_nn, sizeof(float));
	powderRaw = (int64_t*) calloc(pix_nn, sizeof(int64_t));
	powderAssembled = (int64_t*) calloc(image_nn, sizeof(int64_t));
	iceRaw = (int64_t*) calloc(pix_nn, sizeof(int64_t));
	iceAssembled = (int64_t*) calloc(image_nn, sizeof(int64_t));
	waterRaw = (int64_t*) calloc(pix_nn, sizeof(int64_t));
	waterAssembled = (int64_t*) calloc(image_nn, sizeof(int64_t));

	for(long i=0; i<pix_nn; i++) {
		hotpixelmask[i] = 0;
		selfdark[i] = 0;
		powderRaw[i] = 0;
	}
	for(long i=0; i<image_nn; i++) {
		powderAssembled[i] = 0;
	}	

	
	/*
	 *	Set up thread management
	 */
	nActiveThreads = 0;
	threadCounter = 0;
	pthread_mutex_init(&nActiveThreads_mutex, NULL);
	pthread_mutex_init(&hotpixel_mutex, NULL);
	pthread_mutex_init(&selfdark_mutex, NULL);
	pthread_mutex_init(&powdersum1_mutex, NULL);
	pthread_mutex_init(&powdersum2_mutex, NULL);
	pthread_mutex_init(&nhits_mutex, NULL);
	pthread_mutex_init(&framefp_mutex, NULL);
	pthread_mutex_init(&icesumraw_mutex, NULL);
	pthread_mutex_init(&icesumassembled_mutex, NULL);
	pthread_mutex_init(&watersumraw_mutex, NULL);
	pthread_mutex_init(&watersumassembled_mutex, NULL);

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
		savehits = 0;
		hdf5dump = 0;
		saveRaw = 0;
		useAutoHotpixel = 0;
		startFrames = 0;
		powderthresh = 0;
	}
	
	/*
	 *	Other stuff
	 */
	npowder = 0;
	nprocessedframes = 0;
	nhits = 0;
	lastclock = clock()-10;
	gettimeofday(&lasttime, NULL);
	datarate = 1;
	detectorZ = 0;
	runNumber = getRunNumber();
	time(&tstart);
	avgGMD = 0;

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
	printf("Parsing input configuration file:\n",filename);
	printf("\t%s\n",filename);
	
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
	// Processing options
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
		savehits = atoi(value);
	}
	else if (!strcmp(tag, "powdersum")) {
		powdersum = atoi(value);
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
	else if (!strcmp(tag, "icepeakmask")) {
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

	else if (!strcmp(tag, "waterpeakmask")) {
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


	else if (!strcmp(tag, "backgroundpeakmask")) {
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
		printf("%ix%i != %ix%i\n", 8*ROWS, 8*COLS, detector_x.nx, detector_x.ny);
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
	
	
	// Divide array (in m) by pixel size to get pixel location indicies (ijk)
	for(long i=0;i<nn;i++){
		pix_x[i] /= pix_dx;
		pix_y[i] /= pix_dx;
		pix_z[i] /= pix_dx;
	}
	
	
	// Find bounds of image array
	float	xmax = -1e9;
	float	xmin =  1e9;
	float	ymax = -1e9;
	float	ymin =  1e9;
	for(long i=0;i<nn;i++){
		if (pix_x[i] > xmax) xmax = pix_x[i];
		if (pix_x[i] < xmin) xmin = pix_x[i];
		if (pix_y[i] > ymax) ymax = pix_y[i];
		if (pix_y[i] < ymin) ymin = pix_y[i];
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
	printf("\tImage output array will be %i x %i\n",image_nx,image_nx);
	
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
		printf("\tGeometry mismatch: %ix%x != %ix%i\n",temp2d.nx, temp2d.ny, pix_nx, pix_ny);
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
		printf("\tGeometry mismatch: %ix%x != %ix%i\n",temp2d.nx, temp2d.ny, pix_nx, pix_ny);
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
		printf("\tGeometry mismatch: %ix%x != %ix%i\n",temp2d.nx, temp2d.ny, pix_nx, pix_ny);
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
		printf("\tGeometry mismatch: %ix%x != %ix%i\n",temp2d.nx, temp2d.ny, pix_nx, pix_ny);
		printf("\tDefaulting to uniform peak search mask\n");
		return;
	} 
	
	
	// Copy back into array
	for(long i=0;i<pix_nn;i++)
		badpixelmask[i] = (int16_t) temp2d.data[i];
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
	printf("Writing log file: %s\n", logfile);

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
	cleanedfp = fopen (cleanedfile,"w");
	fprintf(cleanedfp, "# Filename, npeaks\n");
	
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
	printf("Writing log file: %s\n", logfile);
	fp = fopen (logfile,"a");
	fprintf(fp, "nFrames: %i,  nHits: %i (%2.2f%%), wallTime: %ihr %imin %isec (%2.1f fps)\n", nprocessedframes, nhits, hitrate, hrs, mins, secs, fps);
	fclose (fp);
	
	
	// Flush frame file buffer
	pthread_mutex_lock(&framefp_mutex);
	fclose(framefp);
	framefp = fopen (framefile,"a");
	fclose(cleanedfp);
	cleanedfp = fopen (cleanedfile,"a");
	pthread_mutex_unlock(&framefp_mutex);
	
}

/*
 *	Write final log file
 */
void cGlobal::writeFinalLog(void){

	
	FILE *fp;
	
	// Logfile name
	printf("Writing log file: %s\n", logfile);
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
	fprintf(fp, "Frames processed: %i\n",nprocessedframes);
	fprintf(fp, "nFrames in powder pattern: %i\n",npowder);
	fprintf(fp, "Number of hits: %i\n",nhits);
	fprintf(fp, "Average hit rate: %2.2f %%\n",hitrate);
	fprintf(fp, "Average data rate: %2.2f fps\n",fps);

	fclose (fp);

	
	// Flush frame file buffer
	pthread_mutex_lock(&framefp_mutex);
	fclose(framefp);
	fclose(cleanedfp);
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
		printf("\tGeometry mismatch: %ix%x != %ix%i\n",temp2d.nx, temp2d.ny, pix_nx, pix_ny);
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
		printf("\tGeometry mismatch: %ix%x != %ix%i\n",temp2d.nx, temp2d.ny, pix_nx, pix_ny);
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
		printf("\tGeometry mismatch: %ix%x != %ix%i\n",temp2d.nx, temp2d.ny, pix_nx, pix_ny);
		printf("\tDefaulting to uniform peak search mask\n");
		return;
	} 
	
	
	// Copy into peakmask array
	for(long i=0;i<pix_nn;i++)
		backgroundfinder.peakmask[i] = (int16_t) temp2d.data[i];
}
