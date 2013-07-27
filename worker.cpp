/*
 *  worker.cpp
 *  cheetah
 *
 *  Created by Anton Barty on 6/2/11.
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
//#include "myana/XtcRun.hh"
#include "release/pdsdata/cspad/ElementHeader.hh"
#include "release/pdsdata/cspad/ElementIterator.hh"
//#include "cspad-gjw/CspadTemp.hh"
//#include "cspad-gjw/CspadGeometry.hh"

#include <string.h>
#include <cmath>
#include <hdf5.h>
#include <fenv.h>
#include <stdlib.h>
#include <iostream>
#include <string>
using std::cout;
using std::endl;
using std::string;

#include "worker.h"
#include "hitfinder.h"
#include "commonmode.h"
#include "background.h"
#include "correlation.h"
#include "arrayclasses.h"
#include "arraydataIO.h"
#include "util.h"


/*
 *	Worker thread function for processing each cspad data frame
 */
void *worker(void *threadarg) {

	/*
	 *	Turn threadarg into a more useful form
	 */
	cGlobal			*global;
	tThreadInfo		*threadInfo;
	cHit 			hit;

	threadInfo = (tThreadInfo*) threadarg; 	
	global = threadInfo->pGlobal;
	
	
	
	/*
	 *	Assemble data from all four quadrants into one large array (rawdata format)
	 */
	threadInfo->corrected_data = (float*) calloc(RAW_DATA_LENGTH,sizeof(float));
	for(int quadrant=0; quadrant<4; quadrant++) {
		long	i,j,ii;
		for(long k=0; k<2*ROWS*8*COLS; k++) {
			i = k % (2*ROWS) + quadrant*(2*ROWS);
			j = k / (2*ROWS);
			ii  = i+(8*ROWS)*j;
			threadInfo->corrected_data[ii] = (float) threadInfo->quad_data[quadrant][k];
		}
	}
	
	
	/*
	 *	Initialize pointers to analysis arrays to NULL
	 */
	threadInfo->angularAvg = NULL;
	threadInfo->angularAvgQ = NULL;
	threadInfo->correlation = NULL;
	threadInfo->image = NULL;
	threadInfo->theta = NULL;
	threadInfo->pix_qx = NULL;
	threadInfo->pix_qy = NULL;
	
	
	/*
     *  Initialize fail check for averaging
     */
	int fail = 0;
	
	
	/*
	 *	Create a unique name for this event
	 */
	nameEvent(threadInfo, global);
		

	/*
	 *	Subtract darkcal image (static electronic offsets)
	 */
	if(global->useDarkcalSubtraction) {
		subtractDarkcal(threadInfo, global);
	}

	
	/*
	 *	Subtract common mode offsets (electronic offsets)
	 */
	if(global->cmModule) {
		cmModuleSubtract(threadInfo, global);
	}
	else if(global->cmSubModule) {
		cmSubModuleSubtract(threadInfo, global);
	}
	
	
	/*
	 *	Apply gain correction
	 */
	if(global->useGaincal) {
		applyGainCorrection(threadInfo, global);
	}
	
	
	/*
	 *	Apply bad pixel map
	 */
	if(global->useBadPixelMask) {
		applyBadPixelMask(threadInfo, global);
	} 
	
	
	/*
	 *	Subtract running photon background
	 */
	if(global->useSubtractPersistentBackground) {
		subtractPersistentBackground(threadInfo, global);
	}
	

	/*
	 *	Identify and remove hot pixels
	 */
	if(global->useAutoHotpixel){
		killHotpixels(threadInfo, global);
	}
	
	
	/*
	 *	Hitfinding
	 */
	hit.standard = 0;
	if(global->hitfinder.use){
		hit.standard = hitfinder(threadInfo, global, &(global->hitfinder));
		hit.standardPeaks = threadInfo->nPeaks;
	}
	
	/*
	 *	Hitfinding - Water
	 */
	hit.water = 0;
	if(global->waterfinder.use){
		hit.water = hitfinder(threadInfo, global, &(global->waterfinder));
		hit.waterPeaks = threadInfo->nPeaks;
	}

	/*
	 *	Hitfinding - Ice
	 */
	hit.ice = 0;
	if(global->icefinder.use){
		hit.ice = hitfinder(threadInfo, global, &(global->icefinder));
		hit.icePeaks = threadInfo->nPeaks;
	}
	 
	/*
	 *	Hitfinding - Background
	 */
	hit.background = 1;
	if(global->backgroundfinder.use){
		hit.background = hitfinder(threadInfo, global, &(global->backgroundfinder));
		hit.backgroundPeaks = threadInfo->nPeaks;
	}
	
	/*
	 *	Hitfinding - List
	 */
	hit.list = 0;
	if(global->listfinder.use){
		hit.list = hitfinder(threadInfo, global, &(global->listfinder));
		hit.listPeaks = threadInfo->nPeaks;
	}
	
	/*
	 *	Hitfinding - Update central hit counter
	 */
	if (hit.standard || hit.water || hit.ice) {
		pthread_mutex_lock(&global->nhits_mutex);
		global->nhits++;
		pthread_mutex_unlock(&global->nhits_mutex);
	}
	
	/*
	 *	Update the running background
	 */
	if(global->useSubtractPersistentBackground) {
		if (global->backgroundfinder.use){	
			updatePersistentBackground(threadInfo, global, hit.background);
		} else {
			updatePersistentBackground(threadInfo, global, hit.standard);
		}
	}

	/*
	 *	Are we still in 'frame digesting' mode?
	 */
	if(threadInfo->threadNum < global->startFrames) {
		printf("r%04u:%i (%3.1f Hz): Digesting initial frames\n", (int)threadInfo->runNumber, (int)threadInfo->threadNum, global->datarate);
		goto cleanup;
        //ATTENTION! goto should not be used at all ( see http://www.cplusplus.com/forum/general/29190/ )
        //should be replaced by a while loop with a break statement... by JF (2011/04/25)
	}
	
	
	/*
	 *	Calculate center correction
	 */
	if (global->calculateCenterCorrectionHit && (hit.standard || hit.water || hit.ice)) {
		calculateCenterCorrection(threadInfo, global, threadInfo->corrected_data, 1);
		if (global->useCenterCorrection) {
			// DOES NOT WORK PROPERLY!!!
			// SCREWS UP THE CENTER CALCULATION (THREADED & NONTHREADED)
			// NEED TO MAKE SURE THREADS DO NOT READ UPDATED VARIABLES VALUES FROM OTHER THREADS WHILE UPDATING ARRAYS/ASSEMBLY/POWDERADD/POWDERSSAVE
			//updatePixelArrays(threadInfo, global);
			//updateImageArrays(global, &hit); // IS NOT THREADSAFE IF POWDERS ARE SAVED
			cout << "\tuseCenterCorrection is currently disabled for hits" << endl;
		}
	}
	
	
	/*
	 *	Calculate scattering angle (theta)
	 */
	if ((global->usePolarizationCorrection || global->useSolidAngleCorrection || global->useCorrelation) && (global->hdf5dump
																			  || global->generateDarkcal
																			  || hit.standard
																		      || hit.water
																			  || hit.ice
																			  || !hit.background )) {
		calculateScatteringAngle(threadInfo, global);
	}
	
	
	/*
	 *	Calculate polarization correction
	 */
	if (global->usePolarizationCorrection && (global->hdf5dump || global->generateDarkcal
															   || hit.standard
															   || hit.water
															   || hit.ice
															   || !hit.background )) {
		calculatePolarizationCorrection(threadInfo, global);
	}
	
	
	/*
	 *	Calculate solid angle correction
	 */
	threadInfo->solidAngle = global->pixelSize*global->pixelSize/(threadInfo->detectorPosition/1000*threadInfo->detectorPosition/1000);
	if (global->useSolidAngleCorrection && (global->hdf5dump || global->generateDarkcal
															 || hit.standard
															 || hit.water
															 || hit.ice
															 || !hit.background )) {
		calculateSolidAngleCorrection(threadInfo, global);
	}
	
	
	/*
	 *	Apply attenuation correction
	 */
	if (global->useAttenuationCorrection > 0) {
		applyAttenuationCorrection(threadInfo, global);
	}
	
	
	/*
     *  Calculate intensity average
     */
	calculateIntensityAvg(threadInfo, global);
	if (global->useIntensityStatistics) {
		pthread_mutex_lock(&global->intensities_mutex);
		if (global->nIntensities >= global->intensityCapacity) global->expandIntensityCapacity();
		if (global->hdf5dump || global->generateDarkcal
							 || hit.standard
							 || hit.water
							 || hit.ice
							 || !hit.background) {
			global->hits[global->nIntensities] = true;
		} else {
			global->hits[global->nIntensities] = false;
		}
		// OBS: - the order of global->intensities will not be exactly the same as the event order, but close enough (within global->nThreads of events)
		//      - if useSolidAngleCorrection and/or usePolarizationCorrection is enabled, they will only be calculated for hits/saved events
		global->intensities[global->nIntensities++] = threadInfo->intensityAvg; 
		if (threadInfo->intensityAvg > global->Imax) global->Imax = threadInfo->intensityAvg;
		if (threadInfo->intensityAvg < global->Imin) global->Imin = threadInfo->intensityAvg;
		pthread_mutex_unlock(&global->intensities_mutex);
	}
	
	
	/*
     *  Calculate angular averages
     */
	if (global->hitAngularAvg && (global->hdf5dump || global->generateDarkcal
												   || (hit.standard && global->hitfinder.savehits) 
												   || (hit.water && global->waterfinder.savehits) 
												   || (hit.ice && global->icefinder.savehits) 
												   || (!hit.background && global->backgroundfinder.savehits) )) {
		calculateAngularAvg(threadInfo, global);
		fail = makeQcalibration(threadInfo, global);
		if (!fail) saveAngularAvg(threadInfo, global);
		else cout << "Failed to calibrate Q for " << threadInfo->eventname << ", angular average NOT saved." << endl;
	}
	
	
	/*
     *  Perform cross-correlation analysis
     */
	if (global->useCorrelation && (global->hdf5dump || global->generateDarkcal
													|| hit.standard
													|| hit.water
													|| hit.ice
													|| !hit.background )) {
		fail = calculatePixelMaps(threadInfo, global);
		if (!fail) correlate(threadInfo, global, &hit);
		else cout << "Failed to make Q-calibrated pixel maps for " << threadInfo->eventname << ", correlation NOT saved." << endl;
	}
	
	
	/*
     *  Calculate single-pixel statistics
     */
	if (global->usePixelStatistics && (global->hdf5dump || global->generateDarkcal
														|| (hit.standard && global->hitfinder.savehits) 
														|| (hit.water && global->waterfinder.savehits) 
														|| (hit.ice && global->icefinder.savehits) 
														|| (!hit.background && global->backgroundfinder.savehits) )) {
		// save pixel intensities of hit without Q-scale (1D array)
		//if (!global->useCorrelation) fail = calculatePixelMaps(threadInfo, global);
		//if (!fail) savePixelIntensities(threadInfo, global);
		//else cout << "Failed to calibrate Q for " << threadInfo->eventname << ", pixel intensities NOT saved." << endl;
		savePixelIntensities(threadInfo, global);
	}
	
	
	/*
	 *	Assemble quadrants into a 'realistic' 2D image
	 */
	if (global->hdf5dump || (hit.standard && (global->hitfinder.savehits || global->powdersum)) 
						 || (hit.water && (global->waterfinder.savehits || global->powdersum)) 
						 || (hit.ice && (global->icefinder.savehits || global->powdersum)) 
						 || (!hit.background && (global->backgroundfinder.savehits || global->powdersum))
        				 || (global->listfinder.use && (global->listfinder.savehits || global->powdersum)) ) {
		assemble2Dimage(threadInfo, global);
	}
	
	
	/*
	 *	Add to powder if it's a hit or if we wish to generateDarkcal(member data of global)
	 */
	if (!fail) addToPowder(threadInfo, global, &hit);
	
	
	/*
	 *	Add to correlation sum if it's a hit and if correlation sum is activated
	 */
	if (!fail) addToCorrelation(threadInfo, global, &hit);
	
	
	/*
	 *	If this is a hit, write out to our favourite HDF5 format
	 */
	if(global->hdf5dump) 
		writeHDF5(threadInfo, global, threadInfo->eventname);
	else {
		if(hit.standard && global->hitfinder.savehits) {
			threadInfo->nPeaks = hit.standardPeaks;
			writeHDF5(threadInfo, global, threadInfo->eventname);
		}
		
		if(hit.water && global->waterfinder.savehits) {
			threadInfo->nPeaks = hit.waterPeaks;
			writeHDF5(threadInfo, global, threadInfo->eventname);
		}
		
		if(hit.ice && global->icefinder.savehits) {
			threadInfo->nPeaks = hit.icePeaks;
			writeHDF5(threadInfo, global, threadInfo->eventname);			
		}
		
		if(!hit.background && global->backgroundfinder.savehits) {
			threadInfo->nPeaks = hit.backgroundPeaks;
			writeHDF5(threadInfo, global, threadInfo->eventname);			
		}
		
//		char eventname[1024];
//		if(hit.water && global->waterfinder.savehits) {
//			sprintf(eventname,"%s_waterhit",threadInfo->eventname);
//			writeHDF5(threadInfo, global, eventname);
//		}
//		if(hit.ice && global->icefinder.savehits) {
//			sprintf(eventname,"%s_icehit",threadInfo->eventname);
//			writeHDF5(threadInfo, global, eventname);
//		}
//		if(!hit.background && global->backgroundfinder.savehits) {
//			sprintf(eventname,"%s_background",threadInfo->eventname);
//			writeHDF5(threadInfo, global, eventname);
//		}
		
	}
	
	
	/*
	 *	Write out diagnostics to screen
	 */	
	if (global->useAutoHotpixel) {
		printf("r%04u:%i (%3.1f Hz): Processed (iavg=%4.2f, hot=%i", (int)threadInfo->runNumber, (int)threadInfo->threadNum, global->datarate, threadInfo->intensityAvg, threadInfo->nHot);
		if (global->hitfinder.use) {
			printf("; hit=%i, nat/npeaks=%i", hit.standard, hit.standardPeaks);
		}
		if (global->icefinder.use) {
			printf("; ice=%i, nat/npeaks=%i", hit.ice, hit.icePeaks);
		}
		if (global->waterfinder.use) {
			printf("; water=%i, nat/npeaks=%i", hit.water, hit.waterPeaks);
		}
		if (global->backgroundfinder.use) {
			printf("; background=%i, nat/npeaks=%i", (hit.background) ? 0 : 1, hit.backgroundPeaks);
		}
		printf(")\n");
	} else if (global->hitfinder.use || global->icefinder.use || global->waterfinder.use || global->backgroundfinder.use || global->listfinder.use) {
		printf("r%04u:%i (%3.1f Hz): Processed (iavg=%4.2f", (int)threadInfo->runNumber, (int)threadInfo->threadNum, global->datarate, threadInfo->intensityAvg);
		if (global->hitfinder.use) {
			printf("; hit=%i, nat/npeaks=%i", hit.standard, hit.standardPeaks);
		}
		if (global->icefinder.use) {
			printf("; ice=%i, nat/npeaks=%i", hit.ice, hit.icePeaks);
		}
		if (global->waterfinder.use) {
			printf("; water=%i, nat/npeaks=%i", hit.water, hit.waterPeaks);
		}
		if (global->backgroundfinder.use) {
			printf("; background=%i, nat/npeaks=%i", (hit.background) ? 0 : 1, hit.backgroundPeaks);
		}
		if (global->listfinder.use) {
			printf("; nat/npeaks=%i", hit.listPeaks);
		}
		printf(")\n");
	}
	
	/*
	 *	Write out information on each frame to a log file
	 */
	pthread_mutex_lock(&global->framefp_mutex);
	fprintf(global->framefp, "%i, %i, %s, %f", (int)threadInfo->threadNum, threadInfo->seconds, threadInfo->eventname, threadInfo->intensityAvg);
	if (global->hitfinder.use) {
		fprintf(global->framefp, ", %i", hit.standardPeaks);
		if (hit.standard) fprintf(global->hitfinder.cleanedfp, "r%04u/%s, %f, %i\n", (int)threadInfo->runNumber, threadInfo->eventname, threadInfo->intensityAvg, hit.standardPeaks);
	}
	if (global->icefinder.use) {
		fprintf(global->framefp, ", %i", hit.icePeaks);
		if (hit.ice) fprintf(global->icefinder.cleanedfp, "r%04u/%s, %f, %i\n", (int)threadInfo->runNumber, threadInfo->eventname, threadInfo->intensityAvg, hit.icePeaks);
	}
	if (global->waterfinder.use) {
		fprintf(global->framefp, ", %i", hit.waterPeaks);
		if (hit.water) fprintf(global->waterfinder.cleanedfp, "r%04u/%s, %f, %i\n", (int)threadInfo->runNumber, threadInfo->eventname, threadInfo->intensityAvg, hit.waterPeaks);
	}
	if (global->backgroundfinder.use) {
		fprintf(global->framefp, ", %i", hit.backgroundPeaks);
		if (hit.background) fprintf(global->backgroundfinder.cleanedfp, "r%04u/%s, %f, %i\n", (int)threadInfo->runNumber, threadInfo->eventname, threadInfo->intensityAvg, hit.backgroundPeaks);
	}
	if (global->listfinder.use) {
		fprintf(global->framefp, ", %i", hit.listPeaks);
		fprintf(global->hitfinder.cleanedfp, "r%04u/%s, %f, %i\n", (int)threadInfo->runNumber, threadInfo->eventname, threadInfo->intensityAvg, hit.listPeaks);
	}
	fprintf(global->framefp, "\n");
	pthread_mutex_unlock(&global->framefp_mutex);
	
	
	/*
	 *	Cleanup and exit
	 */
	cleanup:
	// Decrement thread pool counter by one
	pthread_mutex_lock(&global->nActiveThreads_mutex);
	global->nActiveThreads -= 1;
	pthread_mutex_unlock(&global->nActiveThreads_mutex);
	
	// Free memory
	for(int quadrant=0; quadrant<4; quadrant++) 
		free(threadInfo->quad_data[quadrant]);	
	free(threadInfo->corrected_data);
	free(threadInfo->image);
	free(threadInfo->angularAvg);
	free(threadInfo->angularAvgQ);
	free(threadInfo->correlation);
	delete[] threadInfo->theta;
	delete[] threadInfo->pix_qx;
	delete[] threadInfo->pix_qy;
	free(threadInfo);
	
	// Exit thread
	pthread_exit(NULL);
}


/*
 *	Subtract pre-loaded darkcal file
 */
void subtractDarkcal(tThreadInfo *threadInfo, cGlobal *global){
	
	
	// Do darkcal subtraction
	// Watch out for integer wraparound!
	int32_t diff;
	for(long i=0;i<global->pix_nn;i++) {
		diff = (int32_t) threadInfo->corrected_data[i] - (int32_t) global->darkcal[i];	
		if(diff < -32766) diff = -32767;
		if(diff > 32766) diff = 32767;
		threadInfo->corrected_data[i] = (int16_t) diff;
	}
	
}

/*
 *	Apply gain correction
 *	Assumes the gaincal array is appropriately 'prepared' when loaded so that all we do is a multiplication.
 *	All that checking for division by zero (and inverting when required) needs only be done once, right? 
 */
void applyGainCorrection(tThreadInfo *threadInfo, cGlobal *global){
	
	for(long i=0;i<global->pix_nn;i++) 
		threadInfo->corrected_data[i] *= global->gaincal[i];
	
}


/*
 *	Apply bad pixel mask
 *	Assumes that all we have to do here is a multiplication.
 */
void applyBadPixelMask(tThreadInfo *threadInfo, cGlobal *global){
	
	for(long i=0;i<global->pix_nn;i++) 
		threadInfo->corrected_data[i] *= global->badpixelmask[i];
	
}


/*
 *	Identify and kill hot pixels
 */
void killHotpixels(tThreadInfo *threadInfo, cGlobal *global){
	
	int	nhot = 0;

	pthread_mutex_lock(&global->hotpixel_mutex);
	for(long i=0;i<global->pix_nn;i++){
		global->hotpixelmask[i] = ( (global->hotpixMemory-1)*global->hotpixelmask[i] + ((threadInfo->corrected_data[i]>global->hotpixADC)?(1.0):(0.0))) / global->hotpixMemory;

		if(global->hotpixelmask[i] > global->hotpixFreq) {
			threadInfo->corrected_data[i] = 0;
			nhot++;
		}
	}
	pthread_mutex_unlock(&global->hotpixel_mutex);
	threadInfo->nHot = nhot;
}


/*
 *	Calculate array of scattering angle (theta)
 */
void calculateScatteringAngle(tThreadInfo *threadInfo, cGlobal *global) {
	
	// allocate scattering angle array
	if (threadInfo->theta) {
		delete[] threadInfo->theta;
	}
	threadInfo->theta = new double[global->pix_nn];
	
	// calculate scattering angle (theta) for each pixel
	for (int i = 0; i < global->pix_nn; i++) {
		threadInfo->theta[i] = atan(global->pixelSize*global->pix_r[i]*1000/threadInfo->detectorPosition);
	}
}


/*
 *	Calculate Q-calibrated pixel maps
 */
int calculatePixelMaps(tThreadInfo *threadInfo, cGlobal *global) {
	
	// sanity check
	if (threadInfo->detectorPosition > 60 && threadInfo->detectorPosition < 600 && threadInfo->wavelengthA == threadInfo->wavelengthA) {
		// allocate arrays for Q-calibrated pixel maps (needs to be float precision for current CrossCorrelator)
		if (threadInfo->pix_qx) {
			delete[] threadInfo->pix_qx;
		}
		if (threadInfo->pix_qy) {
			delete[] threadInfo->pix_qy;
		}
		threadInfo->pix_qx = new float[global->pix_nn];
		threadInfo->pix_qy = new float[global->pix_nn];
		
		if (global->correlationQScale == 2) { // |q| [Å-1]
			// calculate the magnitude of the q-vector [Å-1] and create the qx/qy pixel maps from that
			for (int i = 0; i < global->pix_nn; i++) {
				double pix_qperp = 4*M_PI*sin(threadInfo->theta[i]/2)/threadInfo->wavelengthA;
				threadInfo->pix_qx[i] = (float) pix_qperp*sin(global->phi[i]); // phi=0 is defined in +Y direction
				threadInfo->pix_qy[i] = (float) pix_qperp*cos(global->phi[i]);
			}			
		} else if (global->correlationQScale == 3) { // |q_perp|
			// calculate the magnitude of the perpendicular q-component [Å-1] and create the qx/qy pixel maps from that
			for (int i = 0; i < global->pix_nn; i++) {
				double pix_qperp = 2*M_PI*sin(threadInfo->theta[i])/threadInfo->wavelengthA;
				threadInfo->pix_qx[i] = (float) pix_qperp*sin(global->phi[i]); // phi=0 is defined in +Y direction
				threadInfo->pix_qy[i] = (float) pix_qperp*cos(global->phi[i]);
			}			
		} else { // pixels
			// use global->pix_x/global->pix_y as pixel maps
			delete[] threadInfo->pix_qx;
			delete[] threadInfo->pix_qy;
			threadInfo->pix_qx = NULL;
			threadInfo->pix_qy = NULL;
		}
		
		return 0;
	} else return 1;
}


/*
 *	Calculate polarization correction
 */
void calculatePolarizationCorrection(tThreadInfo *threadInfo, cGlobal *global) {
	
	// calculate and apply polarization correction to corrected_data (from Hura et al JCP 2000)
	for (int i = 0; i < global->pix_nn; i++) {
		threadInfo->corrected_data[i] /= global->horizontalPolarization*(1 - sin(global->phi[i])*sin(global->phi[i])*sin(threadInfo->theta[i])*sin(threadInfo->theta[i])) + (1 - global->horizontalPolarization)*(1 - cos(global->phi[i])*cos(global->phi[i])*sin(threadInfo->theta[i])*sin(threadInfo->theta[i]));
	}
}


/*
 *	Calculate solid angle correction
 */
void calculateSolidAngleCorrection(tThreadInfo *threadInfo, cGlobal *global) {
	
	if (global->useSolidAngleCorrection == 2) {
		
		// Azimuthally symmetrical correction
		for (int i = 0; i < global->pix_nn; i++) {
			threadInfo->corrected_data[i] /= cos(threadInfo->theta[i])*cos(threadInfo->theta[i])*cos(threadInfo->theta[i]);
		}
		
	} else {
		
		// Rigorous correction from solid angle of a plane triangle
		for (int i = 0; i < global->pix_nn; i++) {
			
			// allocate local arrays
			double corner_coordinates[4][3]; // array of vector coordinates of pixel corners, first index starts from upper left corner and goes around clock-wise, second index determines X=0/Y=1/Z=2 coordinate
			double corner_distances[4]; // array of distances of pixel corners, index starts from upper left corner and goes around clock-wise
			double determinant;
			double denominator;
			double solid_angle[2]; // array of solid angles of the two plane triangles that form the pixel
			double total_solid_angle;
			
			// upper left corner
			corner_coordinates[0][0] = global->pix_x[i]*global->pixelSize + global->pixelSize/2;
			corner_coordinates[0][1] = global->pix_y[i]*global->pixelSize + global->pixelSize/2;
			// upper right corner
			corner_coordinates[1][0] = global->pix_x[i]*global->pixelSize - global->pixelSize/2;
			corner_coordinates[1][1] = global->pix_y[i]*global->pixelSize + global->pixelSize/2;
			// lower right corner
			corner_coordinates[2][0] = global->pix_x[i]*global->pixelSize - global->pixelSize/2;
			corner_coordinates[2][1] = global->pix_y[i]*global->pixelSize - global->pixelSize/2;
			// lower left corner
			corner_coordinates[3][0] = global->pix_x[i]*global->pixelSize + global->pixelSize/2;
			corner_coordinates[3][1] = global->pix_y[i]*global->pixelSize - global->pixelSize/2;
			// assign Z coordinate as detector distance and calculate length of the vectors to the pixel coordinates
			for (int j = 0; j < 4; j++) {
				corner_coordinates[j][2] = threadInfo->detectorPosition/1000;
				corner_distances[j] = sqrt(corner_coordinates[j][0]*corner_coordinates[j][0] + corner_coordinates[j][1]*corner_coordinates[j][1] + corner_coordinates[j][2]*corner_coordinates[j][2]);
			}
			
			// first triangle made up of upper left, upper right, and lower right corner
			// nominator in expression for solid angle of a plane triangle - magnitude of triple product of first 3 corners
			determinant = fabs( corner_coordinates[0][0]*(corner_coordinates[1][1]*corner_coordinates[2][2] - corner_coordinates[1][2]*corner_coordinates[2][1])
							   - corner_coordinates[0][1]*(corner_coordinates[1][0]*corner_coordinates[2][2] - corner_coordinates[1][2]*corner_coordinates[2][0])
							   + corner_coordinates[0][2]*(corner_coordinates[1][0]*corner_coordinates[2][1] - corner_coordinates[1][1]*corner_coordinates[2][0]) );
			denominator = corner_distances[0]*corner_distances[1]*corner_distances[2] + corner_distances[2]*(corner_coordinates[0][0]*corner_coordinates[1][0] + corner_coordinates[0][1]*corner_coordinates[1][1] + corner_coordinates[0][2]*corner_coordinates[1][2])
			+ corner_distances[1]*(corner_coordinates[0][0]*corner_coordinates[2][0] + corner_coordinates[0][1]*corner_coordinates[2][1] + corner_coordinates[0][2]*corner_coordinates[2][2])
			+ corner_distances[0]*(corner_coordinates[1][0]*corner_coordinates[2][0] + corner_coordinates[1][1]*corner_coordinates[2][1] + corner_coordinates[1][2]*corner_coordinates[2][2]);
			solid_angle[0] = atan2(determinant, denominator);
			if (solid_angle[0] < 0)
				solid_angle[0] += M_PI; // If det > 0 and denom < 0 arctan2 returns < 0, so add PI
			
			// second triangle made up of lower right, lower left, and upper left corner
			// nominator in expression for solid angle of a plane triangle - magnitude of triple product of last 3 corners
			determinant = fabs( corner_coordinates[0][0]*(corner_coordinates[3][1]*corner_coordinates[2][2] - corner_coordinates[3][2]*corner_coordinates[2][1])
							   - corner_coordinates[0][1]*(corner_coordinates[3][0]*corner_coordinates[2][2] - corner_coordinates[3][2]*corner_coordinates[2][0])
							   + corner_coordinates[0][2]*(corner_coordinates[3][0]*corner_coordinates[2][1] - corner_coordinates[3][1]*corner_coordinates[2][0]) );
			denominator = corner_distances[2]*corner_distances[3]*corner_distances[0] + corner_distances[2]*(corner_coordinates[0][0]*corner_coordinates[3][0] + corner_coordinates[0][1]*corner_coordinates[3][1] + corner_coordinates[0][2]*corner_coordinates[3][2])
			+ corner_distances[3]*(corner_coordinates[0][0]*corner_coordinates[2][0] + corner_coordinates[0][1]*corner_coordinates[2][1] + corner_coordinates[0][2]*corner_coordinates[2][2])
			+ corner_distances[0]*(corner_coordinates[3][0]*corner_coordinates[2][0] + corner_coordinates[3][1]*corner_coordinates[2][1] + corner_coordinates[3][2]*corner_coordinates[2][2]);
			solid_angle[1] = atan2(determinant, denominator);
			if (solid_angle[1] < 0)
				solid_angle[1] += M_PI; // If det > 0 and denom < 0 arctan2 returns < 0, so add PI
			
			total_solid_angle = 2*(solid_angle[0] + solid_angle[1]);
			
			threadInfo->corrected_data[i] /= total_solid_angle/threadInfo->solidAngle; // remove constant term to only get theta/phi dependent part of solid angle correction for 2D pattern
		}
	}
}


/*
 *	Apply attenuation correction
 */
void applyAttenuationCorrection(tThreadInfo *threadInfo, cGlobal *global){
	
	for(long i=0; i<global->pix_nn; i++) 
		threadInfo->corrected_data[i] *= threadInfo->attenuation;
	
}


/*
 *	Maintain running powder patterns
 */
void addToPowder(tThreadInfo *threadInfo, cGlobal *global, cHit *hit){

    if (global->generateDarkcal){
		// Sum raw format data
		pthread_mutex_lock(&global->powdersumraw_mutex);
		global->npowder += 1;
		for(long i=0; i<global->pix_nn; i++)
			global->powderRaw[i] += threadInfo->corrected_data[i];
		pthread_mutex_unlock(&global->powdersumraw_mutex);
        
		// Sum variance data
		pthread_mutex_lock(&global->powdersumvariance_mutex);
		for(long i=0; i<global->pix_nn; i++)
			global->powderVariance[i] += ((double) threadInfo->corrected_data[i]*threadInfo->corrected_data[i]);
		pthread_mutex_unlock(&global->powdersumvariance_mutex);
	}
	
    //standard hit
	if (hit->standard || global->listfinder.use){
		if (global->powdersum) {
			// Sum raw format data
			pthread_mutex_lock(&global->powdersumraw_mutex);
			for(long i=0; i<global->pix_nn; i++)
				global->powderRaw[i] += threadInfo->corrected_data[i];
			pthread_mutex_unlock(&global->powdersumraw_mutex);
			
			// Sum assembled data		
			pthread_mutex_lock(&global->powdersumassembled_mutex);
			global->npowder += 1;
			for(long i=0; i<global->image_nn; i++)
				if(threadInfo->image[i] > global->powderthresh)
					global->powderAssembled[i] += threadInfo->image[i];
			pthread_mutex_unlock(&global->powdersumassembled_mutex);
		}		
	}
	
	//ice hit
	if (hit->ice){
		if (global->powdersum) {
			// Sum raw format data 	: ice
			pthread_mutex_lock(&global->icesumraw_mutex);
			for(long i=0; i<global->pix_nn; i++)
				global->iceRaw[i] += threadInfo->corrected_data[i];
			pthread_mutex_unlock(&global->icesumraw_mutex);
			
			// Sum assembled data	: ice
			pthread_mutex_lock(&global->icesumassembled_mutex);
			global->nice += 1;
			for(long i=0; i<global->image_nn; i++)
				if(threadInfo->image[i] > global->powderthresh)
					global->iceAssembled[i] += threadInfo->image[i];
			pthread_mutex_unlock(&global->icesumassembled_mutex);			
		}
	}
	
	//water hit
	if (hit->water){
		if (global->powdersum) {
			// Sum raw format data  : water
			pthread_mutex_lock(&global->watersumraw_mutex);
			for(long i=0; i<global->pix_nn; i++)
				global->waterRaw[i] += threadInfo->corrected_data[i];
			pthread_mutex_unlock(&global->watersumraw_mutex);
			
			// Sum assembled data	: water
			pthread_mutex_lock(&global->watersumassembled_mutex);
			global->nwater += 1;
			for(long i=0; i<global->image_nn; i++)
				if(threadInfo->image[i] > global->powderthresh)
					global->waterAssembled[i] += threadInfo->image[i];
			pthread_mutex_unlock(&global->watersumassembled_mutex);
		}
	}
}


/*
 *	Maintain running correlation sum
 */
void addToCorrelation(tThreadInfo *threadInfo, cGlobal *global, cHit *hit) {

	if (global->useCorrelation && global->sumCorrelation && threadInfo->correlation) {
		
		if (hit->standard || global->listfinder.use) {
			// Sum correlation data
			pthread_mutex_lock(&global->powdersumcorrelation_mutex);
			if (!global->powdersum) global->npowder += 1;
			for(long i=0; i<global->correlation_nn; i++)
				global->powderCorrelation[i] += threadInfo->correlation[i];
			pthread_mutex_unlock(&global->powdersumcorrelation_mutex);
		}
		
		if (hit->ice) {
			// Sum correlation data 	: ice
			pthread_mutex_lock(&global->icesumcorrelation_mutex);
			if (!global->powdersum) global->nice += 1;
			for(long i=0; i<global->correlation_nn; i++)
				global->iceCorrelation[i] += threadInfo->correlation[i];
			pthread_mutex_unlock(&global->icesumcorrelation_mutex);
		}
		
		if (hit->water) {
			// Sum correlation data  : water
			pthread_mutex_lock(&global->watersumcorrelation_mutex);
			if (!global->powdersum) global->nwater += 1;
			for(long i=0; i<global->correlation_nn; i++)
				global->waterCorrelation[i] += threadInfo->correlation[i];
			pthread_mutex_unlock(&global->watersumcorrelation_mutex);
		}
	}
}


/*
 *	Interpolate raw (corrected) cspad data into a physical 2D image
 *	using pre-defined pixel mapping (as loaded from .h5 file)
 *  This function is used for powder averages
 */
void assemble2Dimage(double corrected_data[], double *&image, cGlobal *global){
	
	// Allocate temporary arrays for pixel interpolation (needs to be floating point)
	double	*data = (double*) calloc(global->image_nn,sizeof(double));
	double	*weight = (double*) calloc(global->image_nn,sizeof(double));	
	
	// Loop through all pixels and interpolate onto regular grid
	float	x, y;
	double	pixel_value, w;
	long	ix, iy;
	float	fx, fy;
	long	image_index;
	
	for(long i=0;i<global->pix_nn;i++){
		// Pixel location with (0,0) at array element (0,0) in bottom left corner
		x = global->pix_x[i] + global->image_nx/2;
		y = global->pix_y[i] + global->image_nx/2;
		pixel_value = corrected_data[i];
		
		// Split coordinate into integer and fractional parts
		ix = (long) floor(x);
		iy = (long) floor(y);
		fx = x - ix;
		fy = y - iy;
		
		// Interpolate intensity over adjacent 4 pixels using fractional overlap as the weighting factor
		// (0,0)
		if(ix>=0 && iy>=0 && ix<global->image_nx && iy<global->image_nx) {
			w = (double) (1-fx)*(1-fy);
			image_index = ix + global->image_nx*iy;
			data[image_index] += w*pixel_value;
			weight[image_index] += w;
		}
		// (+1,0)
		if((ix+1)>=0 && iy>=0 && (ix+1)<global->image_nx && iy<global->image_nx) {
			w = (double) (fx)*(1-fy);
			image_index = (ix+1) + global->image_nx*iy;
			data[image_index] += w*pixel_value;
			weight[image_index] += w;
		}
		// (0,+1)
		if(ix>=0 && (iy+1)>=0 && ix<global->image_nx && (iy+1)<global->image_nx) {
			w = (double) (1-fx)*(fy);
			image_index = ix + global->image_nx*(iy+1);
			data[image_index] += w*pixel_value;
			weight[image_index] += w;
		}
		// (+1,+1)
		if((ix+1)>=0 && (iy+1)>=0 && (ix+1)<global->image_nx && (iy+1)<global->image_nx) {
			w = (double) (fx)*(fy);
			image_index = (ix+1) + global->image_nx*(iy+1);
			data[image_index] += w*pixel_value;
			weight[image_index] += w;
		}
	}
	
	
	// Reweight pixel interpolation
	for(long i=0; i<global->image_nn; i++){
		if(weight[i] < 0.05)
			data[i] = 0;
		else
			data[i] /= weight[i];
	}
	
	
	// Allocate memory for output image
	if (image)
		free(image);
	image = (double*) calloc(global->image_nn, sizeof(double));
	
	// Copy interpolated image across into image array
	for(long i=0;i<global->image_nn;i++){
		image[i] = data[i];
	}
	
	// Free temporary arrays
	free(data);
	free(weight);
	
}


/*
 *	Interpolate raw (corrected) cspad data into a physical 2D image
 *	using pre-defined pixel mapping (as loaded from .h5 file)
 *  This function is used for individual shots
 */
void assemble2Dimage(tThreadInfo *threadInfo, cGlobal *global){
	
	// Allocate temporary arrays for pixel interpolation (needs to be floating point)
	float	*data = (float*) calloc(global->image_nn,sizeof(float));
	float	*weight = (float*) calloc(global->image_nn,sizeof(float));	
	
	// Loop through all pixels and interpolate onto regular grid
	float	x, y;
	float	pixel_value, w;
	long	ix, iy;
	float	fx, fy;
	long	image_index;
	
	for(long i=0;i<global->pix_nn;i++){
		// Pixel location with (0,0) at array element (0,0) in bottom left corner
		x = global->pix_x[i] + global->image_nx/2;
		y = global->pix_y[i] + global->image_nx/2;
		pixel_value = threadInfo->corrected_data[i];
		
		// Split coordinate into integer and fractional parts
		ix = (long) floor(x);
		iy = (long) floor(y);
		fx = x - ix;
		fy = y - iy;
		
		// Interpolate intensity over adjacent 4 pixels using fractional overlap as the weighting factor
		// (0,0)
		if(ix>=0 && iy>=0 && ix<global->image_nx && iy<global->image_nx) {
			w = (1-fx)*(1-fy);
			image_index = ix + global->image_nx*iy;
			data[image_index] += w*pixel_value;
			weight[image_index] += w;
		}
		// (+1,0)
		if((ix+1)>=0 && iy>=0 && (ix+1)<global->image_nx && iy<global->image_nx) {
			w = (fx)*(1-fy);
			image_index = (ix+1) + global->image_nx*iy;
			data[image_index] += w*pixel_value;
			weight[image_index] += w;
		}
		// (0,+1)
		if(ix>=0 && (iy+1)>=0 && ix<global->image_nx && (iy+1)<global->image_nx) {
			w = (1-fx)*(fy);
			image_index = ix + global->image_nx*(iy+1);
			data[image_index] += w*pixel_value;
			weight[image_index] += w;
		}
		// (+1,+1)
		if((ix+1)>=0 && (iy+1)>=0 && (ix+1)<global->image_nx && (iy+1)<global->image_nx) {
			w = (fx)*(fy);
			image_index = (ix+1) + global->image_nx*(iy+1);
			data[image_index] += w*pixel_value;
			weight[image_index] += w;
		}
	}
	
	
	// Reweight pixel interpolation
	for(long i=0; i<global->image_nn; i++){
		if(weight[i] < 0.05)
			data[i] = 0;
		else
			data[i] /= weight[i];
	}
	
	
	// Allocate memory for output image
	threadInfo->image = (float*) calloc(global->image_nn,sizeof(float));
	
	// Copy interpolated image across into image array
	for(long i=0;i<global->image_nn;i++){
		threadInfo->image[i] = data[i];
	}
	
	// Free temporary arrays
	free(data);
	free(weight);
	
}



/*
 *	Create filename based on date, time and fiducial for this image
 */
void nameEvent(tThreadInfo *info, cGlobal *global){
//	char outfile[1024];
	char buffer1[80];
	char buffer2[80];	
	time_t eventTime = info->seconds;
	
	//setenv("TZ","US/Pacific",1);		// <--- Dangerous (not thread safe!)
	struct tm *timestatic, timelocal;
	timestatic=localtime_r( &eventTime, &timelocal );	
	strftime(buffer1,80,"%Y_%b%d",&timelocal);
	strftime(buffer2,80,"%H%M%S",&timelocal);
	sprintf(info->eventname,"LCLS_%s_r%04u_%s_%x_cspad",buffer1,info->runNumber,buffer2,info->fiducial);
}



/*
 *	Check if the list of hits contains the current event
 */
bool containsEvent(string event, cGlobal *global) {
	if (global->nhits < global->hitlist.size()) {
		return (event == global->hitlist[global->nhits]);
	}
}



/*
 *	Write out processed data to our 'standard' HDF5 format
 */
void writeHDF5(tThreadInfo *info, cGlobal *global, char *eventname){
	/*
	 *	Create filename based on date, time and fiducial for this image
	 */
	char outfile[1024];
	//char buffer1[80];
	//char buffer2[80];	
	//time_t eventTime = info->seconds;

	//setenv("TZ","US/Pacific",1);		// <--- Dangerous (not thread safe!)
	//struct tm *timestatic, timelocal;
	//timestatic=localtime_r( &eventTime, &timelocal );	
	//strftime(buffer1,80,"%Y_%b%d",&timelocal);
	//strftime(buffer2,80,"%H%M%S",&timelocal);
	//sprintf(outfile,"LCLS_%s_r%04u_%s_%x_cspad.h5",buffer1,info->runNumber,buffer2,info->fiducial);

	sprintf(outfile,"%s.h5",eventname);
	//strcpy(outfile, info->eventname);
	DEBUGL1_ONLY printf("r%04u:%i (%2.1f Hz): Writing data to: %s\n", (int)info->runNumber, (int)info->threadNum, global->datarate, outfile);
	
	
	/* 
 	 *  HDF5 variables
	 */
	hid_t		hdf_fileID;
	hid_t		dataspace_id;
	hid_t		dataset_id;
	hid_t		datatype;
	hsize_t 	size[2],max_size[2];
	herr_t		hdf_error;
	hid_t   	gid;
//	char 		fieldname[100]; 
	
	
	/*
	 *	Create the HDF5 file
	 */
	hdf_fileID = H5Fcreate(outfile,  H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	
	
	
	
	/*
	 *	Save image data into '/data' part of HDF5 file
	 */
	gid = H5Gcreate(hdf_fileID, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( gid < 0 ) {
		ERROR("%i: Couldn't create group\n", (int)info->threadNum);
		H5Fclose(hdf_fileID);
		return;
	}
	
	// Assembled image
	size[0] = global->image_nx;	// size[0] = height
	size[1] = global->image_nx;	// size[1] = width
	max_size[0] = global->image_nx;
	max_size[1] = global->image_nx;
	int16_t *buffer1 = (int16_t*) calloc(global->image_nn, sizeof(int16_t));
	for(long i=0; i<global->image_nn; i++){
		buffer1[i] = (int16_t) info->image[i];
	}
	dataspace_id = H5Screate_simple(2, size, max_size);
	dataset_id = H5Dcreate(gid, "data", H5T_STD_I16LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( dataset_id < 0 ) {
		ERROR("%i: Couldn't create dataset\n", (int)info->threadNum);
		H5Fclose(hdf_fileID);
		return;
	}
	hdf_error = H5Dwrite(dataset_id, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer1);
	if ( hdf_error < 0 ) {
		ERROR("%i: Couldn't write data\n", (int)info->threadNum);
		H5Dclose(dataspace_id);
		H5Fclose(hdf_fileID);
		return;
	}
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	free(buffer1);	
	
	// Save raw data?
	if(global->saveRaw) {
		size[0] = 8*COLS;	// size[0] = height
		size[1] = 8*ROWS;	// size[1] = width
		max_size[0] = 8*COLS;
		max_size[1] = 8*ROWS;
		int16_t *buffer2 = (int16_t*) calloc(RAW_DATA_LENGTH, sizeof(int16_t));
		for(long i=0; i<RAW_DATA_LENGTH; i++){
			buffer2[i] = (int16_t) info->corrected_data[i];
		}
		dataspace_id = H5Screate_simple(2, size, max_size);
		dataset_id = H5Dcreate(gid, "rawdata", H5T_STD_I16LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if ( dataset_id < 0 ) {
			ERROR("%i: Couldn't create dataset\n", (int)info->threadNum);
			H5Fclose(hdf_fileID);
			return;
		}
		hdf_error = H5Dwrite(dataset_id, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer2);
		if ( hdf_error < 0 ) {
			ERROR("%i: Couldn't write data\n", (int)info->threadNum);
			H5Dclose(dataspace_id);
			H5Fclose(hdf_fileID);
			return;
		}
		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		free(buffer2);
	}

	
	// Done with this group
	H5Gclose(gid);
	
	
//	double		phaseCavityTime1;
//	double		phaseCavityTime2;
//	double		phaseCavityCharge1;
//	double		phaseCavityCharge2;
	
	/*
	 *	Write LCLS event information
	 */
	gid = H5Gcreate1(hdf_fileID,"LCLS",0);
	size[0] = 1;
	dataspace_id = H5Screate_simple( 1, size, NULL );
	//dataspace_id = H5Screate(H5S_SCALAR);
	
	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/machineTime", H5T_NATIVE_INT32, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->seconds );
	H5Dclose(dataset_id);
	
	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/fiducial", H5T_NATIVE_INT32, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->fiducial );
	H5Dclose(dataset_id);
		
	// Electron beam data
	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/ebeamCharge", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->fEbeamCharge );
	H5Dclose(dataset_id);

	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/ebeamL3Energy", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->fEbeamL3Energy );
	H5Dclose(dataset_id);
	
	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/ebeamPkCurrBC2", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->fEbeamPkCurrBC2 );
	H5Dclose(dataset_id);

	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/ebeamLTUPosX", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->fEbeamLTUPosX );
	H5Dclose(dataset_id);

	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/ebeamLTUPosY", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->fEbeamLTUPosY );
	H5Dclose(dataset_id);
	
	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/ebeamLTUAngX", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->fEbeamLTUAngX );
	H5Dclose(dataset_id);

	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/ebeamLTUAngY", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->fEbeamLTUAngY );
	H5Dclose(dataset_id);

	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/phaseCavityTime1", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->phaseCavityTime1 );
	H5Dclose(dataset_id);
	
	// Phase cavity information
	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/phaseCavityTime2", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->phaseCavityTime2 );
	H5Dclose(dataset_id);

	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/phaseCavityCharge1", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->phaseCavityCharge1 );
	H5Dclose(dataset_id);

	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/phaseCavityCharge2", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->phaseCavityCharge2 );
	H5Dclose(dataset_id);
	
	// Calculated photon energy
	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/photon_energy_eV", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->photonEnergyeV);
	H5Dclose(dataset_id);
	
	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/photon_wavelength_A", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->wavelengthA);
	H5Dclose(dataset_id);
	
	
	// Gas detector values
	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/f_11_ENRC", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->gmd11 );
	H5Dclose(dataset_id);

	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/f_12_ENRC", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->gmd12 );
	H5Dclose(dataset_id);

	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/f_21_ENRC", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->gmd21 );
	H5Dclose(dataset_id);

	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/f_22_ENRC", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->gmd22 );	
	H5Dclose(dataset_id);

	
	// Motor positions
	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/detectorPosition", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->detectorPosition );	
	H5Dclose(dataset_id);
	
	
	// Attenuation
	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/attenuation", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->attenuation );
	H5Dclose(dataset_id);
	
	
	// Solid Angle
	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/solidAngle", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->solidAngle );
	H5Dclose(dataset_id);
	
	
	// IntensityAvg
	dataset_id = H5Dcreate1(hdf_fileID, "/LCLS/intensityAvg", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &info->intensityAvg );
	H5Dclose(dataset_id);
	
	
	// Finished with scalar dataset ID
	H5Sclose(dataspace_id);
	
	
	// Time in human readable format
	// Writing strings in HDF5 is a little tricky --> this could be improved!
	char* timestr;
	time_t eventTime = info->seconds;
	timestr = ctime(&eventTime);
	dataspace_id = H5Screate(H5S_SCALAR);
	datatype = H5Tcopy(H5T_C_S1);  
	H5Tset_size(datatype,strlen(timestr)+1);
	dataset_id = H5Dcreate1(hdf_fileID, "LCLS/eventTimeString", datatype, dataspace_id, H5P_DEFAULT);
	H5Dwrite(dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, timestr );
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	hdf_error = H5Lcreate_soft( "/LCLS/eventTimeString", hdf_fileID, "/LCLS/eventTime",0,0);
	
	
	// Close group and flush buffers
	H5Gclose(gid);
	H5Fflush(hdf_fileID,H5F_SCOPE_LOCAL);

	
	/*
	 *	Clean up stale HDF5 links
	 *		(thanks Tom/Filipe)
	 */
	int n_ids;
	hid_t ids[256];
	n_ids = (int)H5Fget_obj_ids(hdf_fileID, H5F_OBJ_ALL, 256, ids);
	for ( int i=0; i<n_ids; i++ ) {
		hid_t id;
		H5I_type_t type;
		id = ids[i];
		type = H5Iget_type(id);
		if ( type == H5I_GROUP ) H5Gclose(id);
		if ( type == H5I_DATASET ) H5Dclose(id);
		if ( type == H5I_DATATYPE ) H5Tclose(id);
		if ( type == H5I_DATASPACE ) H5Sclose(id);
		if ( type == H5I_ATTR ) H5Aclose(id);
	}
	
	H5Fclose(hdf_fileID); 
}





/*
 *	Write test data to a simple HDF5 file
 */

void writeSimpleHDF5(const char *filename, const void *data, int width, int height, int type) 
{
	hid_t fh, gh, sh, dh;	/* File, group, dataspace and data handles */
	herr_t r;
	hsize_t size[2];
	hsize_t max_size[2];
	
	fh = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't create file: %s\n", filename);
	}
	
	gh = H5Gcreate(fh, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( gh < 0 ) {
		ERROR("Couldn't create group\n");
		H5Fclose(fh);
	}
	
	size[0] = height;
	size[1] = width;
	max_size[0] = height;
	max_size[1] = width;
	sh = H5Screate_simple(2, size, max_size);
	
	dh = H5Dcreate(gh, "data", type, sh,
	               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset\n");
		H5Fclose(fh);
	}
	
	/* Muppet check */
	H5Sget_simple_extent_dims(sh, size, max_size);
	
	r = H5Dwrite(dh, type, H5S_ALL,
	             H5S_ALL, H5P_DEFAULT, data);
	if ( r < 0 ) {
		ERROR("Couldn't write data\n");
		H5Dclose(dh);
		H5Fclose(fh);
	}
	
	H5Gclose(gh);
	H5Dclose(dh);
	H5Fclose(fh);
}

void writeSimpleHDF5(const char *filename, const void *data, int width, int height, int depth, int type) 
{
	hid_t fh, gh, sh, dh;	/* File, group, dataspace and data handles */
	herr_t r;
	hsize_t size[3];
	hsize_t max_size[3];
	
	fh = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't create file: %s\n", filename);
	}
	
	gh = H5Gcreate(fh, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( gh < 0 ) {
		ERROR("Couldn't create group\n");
		H5Fclose(fh);
	}
	
	size[0] = height;
	size[1] = width;
	size[2] = depth;
	max_size[0] = height;
	max_size[1] = width;
	max_size[2] = depth;
	sh = H5Screate_simple(3, size, max_size);
	
	dh = H5Dcreate(gh, "data", type, sh,
	               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ( dh < 0 ) {
		ERROR("Couldn't create dataset\n");
		H5Fclose(fh);
	}
	
	/* Muppet check */
	H5Sget_simple_extent_dims(sh, size, max_size);
	
	r = H5Dwrite(dh, type, H5S_ALL,
	             H5S_ALL, H5P_DEFAULT, data);
	if ( r < 0 ) {
		ERROR("Couldn't write data\n");
		H5Dclose(dh);
		H5Fclose(fh);
	}
	
	H5Gclose(gh);
	H5Dclose(dh);
	H5Fclose(fh);
}


void flushHDF5() {
	herr_t fail;
	
	// Flush all HDF5 data to disk, close all open identifies, and clean up memory
	cout << "Flushing HDF5 data and memory..." << endl;
	fail = H5close();
	if (fail < 0) cout << "\tError flushing HDF5 data and memory" << endl;
	else cout << "\tFlushing finished successfully!" << endl;
}


void saveRunningSums(cGlobal *global) {
	char	filename[1024];

	if(global->generateDarkcal) {
		/*
		 *	Compute and save darkcal
		 */
		printf("Processing darkcal\n");
		sprintf(filename,"r%04u-darkcal.h5",global->runNumber);
		float *buffer3 = (float*) calloc(global->pix_nn, sizeof(float));
		pthread_mutex_lock(&global->powdersumraw_mutex);
		for(long i=0; i<global->pix_nn; i++)
			buffer3[i] = (float) (global->powderRaw[i]/global->npowder);
		pthread_mutex_unlock(&global->powdersumraw_mutex);
		//for(long i=0; i<global->pix_nn; i++)
		//	if (buffer3[i] < 0) buffer3[i] = 0;
		printf("Saving darkcal to file\n");
		writeSimpleHDF5(filename, buffer3, (int)global->pix_nx, (int)global->pix_ny, H5T_NATIVE_FLOAT);	
		
		/*
		 *	Compute and save variance of darkcal
		 */
		printf("Saving variance of darkcal to file\n");
		sprintf(filename,"r%04u-darkcal_variance.h5",global->runNumber);
		float *buffer12 = (float*) calloc(global->pix_nn, sizeof(float));
		pthread_mutex_lock(&global->powdersumvariance_mutex);
		for(long i=0; i<global->pix_nn; i++)
			buffer12[i] = (float) (global->powderVariance[i]/global->npowder - global->powderRaw[i]*global->powderRaw[i]/(global->npowder*global->npowder));
		pthread_mutex_unlock(&global->powdersumvariance_mutex);
		writeSimpleHDF5(filename, buffer12, (int)global->pix_nx, (int)global->pix_ny, H5T_NATIVE_FLOAT);		
		
		DEBUGL1_ONLY {
			cout << "Pixel statistics:" << endl;
			for(int i=0; i<global->nPixels; i++) {
				cout << "\tPIXEL #" << i << " = " << buffer12[global->pixels[i]] << " ADUs" << endl;
			}
			
			// VERIFY OUTPUT IN ASSEMBLED FORMAT
			array1D<double> *oneX = new array1D<double>(global->pix_x, global->pix_nn);
			array1D<double> *oneY = new array1D<double>(global->pix_y, global->pix_nn);
			array1D<double> *oneDark = new array1D<double>(buffer3, global->pix_nn);		
			array1D<double> *oneDarkVariance = new array1D<double>(buffer12, global->pix_nn);		
			arraydataIO *io = new arraydataIO();
			array2D<double> *two = 0;
			
			/*
			 *	Save assembled darkcal
			 */			
			printf("Saving assembled darkcal image to file\n");
			sprintf(filename,"r%04u-AssembledSum.h5",global->runNumber);
			ns_cspad_util::createAssembledImageCSPAD( oneDark, oneX, oneY, two );		
			io->writeToFile( filename, two );
			
			/*
			 *	Save assembled variance of darkcal
			 */
			printf("Saving assembled variance image to file\n");
			sprintf(filename,"r%04u-AssembledVariance.h5",global->runNumber);
			ns_cspad_util::createAssembledImageCSPAD( oneDarkVariance, oneX, oneY, two );
			io->writeToFile( filename, two );
			
			delete oneX;
			delete oneY;
			delete oneDark;
			delete oneDarkVariance;
			delete two;
			delete io;
		}
		
		free(buffer3);
		free(buffer12);
		
	}

	else {
		if(global->hitfinder.use || global->listfinder.use) {
			if (global->powdersum) {
				
				/*
				 *	Save assembled powder pattern
				 */
				printf("Saving assembled sum data to file\n");
				sprintf(filename,"r%04u-AssembledSum.h5",global->runNumber);
				float *buffer2 = (float*) calloc(global->image_nn, sizeof(float));
				pthread_mutex_lock(&global->powdersumassembled_mutex);
				for(long i=0; i<global->image_nn; i++){
					buffer2[i] = (float) (global->powderAssembled[i]/global->npowder);
				}
				pthread_mutex_unlock(&global->powdersumassembled_mutex);
				writeSimpleHDF5(filename, buffer2, (int)global->image_nx, (int)global->image_nx, H5T_NATIVE_FLOAT);	
				free(buffer2);
				
			}
			
			if (global->powdersum && global->saveRaw) {
				
				/*
				 *	Save powder pattern in raw layout
				 */
				printf("Saving raw sum data to file\n");
				sprintf(filename,"r%04u-RawSum.h5",global->runNumber);
				float *buffer1 = (float*) calloc(global->pix_nn, sizeof(float));
				pthread_mutex_lock(&global->powdersumraw_mutex);
				for(long i=0; i<global->pix_nn; i++)
					buffer1[i] = (float) (global->powderRaw[i]/global->npowder);
				pthread_mutex_unlock(&global->powdersumraw_mutex);
				writeSimpleHDF5(filename, buffer1, (int)global->pix_nx, (int)global->pix_ny, H5T_NATIVE_FLOAT);	
				free(buffer1);
				
			}			
			
			if (global->useCorrelation && global->sumCorrelation) {
				
				/*
				 *	Save correlation sum
				 */
				printf("Saving correlation sum data to file\n");
				sprintf(filename,"r%04u-CorrelationSum.h5",global->runNumber);
				double *buffer8 = (double*) calloc(global->correlation_nn, sizeof(double));
				pthread_mutex_lock(&global->powdersumcorrelation_mutex);
				for(long i=0; i<global->correlation_nn; i++)
					buffer8[i] = global->powderCorrelation[i]/global->npowder;
				pthread_mutex_unlock(&global->powdersumcorrelation_mutex);
				if (global->autoCorrelateOnly) writeSimpleHDF5(filename, buffer8, global->correlationNumDelta, global->correlationNumQ, H5T_NATIVE_DOUBLE);
				else writeSimpleHDF5(filename, buffer8, global->correlationNumDelta, global->correlationNumQ, global->correlationNumQ, H5T_NATIVE_DOUBLE);
				free(buffer8);
				
			}
		}
			
			
		
		
		if(global->icefinder.use) {
			if (global->powdersum) {
				
				/*
				 *	Save assembled powder pattern : ice
				 */
				printf("Saving assembled sum data of ice to file\n");
				sprintf(filename,"r%04u-AssembledSum_ice.h5",global->runNumber);
				float *buffer4 = (float*) calloc(global->image_nn, sizeof(float));
				pthread_mutex_lock(&global->icesumassembled_mutex);
				for(long i=0; i<global->image_nn; i++){
					buffer4[i] = (float) (global->iceAssembled[i]/global->nice);
				}
				pthread_mutex_unlock(&global->icesumassembled_mutex);
				writeSimpleHDF5(filename, buffer4, (int)global->image_nx, (int)global->image_nx, H5T_NATIVE_FLOAT);	
				free(buffer4);
				
			}
			
			if (global->powdersum && global->saveRaw) {
				
				/*
				 *	Save powder pattern in raw layout : ice
				 */
				printf("Saving raw sum data of ice to file\n");
				sprintf(filename,"r%04u-RawSum_ice.h5",global->runNumber);
				float *buffer5 = (float*) calloc(global->pix_nn, sizeof(float));
				pthread_mutex_lock(&global->icesumraw_mutex);
				for(long i=0; i<global->pix_nn; i++)
					buffer5[i] = (float) (global->iceRaw[i]/global->nice);
				pthread_mutex_unlock(&global->icesumraw_mutex);
				writeSimpleHDF5(filename, buffer5, (int)global->pix_nx, (int)global->pix_ny, H5T_NATIVE_FLOAT);	
				free(buffer5);
				
			}
			
			if (global->useCorrelation && global->sumCorrelation) {
				
				/*
				 *	Save correlation sum : ice
				 */
				printf("Saving correlation sum data of ice to file\n");
				sprintf(filename,"r%04u-CorrelationSum_ice.h5",global->runNumber);
				double *buffer9 = (double*) calloc(global->correlation_nn, sizeof(double));
				pthread_mutex_lock(&global->icesumcorrelation_mutex);
				for(long i=0; i<global->correlation_nn; i++)
					buffer9[i] = global->iceCorrelation[i]/global->nice;
				pthread_mutex_unlock(&global->icesumcorrelation_mutex);
				if (global->autoCorrelateOnly) writeSimpleHDF5(filename, buffer9, global->correlationNumDelta, global->correlationNumQ, H5T_NATIVE_DOUBLE);
				else writeSimpleHDF5(filename, buffer9, global->correlationNumDelta, global->correlationNumQ, global->correlationNumQ, H5T_NATIVE_DOUBLE);
				free(buffer9);
				
			}
		}		

		
		
		if(global->waterfinder.use) {
			if (global->powdersum) {
				
				/*
				 *	Save assembled powder pattern : water
				 */
				printf("Saving assembled sum data of water to file\n");
				sprintf(filename,"r%04u-AssembledSum_water.h5",global->runNumber);
				float *buffer6 = (float*) calloc(global->image_nn, sizeof(float));
				pthread_mutex_lock(&global->watersumassembled_mutex);
				for(long i=0; i<global->image_nn; i++){
					buffer6[i] = (float) (global->waterAssembled[i]/global->nwater);
				}
				pthread_mutex_unlock(&global->watersumassembled_mutex);
				writeSimpleHDF5(filename, buffer6, (int)global->image_nx, (int)global->image_nx, H5T_NATIVE_FLOAT);	
				free(buffer6);				
				
			}

			if (global->powdersum && global->saveRaw) {
				
				/*
				 *	Save powder pattern in raw layout : water
				 */
				printf("Saving raw sum data of water to file\n");
				sprintf(filename,"r%04u-RawSum_water.h5",global->runNumber);
				float *buffer7 = (float*) calloc(global->pix_nn, sizeof(float));
				pthread_mutex_lock(&global->watersumraw_mutex);
				for(long i=0; i<global->pix_nn; i++)
					buffer7[i] = (float) (global->waterRaw[i]/global->nwater);
				pthread_mutex_unlock(&global->watersumraw_mutex);
				writeSimpleHDF5(filename, buffer7, (int)global->pix_nx, (int)global->pix_ny, H5T_NATIVE_FLOAT);	
				free(buffer7);				
				
			}
			
			if (global->useCorrelation && global->sumCorrelation) {
				
				/*
				 *	Save correlation sum : water
				 */
				printf("Saving correlation sum data of water to file\n");
				sprintf(filename,"r%04u-CorrelationSum_water.h5",global->runNumber);
				double *buffer10 = (double*) calloc(global->correlation_nn, sizeof(double));
				pthread_mutex_lock(&global->watersumcorrelation_mutex);
				for(long i=0; i<global->correlation_nn; i++)
					buffer10[i] = global->waterCorrelation[i]/global->nwater;
				pthread_mutex_unlock(&global->watersumcorrelation_mutex);
				if (global->autoCorrelateOnly) writeSimpleHDF5(filename, buffer10, global->correlationNumDelta, global->correlationNumQ, H5T_NATIVE_DOUBLE);
				else writeSimpleHDF5(filename, buffer10, global->correlationNumDelta, global->correlationNumQ, global->correlationNumQ, H5T_NATIVE_DOUBLE);
				free(buffer10);
				
			}
		}
		
		
	}

	
}


void calculatePowderAngularAvg(cGlobal *global) {
	
	// calculating angular average from the powder pattern 
	DEBUGL1_ONLY printf("calculating angular average from powder pattern...\n");
	
	// allocate local counter arrays
	unsigned *counterp = NULL;
	unsigned *counteri = NULL;
	unsigned *counterw = NULL;
	if (global->hitfinder.use || global->listfinder.use) {
		counterp = (unsigned*) calloc(global->angularAvg_nn, sizeof(unsigned));
	}
	if (global->icefinder.use) {
		counteri = (unsigned*) calloc(global->angularAvg_nn, sizeof(unsigned));
	}
	if (global->waterfinder.use) {
		counterw = (unsigned*) calloc(global->angularAvg_nn, sizeof(unsigned));
	}
		
	// angular average for each |q|
	DEBUGL1_ONLY printf("# of steps: %d\n",global->angularAvg_nn);
	DEBUGL2_ONLY printf("average SAXS intensity:\n");
	
	for (int i=0; i<global->pix_nn; i++) {
		if ( global->angularAvg_i[i] < global->angularAvg_nn && global->angularAvg_i[i] >= 0 && (!global->useBadPixelMask || global->badpixelmask[i]) ) {
			if (global->hitfinder.use || global->listfinder.use) {
				global->powderAverage[global->angularAvg_i[i]] += global->powderRaw[i];
				counterp[global->angularAvg_i[i]]++;
			}
			if (global->icefinder.use) {
				global->iceAverage[global->angularAvg_i[i]] += global->iceRaw[i];
				counteri[global->angularAvg_i[i]]++;
			}
			if (global->waterfinder.use) {
				global->waterAverage[global->angularAvg_i[i]] += global->waterRaw[i];
				counterw[global->angularAvg_i[i]]++;
			}
		}
	}
	
	if (global->useSolidAngleCorrection) {
		double solidAngle = global->pixelSize*global->pixelSize/(global->detectorZ/1000*global->detectorZ/1000);
		for (int i=0; i<global->angularAvg_nn; i++) {
			if (global->hitfinder.use || global->listfinder.use) {
				if (counterp[i]) global->powderAverage[i] /= counterp[i]*solidAngle;
			}
			if (global->icefinder.use) {
				if (counteri[i]) global->iceAverage[i] /= counteri[i]*solidAngle;
			}
			if (global->waterfinder.use) {
				if (counterw[i]) global->waterAverage[i] /= counterw[i]*solidAngle;
			}
			
			DEBUGL2_ONLY {
				if (global->hitfinder.use || global->listfinder.use) cout << "Q: " << global->angularAvgQ[i] << ",   \t# pixels: " << counterp[i] << ",\tI(powder): " << global->powderAverage[i] << endl;
				if (global->icefinder.use) cout << "Q: " << global->angularAvgQ[i] << ",   \t# pixels: " << counteri[i] << ",\tI(water): " << global->iceAverage[i] << endl;
				if (global->waterfinder.use) cout << "Q: " << global->angularAvgQ[i] << ",   \t# pixels: " << counterw[i] << ",\tI(ice): " << global->waterAverage[i] << endl;
			}
		}
	} else {
		for (int i=0; i<global->angularAvg_nn; i++) {
			if (global->hitfinder.use || global->listfinder.use) {
				if (counterp[i]) global->powderAverage[i] /= counterp[i];
			}
			if (global->icefinder.use) {
				if (counteri[i]) global->iceAverage[i] /= counteri[i];
			}
			if (global->waterfinder.use) {
				if (counterw[i]) global->waterAverage[i] /= counterw[i];
			}
			
			DEBUGL2_ONLY {
				if (global->hitfinder.use || global->listfinder.use) cout << "Q: " << global->angularAvgQ[i] << ",   \t# pixels: " << counterp[i] << ",\tI(powder): " << global->powderAverage[i] << endl;
				if (global->icefinder.use) cout << "Q: " << global->angularAvgQ[i] << ",   \t# pixels: " << counteri[i] << ",\tI(water): " << global->iceAverage[i] << endl;
				if (global->waterfinder.use) cout << "Q: " << global->angularAvgQ[i] << ",   \t# pixels: " << counterw[i] << ",\tI(ice): " << global->waterAverage[i] << endl;
			}
		}
	}
	
	// free memory of local variables
	free(counterp);
	free(counteri);
	free(counterw);
	
}


void savePowderAngularAvg(cGlobal *global) {
	
	if (global->powdersum && global->powderAngularAvg) {
		
		char	filename[1024];		
		double *buffer = (double*) calloc(2*global->angularAvg_nn, sizeof(double));
		for (long i=0; i<global->angularAvg_nn; i++) {
			if (global->useEnergyCalibration) buffer[i] = global->angularAvgQcal[i];
			else buffer[i] = global->angularAvgQ[i];
		}
		
		/*
		 *	Save angular average of powder pattern
		 */
		if (global->hitfinder.use || global->listfinder.use) {
			printf("Saving angular average of powder data to file\n");
			sprintf(filename,"r%04u-angavg.h5",global->runNumber);
			for(long i=0; i<global->angularAvg_nn; i++)
				buffer[global->angularAvg_nn+i] = global->powderAverage[i]/global->npowder;
			writeSimpleHDF5(filename, buffer, global->angularAvg_nn, 2, H5T_NATIVE_DOUBLE);
		}
		
		/*
		 *	Save angualar average of ice powder pattern
		 */
		if (global->icefinder.use) {
			printf("Saving angular average of powder ice data to file\n");
			sprintf(filename,"r%04u-angavg_ice.h5",global->runNumber);
			for(long i=0; i<global->angularAvg_nn; i++)
				buffer[global->angularAvg_nn+i] = global->iceAverage[i]/global->nice;
			writeSimpleHDF5(filename, buffer, global->angularAvg_nn, 2, H5T_NATIVE_DOUBLE);
		}
		
		/*
		 *	Save angular average of water powder pattern
		 */
		if (global->waterfinder.use) {
			printf("Saving angular average of powder water data to file\n");
			sprintf(filename,"r%04u-angavg_water.h5",global->runNumber);
			for(long i=0; i<global->angularAvg_nn; i++)
				buffer[global->angularAvg_nn+i] = global->waterAverage[i]/global->nwater;
			writeSimpleHDF5(filename, buffer, global->angularAvg_nn, 2, H5T_NATIVE_DOUBLE);
		}
		
		free(buffer);
		
	}			

}


void calculateQuadAngularAvg(cGlobal *global, int quad) {
	
	// calculating angular average from the powder pattern on a single quad
	
	// allocate local counter arrays
	unsigned *counterp = NULL;
	unsigned *counteri = NULL;
	unsigned *counterw = NULL;
	if (global->hitfinder.use || global->listfinder.use) {
		counterp = (unsigned*) calloc(global->angularAvg_nn, sizeof(unsigned));
	}
	if (global->icefinder.use) {
		counteri = (unsigned*) calloc(global->angularAvg_nn, sizeof(unsigned));
	}
	if (global->waterfinder.use) {
		counterw = (unsigned*) calloc(global->angularAvg_nn, sizeof(unsigned));
	}
	
	// angular average for each |q| for the pixels in the quad
	for(int mi=0; mi<2; mi++){ // mi decides what col in 2x8 matrix of raw data ASICs for each quad
		for(int mj=0; mj<8; mj++){ // mj decides what row in 2x8 matrix of raw data ASICs for each quad
			for(int i=0; i<ROWS; i++){
				for(int j=0; j<COLS; j++){
					long index = (j + mj*COLS) * (8*ROWS); // [row, 0] in 1552x1480 raw data format
					index += i + mi*ROWS + quad*2*ROWS; // [row, col] in 1552x1480 raw data format
					if ( global->angularAvg_i[index] < global->angularAvg_nn && global->angularAvg_i[index] >= 0 && (!global->useBadPixelMask || global->badpixelmask[index]) ) {
						if (global->hitfinder.use || global->listfinder.use) {
							global->powderAverage[global->angularAvg_i[index]] += global->powderRaw[index];
							counterp[global->angularAvg_i[index]]++;
						}
						if (global->icefinder.use) {
							global->iceAverage[global->angularAvg_i[index]] += global->iceRaw[index];
							counteri[global->angularAvg_i[index]]++;
						}
						if (global->waterfinder.use) {
							global->waterAverage[global->angularAvg_i[index]] += global->waterRaw[index];
							counterw[global->angularAvg_i[index]]++;
						}
					}
				}
			}
		}
	}
	
	if (global->useSolidAngleCorrection) {
		double solidAngle = global->pixelSize*global->pixelSize/(global->detectorZ/1000*global->detectorZ/1000);
		for (int i=0; i<global->angularAvg_nn; i++) {
			if (global->hitfinder.use || global->listfinder.use) {
				if (counterp[i]) global->powderAverage[i] /= counterp[i]*solidAngle;
			}
			if (global->icefinder.use) {
				if (counteri[i]) global->iceAverage[i] /= counteri[i]*solidAngle;
			}
			if (global->waterfinder.use) {
				if (counterw[i]) global->waterAverage[i] /= counterw[i]*solidAngle;
			}			
		}
	} else {
		for (int i=0; i<global->angularAvg_nn; i++) {
			if (global->hitfinder.use || global->listfinder.use) {
				if (counterp[i]) global->powderAverage[i] /= counterp[i];
			}
			if (global->icefinder.use) {
				if (counteri[i]) global->iceAverage[i] /= counteri[i];
			}
			if (global->waterfinder.use) {
				if (counterw[i]) global->waterAverage[i] /= counterw[i];
			}
		}
	}
	
	// free memory of local variables
	free(counterp);
	free(counteri);
	free(counterw);
	
}


void calculateAngularAvg(tThreadInfo *threadInfo, cGlobal *global) {
	
	// allocate local & threadInfo arrays
	unsigned *counter = (unsigned*) calloc(global->angularAvg_nn, sizeof(unsigned));
	if (threadInfo->angularAvg) {
		free(threadInfo->angularAvg);
	}
	threadInfo->angularAvg = (double*) calloc(global->angularAvg_nn, sizeof(double));
	
	// angular average for each |q|
	for (int i=0; i<global->pix_nn; i++) {
		if ( global->angularAvg_i[i] < global->angularAvg_nn && global->angularAvg_i[i] >= 0 && (!global->useBadPixelMask || global->badpixelmask[i]) ) {
			threadInfo->angularAvg[global->angularAvg_i[i]] += threadInfo->corrected_data[i];
			counter[global->angularAvg_i[i]]++;
		}
	}
	
	if (global->useSolidAngleCorrection) {
		for (int i=0; i<global->angularAvg_nn; i++) {
			if (counter[i]) threadInfo->angularAvg[i] /= counter[i]*threadInfo->solidAngle;
		}
	} else {
		for (int i=0; i<global->angularAvg_nn; i++) {
			if (counter[i]) threadInfo->angularAvg[i] /= counter[i];
		}
	}
	
	// free memory of local variables
	free(counter);
	
}


void saveAngularAvg(tThreadInfo *threadInfo, cGlobal *global) {
	
	// save angular average of hit with Q-scale (2D array)
	if (global->hitAngularAvg) {
		
		char	filename[1024];
		double *buffer = (double*) calloc(2*global->angularAvg_nn, sizeof(double));
		sprintf(filename,"%s-angavg.h5",threadInfo->eventname);
		
		for(long i=0; i<global->angularAvg_nn; i++) {
			//buffer[i] = (float) global->angularAvgQ[i];
			buffer[i] = threadInfo->angularAvgQ[i];
			buffer[global->angularAvg_nn+i] = threadInfo->angularAvg[i];
		}
		writeSimpleHDF5(filename, buffer, global->angularAvg_nn, 2, H5T_NATIVE_DOUBLE);
		
		free(buffer);
		
	}
	
}


void calculateIntensityAvg(tThreadInfo *threadInfo, cGlobal *global) {
	
	threadInfo->intensityAvg = 0;
	
	for (long i=0; i<global->pix_nn; i++)
		threadInfo->intensityAvg += threadInfo->corrected_data[i];
	
	threadInfo->intensityAvg /= global->pix_nn;
	
}


void saveIntensities(cGlobal *global) {
	
	if (global->useIntensityStatistics) {
		
		char	filename[1024];
		double *buffer = (double*) calloc(global->nIntensities, sizeof(double));
		printf("Saving average intensities to file\n");
		sprintf(filename,"r%04u-intensities.h5",global->runNumber);
		for(long i=0; i<global->nIntensities; i++)
			buffer[i] = global->intensities[i];
		writeSimpleHDF5(filename, buffer, global->nIntensities, 1, H5T_NATIVE_DOUBLE);
		free(buffer);
		
	}
	
}


void makeIntensityHistograms(cGlobal *global) {
	
	double deltaI = 0.1; // 0.1 ADU steps
	unsigned Ibins = (unsigned) ceil((global->Imax-global->Imin)/deltaI)+1;
	unsigned nHits = 0;
	
	if (Ibins > 1) {
		global->Ihist = (unsigned*) calloc(Ibins, sizeof(unsigned));
		global->IHhist = (unsigned*) calloc(Ibins, sizeof(unsigned));
		for (long i=0; i<global->nIntensities; i++) {
			global->Ihist[int(round((global->intensities[i]-global->Imin)/deltaI))]++;
			if (global->hits[i]) {
				global->IHhist[int(round((global->intensities[i]-global->Imin)/deltaI))]++;
				global->Imean += global->intensities[i];
				nHits++;
			}
		}
		global->Imean /= nHits;
		char	filename[1024];
		float *buffer = (float*) calloc(4*Ibins, sizeof(float));
		printf("Saving histograms of average intensities to file\n");
		printf("\tMean average intensity of hits: %f ADU\n",global->Imean);
		sprintf(filename,"r%04u-intensity_histograms.h5",global->runNumber);
		for(int i=0; i<Ibins; i++) {
			buffer[i] = (float) global->Ihist[i];
			buffer[Ibins+i] = (float) global->IHhist[i];
			buffer[2*Ibins+i] = (float) global->Imin+i*deltaI;
		}
		
		writeSimpleHDF5(filename, buffer, Ibins, 3, H5T_NATIVE_FLOAT);
		free(buffer);
		free(global->Ihist);
		free(global->IHhist);
		global->Ihist = NULL;
		global->IHhist = NULL;
	} else {
		// disable intensity if not enough bins
		global->useIntensityStatistics = 0;
	}
	
}


void savePixelIntensities(tThreadInfo *threadInfo, cGlobal *global) {
	
	// save pixel intensities of hit without Q-scale (1D array)
	if (global->usePixelStatistics) {
		
		char	filename[1024];
		double *buffer = (double*) calloc(global->nPixels, sizeof(double));
		sprintf(filename,"%s-pixels.h5",threadInfo->eventname);
		
		for(long i=0; i<global->nPixels; i++) {
			buffer[i] = (double) threadInfo->corrected_data[global->pixels[i]];
		}
		writeSimpleHDF5(filename, buffer, global->nPixels, 1, H5T_NATIVE_DOUBLE);
		
		free(buffer);
		
	}
	
}


void saveEnergies(cGlobal *global) {
	
	if (global->useEnergyCalibration) {
		
		char	filename[1024];
		float *buffer = (float*) calloc(2*global->nEnergies, sizeof(float));
		printf("Saving energies and wavelengths to file\n");
		sprintf(filename,"r%04u-energies.h5",global->runNumber);
		for(long i=0; i<global->nEnergies; i++)
			buffer[i] = (float) global->energies[i];
		for(long i=0; i<global->nEnergies; i++)
			buffer[global->nEnergies+i] = (float) global->wavelengths[i];
		writeSimpleHDF5(filename, buffer, global->nEnergies, 2, H5T_NATIVE_FLOAT);
		free(buffer);
		
	}
	
}


void makeEnergyHistograms(cGlobal *global) {

	double deltaE = 1; // 1 eV steps
	//double deltaL = 14e-6; // Equivalence of 1 eV steps at 9385 eV
	unsigned Ebins = (unsigned) ceil((global->Emax-global->Emin)/deltaE)+1;
	double deltaL = (global->Lmax-global->Lmin)/(Ebins-1);
	
	if (Ebins > 1) {
		global->Ehist = (unsigned*) calloc(Ebins, sizeof(unsigned));
		global->Lhist = (unsigned*) calloc(Ebins, sizeof(unsigned));
		for (long i=0; i<global->nEnergies; i++) {
			global->Ehist[int(round((global->energies[i]-global->Emin)/deltaE))]++;
			global->Lhist[int(round((global->wavelengths[i]-global->Lmin)/deltaL))]++;
			global->Emean += global->energies[i];
			global->Lmean += global->wavelengths[i];
		}
		global->Emean /= global->nEnergies;
		global->Lmean /= global->nEnergies;
		char	filename[1024];
		float *buffer = (float*) calloc(4*Ebins, sizeof(float));
		printf("Saving histograms of energies and wavelengths to file\n");
		printf("\tMean photon energy: %f eV\n",global->Emean);
		printf("\tMean wavelength: %f A\n",global->Lmean);
		sprintf(filename,"r%04u-energy_histograms.h5",global->runNumber);
		for(int i=0; i<Ebins; i++) {
			buffer[i] = (float) global->Ehist[i];
			buffer[Ebins+i] = (float) global->Emin+i*deltaE;
			buffer[2*Ebins+i] = (float) global->Lhist[i];
			buffer[3*Ebins+i] = (float) global->Lmin+i*deltaL;
		}
		
		writeSimpleHDF5(filename, buffer, Ebins, 4, H5T_NATIVE_FLOAT);
		free(buffer);
		free(global->Ehist);
		free(global->Lhist);
		global->Ehist = NULL;
		global->Lhist = NULL;
	} else {
		// disable energy calibration if not enough bins
		global->useEnergyCalibration = 0;
	}
	
}


void makePowderQcalibration(cGlobal *global) {
	
	// allocate arrays
	if (global->angularAvgQcal) {
		delete[] global->angularAvgQcal;
	}
	global->angularAvgQcal = new double[global->angularAvg_nn];
	
	// calculate Q-calibration (using the detector position from the last event)
	for (int i=0; i<global->angularAvg_nn; i++) {
		global->angularAvgQcal[i] = 4*M_PI*sin(atan(global->pixelSize*global->angularAvgQ[i]*1000/global->detectorZ)/2)/global->Lmean;
	}
	
	if (!global->powdersum || !global->powderAngularAvg) {
		// write data to buffer only if it is not saved to powder angular average
		double *buffer = (double*) calloc(global->angularAvg_nn, sizeof(double));
		char	filename[1024];
		printf("Saving Q-calibration to file\n");
		sprintf(filename,"r%04u-angavg_Q.h5",global->runNumber);
		for(int i=0; i<global->angularAvg_nn; i++) {
			buffer[i] = global->angularAvgQcal[i];
		}	
		// write buffer to HDF5
		writeSimpleHDF5(filename, buffer, global->angularAvg_nn, 1, H5T_NATIVE_DOUBLE);
		// free local arrays
		free(buffer);
	}
	
}


int makeQcalibration(tThreadInfo *threadInfo, cGlobal *global) {
	
	// sanity check
	if (threadInfo->detectorPosition > 60 && threadInfo->detectorPosition < 600 && threadInfo->wavelengthA == threadInfo->wavelengthA) {
		// allocate arrays
		if (threadInfo->angularAvgQ) {
			free(threadInfo->angularAvgQ);
		}
		threadInfo->angularAvgQ = (double*) calloc(global->angularAvg_nn, sizeof(double));
		// calculate Q-calibration (using the detector position and energy from the current event)
		for (int i=0; i<global->angularAvg_nn; i++) {
			threadInfo->angularAvgQ[i] = 4*M_PI*sin(atan(global->pixelSize*global->angularAvgQ[i]*1000/threadInfo->detectorPosition)/2)/threadInfo->wavelengthA;
		}		
		return 0;
	} else return 1;
	
}


void calculateCenterCorrection(cGlobal *global, double intensities[], double normalization) {
	
	// calculating center correction from an array of doubles with intensities in raw format (index number matches pix_x and pix_y)
	printf("calculating center correction...\n");
	
	// calculate size of arrays
	DEBUGL2_ONLY cout << "MinR: " << global->centerCorrectionMinR << ", MaxR: " << global->centerCorrectionMaxR << ", DeltaR: " << global->centerCorrectionDeltaR << endl;
	int nR = (int) ceil((global->centerCorrectionMaxR - global->centerCorrectionMinR)/global->centerCorrectionDeltaR) + 1;
	int nC = (int) ceil((2*global->centerCorrectionMaxC)/global->centerCorrectionDeltaC) + 1;
	DEBUGL2_ONLY cout << "nR: " << nR << ", nC: " << nC << endl;
	
	// allocate hough array as 3D array (in doubles) initialized to zero
	array3D<double> *hough = new array3D<double>(nR, nC, nC);
	
	// calculate hough array
	for (int n=0; n<global->pix_nn; n++) {
        if (intensities[n]/normalization > global->centerCorrectionThreshold) {
			float x = global->pix_x[n];
			float y = global->pix_y[n];
			for (int na=0; na<nC; na++) {
				double a = na*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
				for (int nb=0; nb<nC; nb++) {
					double b = nb*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
					double r = sqrt(((double) x-a)*((double) x-a) + ((double) y-b)*((double) y-b));
					if (r>global->centerCorrectionMinR && r<global->centerCorrectionMaxR) {
						int nr = (int) round((r-global->centerCorrectionMinR)/global->centerCorrectionDeltaR);
						hough->set(nr, na, nb, hough->get(nr, na, nb)+1);
					}
				}
			}
		}
	}
	
	// find maximum in hough array
	double max = 0;
	int imax, jmax, kmax;
	for(int i=0; i<nR; i++) {
        for(int j=0; j<nC; j++) {
			for(int k=0; k<nC; k++) {
				if(hough->get(i, j, k) > max) {
					max = hough->get(i, j, k);
					imax = i;
					jmax = j;
					kmax = k;
				}
			}
		}
	}
	DEBUGL2_ONLY cout << "imax: " << imax << ", jmax: " << jmax << ", kmax: " << kmax << ", houghmax: " << max << endl;
	
	// assign new center to global variables
	pthread_mutex_lock(&global->pixelcenter_mutex);
	global->pixelCenterX = (float) jmax*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
	global->pixelCenterY = (float) kmax*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
	cout << "\tCorrected center in x: " << global->pixelCenterX << endl;
	cout << "\tCorrected center in y: " << global->pixelCenterY << endl;
	pthread_mutex_unlock(&global->pixelcenter_mutex);
	
	// cleanup of allocated memory
	delete hough;
}

void calculateCenterCorrection(tThreadInfo *info, cGlobal *global, float intensities[], float normalization) {
	
	// calculating center correction from an array of floats with intensities in raw format (index number matches pix_x and pix_y)
	DEBUGL1_ONLY printf("calculating center correction...\n");
	
	// calculate size of arrays
	int nR = (int) ceil((global->centerCorrectionMaxR - global->centerCorrectionMinR)/global->centerCorrectionDeltaR) + 1;
	int nC = (int) ceil((2*global->centerCorrectionMaxC)/global->centerCorrectionDeltaC) + 1;
	
	// allocate hough array as 3D array (in doubles) initialized to zero
	array3D<double> *hough = new array3D<double>(nR, nC, nC);
	
	// calculate hough array
	for (int n=0; n<global->pix_nn; n++) {
        if (intensities[n]/normalization > global->centerCorrectionThreshold) {
			float x = global->pix_x[n];
			float y = global->pix_y[n];
			for (int na=0; na<nC; na++) {
				double a = na*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
				for (int nb=0; nb<nC; nb++) {
					double b = nb*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
					double r = sqrt(((double) x-a)*((double) x-a) + ((double) y-b)*((double) y-b));
					if (r>global->centerCorrectionMinR && r<global->centerCorrectionMaxR) {
						int nr = (int) round((r-global->centerCorrectionMinR)/global->centerCorrectionDeltaR);
						hough->set(nr, na, nb, hough->get(nr, na, nb)+1);
					}
				}
			}
		}
	}
	
	// find maximum in hough array
	double max = 0;
	int imax, jmax, kmax;
	for(int i=0; i<nR; i++) {
        for(int j=0; j<nC; j++) {
			for(int k=0; k<nC; k++) {
				if(hough->get(i, j, k) > max) {
					max = hough->get(i, j, k);
					imax = i;
					jmax = j;
					kmax = k;
				}
			}
		}
	}
	
	// assign new center to variables in threadInfo
	
	info->pixelCenterX = (float) jmax*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
	info->pixelCenterY = (float) kmax*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
	DEBUGL1_ONLY cout << "\tCorrected center in x: " << info->pixelCenterX << endl;
	DEBUGL1_ONLY cout << "\tCorrected center in y: " << info->pixelCenterY << endl;
	
	// cleanup of allocated memory
	delete hough;
}

void calculateCenterCorrectionQuad(cGlobal *global, double intensities[], double normalization) {
	
	// calculating center correction from an array of doubles with intensities in raw format (index number matches pix_x and pix_y)
	printf("calculating center correction for each quad...\n");
	
	// calculate size of arrays
	DEBUGL2_ONLY cout << "MinR: " << global->centerCorrectionMinR << ", MaxR: " << global->centerCorrectionMaxR << ", DeltaR: " << global->centerCorrectionDeltaR << endl;
	int nR = (int) ceil((global->centerCorrectionMaxR - global->centerCorrectionMinR)/global->centerCorrectionDeltaR) + 1;
	int nC = (int) ceil((2*global->centerCorrectionMaxC)/global->centerCorrectionDeltaC) + 1;
	DEBUGL2_ONLY cout << "nR: " << nR << ", nC: " << nC << endl;
	
	double radius[4];
	double dx[4];
	double dy[4];
	// loop over quads
	for (int quad=0; quad<4; quad++) {
		
		// declare hough array as 3D array (in doubles) initialized to zero
		array3D<double> hough = array3D<double>(nR, nC, nC);
		
		// loop over all pixels in each quad
		for(int mi=0; mi<2; mi++) { // mi decides what col in 2x8 matrix of raw data ASICs for each quad
			for(int mj=0; mj<8; mj++) { // mj decides what row in 2x8 matrix of raw data ASICs for each quad
				for(int i=0; i<ROWS; i++) {
					for(int j=0; j<COLS; j++) {
						
						long n = (j + mj*COLS) * (8*ROWS); // [row, 0] in 1552x1480 raw data format
						n += i + mi*ROWS + quad*2*ROWS; // [row, col] in 1552x1480 raw data format
						
						// calculate hough array
						if (intensities[n]/normalization > global->centerCorrectionThreshold) {
							float x = global->pix_x[n];
							float y = global->pix_y[n];
							for (int na=0; na<nC; na++) {
								double a = na*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
								for (int nb=0; nb<nC; nb++) {
									double b = nb*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
									double r = sqrt(((double) x-a)*((double) x-a) + ((double) y-b)*((double) y-b));
									if (r>global->centerCorrectionMinR && r<global->centerCorrectionMaxR) {
										int nr = (int) round((r-global->centerCorrectionMinR)/global->centerCorrectionDeltaR);
										hough.set(nr, na, nb, hough.get(nr, na, nb)+1);
									}
								}
							}
						}
					}
				}
			}
		}
		
		// find maximum in hough array
		double max = 0;
		int imax, jmax, kmax;
		for(int i=0; i<nR; i++) {
			for(int j=0; j<nC; j++) {
				for(int k=0; k<nC; k++) {
					if(hough.get(i, j, k) > max) {
						max = hough.get(i, j, k);
						imax = i;
						jmax = j;
						kmax = k;
					}
				}
			}
		}
		DEBUGL2_ONLY cout << "imax: " << imax << ", jmax: " << jmax << ", kmax: " << kmax << ", houghmax: " << max << endl;
		
		// assign new center to local variables
		radius[quad] = imax*global->centerCorrectionDeltaR + global->centerCorrectionMinR;
		dx[quad] = jmax*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
		dy[quad] = kmax*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
		
	}
	
	if (global->calculateCenterCorrectionQuad == 1) {
		// calculate quad shift for global pixel arrays with unified center		
		for (int quad=0; quad<4; quad++) {		
			global->quad_dx[quad] = (float) -dx[quad];
			global->quad_dy[quad] = (float) -dy[quad];
			DEBUGL1_ONLY cout << "\tQuad" << quad << " radius = " << radius[quad] << endl;
		}
		
	} else {
		// calculate max radius from hough transform
		double rmax = 0;
		for (int quad=0; quad<4; quad++) {
			if (radius[quad] > rmax)
				rmax = radius[quad];
		}
		
		// calculate quad shift for global pixel arrays with unified radius
		// Q0 (+X, +Y)
		global->quad_dx[0] = -dx[0] + (rmax - radius[0])/sqrt(2.0);
		global->quad_dy[0] = -dy[0] + (rmax - radius[0])/sqrt(2.0);
		// Q1 (-X, +Y)
		global->quad_dx[1] = -dx[1] - (rmax - radius[1])/sqrt(2.0);
		global->quad_dy[1] = -dy[1] + (rmax - radius[1])/sqrt(2.0);
		// Q2 (-X, -Y)
		global->quad_dx[2] = -dx[2] - (rmax - radius[2])/sqrt(2.0);
		global->quad_dy[2] = -dy[2] - (rmax - radius[2])/sqrt(2.0);
		// Q3 (+X, -Y)
		global->quad_dx[3] = -dx[3] + (rmax - radius[3])/sqrt(2.0);
		global->quad_dy[3] = -dy[3] - (rmax - radius[3])/sqrt(2.0);
	}
	
	// shift global pixel arrays to the correct position
	global->shiftQuads(global->pix_x, global->quad_dx, global->pix_y, global->quad_dy);
}

void shiftPixelCenter(cGlobal *global) {
	
	// shift center in pixel arrays
	pthread_mutex_lock(&global->pixelcenter_mutex);
	for (int i=0; i<global->pix_nn; i++) {
		global->pix_x[i] -= global->pixelCenterX;
		global->pix_y[i] -= global->pixelCenterY;
	}
	pthread_mutex_unlock(&global->pixelcenter_mutex);
	
}

void updatePixelArrays(cGlobal *global) {
	
	/*
	 *	Find bounds of image array and update global radius variables
	 */
	float	xmax = -1e9;
	float	xmin =  1e9;
	float	ymax = -1e9;
	float	ymin =  1e9;
	float	rmax = 0;
	pthread_mutex_lock(&global->pixelcenter_mutex);
	for (int i=0; i<global->pix_nn; i++) {
		if (global->pix_x[i] > xmax) xmax = global->pix_x[i];
		if (global->pix_x[i] < xmin) xmin = global->pix_x[i];
		if (global->pix_y[i] > ymax) ymax = global->pix_y[i];
		if (global->pix_y[i] < ymin) ymin = global->pix_y[i];
		global->pix_r[i] = sqrt(((double) global->pix_x[i])*global->pix_x[i] + ((double) global->pix_y[i])*global->pix_y[i]);		
		if (global->pix_r[i] > rmax) rmax = global->pix_r[i];
	}
	global->pix_xmax = xmax;
	global->pix_xmin = xmin;
	global->pix_ymax = ymax;
	global->pix_ymin = ymin;
	global->pix_rmax = rmax;
	
	/*
	 *	Update global angle variables
	 */
	if (global->usePolarizationCorrection || global->useCorrelation) {
		if (global->phi)
			delete[] global->phi;
		global->phi = new double[global->pix_nn];
		// calculate azimuthal angle (phi) for each pixel
		for (int i = 0; i < global->pix_nn; i++) {
			double phii = atan2(global->pix_x[i], global->pix_y[i]);
			if (phii < 0) { // make sure the angle is between 0 and 2PI
				phii += 2*M_PI;
			}
			// assign phi to each pixel
			global->phi[i] = phii;
		}
		
	}
	pthread_mutex_unlock(&global->pixelcenter_mutex);	
	
}

void updatePixelArrays(tThreadInfo *info, cGlobal *global) {
	// update pixel arrays
	
	pthread_mutex_lock(&global->pixelcenter_mutex);
	float deltaX = info->pixelCenterX - global->pixelCenterX;
	float deltaY = info->pixelCenterY - global->pixelCenterY;	
	for (int i=0; i<global->pix_nn; i++) {
		global->pix_x[i] -= deltaX;
		global->pix_y[i] -= deltaY;
		float rtemp = sqrt(global->pix_x[i]*global->pix_x[i] + global->pix_y[i]*global->pix_y[i]);
		if (rtemp > global->pix_rmax) global->pix_rmax = rtemp;
	}
	global->pix_xmax -= deltaX;
	global->pix_xmin -= deltaX;
	global->pix_ymax -= deltaY;
	global->pix_ymin -= deltaY;
	global->pixelCenterX = info->pixelCenterX;
	global->pixelCenterY = info->pixelCenterY;
	pthread_mutex_unlock(&global->pixelcenter_mutex);
	
}

void updateImageArrays(cGlobal *global) {
	// update image arrays after powders have been updated (used for powder sums)
	DEBUGL1_ONLY printf("Updating detector configuration:\n");
	
	int deltax = (int) round(global->pixelCenterX);
	int deltay = (int) round(global->pixelCenterY);
	float xmax = global->pix_xmax;
	float xmin = global->pix_xmin;
	float ymax = global->pix_ymax;
	float ymin = global->pix_ymin;
	
	fesetround(1);
	xmax = lrint(xmax);
	xmin = lrint(xmin);
	ymax = lrint(ymax);
	ymin = lrint(ymin);
	DEBUGL1_ONLY printf("\tImage bounds:\n");
	DEBUGL1_ONLY printf("\tx range %f to %f\n",xmin,xmax);
	DEBUGL1_ONLY printf("\ty range %f to %f\n",ymin,ymax);
	
	// calculate size of output image array
	float max = xmax;
	if(ymax > max) max = ymax;
	if(fabs(xmin) > max) max = fabs(xmin);
	if(fabs(ymin) > max) max = fabs(ymin);
	long image_nx = 2*(unsigned)max;
	long image_nn = image_nx*image_nx;	
	DEBUGL1_ONLY printf("\tUpdated image output array will be %i x %i\n",(int)image_nx,(int)image_nx);
		
	// update global variables
	pthread_mutex_lock(&global->image_mutex);
	global->image_nx = image_nx;
	global->image_nn = image_nn;
	
	// update image arrays from raw data
	if (global->generateDarkcal) {
		// Update assembled data
		pthread_mutex_lock(&global->powdersumassembled_mutex);
		assemble2Dimage(global->powderRaw, global->powderAssembled, global);
		pthread_mutex_unlock(&global->powdersumassembled_mutex);
	}
	
	if ((global->hitfinder.use || global->listfinder.use) && global->powdersum) {
		// Update assembled data
		pthread_mutex_lock(&global->powdersumassembled_mutex);
		assemble2Dimage(global->powderRaw, global->powderAssembled, global);
		pthread_mutex_unlock(&global->powdersumassembled_mutex);
	}
	
	if (global->icefinder.use && global->powdersum) {
		// Update assembled data	: ice
		pthread_mutex_lock(&global->icesumassembled_mutex);
		assemble2Dimage(global->iceRaw, global->iceAssembled, global);
		pthread_mutex_unlock(&global->icesumassembled_mutex);			
	}
	
	
	if (global->waterfinder.use && global->powdersum) {
		// Update assembled data	: water
		pthread_mutex_lock(&global->watersumassembled_mutex);
		assemble2Dimage(global->waterRaw, global->waterAssembled, global);
		pthread_mutex_unlock(&global->watersumassembled_mutex);
	}
	
	pthread_mutex_unlock(&global->image_mutex);
	
}

void updateImageArrays(cGlobal *global, cHit *hit) {
	// update image arrays before powders have been updated (used for each hit)
	DEBUGL1_ONLY printf("Updating detector configuration:\n");
	
	float xmax = global->pix_xmax;
	float xmin = global->pix_xmin;
	float ymax = global->pix_ymax;
	float ymin = global->pix_ymin;
	
	fesetround(1);
	xmax = lrint(xmax);
	xmin = lrint(xmin);
	ymax = lrint(ymax);
	ymin = lrint(ymin);
	DEBUGL1_ONLY printf("\tImage bounds:\n");
	DEBUGL1_ONLY printf("\tx range %f to %f\n",xmin,xmax);
	DEBUGL1_ONLY printf("\ty range %f to %f\n",ymin,ymax);
	
	// calculate size of output image array
	float max = xmax;
	if(ymax > max) max = ymax;
	if(fabs(xmin) > max) max = fabs(xmin);
	if(fabs(ymin) > max) max = fabs(ymin);
	long image_nx = 2*(unsigned)max;
	long image_nn = image_nx*image_nx;
	
	// reallocate assembled arrays if necessary
	if (image_nx > global->image_nx) {
		
		DEBUGL1_ONLY printf("\tUpdated image output array will be %i x %i\n",(int)image_nx,(int)image_nx);
		
		// for the assembled powder, we have to reallocate the already saved image
		long deltanx = image_nx-global->image_nx;
		double *buffer;
		
		if (global->generateDarkcal) {
			// Allocate temporary buffer
			buffer = (double*) calloc(image_nn, sizeof(double));
			// Update assembled data
			pthread_mutex_lock(&global->powdersumassembled_mutex);
			for(long i=0; i<global->image_nn; i++)
				buffer[image_nx*deltanx/2+deltanx/2+deltanx*int(i/global->image_nx)+i] += global->powderAssembled[i];
			// Free old memory and update pointers
			if (global->powderAssembled)
				free(global->powderAssembled);
			global->powderAssembled = buffer;
			buffer = NULL;
			pthread_mutex_unlock(&global->powdersumassembled_mutex);
		}
		
		if (hit->standard && global->powdersum){
			// Allocate temporary buffer
			buffer = (double*) calloc(image_nn, sizeof(double));
			// Update assembled data
			pthread_mutex_lock(&global->powdersumassembled_mutex);
			for(long i=0; i<global->image_nn; i++)
				buffer[image_nx*deltanx/2+deltanx/2+deltanx*int(i/global->image_nx)+i] = global->powderAssembled[i];
			// Free old memory and update pointers
			if (global->powderAssembled)
				free(global->powderAssembled);
			global->powderAssembled = buffer;
			buffer = NULL;
			pthread_mutex_unlock(&global->powdersumassembled_mutex);
		}
	
		if (hit->ice && global->powdersum) {
			// Allocate temporary buffer
			buffer = (double*) calloc(image_nn, sizeof(double));
			// Update assembled data	: ice
			pthread_mutex_lock(&global->icesumassembled_mutex);
			for(long i=0; i<global->image_nn; i++)
				buffer[image_nx*deltanx/2+deltanx/2+deltanx*int(i/global->image_nx)+i] = global->iceAssembled[i];
			// Free old memory and update pointers
			if (global->iceAssembled)
				free(global->iceAssembled);
			global->iceAssembled = buffer;
			buffer = NULL;
			pthread_mutex_unlock(&global->icesumassembled_mutex);			
		}


		if (hit->water && global->powdersum) {
			// Allocate temporary buffer
			buffer = (double*) calloc(image_nn, sizeof(double));
			// Update assembled data	: water
			pthread_mutex_lock(&global->watersumassembled_mutex);
			for(long i=0; i<global->image_nn; i++)
				buffer[image_nx*deltanx/2+deltanx/2+deltanx*int(i/global->image_nx)+i] = global->waterAssembled[i];
			// Free old memory and update pointers
			if (global->waterAssembled)
				free(global->waterAssembled);
			global->waterAssembled = buffer;
			buffer = NULL;
			pthread_mutex_unlock(&global->watersumassembled_mutex);
		}
		
		// all that needs to be done for assemble2DImage to work properly is to update image_nx/nn
		// OBS: THIS MUST BE DONE AFTER THE REALLOCATION OF THE ARRAYS FOR THE CODE TO WORK
		pthread_mutex_lock(&global->image_mutex);
		global->image_nx = image_nx;
		global->image_nn = image_nn;
		pthread_mutex_unlock(&global->image_mutex);
		
	} else DEBUGL1_ONLY cout << "\tSize of image output array is unchanged." << endl;
}

void updateAngularAvgArrays(cGlobal *global) {
	// update size of angular average arrays before they are filled
	
	/*
	 *	Update global angular average variables
	 */
	if ((global->powdersum && (global->powderAngularAvg || global->refineMetrology)) || global->hitAngularAvg) {
		if (global->angularAvgStopQ == 0 || global->angularAvgStopQ < global->angularAvgStartQ) {
			global->angularAvg_nn = (unsigned) ceil(global->pix_rmax/global->angularAvgDeltaQ)+1;
			global->angularAvgStartQ = 0;
		} else {
			global->angularAvg_nn = (unsigned) ceil((global->angularAvgStopQ-global->angularAvgStartQ)/global->angularAvgDeltaQ)+1;
		}
		if (global->hitfinder.use || global->listfinder.use) {
			if (global->powderAverage)
				free(global->powderAverage);
			global->powderAverage = (double*) calloc(global->angularAvg_nn, sizeof(double));
		}
		if (global->icefinder.use) {
			if (global->iceAverage)
				free(global->iceAverage);
			global->iceAverage = (double*) calloc(global->angularAvg_nn, sizeof(double));
		}
		if (global->waterfinder.use) {
			if (global->waterAverage)
				free(global->waterAverage);
			global->waterAverage = (double*) calloc(global->angularAvg_nn, sizeof(double));
		}
		if (global->angularAvgQ)
			delete[] global->angularAvgQ;
		global->angularAvgQ = new double[global->angularAvg_nn];
		for (int i=0; i<global->angularAvg_nn; i++) {
			global->angularAvgQ[i] = global->angularAvgStartQ + i*global->angularAvgDeltaQ;
		}
		// calculate index for each pixel with correct bin lengths in angular average array
		if (global->angularAvg_i)
			delete[] global->angularAvg_i;
		global->angularAvg_i = new int[global->pix_nn];
		for (int i=0; i<global->pix_nn; i++) {
			global->angularAvg_i[i] = (int) round( (global->pix_r[i] - global->angularAvgStartQ) / global->angularAvgDeltaQ );
		}	
	}	
}


void translateQuads(cGlobal *global) {
	
	float dx, dy;
	int refinementNumC = int(2*round(global->refinementMaxC/global->refinementDeltaC) + 1);
	double imax[4];
	double rmax[4];
	for (int i=0; i<4; i++) {
		imax[i] = 0;
		rmax[i] = 0;		
	}
	
	float *pix_x = global->pix_x; // let original pix_x be stored through local pointer
	float *pix_y = global->pix_y; // let original pix_y be stored through local pointer
	global->pix_x = (float *) calloc(global->pix_nn, sizeof(float));
	global->pix_y = (float *) calloc(global->pix_nn, sizeof(float));
	for (int i=0; i<global->pix_nn; i++) {
		global->pix_x[i] = pix_x[i];
		global->pix_y[i] = pix_y[i];
	}
	
	// for each quad
	for (int quad=0; quad<4; quad++) {
		
		// refine x
		for (int x=0; x<refinementNumC; x++) {
			dx = x*global->refinementDeltaC - global->refinementMaxC;
			DEBUGL2_ONLY cout << "dx: " << dx << endl;
			
			// refine y
			for (int y=0; y<refinementNumC; y++) {
				dy = y*global->refinementDeltaC - global->refinementMaxC; 
				DEBUGL2_ONLY cout << "\tdy: " << dy << endl;
				
				// shift global->pix_x and global->pix_y
				for(int mi=0; mi<2; mi++){ // mi decides what col in 2x8 matrix of raw data ASICs for each quad
					for(int mj=0; mj<8; mj++){ // mj decides what row in 2x8 matrix of raw data ASICs for each quad
						for(int i=0; i<ROWS; i++){
							for(int j=0; j<COLS; j++){
								long index = (j + mj*COLS) * (8*ROWS); // [row, 0] in 1552x1480 raw data format
								index += i + mi*ROWS + quad*2*ROWS; // [row, col] in 1552x1480 raw data format
								global->pix_x[index] += dx;
								global->pix_y[index] += dy;
							}
						}
					}
				}
				
				// update pixel distances, pixel indices, and recalculate angular average
				for(int i=0; i<global->pix_nn; i++) {
					global->pix_r[i] = sqrt(((double) global->pix_x[i])*global->pix_x[i] + ((double) global->pix_y[i])*global->pix_y[i]);
				}
				for (int i=0; i<global->pix_nn; i++) {
					global->angularAvg_i[i] = (int) round( (global->pix_r[i] - global->angularAvgStartQ) / global->angularAvgDeltaQ );
				}
				if (global->refineMetrology == 2) calculatePowderAngularAvg(global);
				else calculateQuadAngularAvg(global, quad);
				
				// evaluate maximum intensity
				double imaxtemp = 0;
				double rmaxtemp = 0;
				for (int i=0; i<global->angularAvg_nn; i++) {
					if (global->powderAverage[i] > imaxtemp) {
						imaxtemp = global->powderAverage[i];
						rmaxtemp = i*global->angularAvgDeltaQ + global->angularAvgStartQ;
					}
					global->powderAverage[i] = 0;
				}
				DEBUGL2_ONLY cout << "\tQ" << quad << ", Imax: " << imaxtemp<< ", rmax: " << rmaxtemp << endl;
				
				// update best dx/dy based on maximum intensity
				if (imaxtemp > imax[quad]) {
					imax[quad] = imaxtemp;
					rmax[quad] = rmaxtemp;
					global->quad_dx[quad] = dx;
					global->quad_dy[quad] = dy;
					cout << "\tQ" << quad << ": (" << dx << "," << dy << ")";
					DEBUGL1_ONLY cout << " at rmax = " << rmaxtemp << " pixels (Imax = " << imaxtemp << ")";
					cout << endl;
				}
				
				// restore global arrays
				for(int mi=0; mi<2; mi++){ // mi decides what col in 2x8 matrix of raw data ASICs for each quad
					for(int mj=0; mj<8; mj++){ // mj decides what row in 2x8 matrix of raw data ASICs for each quad
						for(int i=0; i<ROWS; i++){
							for(int j=0; j<COLS; j++){
								int index = (j + mj*COLS) * (8*ROWS);
								index += i + mi*ROWS + quad*2*ROWS;
								global->pix_x[index] = pix_x[index];
								global->pix_y[index] = pix_y[index];
							}
						}
					}
				}
				
			}
			
		}
		
		// shift global arrays by optimized position from the current quadrant to use for optimization of subsequent quadrants
		for(int mi=0; mi<2; mi++){ // mi decides what col in 2x8 matrix of raw data ASICs for each quad
			for(int mj=0; mj<8; mj++){ // mj decides what row in 2x8 matrix of raw data ASICs for each quad
				for(int i=0; i<ROWS; i++){
					for(int j=0; j<COLS; j++){
						int index = (j + mj*COLS) * (8*ROWS);
						index += i + mi*ROWS + quad*2*ROWS;
						global->pix_x[index] += global->quad_dx[quad];
						global->pix_y[index] += global->quad_dy[quad];
					}
				}
			}
		}
		
	}
	
	// clear global arrays
	if (global->pix_x)
		free(global->pix_x);
	if (global->pix_y)
		free(global->pix_y);
	
	// shift pix_x and pix_y to the correct position and update global arrays
	global->shiftQuads(pix_x, global->quad_dx, pix_y, global->quad_dy);
	global->pix_x = pix_x;
	global->pix_y = pix_y;
}

void rotateQuads(cGlobal *global) {
	DEBUGL1_ONLY cout << "rotateQuads" << endl;
}

void optimizeGap(cGlobal *global) {
	DEBUGL1_ONLY cout << "optimizeGap" << endl;
}
