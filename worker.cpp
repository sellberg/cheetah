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
#include "myana/XtcRun.hh"
#include "release/pdsdata/cspad/ElementHeader.hh"
#include "release/pdsdata/cspad/ElementIterator.hh"
#include "cspad-gjw/CspadTemp.hh"
#include "cspad-gjw/CspadGeometry.hh"

#include <string.h>
#include <cmath>
#include <hdf5.h>
#include <stdlib.h>
#include <iostream>
using std::cout;
using std::endl;

#include "worker.h"
#include "hitfinder.h"
#include "commonmode.h"
#include "background.h"
#include "correlation.h"

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
	threadInfo->raw_data = (uint16_t*) calloc(RAW_DATA_LENGTH,sizeof(uint16_t));
	for(int quadrant=0; quadrant<4; quadrant++) {
		long	i,j,ii;
		for(long k=0; k<2*ROWS*8*COLS; k++) {
			i = k % (2*ROWS) + quadrant*(2*ROWS);
			j = k / (2*ROWS);
			ii  = i+(8*ROWS)*j;
			threadInfo->raw_data[ii] = threadInfo->quad_data[quadrant][k];
		}
	}
	threadInfo->corrected_data = (float*) calloc(RAW_DATA_LENGTH,sizeof(float));
	for(long i=0;i<global->pix_nn;i++)
		threadInfo->corrected_data[i] = threadInfo->raw_data[i];
	
	
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
	}
	
	/*
	 *	Hitfinding - Water
	 */
	hit.water = 0;
	if(global->waterfinder.use){
		hit.water = hitfinder(threadInfo, global, &(global->waterfinder));
	}

	/*
	 *	Hitfinding - Ice
	 */
	hit.ice = 0;
	if(global->icefinder.use){
		hit.ice = hitfinder(threadInfo, global, &(global->icefinder));
	}
	 
	/*
	 *	Hitfinding - Background
	 */
	hit.background = 0;
	if(global->backgroundfinder.use){
		hit.background = hitfinder(threadInfo, global, &(global->backgroundfinder));
	}
	
	/*
	 *	Update the running background
	 */
	if (global->backgroundfinder.use){	
		updatePersistentBackground(threadInfo, global, hit.background);
	} else {
		updatePersistentBackground(threadInfo, global, hit.standard);
	}

	/*
	 *	Are we still in 'frame digesting' mode?
	 */
	if(threadInfo->threadNum < global->startFrames) {
		printf("r%04u:%i (%3.1fHz): Digesting initial frames\n", (int)global->runNumber, (int)threadInfo->threadNum, global->datarate);
		threadInfo->image = NULL;
		goto cleanup;
        //ATTENTION! goto should not be used at all ( see http://www.cplusplus.com/forum/general/29190/ )
        //should be replaced by a while loop with a break statement... by JF (2011/04/25)
	}
	
	
	/*
	 *	Calculate center correction
	 */
	if (global->calculateCenterCorrectionHit && (hit.standard || hit.water || hit.ice)) {
		calculateCenterCorrection(global, corrected_data);
		if (global->useCenterCorrection) {
			updatePixelArrays(global);
			updateImageArrays(global);
		}
	}
	
	
	/*
	 *	Apply attenuation correction
	 */
	if (global->useAttenuationCorrection > 0) {
		applyAttenuationCorrection(threadInfo, global);
	}
	
	
	/*
	 *	Assemble quadrants into a 'realistic' 2D image
	 */
	assemble2Dimage(threadInfo, global);
	
	
	/*
	 *	Add to powder if it's a hit or if we wish to generateDarkcal(member data of global)
	 */
	addToPowder(threadInfo, global, &hit);
	
	
	/*
     *  Perform cross-correlation analysis
     */
	if (global->useCorrelation && (global->hdf5dump || (hit.standard && global->hitfinder.savehits) 
								   || (hit.water && global->waterfinder.savehits) 
								   || (hit.ice && global->icefinder.savehits) 
								   || (!hit.background && global->backgroundfinder.savehits) )) {
		correlate(threadInfo, global);
	}
	
	
	/*
	 *	If this is a hit, write out to our favourite HDF5 format
	 */
	if(global->hdf5dump) 
		writeHDF5(threadInfo, global, threadInfo->eventname, global->hitfinder.cleanedfp);
	else {
		if(hit.standard && global->hitfinder.savehits)
			writeHDF5(threadInfo, global, threadInfo->eventname, global->hitfinder.cleanedfp);

		char eventname[1024];
		if(hit.water && global->waterfinder.savehits) {
			sprintf(eventname,"%s_waterhit",threadInfo->eventname);
			writeHDF5(threadInfo, global, eventname, global->waterfinder.cleanedfp);
		}
		if(hit.ice && global->icefinder.savehits) {
			sprintf(eventname,"%s_icehit",threadInfo->eventname);
			writeHDF5(threadInfo, global, eventname, global->icefinder.cleanedfp);
		}
		if(!hit.background && global->backgroundfinder.savehits) {
			sprintf(eventname,"%s_background",threadInfo->eventname);
			writeHDF5(threadInfo, global, eventname, global->backgroundfinder.cleanedfp);
		}
	}
	printf("r%04u:%i (%3.1fHz): Processed (npeaks=%i)\n", (int)global->runNumber, (int)threadInfo->threadNum,global->datarate, threadInfo->nPeaks);


	/*
	 *	Write out information on each frame to a log file
	 */
	pthread_mutex_lock(&global->framefp_mutex);
	fprintf(global->framefp, "%i, %i, %s, %i\n", (int)threadInfo->threadNum, threadInfo->seconds, threadInfo->eventname, threadInfo->nPeaks);
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
	free(threadInfo->raw_data);
	free(threadInfo->corrected_data);
	free(threadInfo->image);
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
        
        
		// Sum assembled data
		pthread_mutex_lock(&global->powdersumassembled_mutex);
		for(long i=0; i<global->image_nn; i++)
            global->powderAssembled[i] += threadInfo->image[i];
		pthread_mutex_unlock(&global->powdersumassembled_mutex);
	}

    
	if (hit->standard){
		// Sum raw format data
		if (global->powdersum && global->saveRaw) {
			pthread_mutex_lock(&global->powdersumraw_mutex);
			for(long i=0; i<global->pix_nn; i++)
				global->powderRaw[i] += threadInfo->corrected_data[i];
			pthread_mutex_unlock(&global->powdersumraw_mutex);
		}
		
		// Sum assembled data		
		if (global->powdersum) {
			pthread_mutex_lock(&global->powdersumassembled_mutex);
			global->npowder += 1;
			for(long i=0; i<global->image_nn; i++)
				if(threadInfo->image[i] > global->powderthresh)
					global->powderAssembled[i] += threadInfo->image[i];
			pthread_mutex_unlock(&global->powdersumassembled_mutex);
		}
	}

	if (hit->ice){
		// Sum raw format data 	: ice
		if (global->powdersum && global->saveRaw) {
			pthread_mutex_lock(&global->icesumraw_mutex);
			for(long i=0; i<global->pix_nn; i++)
				global->iceRaw[i] += threadInfo->corrected_data[i];
			pthread_mutex_unlock(&global->icesumraw_mutex);			
		}
	
		// Sum assembled data	: ice
		if (global->powdersum) {
			pthread_mutex_lock(&global->icesumassembled_mutex);
			global->nice += 1;
			for(long i=0; i<global->image_nn; i++)
				if(threadInfo->image[i] > global->powderthresh)
					global->iceAssembled[i] += threadInfo->image[i];
			pthread_mutex_unlock(&global->icesumassembled_mutex);			
		}
	}

	if (hit->water){
		// Sum raw format data  : water
		if (global->powdersum && global->saveRaw) {
			pthread_mutex_lock(&global->watersumraw_mutex);
			for(long i=0; i<global->pix_nn; i++)
				global->waterRaw[i] += threadInfo->corrected_data[i];
			pthread_mutex_unlock(&global->watersumraw_mutex);
		}
	
		// Sum assembled data	: water
		if (global->powdersum) {
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
 *	Interpolate raw (corrected) cspad data into a physical 2D image
 *	using pre-defined pixel mapping (as loaded from .h5 file)
 */
void assemble2Dimage(tThreadInfo *threadInfo, cGlobal *global){
	
	
	// Allocate temporary arrays for pixel interpolation (needs to be floating point)
	float	*data = (float*) calloc(global->image_nn,sizeof(float));
	float	*weight = (float*) calloc(global->image_nn,sizeof(float));
	for(long i=0; i<global->image_nn; i++){
		data[i] = 0;
		weight[i]= 0;
	}
	
	
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



void nameEvent(tThreadInfo *info, cGlobal *global){
	/*
	 *	Create filename based on date, time and fiducial for this image
	 */
//	char outfile[1024];
	char buffer1[80];
	char buffer2[80];	
	time_t eventTime = info->seconds;
	
	//setenv("TZ","US/Pacific",1);		// <--- Dangerous (not thread safe!)
	struct tm *timestatic, timelocal;
	timestatic=localtime_r( &eventTime, &timelocal );	
	strftime(buffer1,80,"%Y_%b%d",&timelocal);
	strftime(buffer2,80,"%H%M%S",&timelocal);
	sprintf(info->eventname,"LCLS_%s_r%04u_%s_%x_cspad",buffer1,global->runNumber,buffer2,info->fiducial);
}
	
	
/*
 *	Write out processed data to our 'standard' HDF5 format
 */
void writeHDF5(tThreadInfo *info, cGlobal *global, char *eventname, FILE* hitfp){
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
	//sprintf(outfile,"LCLS_%s_r%04u_%s_%x_cspad.h5",buffer1,global->runNumber,buffer2,info->fiducial);

	sprintf(outfile,"%s.h5",eventname);
	//strcpy(outfile, info->eventname);
	printf("r%04u:%i (%2.1f Hz): Writing data to: %s\n", (int)global->runNumber, (int)info->threadNum, global->datarate, outfile);

	pthread_mutex_lock(&global->framefp_mutex);
	fprintf(hitfp, "r%04u/%s, %i\n", (int)global->runNumber, info->eventname, info->nPeaks);
	pthread_mutex_unlock(&global->framefp_mutex);
	
		
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


void saveRunningSums(cGlobal *global) {
	char	filename[1024];

	/*
	 *	Compute and save darkcal
	 */
	if(global->generateDarkcal) {
		printf("Processing darkcal\n");
		sprintf(filename,"r%04u-darkcal.h5",global->runNumber);
		int16_t *buffer3 = (int16_t*) calloc(global->pix_nn, sizeof(int16_t));
		pthread_mutex_lock(&global->powdersumraw_mutex);
		for(long i=0; i<global->pix_nn; i++)
			buffer3[i] = (int16_t) (global->powderRaw[i]/(float)global->npowder);
		pthread_mutex_unlock(&global->powdersumraw_mutex);
		//for(long i=0; i<global->pix_nn; i++)
		//	if (buffer3[i] < 0) buffer3[i] = 0;
		printf("Saving darkcal to file\n");
		writeSimpleHDF5(filename, buffer3, (int)global->pix_nx, (int)global->pix_ny, H5T_STD_I16LE);	
		free(buffer3);
	}

	else {
		if(global->hitfinder.use) {
			if (global->powdersum) {
				
				/*
				 *	Save assembled powder pattern
				 */
				printf("Saving assembled sum data to file\n");
				sprintf(filename,"r%04u-AssembledSum.h5",global->runNumber);
				float *buffer2 = (float*) calloc(global->image_nn, sizeof(float));
				pthread_mutex_lock(&global->powdersumassembled_mutex);
				for(long i=0; i<global->image_nn; i++){
					buffer2[i] = (float) global->powderAssembled[i]/global->npowder;
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
					buffer1[i] = (float) global->powderRaw[i]/global->npowder;
				pthread_mutex_unlock(&global->powdersumraw_mutex);
				writeSimpleHDF5(filename, buffer1, (int)global->pix_nx, (int)global->pix_ny, H5T_NATIVE_FLOAT);	
				free(buffer1);
				
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
					buffer4[i] = (float) global->iceAssembled[i]/global->nice;
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
					buffer5[i] = (float) global->iceRaw[i]/global->nice;
				pthread_mutex_unlock(&global->icesumraw_mutex);
				writeSimpleHDF5(filename, buffer5, (int)global->pix_nx, (int)global->pix_ny, H5T_NATIVE_FLOAT);	
				free(buffer5);
				
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
					buffer6[i] = (float) global->waterAssembled[i]/global->nwater;
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
					buffer7[i] = (float) global->waterRaw[i]/global->nwater;
				pthread_mutex_unlock(&global->watersumraw_mutex);
				writeSimpleHDF5(filename, buffer7, (int)global->pix_nx, (int)global->pix_ny, H5T_NATIVE_FLOAT);	
				free(buffer7);				
				
			}
		}
		
		
	}

	
}


void calculatePowderSAXS(cGlobal *global) {
	
	// calculating SAXS average from the powder pattern 
	DEBUGL1_ONLY printf("calculating angular average from powder pattern...\n");
	
	// allocate local arrays
	double *pix_r = new double[global->pix_nn];
	
	// calculate |q| for each pixel and bin lengths with correct resolution
	for (int i=0; i<global->pix_nn; i++) {
		pix_r[i] = (double) round( sqrt(((double) global->pix_x[i])*global->pix_x[i] + ((double) global->pix_y[i])*global->pix_y[i]) / global->deltaqSAXS) * global->deltaqSAXS;
	}
	
	// angular average for each |q|
	DEBUGL1_ONLY printf("# of steps: %d\n",global->powder_nn);
	DEBUGL2_ONLY printf("average SAXS intensity:\n");
	
	for (int i=0; i<global->powder_nn; i++) {
		global->powderQ[i] = i*global->deltaqSAXS;
		int counterp = 0; // reset counter
		int counteri = 0; // reset counter
		int counterw = 0; // reset counter
		for (int j=0; j<global->pix_nn; j++) {
			if ( pix_r[j] == global->powderQ[i] ) {
				if (global->hitfinder.use) {
					global->powderAverage[i] += global->powderRaw[j];
					counterp++;
				}
				if (global->icefinder.use) {
					global->iceAverage[i] += global->iceRaw[j];
					counteri++;
				}
				if (global->waterfinder.use) {
					global->waterAverage[i] += global->waterRaw[j];
					counterw++;
				}
			}
		}
		if (global->hitfinder.use) global->powderAverage[i] /= counterp;
		if (global->icefinder.use) global->iceAverage[i] /= counteri;
		if (global->waterfinder.use) global->waterAverage[i] /= counterw;
		
		DEBUGL2_ONLY {
			if (global->hitfinder.use) cout << "Q: " << global->powderQ[i] << ",   \t# pixels: " << counterp << ",\tI(powder): " << global->powderAverage[i] << endl;
			if (global->icefinder.use) cout << "Q: " << global->powderQ[i] << ",   \t# pixels: " << counteri << ",\tI(water): " << global->iceAverage[i] << endl;
			if (global->waterfinder.use) cout << "Q: " << global->powderQ[i] << ",   \t# pixels: " << counterw << ",\tI(ice): " << global->waterAverage[i] << endl;
		}
	}
	
	// free memory of local variables
	delete[] pix_r;
	
}

void calculateCenterCorrection(cGlobal *global, double *intensities, double normalization) {
	
	// calculating center correction from an array of doubles with intensities in raw format (index number matches pix_x and pix_y)
	printf("calculating center correction...\n");
	
	// calculate size of arrays
	int nR = (int) ceil((global->centerCorrectionMaxR-global->centerCorrectionMaxR)/global->centerCorrectionDeltaR) + 1;
	int nC = (int) ceil((2*global->centerCorrectionMaxC)/global->centerCorrectionDeltaC) + 1;
	
	// allocate hough array and initialize to zero
	int *hough = new int[nR][nC][nC];
	for (int i=0; i<nR; i++) {
		for (int j=0; j<nC; j++) {
			for (int k=0; k<nC; k++) {
				hough[i][j][k] = 0;
			}
		}
	}
	
	// calculate hough array
	for (int n=0; n<global->pix_nn; n++) {
        if (intensities[n]/normalization > global->centerCorrectionThreshold) {
			float x = global->pix_x[n];
			float y = global->pix_y[n];
			for (int na=0; a<nC; na++) {
				double a = na*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
				for (int nb=0; nb<nC; nb++) {
					double b = nb*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
					double r = sqrt(((double) x-a)*((double) x-a) + ((double) y-b)*((double) y-b));
					if (r>global->centerCorrectionMinR && r<global->centerCorrectionMaxR) {
						int nr = (int) round((r-global->centerCorrectionMinR)/global->centerCorrectionDeltaR);
						hough[nr][na][nb] += 1;
					}
				}
			}
		}
	}
	
	// find maximum in hough array
	int max = 0;
	int imax, jmax, kmax;
	for(int i=0; i<nR; i++) {
        for(int j=0; j<nC; j++) {
			for(int k=0; k<nC; k++) {
				if(hough[i][j][k] > max) {
					max = hough[i][j][k];
					imax = i;
					jmax = j;
					kmax = k;
				}
			}
		}
	}
	
	// assign new center to global variables
	global->pixelCenterX = (float) jmax*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
	global->pixelCenterY = (float) kmax*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
	
	cout << "\tCorrected center in x: " << global->pixelCenterX << endl;
	cout << "\tCorrected center in y: " << global->pixelCenterY << endl;
	
	// cleanup of allocated memory
	delete[] hough;
}

void calculateCenterCorrection(cGlobal *global, float *intensities) {
	
	// calculating center correction from an array of floats with intensities in raw format (index number matches pix_x and pix_y)
	DEBUGL1_ONLY printf("calculating center correction...\n");
	
	// calculate size of arrays
	int nR = (int) ceil((global->centerCorrectionMaxR-global->centerCorrectionMaxR)/global->centerCorrectionDeltaR) + 1;
	int nC = (int) ceil((2*global->centerCorrectionMaxC)/global->centerCorrectionDeltaC) + 1;
	
	// allocate hough array and initialize to zero
	int *hough = new int[nR][nC][nC];
	for (int i=0; i<nR; i++) {
		for (int j=0; j<nC; j++) {
			for (int k=0; k<nC; k++) {
				hough[i][j][k] = 0;
			}
		}
	}
	
	// calculate hough array
	for (int n=0; n<global->pix_nn; n++) {
        if (intensities[n] > global->centerCorrectionThreshold) {
			float x = global->pix_x[n];
			float y = global->pix_y[n];
			for (int na=0; a<nC; na++) {
				double a = na*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
				for (int nb=0; nb<nC; nb++) {
					double b = nb*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
					double r = sqrt(((double) x-a)*((double) x-a) + ((double) y-b)*((double) y-b));
					if (r>global->centerCorrectionMinR && r<global->centerCorrectionMaxR) {
						int nr = (int) round((r-global->centerCorrectionMinR)/global->centerCorrectionDeltaR);
						hough[nr][na][nb] += 1;
					}
				}
			}
		}
	}
	
	// find maximum in hough array
	int max = 0;
	int imax, jmax, kmax;
	for(int i=0; i<nR; i++) {
        for(int j=0; j<nC; j++) {
			for(int k=0; k<nC; k++) {
				if(hough[i][j][k] > max) {
					max = hough[i][j][k];
					imax = i;
					jmax = j;
					kmax = k;
				}
			}
		}
	}
	
	// assign new center to global variables
	// MUST BE MADE THREAD SAFE WITH ITS OWN MUTEX LOCK/UNLOCK
	global->pixelCenterX = (float) jmax*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
	global->pixelCenterY = (float) kmax*global->centerCorrectionDeltaC - global->centerCorrectionMaxC;
	
	DEBUGL1_ONLY cout << "\tCorrected center in x: " << global->pixelCenterX << endl;
	DEBUGL1_ONLY cout << "\tCorrected center in y: " << global->pixelCenterY << endl;
	
	// cleanup of allocated memory
	delete[] hough;
}

void updatePixelArrays(cGlobal *global) {
	// update pixel arrays
	
	//MUTEX LOCK
	for (int i=0; i<global->pix_nn; i++) {
		global->pix_x[i] -= global->pixelCenterX;
		global->pix_y[i] -= global->pixelCenterY;
		float rtemp = sqrt(pix_x[i]*pix_x[i] + pix_y[i]*pix_y[i]);
		if (rtemp > global->pix_rmax) global->rmax = rtemp;
	}
	global->pix_xmax -= global->pixelCenterX;
	global->pix_xmin -= global->pixelCenterX;
	global->pix_ymax -= global->pixelCenterY;
	global->pix_ymin -= global->pixelCenterY;
	//MUTEX UNLOCK
	
}

void updateImageArrays(cGlobal *global) {
	// update size of image array
	DEBUGL1_ONLY printf("Updating detector configuration:\n");
	
	float xmax = global->pix_xmax;
	float xmin = global->pix_xmin;
	float ymax = global->pix_ymax;
	float ymin = global->pix_ymin;
	
	// NEED HEADER FOR lrint???
	fesetround(1);
	xmax = lrint(xmax);
	xmin = lrint(xmin);
	ymax = lrint(ymax);
	ymin = lrint(ymin);
	DEBUGL1_ONLY printf("\tImage bounds:\n");
	DEBUGL1_ONLY printf("\tx range %f to %f\n",xmin,xmax);
	DEBUGL1_ONLY printf("\ty range %f to %f\n",ymin,ymax);
	
	// How big must the output image be?
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
		int deltanx = image_nx-global->image_nx;
		double *buffer;
		
		if (global->generateDarkcal){
			// Allocate temporary buffer
			buffer = (double*) calloc(image_nn, sizeof(double));
			// Update assembled data
			pthread_mutex_lock(&global->powdersumassembled_mutex);
			for(long i=0; i<global->image_nn; i++)
				buffer[image_nx*deltanx/2+deltanx/2+deltanx*floor(i/global->image_nx)+i] += global->powderAssembled[i];
			// Free old memory and update pointers
			free(global->powderAssembled);
			global->powderAssembled = buffer;
			buffer = NULL;
			DEBUGL2_ONLY cout << "last index of new array: " << global->powderAssembled[image_nx-1] << endl;
			pthread_mutex_unlock(&global->powdersumassembled_mutex);
		}
		
		if (hit->standard && global->powdersum){
			// Allocate temporary buffer
			buffer = (double*) calloc(image_nn, sizeof(double));
			// Update assembled data
			pthread_mutex_lock(&global->powdersumassembled_mutex);
			for(long i=0; i<global->image_nn; i++)
				buffer[image_nx*deltanx/2+deltanx/2+deltanx*floor(i/global->image_nx)+i] = global->powderAssembled[i];
			// Free old memory and update pointers
			free(global->powderAssembled);
			global->powderAssembled = buffer;
			buffer = NULL;
			DEBUGL2_ONLY cout << "last index of new array: " << global->powderAssembled[image_nx-1] << endl;
			pthread_mutex_unlock(&global->powdersumassembled_mutex);
		}
	
		if (hit->ice && global->powdersum){
			// Allocate temporary buffer
			buffer = (double*) calloc(image_nn, sizeof(double));
			// Update assembled data	: ice
			pthread_mutex_lock(&global->icesumassembled_mutex);
			for(long i=0; i<global->image_nn; i++)
				buffer[image_nx*deltanx/2+deltanx/2+deltanx*floor(i/global->image_nx)+i] = global->iceAssembled[i];
			// Free old memory and update pointers
			free(global->iceAssembled);
			global->iceAssembled = buffer;
			buffer = NULL;
			DEBUGL2_ONLY cout << "last index of new array: " << global->iceAssembled[image_nx-1] << endl;
			pthread_mutex_unlock(&global->icesumassembled_mutex);			
		}


		if (hit->water && global->powdersum){
			// Allocate temporary buffer
			buffer = (double*) calloc(image_nn, sizeof(double));
			// Update assembled data	: water
			pthread_mutex_lock(&global->watersumassembled_mutex);
			for(long i=0; i<global->image_nn; i++)
				buffer[image_nx*deltanx/2+deltanx/2+deltanx*floor(i/global->image_nx)+i] = global->waterAssembled[i];
			// Free old memory and update pointers
			free(global->waterAssembled);
			global->waterAssembled = buffer;
			buffer = NULL;
			DEBUGL2_ONLY cout << "last index of new array: " << global->waterAssembled[image_nx-1] << endl;
			pthread_mutex_unlock(&global->watersumassembled_mutex);
		}
		
		// all that needs to be done for assemble2DImage to work properly is to update image_nx/nn
		// OBS: THIS MUST BE DONE AFTER THE REALLOCATION OF THE ARRAYS FOR THE CODE TO WORK
		// MUTEX LOCK
		global->image_nx = image_nx;
		global->image_nn = image_nn;
		// MUTEX UNLOCK
	} else DEBUGL1_ONLY cout << "\tSize of image output array is unchanged." << endl;
}