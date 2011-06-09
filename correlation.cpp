/*
 *  correlation.cpp
 *  cheetah
 *
 *  Created by Jan Feldkamp on 25/April/11.
 *  Copyright 2011 SLAC. All rights reserved.
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


//#include <stdio.h>
//#include <string.h>
//#include <pthread.h>
//#include <stdlib.h>

#include <pthread.h>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <cmath>

#include "setup.h"
#include "worker.h"
#include "correlation.h"
#include "crosscorrelator.h"





//----------------------------------------------------------correlate
// apply angular cross correlation
//-------------------------------------------------------------------
void correlate(tThreadInfo *threadInfo, cGlobal *global) {
    
    DEBUGL1_ONLY cout << "CORRELATING... in thread #" << threadInfo->threadNum << "." << endl;
	
    //create cross correlator object that takes care of the computations
	CrossCorrelator *cc = new CrossCorrelator( threadInfo->corrected_data, global->pix_x, global->pix_y, RAW_DATA_LENGTH, global->correlationNumQ, global->correlationNumPhi );

    DEBUGL1_ONLY cc->setDebug(1);                           //turn on debug level inside the CrossCorrelator, if needed
    DEBUGL2_ONLY cc->setDebug(2);                           //turn on debug level inside the CrossCorrelator, if needed

	//--------------------------------------------------------------------------------------------alg1
	if (global->useCorrelation == 1) {					
		
		DEBUGL1_ONLY cout << "XCCA regular" << endl;
		
		cc->calculatePolarCoordinates(global->correlationStartQ, global->correlationStopQ);
		cc->calculateSAXS();
		cc->calculateXCCA();	
		
		//writeSAXS(threadInfo, global, cc, threadInfo->eventname); // writes SAXS only to binary
		writeXCCA(threadInfo, global, cc, threadInfo->eventname); // writes XCCA+SAXS to binary
		
		
	//--------------------------------------------------------------------------------------------alg2
	} else if (global->useCorrelation == 2) {					
		DEBUGL2_ONLY cout << "XCCA fast" << endl;

        array2D *polar = new array2D;
        array2D *corr = new array2D;
        
        //cc->createLookupTable( 100, 100 );		// lookup table should in general not be created here (i.e., shot-by-shot)
		cc->setLookupTable( global->correlationLUT, global->correlationLUTdim1, global->correlationLUTdim2 );

		//transform data to polar coordinates as determined by the cheetah ini file	(in detector pixels)
		//to the q-calibrated values the cross-correlator expects
        double start_q = global->correlationStartQ;
        double stop_q = global->correlationStopQ;
        int number_q = global->correlationNumQ;
        double start_phi = global->correlationStartPhi;
        double stop_phi = global->correlationStopPhi;
        int number_phi = global->correlationNumPhi;
		
		//cout << "thread #" << threadInfo->threadNum << " locking mutex" << endl;
		pthread_mutex_lock(&global->correlation_mutex);	
		
		try{
			//calculate polar coordinates: returns an array2D in the first polar argument
			int polar_fail = cc->calculatePolarCoordinates_FAST(polar, number_q, start_q, stop_q, number_phi, start_phi, stop_phi);
			if (polar_fail){
				cout << "ERROR in correlate! Could not calculate polar coordinates!" << endl;
			}

			//perform the correlation
			int XCCA_fail = cc->calculateXCCA_FAST( polar, corr );
			if (XCCA_fail){
				cout << "ERROR in correlate! Could not calculate XCCA!" << endl;
			}
			
			//write output
			//writeXCCA(threadInfo, global, cc, threadInfo->eventname);        //jas: writeXCCA doesn't currently work for alg 2 because it does not use array3D *crossCorrelation to store data.
		} catch(...) {
			cerr << "Exception caught in correlate! Thread #" << threadInfo->threadNum << endl;
			cerr << "Aborting..." << endl;
			delete polar;
			delete corr;
			delete cc;
			exit(1);
		}
				
		//cout << "thread #" << threadInfo->threadNum << " unlocking mutex" << endl;
		pthread_mutex_unlock(&global->correlation_mutex);
		
		//clean up
        delete polar;
        delete corr;
		
		DEBUGL2_ONLY cout << "XCCA done." << endl;

	//--------------------------------------------------------------------------------------------default case
	} else {
		cerr << "Error in correlate()! Correlation algorithm " << global->useCorrelation << " not known." << endl;
	}
	
	delete cc;
}



//----------------------------------------------------------------------------writeSAXS
// write SAXS intensity to binary
//-------------------------------------------------------------------------------------
void writeSAXS(tThreadInfo *info, cGlobal *global, CrossCorrelator *cc, char *eventname) {	
	
	DEBUGL1_ONLY printf("writing SAXS to file...\n");
	FILE *filePointerWrite;
	char outfile[1024];
	double samplingLengthD = (double) cc->samplingLength(); // save everything as doubles
	double *buffer;
	buffer = (double*) calloc(cc->samplingLength(), sizeof(double));
	
	sprintf(outfile,"%s-saxs.bin",eventname);
	printf("r%04u:%i (%2.1f Hz): Writing data to: %s\n", (int)global->runNumber, (int)info->threadNum, global->datarate, outfile);
	
	filePointerWrite = fopen(outfile,"w+");

	// angular averages
	for (int i=0; i<cc->samplingLength(); i++) {
		buffer[i] = cc->getIave(i);
	}
	fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite); // saving dimensions of array before the actual data
	fwrite(&buffer[0],sizeof(double),cc->samplingLength(),filePointerWrite);
	
	// q binning
	for (int i=0; i<cc->samplingLength(); i++) {
		buffer[i] = cc->getQave(i);
	}
	fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
	fwrite(&buffer[0],sizeof(double),cc->samplingLength(),filePointerWrite);
	
	fclose(filePointerWrite);
	free(buffer);
	
	DEBUGL1_ONLY cout << "writeSAXS done" << endl;
}



//----------------------------------------------------------------------------writeXCCA
// write cross-correlation to binary
//-------------------------------------------------------------------------------------
void writeXCCA(tThreadInfo *info, cGlobal *global, CrossCorrelator *cc, char *eventname) {
	
	DEBUGL1_ONLY cout << "writing XCCA to file..." << endl;
	FILE *filePointerWrite;
	char outfile[1024];
	double samplingLengthD = (double) cc->samplingLength(); // save everything as doubles
	double samplingLagD = (double) cc->samplingLag();
	double samplingAngleD = (double) cc->samplingAngle();
	double *buffer;
	buffer = (double*) calloc(cc->samplingLength()*cc->samplingLength()*cc->samplingLag(), sizeof(double));
	info->correlation = (double*) calloc(global->correlation_nn, sizeof(double));
	
	if (global->autoCorrelationOnly)
		sprintf(outfile,"%s-xaca.bin",eventname);
	else {
		sprintf(outfile,"%s-xcca.bin",eventname);
	}
	printf("r%04u:%i (%2.1f Hz): Writing data to: %s\n", (int)global->runNumber, (int)info->threadNum, global->datarate, outfile);
	
	filePointerWrite = fopen(outfile,"w+");
	
	// angular averages
	for (int i=0; i<cc->samplingLength(); i++) {
		buffer[i] = cc->getIave(i);
	}
	fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite); // saving dimensions of array before the actual data
	fwrite(&buffer[0],sizeof(double),cc->samplingLength(),filePointerWrite);
	
	// q binning
	for (int i=0; i<cc->samplingLength(); i++) {
		buffer[i] = cc->getQave(i);
	}
	fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
	fwrite(&buffer[0],sizeof(double),cc->samplingLength(),filePointerWrite);
	
	// angle binning
	for (int i=0; i<cc->samplingAngle(); i++) {
		buffer[i] = cc->getPhiave(i);
	}
	fwrite(&samplingAngleD,sizeof(double),1,filePointerWrite);
	fwrite(&buffer[0],sizeof(double),cc->samplingAngle(),filePointerWrite);
	
	// cross-correlation
	if (global->autoCorrelationOnly) {

		// autocorrelation only (q1=q2)
		for (int i=0; i<cc->samplingLength(); i++) {
			for (int k=0; k<cc->samplingLag(); k++) {
				info->correlation[i*cc->samplingLag()+k] = cc->getCrossCorrelation(i,i,k);
				
			}
		}
		fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
		fwrite(&samplingLagD,sizeof(double),1,filePointerWrite);
		
	} else {
		
		// full version
		for (int i=0; i<cc->samplingLength(); i++) {
			for (int j=0; j<cc->samplingLength(); j++) {
				for (int k=0; k<cc->samplingLag(); k++) {
					info->correlation[i*cc->samplingLength()*cc->samplingLag()+j*cc->samplingLag()+k] = cc->getCrossCorrelation(i,j,k);
				}
			}
		}
		fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
		fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
		fwrite(&samplingLagD,sizeof(double),1,filePointerWrite);
		
	}
	fwrite(&info->correlation[0],sizeof(double),global->correlation_nn,filePointerWrite);
	
	fclose(filePointerWrite);
	free(buffer);
	
	DEBUGL1_ONLY cout << "writeXCCA done" << endl;
}

