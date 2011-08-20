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

#include <sstream>

#include <cmath>

#include "setup.h"
#include "worker.h"
#include "correlation.h"
#include "crosscorrelator.h"

#include "arrayclasses.h"
#include "arraydataIO.h"





//----------------------------------------------------------correlate
// apply angular cross correlation
//-------------------------------------------------------------------
void correlate(tThreadInfo *threadInfo, cGlobal *global) {
    
    DEBUGL1_ONLY cout << "CORRELATING... in thread #" << threadInfo->threadNum << "." << endl;
	
    //create cross correlator object that takes care of the computations
	//the arguments that are passed to the constructor determine 2D/3D calculations with/without mask
	CrossCorrelator *cc = NULL;
	if (global->autoCorrelateOnly) {
		if (global->useBadPixelMask) 
			cc = new CrossCorrelator( //auto-correlation 2D case, with mask
						threadInfo->corrected_data, global->pix_x, global->pix_y, RAW_DATA_LENGTH, 
						global->correlationNumPhi, global->correlationNumQ, 0, global->badpixelmask );
		else 
			cc = new CrossCorrelator( //auto-correlation 2D case, no mask
						threadInfo->corrected_data, global->pix_x, global->pix_y, RAW_DATA_LENGTH, 
						global->correlationNumPhi, global->correlationNumQ );
	} else {
		if (global->useBadPixelMask) 
			cc = new CrossCorrelator( //full cross-correlation 3D case, with mask
						threadInfo->corrected_data, global->pix_x, global->pix_y, RAW_DATA_LENGTH, 
						global->correlationNumPhi, global->correlationNumQ, global->correlationNumQ, global->badpixelmask );
		else 
			cc = new CrossCorrelator( //full cross-correlation 3D case, no mask
						threadInfo->corrected_data, global->pix_x, global->pix_y, RAW_DATA_LENGTH, 
						global->correlationNumQ, global->correlationNumQ, global->correlationNumPhi );	
	}
	
	//turn on debug level inside the CrossCorrelator, if needed
    DEBUGL1_ONLY cc->setDebug(1); 
	DEBUGL2_ONLY cc->setDebug(2);
	
	//--------------------------------------------------------------------------------------------alg1
	if (global->useCorrelation == 1) {					
		
		DEBUGL1_ONLY cout << "XCCA regular" << endl;
		
		cc->calculatePolarCoordinates(global->correlationStartQ, global->correlationStopQ);
		cc->calculateSAXS();
		cc->calculateXCCA();
		
		//writeSAXS(threadInfo, global, cc, threadInfo->eventname); // writes SAXS only to binary
		writeXCCA(threadInfo, global, cc, threadInfo->eventname); // writes XCCA/XACA+SAXS to binary
		
		
	//--------------------------------------------------------------------------------------------alg2
	} else if (global->useCorrelation == 2) {					
		DEBUGL2_ONLY cout << "XCCA fast" << endl;

        //cc->createLookupTable( 100, 100 );		// lookup table should in general not be created here (i.e., shot-by-shot)
		cc->setLookupTable( global->correlationLUT, global->correlationLUTdim1, global->correlationLUTdim2 );

		//transform data to polar coordinates as determined by the cheetah ini file	(in detector pixels)
		//to the q-calibrated values the cross-correlator expects
		int polar_fail = cc->calculatePolarCoordinates_FAST(global->correlationStartQ, global->correlationStopQ);
		if (polar_fail){
			cout << "ERROR in correlate! Could not calculate polar coordinates!" << endl;
		}

		//perform the correlation
		pthread_mutex_lock(&global->correlationFFT_mutex);			// need to protect the FFTW at the core, not thread-safe!!
		int XCCA_fail = cc->calculateXCCA_FAST();
		if (XCCA_fail){
			cout << "ERROR in correlate! Could not calculate XCCA!" << endl;
		}
		pthread_mutex_unlock(&global->correlationFFT_mutex);
		
		//write output
		pthread_mutex_lock(&global->correlation_mutex);
		writeXCCA(threadInfo, global, cc, threadInfo->eventname); // writes XCCA+SAXS to binary
		pthread_mutex_unlock(&global->correlation_mutex);
		
		//writing Tiff is not explicitly thread-safe... is this a problem for the cheetah?
		std::ostringstream name_osst;
		name_osst << threadInfo->eventname << "-xaca.tif";
		arraydataIO *io = new arraydataIO;
		io->writeToTiff( name_osst.str(), cc->corr(), 1 );            //dump scaled output from correlation
		delete io;
		
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
	
	DEBUGL1_ONLY cout << "writing SAXS to file..." << std::flush;
	FILE *filePointerWrite;
	char outfile[1024];
	double nQD = (double) cc->nQ(); // save everything as doubles
	double *buffer;
	buffer = (double*) calloc(cc->nQ(), sizeof(double));
	
	sprintf(outfile,"%s-saxs.bin",eventname);
	DEBUGL1_ONLY printf("r%04u:%i (%2.1f Hz): Writing data to: %s\n", (int)global->runNumber, (int)info->threadNum, global->datarate, outfile);
	
	filePointerWrite = fopen(outfile,"w+");

	// angular averages
	for (int i=0; i<cc->nQ(); i++) {
		buffer[i] = cc->getIave(i);
	}
	fwrite(&nQD,sizeof(double),1,filePointerWrite); // saving dimensions of array before the actual data
	fwrite(&buffer[0],sizeof(double),cc->nQ(),filePointerWrite);
	
	// q binning
	for (int i=0; i<cc->nQ(); i++) {
		buffer[i] = cc->getQave(i);
	}
	fwrite(&nQD,sizeof(double),1,filePointerWrite);
	fwrite(&buffer[0],sizeof(double),cc->nQ(),filePointerWrite);
	
	fclose(filePointerWrite);
	free(buffer);
	
	DEBUGL1_ONLY cout << "writeSAXS done" << endl;
}



//----------------------------------------------------------------------------writeXCCA
// write cross-correlation to binary
//-------------------------------------------------------------------------------------
void writeXCCA(tThreadInfo *info, cGlobal *global, CrossCorrelator *cc, char *eventname) {
	
	DEBUGL1_ONLY cout << "writing XCCA to file..." << std::flush;
	FILE *filePointerWrite;
	char outfile[1024];
	double nQD = (double) cc->nQ(); // save everything as doubles
	double nLagD = (double) cc->nLag();
	double nPhiD = (double) cc->nPhi();
	double *buffer;
	buffer = (double*) calloc(cc->nQ()*cc->nQ()*cc->nLag(), sizeof(double));
	if (global->sumCorrelation) info->correlation = (double*) calloc(global->correlation_nn, sizeof(double));
	
	if (global->autoCorrelateOnly){
		sprintf(outfile,"%s-xaca.bin",eventname);
	} else {
		sprintf(outfile,"%s-xcca.bin",eventname);
	}
	DEBUGL1_ONLY printf("r%04u:%i (%2.1f Hz): Writing data to: %s\n", (int)global->runNumber, (int)info->threadNum, global->datarate, outfile);
	
	filePointerWrite = fopen(outfile,"w+");
	
	// angular averages
	for (int i=0; i<cc->nQ(); i++) {
		buffer[i] = cc->getIave(i);
	}
	fwrite(&nQD,sizeof(double),1,filePointerWrite); // saving dimensions of array before the actual data
	fwrite(&buffer[0],sizeof(double),cc->nQ(),filePointerWrite);
	
	// q binning
	for (int i=0; i<cc->nQ(); i++) {
		buffer[i] = cc->getQave(i);
	}
	fwrite(&nQD,sizeof(double),1,filePointerWrite);
	fwrite(&buffer[0],sizeof(double),cc->nQ(),filePointerWrite);
	
	// angle binning
	for (int i=0; i<cc->nPhi(); i++) {
		buffer[i] = cc->getPhiave(i);
	}
	fwrite(&nPhiD,sizeof(double),1,filePointerWrite);
	fwrite(&buffer[0],sizeof(double),cc->nPhi(),filePointerWrite);
	
	// cross-correlation
	if (global->useCorrelation) {
		if (global->autoCorrelateOnly) {
			// autocorrelation only (q1=q2)
			for (int i=0; i<cc->nQ(); i++) {
				for (int k=0; k<cc->nLag(); k++) {
					if (global->sumCorrelation)
						info->correlation[i*cc->nLag()+k] = cc->getCrossCorrelation(i,k);
					else buffer[i*cc->nLag()+k] = cc->getCrossCorrelation(i,k);
				}
			}
			fwrite(&nQD,sizeof(double),1,filePointerWrite);
			fwrite(&nLagD,sizeof(double),1,filePointerWrite);
			
		} else {
			// full version
			for (int i=0; i<cc->nQ(); i++) {
				for (int j=0; j<cc->nQ(); j++) {
					for (int k=0; k<cc->nLag(); k++) {
						if (global->sumCorrelation)
							info->correlation[i*cc->nQ()*cc->nLag()+j*cc->nLag()+k] = cc->getCrossCorrelation(i,j,k);
						else buffer[i*cc->nQ()*cc->nLag()+j*cc->nLag()+k] = cc->getCrossCorrelation(i,j,k);
					}
				}
			}
			fwrite(&nQD,sizeof(double),1,filePointerWrite);
			fwrite(&nQD,sizeof(double),1,filePointerWrite);
			fwrite(&nLagD,sizeof(double),1,filePointerWrite);
			
		}
		if (global->sumCorrelation)
			fwrite(&info->correlation[0],sizeof(double),global->correlation_nn,filePointerWrite);
		else fwrite(&buffer[0],sizeof(double),global->correlation_nn,filePointerWrite); 
	}
	
	fclose(filePointerWrite);
	free(buffer);
	
	DEBUGL1_ONLY cout << "writeXCCA done" << endl;
}

