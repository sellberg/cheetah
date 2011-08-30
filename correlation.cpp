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
#include <string>
using std::string;

#include <cmath>

#include "correlation.h"
#include "arrayclasses.h"
#include "arraydataIO.h"





//----------------------------------------------------------correlate
// apply angular cross correlation
//-------------------------------------------------------------------
void correlate(tThreadInfo *threadInfo, cGlobal *global) {
    
    DEBUGL1_ONLY cout << "CORRELATING... in thread #" << threadInfo->threadNum << "." << endl;

	//prepare some things for output
	bool tif_out = true;
	bool bin_out = true;
	bool h5_out = true;
	
	arraydataIO *io = new arraydataIO;
	std::ostringstream osst;
	osst << threadInfo->eventname;
	string eventname_str = osst.str();
	
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
		DEBUGL1_ONLY cout << "XCCA regular (algorithm 1)" << endl;
		
		cc->calculatePolarCoordinates(global->correlationStartQ, global->correlationStopQ);
		cc->calculateSAXS();
		cc->calculateXCCA();
		
	//--------------------------------------------------------------------------------------------alg2
	} else if (global->useCorrelation == 2) {					
		DEBUGL1_ONLY cout << "XCCA fast (algorithm 2)" << endl;
		cc->setLookupTable( global->correlationLUT, global->correlationLUTdim1, global->correlationLUTdim2 );
		cc->calculatePolarCoordinates_FAST(global->correlationStartQ, global->correlationStopQ);
		
		// need to protect the FFTW at the core with a mutex, not thread-safe!!
		pthread_mutex_lock(&global->correlationFFT_mutex);
		cc->calculateXCCA_FAST();
		pthread_mutex_unlock(&global->correlationFFT_mutex);
		
		io->writeToTiff( eventname_str+"-polar.tif", cc->polar(), 1, 1 );		// 0: unscaled, 1: scaled

	} else {
		cerr << "Error in correlate()! Correlation algorithm " << global->useCorrelation << " not known." << endl;
	}
	
	
	//--------------------------------------------------------------------------------------------output for this shot
	// by JF: not sure if the mutex is really necessary... maybe not
	pthread_mutex_lock(&global->correlation_mutex);
	
	//bin output
	if (bin_out){
		writeXCCA(threadInfo, global, cc, threadInfo->eventname); 			// writes XCCA+SAXS to binary
	}
	//TIFF image output
	if (tif_out){
		if (global->autoCorrelateOnly){
			io->writeToTiff( eventname_str+"-xaca.tif", cc->autoCorr(), 1 );			// 0: unscaled, 1: scaled
		}else{
			cerr << "WARNING. No tiff output for 3D cross-correlation case implemented, yet!" << endl;
			//one possibility would be to write a stack of tiffs, one for each of the outer q values
		}
	}
	//HDF5 output
	if (h5_out){
		if (global->autoCorrelateOnly){
			io->writeToHDF5( osst.str()+"-xaca.h5", cc->autoCorr() );
		}else{
			cerr << "WARNING. No HDF5 output for 3D cross-correlation case implemented, yet!" << endl;
		}
	}
	
	pthread_mutex_unlock(&global->correlation_mutex);
	
	delete io;
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
		buffer[i] = cc->getIavg(i);
	}
	fwrite(&nQD,sizeof(double),1,filePointerWrite); // saving dimensions of array before the actual data
	fwrite(&buffer[0],sizeof(double),cc->nQ(),filePointerWrite);
	
	// q binning
	for (int i=0; i<cc->nQ(); i++) {
		buffer[i] = cc->getQavg(i);
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
	
	FILE *filePointerWrite;
	char outfile[1024];
	
	double nQD = (double) cc->nQ(); // save everything as doubles
	double nLagD = (double) cc->nLag();
	double nPhiD = (double) cc->nPhi();
	
	const int nQ = cc->nQ();
	const int nPhi = cc->nPhi();
	const int nLag = cc->nLag();
	
	double *buffer;
	buffer = (double*) calloc(nQ*nQ*nLag, sizeof(double));
	if (global->sumCorrelation) info->correlation = (double*) calloc(global->correlation_nn, sizeof(double));
	
	if (global->autoCorrelateOnly){
		sprintf(outfile,"%s-xaca.bin",eventname);
	} else {
		sprintf(outfile,"%s-xcca.bin",eventname);
	}
	DEBUGL1_ONLY printf("r%04u:%i (%2.1f Hz): Writing data to: %s\n", (int)global->runNumber, (int)info->threadNum, global->datarate, outfile);
	
	filePointerWrite = fopen(outfile,"w+");
	
	// angular averages
	for (int i=0; i<nQ; i++) {
		buffer[i] = cc->getIavg(i);
	}
	fwrite(&nQD,sizeof(double),1,filePointerWrite); // saving dimensions of array before the actual data
	fwrite(&buffer[0],sizeof(double),nQ,filePointerWrite);
	
	// q binning
	for (int i=0; i<nQ; i++) {
		buffer[i] = cc->getQavg(i);
	}
	fwrite(&nQD,sizeof(double),1,filePointerWrite);
	fwrite(&buffer[0],sizeof(double),nQ,filePointerWrite);
	
	// angle binning
	for (int i=0; i<nPhi; i++) {
		buffer[i] = cc->getPhiavg(i);
	}
	fwrite(&nPhiD,sizeof(double),1,filePointerWrite);
	fwrite(&buffer[0],sizeof(double),nPhi,filePointerWrite);
	
	// cross-correlation
	if (global->useCorrelation) {
		if (global->autoCorrelateOnly) {
			// autocorrelation only (q1=q2)
			for (int i=0; i < nQ; i++) {
				for (int k=0; k < nLag; k++) {
					if (global->sumCorrelation)
						info->correlation[i*nLag + k] = cc->getAutoCorrelation(i,k);
					else 
						buffer[i*nLag + k] = cc->getAutoCorrelation(i,k);
				}
			}
			fwrite(&nQD,sizeof(double),1,filePointerWrite);
			fwrite(&nLagD,sizeof(double),1,filePointerWrite);
			
		} else {
			// full version
			for (int i=0; i < nQ; i++) {
				for (int j=0; j < nQ; j++) {
					for (int k=0; k < nLag; k++) {
						if (global->sumCorrelation)
							info->correlation[i*nQ*nLag + j*nLag + k] = cc->getCrossCorrelation(i,j,k);
						else 
							buffer[i*nQ*nLag + j*nLag + k] = cc->getCrossCorrelation(i,j,k);
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
}

