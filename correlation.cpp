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
//#include <math.h>
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

    int arraylength = RAW_DATA_LENGTH;
	
    //jas: calculate center of CArray and shift qx, qy accordingly
	double x0 = centerX(global->pix_x);
    double y0 = centerY(global->pix_y);
    for (int i=0; i<RAW_DATA_LENGTH; i++) {                 //previously function shiftCenter() in crosscorrelator
        global->pix_x[i] = global->pix_x[i] - x0;
        global->pix_y[i] = global->pix_y[i] - y0;
    }
    
    //create cross correlator object that takes care of the computations
	CrossCorrelator *cc = new CrossCorrelator( threadInfo->corrected_data, global->pix_x, global->pix_y, RAW_DATA_LENGTH );

    DEBUGL1_ONLY cc->setDebug(1);                           //turn on debug level inside the CrossCorrelator, if needed
    DEBUGL2_ONLY cc->setDebug(2);                           //turn on debug level inside the CrossCorrelator, if needed

	if (global->useCorrelation == 1) {
		
		DEBUGL1_ONLY cout << "XCCA regular" << endl;
		
		cc->calculatePolarCoordinates();
		cc->calculateSAXS();
		cc->calculateXCCA();	
		
		//writeSAXS(threadInfo, global, cc, threadInfo->eventname);
		writeXCCA(threadInfo, global, cc, threadInfo->eventname);
		
	} else if (global->useCorrelation == 2) {
		
		DEBUGL1_ONLY cout << "XCCA fast" << endl;

        array2D *polar = new array2D;
        array2D *corr = new array2D;
        
        cc->createLookupTable();
        
        double start_q = 5*cc->deltaq();
        double stop_q = cc->qmax();
        double number_q = 20;
        double start_phi = 0;
        double stop_phi = 360;
        double number_phi = 128;
        cc->calculatePolarCoordinates_FAST(polar, start_q, stop_q, number_q, start_phi, stop_phi, number_phi);

        cc->calculateXCCA_FAST( polar, corr );
        
        delete polar;
        delete corr;

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
	
	DEBUGL1_ONLY printf("writing XCCA to file...\n");
	FILE *filePointerWrite;
	char outfile[1024];
	double samplingLengthD = (double) cc->samplingLength(); // save everything as doubles
	double samplingLagD = (double) cc->samplingLag();
	double samplingAngleD = (double) cc->samplingAngle();
	double *buffer;
	buffer = (double*) calloc(cc->samplingLength()*cc->samplingLength()*cc->samplingLag(), sizeof(double));
	
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
				buffer[i*cc->samplingLag()+k] = cc->getCrossCorrelation(i,i,k);
			}
		}
		fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
		fwrite(&samplingLagD,sizeof(double),1,filePointerWrite);
		fwrite(&buffer[0],sizeof(double),cc->samplingLength()*cc->samplingLag(),filePointerWrite);
		
	} else {
		
		// full version
		for (int i=0; i<cc->samplingLength(); i++) {
			for (int j=0; j<cc->samplingLength(); j++) {
				for (int k=0; k<cc->samplingLag(); k++) {
					buffer[i*cc->samplingLength()*cc->samplingLag()+j*cc->samplingLag()+k] = cc->getCrossCorrelation(i,j,k);
				}
			}
		}
		fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
		fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
		fwrite(&samplingLagD,sizeof(double),1,filePointerWrite);
		fwrite(&buffer[0],sizeof(double),cc->samplingLength()*cc->samplingLength()*cc->samplingLag(),filePointerWrite);
		
	}
	
	fclose(filePointerWrite);
	free(buffer);
	
	DEBUGL1_ONLY cout << "writeXCCA done" << endl;
}



//----------------------------------------------------------centerX
//find center of scattering image in x-direction
//----------------------------------------------------------
double centerX( float *qx ) {
    double center = 0;
	int quads = 4;
	// Loop over quads and pick out closest pixel to center
	for (int i=0; i<quads; i++) {
		center += (double) qx[8*ROWS*(2*COLS-1)+i*2*ROWS];
	}
	cout << "corrected center in X: " << center/quads << endl;
	return center/quads;
}


//----------------------------------------------------------centerY
//find center of scattering image in x-direction
//----------------------------------------------------------
double centerY( float *qy ) {
    double center = 0;
	int quads = 4;
	// Loop over quads and pick out closest pixel to center
	for (int i=0; i<quads; i++) {
		center += (double) qy[8*ROWS*(2*COLS-1)+i*2*ROWS];
	}
	cout << "corrected center in Y: " << center/quads << endl;
	return center/quads;
}



