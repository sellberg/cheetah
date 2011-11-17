/*
 *  commonmode.cpp
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


#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <stdlib.h>
#include <hdf5.h>
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include "commonmode.h"
#include "worker.h"
#include "peakdetect.h"
#include "point.h"


/*
 *	Subtract common mode on each 2x1 module
 *	This is done in a very slow way now - speed up later once we know it works!
 */
void cmModuleSubtract(tThreadInfo *threadInfo, cGlobal *global){

	DEBUGL2_ONLY printf("cmModuleSubtract\n");
	
	long		e;
	long		counter;
//	uint16_t	value;
	long		median;
	char		filename[1024];
	
	// Create histogram array
	int			noffset = 65535;
	long		nhist = 65536+noffset;
	uint16_t	*histograms;
	//histogram = (uint16_t*) calloc(nhist, sizeof(uint16_t));
	histograms = (uint16_t*) calloc(64*nhist, sizeof(uint16_t)); //histograms = (uint16_t*) calloc(32*nhist, sizeof(uint16_t));
	
	// Peak detection parameters
	int			peakfinderStart = global->cmStart;
	int			peakfinderStop = global->cmStop;
	float		peakfinderDelta = global->cmDelta;
	int			nPeakfinder = peakfinderStop-peakfinderStart+1;
	
	// Loop over each ASIC (8x8 array)
	for(long mi=0; mi<8; mi++){ //for(long mi=0; mi<4; mi++){
		for(long mj=0; mj<8; mj++){

			// Zero histogram
			//memset(histogram, 0, nhist*sizeof(uint16_t));
			
			// Loop over pixels within an ASIC
			for(int i=0; i<ROWS; i++){ //for(long i=0; i<2*ROWS; i++){
				for(int j=0; j<COLS; j++){
					e = (j + mj*COLS) * (8*ROWS);
					e += i + mi*ROWS; //e += i + mi*2*ROWS;
					//histogram[int(round(threadInfo->corrected_data[e])+noffset)] += 1;
					histograms[int(round(threadInfo->corrected_data[e])+noffset)+mj*nhist+mi*8*nhist]++;
				}
			}
			
			if (global->cmModule == 1) {
				// NEW COMMON-MODE ALGORITHM BASED ON PEAK DETECTION
				if (peakfinderStart < peakfinderStop && noffset+peakfinderStart >= 0 && noffset+peakfinderStop < nhist && peakfinderDelta > 0) { // sanity check
					uint16_t *peakfinderHist = (uint16_t*) calloc(nPeakfinder, sizeof(uint16_t));
					int *peakfinderHistX = (int*) calloc(nPeakfinder, sizeof(int));
					for (long i=0; i<nPeakfinder; i++) {
						peakfinderHist[i] = histograms[i+noffset+peakfinderStart+mj*nhist+mi*8*nhist];
						peakfinderHistX[i] = peakfinderStart+i;
					}
					
					PeakDetect peakfinder(peakfinderHistX, peakfinderHist, nPeakfinder);
					peakfinder.findAll(peakfinderDelta);
					
					Point *min, *max;
					if (peakfinder.maxima->size() > 0) {
						for (int k=0; k<peakfinder.maxima->size(); k++) {
							min = peakfinder.minima->get(k);
							max = peakfinder.maxima->get(k);
							DEBUGL1_ONLY cout << "max->getX()-min->getX() = " << max->getX()-min->getX() << ", max->getY()-peakfinderHist[max->getX()-peakfinderStart-1] = " << max->getY()-peakfinderHist[max->getX()-peakfinderStart-1] << ", max->getY()-peakfinderHist[max->getX()-peakfinderStart+1] = " << max->getY()-peakfinderHist[max->getX()-peakfinderStart+1] << endl;
							if (max->getX()-min->getX() > 4) { //if (max->getX()-min->getX() > 2 && max->getY()-peakfinderHist[max->getX()-peakfinderStart-1] < peakfinderDelta && max->getY()-peakfinderHist[max->getX()-peakfinderStart+1] < peakfinderDelta) {
								int commonmode = max->getX();
								for (int i=0; i<ROWS; i++) { //for (int i=0; i<2*ROWS; i++) {
									for (int j=0; j<COLS; j++) {
										e = (j + mj*COLS) * (8*ROWS);
										e += i + mi*ROWS; //e += i + mi*2*ROWS;
										threadInfo->corrected_data[e] -= commonmode;
									}
								}
								DEBUGL1_ONLY printf("r%04u:%i ", (int)threadInfo->runNumber, (int)threadInfo->threadNum);
								DEBUGL1_ONLY cout << "Commonmode (Q" << mi/2 << ", S" << mi%2+2*mj << "): " << commonmode << endl; //DEBUGL1_ONLY cout << "Commonmode (Q" << mi << ", S" << mj << "): " << commonmode << endl;
								break;
							} else if (k == peakfinder.maxima->size()-1) {
								// June data
								//cout << "Commonmode (Q" << mi << ", S" << mj << "): N/A" << endl;
								// Feb data (ASICs missing)
								if ((mi != 0 || mj != 5) && (mi != 1 || mj != 5) && (mi != 5 || mj != 3) && (mi != 6 || mj != 4) && (mi != 7 || mj != 4) && (mi != 6 || mj != 6) && (mi != 7 || mj != 6)) { //if ((mi != 0 || mj != 5) && (mi != 3 || mj != 4) && (mi != 3 || mj != 6)) {
									printf("r%04u:%i ", (int)threadInfo->runNumber, (int)threadInfo->threadNum);
									cout << "Commonmode (Q" << mi/2 << ", S" << mi%2+2*mj << "): N/A" << endl; //cout << "Commonmode (Q" << mi << ", S" << mj << "): N/A" << endl;
								}
							}
						}
					} else {
						printf("r%04u:%i ", (int)threadInfo->runNumber, (int)threadInfo->threadNum);
						cout << "Commonmode (Q" << mi/2 << ", S" << mi%2+2*mj << "): N/A (no maxima)" << endl; //cout << "Commonmode (Q" << mi << ", S" << mj << "): N/A (no maxima)" << endl;
					}
					
					free(peakfinderHist);
					free(peakfinderHistX);
					
				} else {
					cerr << "ERROR in cmModuleSubtract: Input parameters are out of range." << endl;
					cout << "\tpeakfinderStart: " << peakfinderStart << endl;
					cout << "\tpeakfinderStop: " << peakfinderStop << endl;
					cout << "\tpeakfinderDelta: " << peakfinderDelta << endl;
				}				
			} else if (global->cmModule == 2) {
			// OLD COMMON-MODE ALGORITHM BASED ON MEDIAN
				// Find median value
				counter = 0;
				for(long i=0; i<nhist; i++){
					//counter += histogram[i];
					counter += histograms[i+mj*nhist+mi*8*nhist];
					if(counter > (global->cmFloor*ROWS*COLS)) { //if(counter > (global->cmFloor*2*ROWS*COLS)) {
						median = i-noffset;
						
						if (global->cmSaveHistograms) {
							if (histograms[mj*nhist+mi*8*nhist] == 0) histograms[mj*nhist+mi*8*nhist] = i;
							else {
								cout << "1st element of histogram non-zero! Save median to 11th element." << endl;
								if (histograms[10+mj*nhist+mi*8*nhist] == 0) histograms[10+mj*nhist+mi*8*nhist] = i;
								else {
									cout << "11th element of histogram non-zero! Save median to 101st element." << endl;
									if (histograms[100+mj*nhist+mi*8*nhist] == 0) histograms[100+mj*nhist+mi*8*nhist] = i;
									else cout << "101st element of histogram non-zero! Aborting save of median..." << endl;
								}							
							}
						}
						
						break;
					}
					
				}
				
				
				// Ignore common mode for ASICs without wires (only Feb run)
				if ((mi == 1 && mj == 6) || (mi == 2 && mj == 5) || (mi == 3 && mj == 5) || (mi == 4 && mj == 5) || (mi == 4 && mj == 6)) {
					median = 0;
				}
				
				// Subtract median value
				for(long i=0; i<ROWS; i++){ //for(long i=0; i<2*ROWS; i++){
					for(long j=0; j<COLS; j++){
						e = (j + mj*COLS) * (8*ROWS);
						e += i + mi*ROWS; //e += i + mi*2*ROWS;
						threadInfo->corrected_data[e] -= median;
					}
				}
			} else {
				cerr << "WARNING in cmModuleSubtract: No such common-mode algorithm exists, common-mode correction disabled." << endl;
			}
		}
	}
	
	//free(histogram);
	if (global->cmSaveHistograms) {
		sprintf(filename,"%s-hist.h5",threadInfo->eventname);
		//writeSimpleHDF5(filename, histograms, (int)nhist, 8, 4, H5T_STD_U16LE);
		// save as 2D array with each ASIC sorted by quad
		writeSimpleHDF5(filename, histograms, (int)nhist, 64, H5T_STD_U16LE); //writeSimpleHDF5(filename, histograms, (int)nhist, 32, H5T_STD_U16LE);
	}
	free(histograms);
}

/*
 *	Subtract common mode on each module
 *	This is done in a very slow way now - speed up later once we know it works!
 */
void cmSubModuleSubtract(tThreadInfo *threadInfo, cGlobal *global){
	
	// ROWS = 194;
	// COLS = 185;
	
	long		e;
	long		ii,jj;
	long		counter;
//	uint16_t	value;
	uint16_t	median;
	
	// Create histogram array
	int			nhist = 131071;
	uint16_t	*histogram;
	histogram = (uint16_t*) calloc(nhist, sizeof(uint16_t));
	
	// Subunits
	long		nn = global->cmSubModule;		// Multiple of 2 please!
	if (nn == 1 || nn % 2) {
		cout << "subtractCMSubModule must be a multiple of 2, aborting common-mode correction." << endl;
		return;
	}
	
	// Loop over whole modules (8x8 array)
	for(long mi=0; mi<8; mi++){
		for(long mj=0; mj<8; mj++){
			
			
			// Loop over sub-modules
			for(long smi=0; smi<ROWS; smi+=ROWS/nn){
				for(long smj=0; smj<COLS; smj+=COLS/nn){
				
					// Zero histogram
					memset(histogram, 0, nhist*sizeof(uint16_t));
				
				
					// Loop over pixels within this subregion
					for(long i=0; i<ROWS/nn && (i+smi)<ROWS; i++){
						for(long j=0; j<COLS/nn && (j+smj)<COLS; j++){
							jj = smj + j + mj*COLS;
							ii = smi + i + mi*ROWS;
							e = ii + jj*8*ROWS;
							histogram[int(round(threadInfo->corrected_data[e]))+65535] += 1;
						}
					}
					
					// Find median value
					counter = 0;
					for(long i=0; i<nhist; i++){
						counter += histogram[i];
						if(counter > (global->cmFloor*ROWS*COLS/(nn*nn))) {
							median = i-65535;
							break;
						}
					}
				
					// Subtract median value
					for(long i=0; i<ROWS/nn && (i+smi)<ROWS; i++){
						for(long j=0; j<COLS/nn && (j+smj)<COLS; j++){
							jj = smj + j + mj*COLS;
							ii = smi + i + mi*ROWS;
							e = ii + jj*8*ROWS;
							threadInfo->corrected_data[e] -= median;
						}
					}
				}
			}
		}
	}
	free(histogram);
}


