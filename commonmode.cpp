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

#include "commonmode.h"
#include "worker.h"


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
	long		nhist = 131071;
	uint16_t	*histogram, *histograms;
	//histogram = (uint16_t*) calloc(nhist, sizeof(uint16_t));
	histograms = (uint16_t*) calloc(32*nhist, sizeof(uint16_t));
	
	// Loop over 2x1 modules (4x8 array)
	for(long mi=0; mi<4; mi++){ //for(long mi=0; mi<8; mi++){
		for(long mj=0; mj<8; mj++){

			// Zero histogram
			//memset(histogram, 0, nhist*sizeof(uint16_t));
			
			// Loop over pixels within a module
			for(long i=0; i<2*ROWS; i++){ //for(long i=0; i<ROWS; i++){
				for(long j=0; j<COLS; j++){
					e = (j + mj*COLS) * (8*ROWS);
					e += i + mi*2*ROWS; //e += i + mi*ROWS;
					//histogram[int(round(threadInfo->corrected_data[e])+65535)] += 1;
					histograms[int(round(threadInfo->corrected_data[e])+65535)+mj*nhist+mi*8*nhist]++;
				}
			}
			
			// Find median value
			counter = 0;
			for(long i=0; i<nhist; i++){
				//counter += histogram[i];
				counter += histograms[i+mj*nhist+mi*8*nhist];
				if(counter > (global->cmFloor*2*ROWS*COLS)) { //if(counter > (global->cmFloor*ROWS*COLS)) {
					median = i-65535;
					
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
//			if ((mi == 1 && mj == 6) || (mi == 2 && mj == 5) || (mi == 3 && mj == 5) || (mi == 4 && mj == 5) || (mi == 4 && mj == 6)) {
//				median = 0;
//			}
			
			// Subtract median value
			for(long i=0; i<2*ROWS; i++){
				for(long j=0; j<COLS; j++){
					e = (j + mj*COLS) * (8*ROWS);
					e += i + mi*2*ROWS; //e += i + mi*ROWS;
					threadInfo->corrected_data[e] -= median;
					
				}
			}
		}
	}
	//free(histogram);
	if (global->cmSaveHistograms) {
		sprintf(filename,"%s_cm-hist.h5",threadInfo->eventname);
		writeSimpleHDF5(filename, histograms, (int)nhist, 8, 4, H5T_STD_U16LE);
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


