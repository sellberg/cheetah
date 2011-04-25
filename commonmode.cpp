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
#include <math.h>
#include <stdlib.h>


#include "commonmode.h"


/*
 *	Subtract common mode on each module
 *	This is done in a very slow way now - speed up later once we know it works!
 */
void cmModuleSubtract(tThreadInfo *threadInfo, cGlobal *global){

	DEBUGL2_ONLY printf("cmModuleSubtract\n");
	
	long		e;
	long		counter;
	uint16_t	value;
	uint16_t	median;
	
	// Create histogram array
	int			nhist = 65535;
	uint16_t	*histogram;
	histogram = (uint16_t*) calloc(nhist, sizeof(uint16_t));
	
	// Loop over modules (8x8 array)
	for(long mi=0; mi<8; mi++){
		for(long mj=0; mj<8; mj++){

			// Zero histogram
			memset(histogram, 0, nhist*sizeof(uint16_t));
			
			
			// Loop over pixels within a module
			for(long i=0; i<ROWS; i++){
				for(long j=0; j<COLS; j++){
					e = (j + mj*COLS) * (8*ROWS);
					e += i + mi*ROWS;
					histogram[threadInfo->corrected_data[e]] += 1;
				}
			}
			
			// Find median value
			counter = 0;
			for(long i=0; i<nhist; i++){
				counter += histogram[i];
				if(counter > (global->cmFloor*ROWS*COLS)) {
					median = i;
					break;
				}
			}
			//DEBUGL2_ONLY printf("Median of module (%i,%i) = %i\n",mi,mj,median);

			// Subtract median value
			for(long i=0; i<ROWS; i++){
				for(long j=0; j<COLS; j++){
					e = (j + mj*COLS) * (8*ROWS);
					e += i + mi*ROWS;
					threadInfo->corrected_data[e] -= median;

					// Zero checking only needed if corrected data is uint16
					//value = threadInfo->corrected_data[e];
					//if(value > median)
					//	threadInfo->corrected_data[e] -= median;
					//else
					//	threadInfo->corrected_data[e] = 0;
				}
			}
		}
	}
	free(histogram);
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
	uint16_t	value;
	uint16_t	median;
	
	// Create histogram array
	int			nhist = 65535;
	uint16_t	*histogram;
	histogram = (uint16_t*) calloc(nhist, sizeof(uint16_t));
	
	// Subunits
	long	nn=global->cmSubModule;		// Multiple of 2 please!
	if(nn < 2 )
		return;
	
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
							histogram[threadInfo->corrected_data[e]] += 1;
						}
					}
					
					// Find median value
					counter = 0;
					for(long i=0; i<nhist; i++){
						counter += histogram[i];
						if(counter > (0.25*ROWS*COLS/(nn*nn))) {
							median = i;
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


