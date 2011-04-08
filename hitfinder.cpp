/*
 *  hitfinder.cpp
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
#include <math.h>
#include <stdlib.h>

#include "setup.h"
#include "worker.h"
#include "hitfinder.h"


/*
 *	A basic hitfinder
 */
int  hitfinder(tThreadInfo *threadInfo, cGlobal *global, cHitfinder *hitf){

	long	nat, lastnat;
	long	counter;
	int		hit=0;
	long	ii,jj,nn;

	nat = 0;
	counter = 0;

		/*
	 *	Use a data buffer so we can zero out pixels already counted
	 */
	int16_t *temp = (int16_t*) calloc(global->pix_nn, sizeof(int16_t));
	memcpy(temp, threadInfo->corrected_data, global->pix_nn*sizeof(int16_t));
	
	
	/*
	 *	Apply peak search mask 
	 *	(multiply data by 0 to ignore regions)
	 */
	if(hitf->UsePeakmask) {
		for(long i=0;i<global->pix_nn;i++){
			temp[i] *= hitf->peakmask[i]; 
		}
	}
	
	
	/*
	 *	Use one of various hitfinder algorithms
	 */
	switch(hitf->Algorithm) {
		
		case 1 :	// Simply count the number of pixels above ADC threshold (very basic)
			for(long i=0;i<global->pix_nn;i++){
				if(temp[i] > hitf->ADC){
					nat++;
				}
			}
			if(nat >= hitf->NAT)
				hit = 1;
			break;

	
		case 2 :	//	Count clusters of pixels above threshold
			for(long j=1; j<8*COLS-1; j++){
				for(long i=1; i<8*ROWS-1; i++) {
					nn = 0;
					ii = i+(8*ROWS)*j;
					if(temp[i+(8*ROWS)*j] > hitf->ADC) {
						nn += 1;
						if(temp[i+1+(8*ROWS)*j] > hitf->ADC) nn++;
						if(temp[i-1+(8*ROWS)*j] > hitf->ADC) nn++;
						if(temp[i+(8*ROWS)*(j+1)] > hitf->ADC) nn++;
						if(temp[i+1+(8*ROWS)*(j+1)] > hitf->ADC) nn++;
						if(temp[i-1+(8*ROWS)*(j+1)] > hitf->ADC) nn++;
						if(temp[i+(8*ROWS)*(j-1)] > hitf->ADC) nn++;
						if(temp[i+1+(8*ROWS)*(j-1)] > hitf->ADC) nn++;
						if(temp[i-1+(8*ROWS)*(j-1)] > hitf->ADC) nn++;
					}
					if(nn >= hitf->Cluster) {
						nat++;
						temp[i+(8*ROWS)*j] = 0;
						temp[i+1+(8*ROWS)*j] = 0;
						temp[i-1+(8*ROWS)*j] = 0;
						temp[i+(8*ROWS)*(j+1)] = 0;
						temp[i+1+(8*ROWS)*(j+1)] = 0;
						temp[i-1+(8*ROWS)*(j+1)] = 0;
						temp[i+(8*ROWS)*(j-1)] = 0;
						temp[i+1+(8*ROWS)*(j-1)] = 0;
						temp[i-1+(8*ROWS)*(j-1)] = 0;
					}
				}
			}
			threadInfo->nPeaks = nat;
			if(nat >= hitf->MinPixCount)
				hit = 1;
			break;


	
	
		case 3 : 	// Real peak counter
		default:
			int search_x[] = {-1,0,1,-1,1,-1,0,1};
			int search_y[] = {-1,-1,-1,0,0,1,1,1};
			int	search_n = 8;
			long e;
			long *inx = (long *) calloc(global->pix_nn, sizeof(long));
			long *iny = (long *) calloc(global->pix_nn, sizeof(long));
			// Loop over modules (8x8 array)
			for(long mj=0; mj<8; mj++){
				for(long mi=0; mi<8; mi++){
					
					// Loop over pixels within a module
					for(long j=1; j<COLS-1; j++){
						for(long i=1; i<ROWS-1; i++){

							e = (j+mj*COLS)*global->pix_nx;
							e += i+mi*ROWS;

							//if(e >= global->pix_nn)
							//	printf("Array bounds error: e=%i\n");
							
							if(temp[e] > hitf->ADC){
								// This might be the start of a peak - start searching
								inx[0] = i;
								iny[0] = j;
								nat = 1;
								
								// Keep looping until the pixel count within this peak does not change (!)
								do {
									lastnat = nat;
									// Loop through points known to be within this peak
									for(long p=0; p<nat; p++){
										// Loop through search pattern
										for(long k=0; k<search_n; k++){
											// Array bounds check
											if((inx[p]+search_x[k]) < 0)
												continue;
											if((inx[p]+search_x[k]) >= ROWS)
												continue;
											if((iny[p]+search_y[k]) < 0)
												continue;
											if((iny[p]+search_y[k]) >= COLS)
												continue;
											
											// Neighbour point 
											e = (iny[p]+search_y[k]+mj*COLS)*global->pix_nx;
											e += inx[p]+search_x[k]+mi*ROWS;
											
											//if(e < 0 || e >= global->pix_nn){
											//	printf("Array bounds error: e=%i\n",e);
											//	continue;
											//}
											
											// Above threshold?
											if(temp[e] > hitf->ADC){
												//if(nat < 0 || nat >= global->pix_nn) {
												//	printf("Array bounds error: nat=%i\n",nat);
												//	break
												//}
												temp[e] = 0;
												inx[nat] = inx[p]+search_x[k];
												iny[nat] = iny[p]+search_y[k];
												nat++; 
											}
										}
									}
								} while(lastnat != nat);
								
								// Peak or junk?
								if(nat>=hitf->MinPixCount && nat<=hitf->MaxPixCount) {
									counter ++;
								}
							}
						}
					}
				}
			}	
			// Hit?
			threadInfo->nPeaks = counter;
			if(counter >= hitf->Npeaks && counter <= hitf->NpeaksMax)
				hit = 1;
			
			free(inx);
			free(iny);
			break;
	}
		
	
	
	// Update central hit counter
	if(hit) {
		pthread_mutex_lock(&global->nhits_mutex);
		global->nhits++;
		pthread_mutex_unlock(&global->nhits_mutex);
	}
	
	free(temp);
	return(hit);
}
