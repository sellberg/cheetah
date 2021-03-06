/*
 *  background.cpp
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

#include "background.h"


/*
 *	Update self generated darkcal file
 */
void updatePersistentBackground(tThreadInfo *threadInfo, cGlobal *global, int hit){
	float	gmd;
	// Add current (uncorrected) image to self darkcal
	if (!hit) {	
		pthread_mutex_lock(&global->selfdark_mutex);
		for(long i=0;i<global->pix_nn;i++){
			global->selfdark[i] = ( threadInfo->corrected_data[i] + (global->bgMemory-1)*global->selfdark[i]) / global->bgMemory;
		}
		gmd = (threadInfo->gmd21+threadInfo->gmd22)/2;
		global->avgGMD = ( gmd + (global->bgMemory-1)*global->avgGMD) / global->bgMemory;
		pthread_mutex_unlock(&global->selfdark_mutex);
	}

}


/*
 *	Subtract self generated darkcal file
 */
void subtractPersistentBackground(tThreadInfo *threadInfo, cGlobal *global){
	
	float	top = 0;
	float	s1 = 0;
	float	s2 = 0;
	float	v1, v2;
	float	factor;
	
	
	
	// Find appropriate scaling factor 
	if(global->scaleBackground) {
		for(long i=0;i<global->pix_nn;i++){
			//v1 = pow(global->selfdark[i], 0.25);
			//v2 = pow(threadInfo->corrected_data[i], 0.25);
			v1 = global->selfdark[i];
			v2 = threadInfo->corrected_data[i];
			if(v2 > global->hitfinder.ADC)
				continue;
			
			// Simple inner product gives cos(theta), which is always less than zero
			// Want ( (a.b)/|b| ) * (b/|b|)
			top += v1*v2;
			s1 += v1*v1;
			s2 += v2*v2;
		}
		factor = top/s1;
	}
	else 
		factor=1;
	
	
	// Do the weighted subtraction
	// Watch out for integer wraparound! Not a problem anymore since corrected_data is float
	float diff;
	for(long i=0;i<global->pix_nn;i++) {
		diff = threadInfo->corrected_data[i] - (factor*global->selfdark[i]);	
		//if(diff < -32766) diff = -32767;
		//if(diff > 32766) diff = 32767;
		threadInfo->corrected_data[i] = diff;
	}	
}

