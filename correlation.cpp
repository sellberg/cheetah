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
    
    int alg = 0;        //select the correlation algorithm... 
                        //should be governed by a setting in the ini file at some point
                        
    cout << "CORRELATING... using alg " << alg << " in thread #" << threadInfo->threadNum << "." << endl;
    
    //jas: calculate center of CArray and shift qx, qy accordingly
	double centerX = centerXCArray(global->pix_x);
    double centerY = centerYCArray(global->pix_y);
    for (int i=0; i<RAW_DATA_LENGTH; i++) {                 //previously function shiftCenter() in crosscorrelator
        global->pix_x[i] = global->pix_x[i] - centerX;
        global->pix_y[i] = global->pix_y[i] - centerY;
    }
    
    //create cross correlator object that takes care of the computations
	CrossCorrelator *cc = new CrossCorrelator( threadInfo->corrected_data, global->pix_x, global->pix_y, RAW_DATA_LENGTH );
    
    
    switch (alg) {
        case 0:{
                cout << "XCCA regular" << endl;
                cc->calculatePolarCoordinates();
                cc->calculateSAXS();
                cc->calculateXCCA();	
                
                cc->writeSAXS();
                //cc->writeXCCA();
            }
            break;
        case 1:{
                cout << "XCCA FAST" << endl;
                array2D *polar = new array2D;
                array2D *corr = new array2D;
                cc->calculatePolarCoordinates_FAST( polar );
                cc->calculateXCCA_FAST( polar, corr );
                delete polar;
                delete corr;
            }
            break;
        default:
            cerr << "Error in correlate()! Correlation algorithm " << alg << " not known." << endl;
    }

    delete cc;
}




//----------------------------------------------------------centerXCArray
//find center of scattering image in x-direction
//----------------------------------------------------------
double centerXCArray( float *qxCArray ) {
    double center = 0;
	int quads = 4;
	// Loop over quads and pick out closest pixel to center
	for (int i=0; i<quads; i++) {
		center += (double) qxCArray[8*ROWS*(2*COLS-1)+i*2*ROWS];
	}
	cout << "new Center in X: " << center/quads << endl;
	return center/quads;
}


//----------------------------------------------------------centerXCArray
//find center of scattering image in x-direction
//----------------------------------------------------------
double centerYCArray( float *qyCArray ) {
    double center = 0;
	int quads = 4;
	// Loop over quads and pick out closest pixel to center
	for (int i=0; i<quads; i++) {
		center += (double) qyCArray[8*ROWS*(2*COLS-1)+i*2*ROWS];
	}
	cout << "new Center in Y: " << center/quads << endl;
	return center/quads;
}



///////////////OLD STUFF////////////////

/*
//----------------------------------------------------------shiftCenter
//
//----------------------------------------------------------
void shiftCenter() {
	for (int i=0; i<arraySize(); i++) {
		qx->set(i, qx->get(i)-centerX());
		qy->set(i, qy->get(i)-centerY());
	}
}
*/
