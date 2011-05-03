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

    //int arraylength = RAW_DATA_LENGTH;
	
	CrossCorrelator *cc = new CrossCorrelator( threadInfo->corrected_data, global->pix_x, global->pix_y );
    
    switch (alg) {
        case 0:{
                cout << "XCCA regular" << endl;
                cc->calculatePolarCoordinates();
                cc->calculateSAXS();
                //cc->calculateXCCA();	
                
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


