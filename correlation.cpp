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
using std::endl;

#include <cmath>

#include "correlation.h"
#include "crosscorrelator.h"


//----------------------------------------------------------correlate
// apply angular cross correlation
//-------------------------------------------------------------------
void correlate(tThreadInfo *threadInfo, cGlobal *global){

    cout << "CORRELATING... in thread #" << threadInfo->threadNum << "." << endl;

    int arraylength = RAW_DATA_LENGTH;    
   	CrossCorrelator *cc = new CrossCorrelator( threadInfo->corrected_data, arraylength );
    

	cc->calculatePolarCoordinates();
	cc->calculateSAXS();
	cc->calculateXCCA();	
	
	cc->writeSAXS();
	cc->writeXCCA();


}


