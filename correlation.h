/*
 *  correlation.h
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

#ifndef _correlation_h
#define _correlation_h

#include "setup.h"
#include "worker.h"
#include "crosscorrelator.h"


/*
 *	Function prototypes
 */
void correlate(tThreadInfo*, cGlobal*);
void writeSAXS(tThreadInfo *info, cGlobal *global, CrossCorrelator *cc, char *eventname);
void writeXCCA(tThreadInfo *info, cGlobal *global, CrossCorrelator *cc, char *eventname);
double centerX( float *qx ); // calculates center in X from CArray of X positions
double centerY( float *qy ); // calculates center in Y from CArray of Y positions


#endif
