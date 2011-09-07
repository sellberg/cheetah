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


//  For the correlation functionality, the cheetah needs the 'giraffe' library
//  available at git://github.com/feldkamp/giraffe.git
//
//  The giraffe library wraps the classes that, previous to Sep 2011 was in cheetah's own
//     arrayclasses.cpp,h
//     arraydataIO.cpp,h
//     crosscorrelator.cpp,h
//     fouriertransformer.cpp,h
//  
//  If the preprocessor flag CORRELATION_ENABLED is not set, cheetah does not need 
//  the giraffe library (or the supporting libraries for fftw and tiff) 
//  but, in turn, it cannot perform correlation calculations, either
//
#ifdef CORRELATION_ENABLED

	#include "crosscorrelator.h"
	
	void correlate(tThreadInfo*, cGlobal*);
	void writeSAXS(tThreadInfo *info, cGlobal *global, CrossCorrelator *cc, char *eventname);
	void writeXCCA(tThreadInfo *info, cGlobal *global, CrossCorrelator *cc, char *eventname);

	#else //no correlation functionality --> define a dummy version of correlate that does nothing

	void correlate(tThreadInfo*, cGlobal*);

	#endif
#endif
