//
//  arraydataIO.h
//  xcca_commandline
//
//  Created by Feldkamp on 7/3/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//
//
//	handles input/output to/from arraydata objects
//
//	INFO: this class relies on other libraries to provide the io for different formats
//	INFO: the dependence can be controlled by the ARRAYDATAIO_xxx preprocessor macros 
//	INFO: if none of the i/o is actually needed, you can undefine all of the macros
//

#ifndef _arraydataIO_h
#define _arraydataIO_h

#include <string>
#include "arrayclasses.h"

class arraydataIO{


//#define ARRAYDATAIO_EDF			//use EDF
#undef ARRAYDATAIO_EDF				//do not use EDF

#define ARRAYDATAIO_TIFF			//use TIFF
//#undef ARRAYDATAIO_TIFF				//do not use TIFF

#define ARRAYDATAIO_HDF5			//use HDF5
//#undef ARRAYDATAIO_HDF5			//do not use HDF5

public:
	arraydataIO();
	~arraydataIO();

	//------------------------------------------------------------------------------EDF
	int readFromEDF( std::string filename, array2D *&dest ) const;										// works
	int writeToEDF( std::string filename, array2D *src ) const;											// works
	
	//------------------------------------------------------------------------------Tiff	
	int readFromTiff( std::string filename, array2D *&dest ) const;										// works, exept for a scaling factor --> needs fix
	
	// scaleFlag = 0 --> do not scale data
	// scaleFlag = 1 --> scale data to 0..65535
	int writeToTiff( std::string filename, array2D *src, int scaleFlag = 0, int verbose = 0 ) const;    // works
	
	
	//------------------------------------------------------------------------------HDF5
	int readFromHDF5( std::string filename, array2D *&dest ) const;										// seems to work, still testing, use with caution
	
	// dataType = 0 --> write as doubles (default) 
	// dataType = 1 --> write as float
	// dataType = 2 --> write as int
	int writeToHDF5( std::string filename, array2D *src, int dataType = 0, int verbose = 0 ) const;		// works
};


#endif
