/*
 *  CrossCorrelator.h
 *  xcca
 *
 *  Created by Feldkamp on 2/17/11.
 *  Last changed on 04/27/11.
 *  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
 *
 */
#ifndef _crosscorrelator_h
#define _crosscorrelator_h

#include "arrayclasses.h"

#include <string>

class CrossCorrelator {

private:
	int p_arraySize;
	double p_qmax;
	double p_deltaq;		// step size (bin length) in q-direction
	double p_deltaphi;		// step size (bin length) in phi direction
	int p_samplingLength;
	int p_samplingAngle;
	int p_samplingLag;

	array1D *data;			//data storage
	array1D *qx;				//pixel x coordinate
	array1D *qy;				//pixel y coordinate
	
	array1D *q;				//magnitude of q vector (1st dimension in correlation)
	array1D *phi;			//angle (2nd dimension in correlation)
	
	array1D *qave;
	array1D *iave;
	array1D *phiave;
	array1D *crossCorrelation;	
	
	
	void updateDependentVariables();
	
public:
	//---------------------------------------------constructors & destructor
    CrossCorrelator( int arraylength=100 );                     //init with 1D data size
    CrossCorrelator( int16_t *dataCArray, int arraylength );    //init with actual data
	~CrossCorrelator();
	
	//---------------------------------------------input/output
	void initFromFile( std::string filename, int type=0 );
    void initWithTestPattern( int type=0 );                     //generate some test
	void printRawData(uint16_t *buffer,long lSize);
	void writeSAXS();
	void writeXCCA();
	void dumpResults( std::string filename );
	
	//---------------------------------------------calculations (Jonas's way)
	void calculatePolarCoordinates();
	void calculateSAXS();
	void calculateXCCA();
	
    //---------------------------------------------alternative approach (Jan's way)
    // these functions have the byname FAST to distinguish them from the ones above 
    // for lack of better name and in hope that it may be fast. We'll see...
    
    // 'calculatePolarCoordinates' returns a 2D pattern in polar coordinates (r vs. phi)
    int calculatePolarCoordinates_FAST(array2D* polar); 
    
    // 'calculateXCCA' returns the autocorrelation function (r = const.)          
    int calculateXCCA_FAST( array2D *polarData, array2D *corr );
    
    // looks up the value closest to xcoord, ycoord in the data
    double lookup( double xcoord, double ycoord );
    
    //compute 1D correlations using FFT, the result is returned in f, respectively
    int correlateFFT( array1D *f, array1D *g );    // result of corr(f,g) -> f
    int autocorrelateFFT( array1D *f );            // result of corr(f,f) -> f
        
	//---------------------------------------------setters & getters
	int arraySize();                                    // returns private variable p_arraySize
	void setArraySize( int arraySize_val );
	int matrixSize();                                   // returns sqrt(p_arraySize)
	void setMatrixSize( int matrixSize_val );
                                                        // matrixSize is just a little helper for now
                                                        // to come up with a q-calibration
                                                        // need to change this soon...
	
	double qmax();
	void setQmax( double qmax_val );
	
	//---------------------------------------------getters for dependent variables
	double deltaq();
	double deltaphi();
	int samplingLength();
	int samplingAngle();
	int samplingLag();
	
};


#endif
