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
	void printRawData(uint16_t *buffer,long lSize);
	void writeSAXS();
	void writeXCCA();
	void dumpResults( std::string filename );
	
	//---------------------------------------------calculations
	void calculatePolarCoordinates();
	void calculateSAXS();
	void calculateXCCA();
	
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
