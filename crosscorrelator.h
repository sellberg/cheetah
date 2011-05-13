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
	double p_qmin;
	double p_deltaq;		// step size (bin length) in q-direction
	double p_deltaphi;		// step size (bin length) in phi direction
	int p_samplingLength;
	int p_samplingAngle;
	int p_samplingLag;
    std::string p_outputdir;  // the output directory if anything is dumped from withing this class (default is working dir)

	array1D *data;			//data storage
	array1D *qx;				//pixel x coordinate
	array1D *qy;				//pixel y coordinate
	
	array1D *q;				//magnitude of q vector (1st dimension in correlation)
	array1D *phi;			//angle (2nd dimension in correlation)
	
	array1D *qave;
	array1D *iave;
	array1D *phiave;
	array3D *crossCorrelation;	
	void updateDependentVariables();
    
    array2D *table;         //lookup table

    //-----some debug features that can turned on via setDebug()
    int p_debug;
    array1D *check1D;
    
public:
	//---------------------------------------------constructors & destructor
    CrossCorrelator( int arraylength=1 );                                                           //init with 1D data size
	CrossCorrelator( int16_t *dataCArray, int arraylength );                                        //init with actual data
    CrossCorrelator( int16_t *dataCArray, float *qxCArray, float *qyCArray, int arraylength );      //init with actual data + centered(!) q calibration
	CrossCorrelator( int16_t *dataCArray, float *qxCArray, float *qyCArray, int arraylength, double qMax, double qMin );
	
	~CrossCorrelator();
	
	//---------------------------------------------input/output
    void initPrivateVariables();
    void initDefaultQ();
	void initFromFile( std::string filename, int type=0 );
    void initWithTestPattern( int sizex, int sizey, int type=0 );                           //generate some test patterns
	void printRawData(uint16_t *buffer,long lSize);
	void dumpResults( std::string filename );
	
	//---------------------------------------------calculations (Jonas's way)
	void calculatePolarCoordinates();
	void calculateSAXS();
	void calculateXCCA();
	
    //---------------------------------------------alternative approach (Jan's way)
    // these functions have the byname _FAST to distinguish them from the ones above 
    // (for lack of a better name and in the hope that they may be fast. We'll see...)
    
    // 'calculatePolarCoordinates' returns a 2D pattern in polar coordinates (r vs. phi)
    int calculatePolarCoordinates_FAST(array2D* polar, double number_q=20, double number_phi=128); 
    int calculatePolarCoordinates_FAST(array2D* polar, double number_q, double start_q, double stop_q,
                                                        double number_phi=128, double start_phi=0, double stop_phi=360 );
    
    // 'calculateXCCA' returns the autocorrelation function (r = const.)          
    int calculateXCCA_FAST( array2D *polarData, array2D *corr );
    
    // looks up the value closest to xcoord, ycoord in the data
    int createLookupTable();
    double lookup( double xcoord, double ycoord );
    
    //compute 1D correlations using FFT, the result is returned in f, respectively
    int correlateFFT( array1D *f, array1D *g );    // result of corr(f,g) -> f
    int autocorrelateFFT( array1D *f );            // result of corr(f,f) -> f
        
	//---------------------------------------------setters & getters
	int arraySize() const;                              // returns private variable p_arraySize
	void setArraySize( int arraySize_val );
	int matrixSize() const;                             // returns sqrt(p_arraySize)
	void setMatrixSize( int matrixSize_val );
                                                        // matrixSize is just a little helper for now
                                                        // to come up with a q-calibration
                                                        // need to change this soon...
	
	void setQmaxmin( double qmax_val, double qmin_val );	
	double qmax() const;
	void setQmax( double qmax_val );
	double qmin() const;
	void setQmin( double qmax_val );
	
	// jas: qmaxCArray() could easily be rewritten to use array1D qx, qy instead of CArrays if preferable
	double qmaxCArray( float *qxCArray, float *qyCArray, int arraylength ); // calculates qmax from CArrays of X/Y positions
    double qxmaxCArray( float *qxCArray, int arraylength ); // calculates qmax from CArray of X positions
	double qymaxCArray( float *qyCArray, int arraylength ); // calculates qmax from CArray of Y positions
	
    void setOutputdir( std::string dir );
    std::string outputdir();

    int debug();                                        // control the amount of (commandline) talking while running, default: 0
    void setDebug( int debuglevel );

	//---------------------------------------------getters for dependent variables
	double deltaq() const;
	double deltaphi() const;
	int samplingLength() const;
	int samplingAngle() const;
	int samplingLag() const;
	
	//---------------------------------------------getters for calculated arrays
	double getQave(unsigned index) const;
	double getPhiave(unsigned index) const;
	double getIave(unsigned index) const;
	double getCrossCorrelation(unsigned index1, unsigned index2, unsigned index3) const;
};





#endif
