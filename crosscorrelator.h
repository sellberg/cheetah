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
	double p_qmin;
	double p_qmax;
	double p_deltaq;		// step size (bin length) in q-direction
	double p_phimin;
	double p_phimax;
	double p_deltaphi;		// step size (bin length) in phi direction
	int p_samplingLength;
	int p_samplingAngle;
	int p_samplingLag;
	bool p_mask;			// enables/disables masking of bad pixels
    std::string p_outputdir;  // the output directory if anything is dumped from withing this class (default is working dir)

	array1D *data;			//data storage
	array1D *qx;			//pixel x coordinate
	array1D *qy;			//pixel y coordinate
	array1D *mask;			//mask used to remove bad pixels
	
	array1D *q;				//magnitude of q vector (1st dimension in correlation)
	array1D *phi;			//angle (2nd dimension in correlation)
	
	array1D *qave;
	array1D *iave;
	array1D *phiave;
	array3D *crossCorrelation;	
	void updateDependentVariables();
    
    array2D *table;         //lookup table
	fftw_plan p_fplan;		//forward FFT plan
	fftw_plan p_bplan;		//backward FFT plan

    //-----some debug features that can turned on via setDebug()
    int p_debug;
//    array1D *check1D;
	float *p_dataCArray; 	//!!!only temporary!!! keep the cheetah c-array to be able to write to it for debugging purposes
    
public:
	//---------------------------------------------constructors & destructor
    CrossCorrelator( int arraylength=1 );                                                           //init with 1D data size
	CrossCorrelator( float *dataCArray, int arraylength );                                        //init with actual data
    CrossCorrelator( float *dataCArray, float *qxCArray, float *qyCArray, int arraylength );      //init with actual data + centered(!) q calibration
	CrossCorrelator( float *dataCArray, float *qxCArray, float *qyCArray, int arraylength, double qMax, double qMin );
	CrossCorrelator( float *dataCArray, float *qxCArray, float *qyCArray, int arraylength, int nq, int nphi );
	CrossCorrelator( float *dataCArray, float *qxCArray, float *qyCArray, int16_t *maskCArray, int arraylength, int nq, int nphi );
	
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
	void calculatePolarCoordinates( double start_q, double stop_q ); // for cheetah.ini control
	void calculateSAXS();
	void calculateXCCA();
	void calculateXACA();
	
    //---------------------------------------------alternative approach (Jan's way)
    // these functions have the byname _FAST to distinguish them from the ones above 
    // (for lack of a better name and in the hope that they may be fast. We'll see...)
    
    // 'calculatePolarCoordinates' returns a 2D pattern in polar coordinates (r vs. phi)
    int calculatePolarCoordinates_FAST(array2D *&polar, int number_q=20, int number_phi=128); 
    int calculatePolarCoordinates_FAST(array2D *&polar, int number_q, double start_q, double stop_q,
                                                        int number_phi=128, double start_phi=0, double stop_phi=360 );

    
    // 'calculateXCCA' returns the autocorrelation function (r = const.)          
    int calculateXCCA_FAST( array2D *&polar, array2D *&corr, int writeToInternalDataStructure=1 );
    
    // looks up the value closest to xcoord, ycoord in the data
    double lookup( double xcoord, double ycoord );
	
	// sets the lookup table necessary for lookup() to work. this needs to be provided externally, (a default one is available)
	void setLookupTable( array2D *LUT );
	void setLookupTable( const int *cLUT, unsigned int LUT_dim1, unsigned int LUT_dim2 );
	
	// sets the two plans needed for forward and backward Fourier Transform --> correlation
	void setFFTWPlans( fftw_plan forwardplan, fftw_plan backwardplan );

	// a lookup table can be created within this object, if none was available externally
	// (for cheetah, performance is better if this is calculated once in advance 
	// and then handed to the CrossCorrelator for each shot)
	int createLookupTable(int Nx, int Ny);
	void calcLUTvariables( int lutNx, int lutNy, double &qx_min, double &qx_stepsize, double &qy_min, double &qy_stepsize );
    
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
	void setQmin( double qmin_val );
	double phimin() const;
	void setPhimin( double phimin_val );
	double phimax() const;
	void setPhimax( double phimax_val );
	
	// jas: qmaxCArray() could easily be rewritten to use array1D qx, qy instead of CArrays if preferable
	double qmax2CArray( float *qxCArray, float *qyCArray, int arraylength ); // calculates qmax from 2 CArrays
    double qmax1CArray( float *qCArray, int arraylength ); // calculates qmax from 1 CArray
	
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
