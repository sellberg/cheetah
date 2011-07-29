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

	//-------------------------------------------------required for both algorithms
	int p_arraySize;
	bool p_mask_enable;			// enables/disables masking of bad pixels
	bool p_xcca_enable;			// enables/disables cross-correlations (if disabled only autocorrelations are calculated)
    std::string p_outputdir;  	// the output directory if anything is dumped from withing this class (default is working dir)

	array1D *p_data;			//data storage
	array1D *p_qx;				//pixel x coordinate
	array1D *p_qy;				//pixel y coordinate
	array1D *p_mask;			//mask used to remove bad pixels

	array2D *p_autoCorrelation;
	array3D *p_crossCorrelation;
	
	int p_debug;		
	
	//-------------------------------------------------required for alg1
	double p_qmin;
	double p_qmax;
	double p_deltaq;			// step size (bin length) in q-direction
	double p_phimin;
	double p_phimax;
	double p_deltaphi;			// step size (bin length) in phi direction
	int p_samplingLength;
	int p_samplingAngle;
	int p_samplingLag;

	array1D *p_q;				//magnitude of q vector (1st dimension in correlation)
	array1D *p_phi;				//angle (2nd dimension in correlation)
	
	array1D *p_qave;
	array1D *p_iave;
	array1D *p_phiave;
	void updateDependentVariables();
    
    
	
	//-------------------------------------------------required for alg2
	array2D *p_polar;
	array2D *p_corr;
	array2D *p_mask_polar;
	array2D *p_mask_corr;
	array2D *p_table;			//lookup table
	fftw_plan p_fplan;			//forward FFT plan
	fftw_plan p_bplan;			//backward FFT plan

	double p_qxmax;
	double p_qymax;	
	double p_qxmin;
	double p_qymin;
	double p_qxdelta;
	double p_qydelta;

	float *p_dataCArray; 	//!!!only temporary!!! keep the cheetah c-array to be able to write to it for debugging purposes
	//array1D *check1D;    

public:
	//---------------------------------------------constructors & destructor
    CrossCorrelator( int arraylength=1 );                                                           //init with 1D data size
	CrossCorrelator( float *dataCArray, int arraylength );                                        //init with actual data
    CrossCorrelator( float *dataCArray, float *qxCArray, float *qyCArray, int arraylength );      //init with actual data + centered(!) q calibration
	CrossCorrelator( float *dataCArray, float *qxCArray, float *qyCArray, int arraylength, double qMax, double qMin );
	CrossCorrelator( float *dataCArray, float *qxCArray, float *qyCArray, int arraylength, int nq, int nphi ); // autocorrelation only
	CrossCorrelator( float *dataCArray, float *qxCArray, float *qyCArray, int arraylength, int nq1, int nq2, int nphi ); // full xcca
	CrossCorrelator( float *dataCArray, float *qxCArray, float *qyCArray, int16_t *maskCArray, int arraylength, int nq, int nphi ); // autocorrelation only
	CrossCorrelator( float *dataCArray, float *qxCArray, float *qyCArray, int16_t *maskCArray, int arraylength, int nq1, int nq2, int nphi ); // full xcca
	CrossCorrelator( array2D *dataArray, array2D *qxArray, array2D *qyArray, int nq, int nphi ); // autocorrelation only
	
	~CrossCorrelator();
	
	//---------------------------------------------input/output
    void initPrivateVariables();
    void initDefaultQ();
	void initFromFile( std::string filename, int type=0 );
    void initWithTestPattern( int sizex, int sizey, int type=0 );                           //generate some test patterns
	void printRawData(uint16_t *buffer,long lSize);
//	void dumpResults( std::string filename );				// (arraydataIO needed for this)
	
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
	int calculatePolarCoordinates_FAST(int number_phi=512, int number_q=20 );
    int calculatePolarCoordinates_FAST(int number_phi, int number_q, double start_q, double stop_q );
	
	int calculatePolarCoordinates_FAST( array1D* image, array2D* polar2D, 
										int number_phi, int number_q, double start_q, double stop_q ); 

    
    // 'calculateXCCA' returns the autocorrelation function (r = const.)          
    int calculateXCCA_FAST();
	
    // looks up the value closest to xcoord, ycoord in the data
    double lookup( double xcoord, double ycoord, array1D *dataArray );
	
	// sets the lookup table necessary for lookup() to work. this needs to be provided externally, (a default one is available)
	array2D *lookupTable();
	void setLookupTable( array2D *LUT );
	void setLookupTable( const int *cLUT, unsigned int LUT_dim1, unsigned int LUT_dim2 );
	
	// sets the two plans needed for forward and backward Fourier Transform --> correlation
	void setFFTWPlans( fftw_plan forwardplan, fftw_plan backwardplan );

	// a lookup table can be created within this object, if none was available externally
	// (for cheetah, performance is better if this is calculated once in advance 
	// and then handed to the CrossCorrelator for each shot)
	int createLookupTable( int Nx, int Ny );
	void calcLUTvariables( int lutNx, int lutNy );
    
	//compute correlations in a 2D polar coordinate matrix
	int autocorrelateFFT(array2D *polar2D, array2D *corr2D);
	int crosscorrelateFFT(array2D *polar2D, array3D *corr3D);
	
	int normalizeToSAXS();
	
    //compute 1D correlations using FFT, the result is returned in f, respectively
//    int correlateFFT( array1D *f, array1D *g );    // result of corr(f,g) -> f
//    int autocorrelateFFT( array1D *f );            // result of corr(f,f) -> f
//	int autocorrelateFFTwithMask( array1D *f, array1D *m );
        
	//---------------------------------------------setters & getters
	array1D *data() const;
	void setData( array1D *data );
	void setData( array2D *data );
	void setData( float *dataCArray, unsigned int size );
	
	array1D *qx() const;
	void setQx( array1D *qx );
	void setQx( array2D *qx );
	void setQx( float *qxCArray, unsigned int size );

	array1D *qy() const;
	void setQy( array1D *qy );
	void setQy( array2D *qy );
	void setQy( float *qyCArray, unsigned int size );

	array1D *mask() const;
	void setMask( array1D *mask );
	void setMask( array2D *mask );	
	void setMask( int16_t *maskCArray, unsigned int size );
	void normalizeMask();
	
	array2D *polar() const;
	array2D *corr() const;
	array2D *mask_polar() const;
	array2D *mask_corr() const;
	
	
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
	double getCrossCorrelation(unsigned index1, unsigned index2) const;
	double getCrossCorrelation(unsigned index1, unsigned index2, unsigned index3) const;
};





#endif



// call the functions below in a reasonable way (less fine-grained control, but less to worry about)
//int performAnalysis_FAST();
//    int calculatePolarCoordinates_FAST(array2D *&polar, int number_q=20, int number_phi=128); 
//    int calculatePolarCoordinates_FAST(array2D *&polar, int number_q, double start_q, double stop_q,
//                                                        int number_phi=128, double start_phi=0, double stop_phi=360 );
//    int calculateXCCA_FAST( array2D *&polar, array2D *&corr, int writeToInternalDataStructure=1 );

