/*
 *  ArrayClasses.h
 *  xcca
 *
 *  Created by Feldkamp on 2/17/11.
 *  Last changed on 04/27/11.
 *  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
 *
 */

#ifndef _arrayclasses_h
#define _arrayclasses_h

#include <string>

#include <fftw3.h>

//tell the compiler that these classes are going to be defined
class array1D;
class array2D;
class array3D;
class FourierTransformer;



//=================================================================================
//
// arraydata (the base class for array1D, array2D, array3D)
//
//=================================================================================
class arraydata {
protected:
	double *p_data;												//pointer to the data array
	unsigned int p_size;	
	
public:
//	arraydata();                                                //default constructor
	arraydata( const unsigned int sizeval = 1 );                     
    arraydata( const int16_t *CArray, const unsigned int size_val );     //initialize with C-array of ints of a given length
	arraydata( const float *CArray, const unsigned int size_val );			//initialize with C-array of floats of a given length
    arraydata( const arraydata &src );                          //copy constructor
    arraydata & operator=(const arraydata & src);               //assignment operator
	~arraydata();

    //helper functions
    void init();
    void copy( const double* src_data, const unsigned int size_val );
    void copy( const int* src_data, const unsigned int size_val );
    void copy( const arraydata& src );
    void copy( const array1D& src );
    void destroy();
    double *data() const;                                         //return pointer to internal raw data array
    
    void zero();												//set all elements to zero
	void ones();												//set all elements to 1
    
	// 'atAbsoluteIndex' functions:
	// the following functions change or return properties
	// assuming a one-dimensional array (consistent with the internal data structure), 
    // no matter what dimension the actual subclass may have
	unsigned int size() const;										//total size of the array	
	double get_atAbsoluteIndex( unsigned int index) const;					// get element value
	void set_atAbsoluteIndex( unsigned int index, double val);			// set element value

	double calcMin() const;
	double calcMax() const;
	
    std::string getASCIIdataAsRow() const;                         //can/should be overridden by subclasses
    std::string getASCIIdataAsColumn() const;
    
	void readFromRawBinary( std::string filename );
	void writeToRawBinary( std::string filename );
    
    //perform some basic math on array
    int multiplyByFactor( double factor );
    int multiplyByArrayElementwise( const arraydata *secondFactor );

//	int verbose() const;
//	void setVerbose( int verbosity );
};




//=================================================================================
//
// array1D
//
//=================================================================================
class array1D : public arraydata{
	
private:
	int p_dim1;
	
public:
//    array1D();                                                  //default constructor
	array1D( unsigned int size_dim1 = 1);                          
    array1D( int16_t *CArray, unsigned int size_val );          //init with C-style array of ints
	array1D( float *CArray, unsigned int size_val );			//init with C-style array of floats
    array1D( array2D* dataTwoD );                               //init with an array2D object
	~array1D();
    
    void copy( const array1D& src );
	
	unsigned int dim1() const;
	void setDim1( unsigned int size_dim1 );
	
	double get( unsigned int i ) const;
	void set( unsigned int i, double value );
    
                                   
    
    std::string getASCIIdata() const;
    int writeToASCII( std::string filename ) const;
};




//=================================================================================
//
// array2D
//
//=================================================================================
class array2D : public arraydata{
	
private:
	int p_dim1;
	int p_dim2;
	
public:
//    array2D();
	array2D( unsigned int size_dim1 = 1, unsigned int size_dim2 = 1 );              //default constructor
    array2D( array1D* dataOneD, unsigned int size_dim1, unsigned int size_dim2);   // use 1D data to initialize
	~array2D();

    void copy( const array2D& src );
        
	unsigned int dim1() const;
	void setDim1( unsigned int size_dim1 );
	unsigned int dim2() const;
	void setDim2( unsigned int size_dim2 );
	
	double get( unsigned int i, unsigned int j ) const;                 //returns single pixel value
	void set( unsigned int i, unsigned int j, double value );
    
    int getRow( int rownum, array1D *&row ) const; 	                    //returns one-dimensional 'row' or 'col'
    int getCol( int colnum, array1D *&col ) const;						//return value is 0 if successful
    void setRow( int rownum, const array1D *row );                              //sets a one-dimensional row
    void setCol( int colnum, const array1D *col );
            
	void readFromHDF5( std::string filename );
	
	int writeToTiff( std::string filename, int scaleFlag = 0 ) const;     //needs libtiff
	int writeToHDF5( std::string filename ) const;
    
    std::string getASCIIdata() const;
	int writeToASCII( std::string filename ) const;
	
	void generateTestPattern( int type );				//for debugging
};




//=================================================================================
//
// array3D
//
//=================================================================================
class array3D : public arraydata{
	
private:
	int p_dim1;
	int p_dim2;
	int p_dim3;
	
public:
//    array3D();
	array3D( unsigned int size_dim1 = 1, unsigned int size_dim2 = 1, unsigned int size_dim3 = 1 );  //default constructor
	~array3D();
    
    void copy( const array3D& src );
	
	unsigned int dim1() const;
	void setDim1( unsigned int size_dim1 );
	unsigned int dim2() const;
	void setDim2( unsigned int size_dim2 );
	unsigned int dim3() const;
	void setDim3( unsigned int size_dim3 );
	

	double get( unsigned int i, unsigned int j, unsigned int k ) const;
	void set( unsigned int i, unsigned int j, unsigned int k, double value );

    std::string getASCIIdata() const;
    int writeToASCII( std::string filename ) const;
};





//=================================================================================
//
// FourierTransformer. helper class. it does the Fourier transform
//
//=================================================================================
class FourierTransformer{

private:
    int verbose;
	fftw_complex *p_in;				//internal complex input array
	fftw_complex *p_out;			//internal complex output array
	fftw_plan p_forward_plan;
	fftw_plan p_backward_plan;
	int p_n;						//size of the input (and output) arrays
	
	void createPlans();
	void destroyPlans();	
	
public:
    FourierTransformer( const array1D *real, const array1D *imag );	//initialize with data
    ~FourierTransformer();	

    //wrapper for the FFTW discrete Fourier Transform
    int transform( int direction=1 );
	int transformWithNewPlans( int direction=1 );

	//after the transform, use these functions to ask for the transformed data
	void getData( array1D *&real, array1D *&imag ) const;		// return within passed arguments
	array1D getReal() const;									// return by copy (may be slower)
	array1D getImag() const;									// return by copy (may be slower)
};


#endif
