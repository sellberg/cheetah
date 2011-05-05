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

//tell the compiler that these classes are going to be defined
class array1D;
class array2D;
class array3D;
class FourierTransformer;

//*********************************************************************************
class arraydata {
protected:
	double *p_data;												//pointer to the data array
	unsigned int p_size;	
	int p_verbose;												//1 = talk, 0 = silent


	
public:
	arraydata();                                                //default constructor
	arraydata( unsigned int sizeval );                     
    arraydata( int16_t *CArray, unsigned int size_val );     //initialize with C-array of ints of a given length
	arraydata( float *CArray, unsigned int size_val );			//initialize with C-array of floats of a given length
    arraydata( const arraydata &src );                          //copy constructor
    arraydata & operator=(const arraydata & src);               //assignment operator
	~arraydata();

    //helper functions
    void init();
    void copy( const double* src_data, int arraysize );
    void copy( const arraydata& src );
    void copy( const array1D& src );
    void destroy();
    double *data() const;                                         //return pointer to internal raw data array
//    void data_copy(double *copy);                               //return a new array pointer with a copy of the data
    
    void zero();												//set all elements to zero
	void ones();												//set all elements to 1
    
	// 'atAbsoluteIndex' functions:
	// the following functions change or return properties
	// assuming a one-dimensional array (consistent with the internal data structure), 
    // no matter what dimension the actual subclass may have
	unsigned int size() const;										//total size of the array	
	double get_atAbsoluteIndex( unsigned int index) const;					// get element value
	void set_atAbsoluteIndex( unsigned int index, double val);			// set element value

	double getMin() const;
	double getMax() const;
	
    std::string getASCIIdataAsRow() const;                         //can/should be overridden by subclasses
    std::string getASCIIdataAsColumn() const;
    
	void readFromRawBinary( std::string filename );
	void writeToRawBinary( std::string filename );
    
    //perform some basic math on array
    int multiplyByFactor( double factor );
    int multiplyByArrayElementwise( arraydata *secondFactor );

	int verbose() const;
	void setVerbose( int verbosity );
};




//*********************************************************************************
class array1D : public arraydata{
	
private:
	int p_dim1;
	
public:
    array1D();                                                  //default constructor
	array1D( unsigned int size_dim1 );                          
    array1D( int16_t *CArray, unsigned int size_val );          //init with C-style array of ints
	array1D( float *CArray, unsigned int size_val );			//init with C-style array of floats
//    array1D( const arraydata &src );                            //copy constructor
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




//*********************************************************************************
class array2D : public arraydata{
	
private:
	int p_dim1;
	int p_dim2;
	
public:
    array2D();
	array2D( unsigned int size_dim1, unsigned int size_dim2 );              //default constructor
    array2D( array1D* dataOneD, unsigned int size_dim1, unsigned int size_dim2);   // use 1D data to initialize
	~array2D();

    void copy( const array2D& src );
        
	unsigned int dim1() const;
	void setDim1( unsigned int size_dim1 );
	unsigned int dim2() const;
	void setDim2( unsigned int size_dim2 );
	
	double get( unsigned int i, unsigned int j ) const;                   //returns single pixel value
	void set( unsigned int i, unsigned int j, double value );
    
    void getRow( int rownum, array1D *row ) const;                        //returns one-dimensional row
    void getCol( int colnum, array1D *col ) const;
        
	void readFromHDF5( std::string filename );
	
	int writeToTiff( std::string filename, int scaleFlag = 0 ) const;     //needs libtiff
	int writeToHDF5( std::string filename ) const;
    
    std::string getASCIIdata() const;
	int writeToASCII( std::string filename ) const;
	
	void generateTestPattern( int type );				//for debugging
};



//*********************************************************************************
class array3D : public arraydata{
	
private:
	int p_dim1;
	int p_dim2;
	int p_dim3;
	
public:
    array3D();
	array3D( unsigned int size_dim1, unsigned int size_dim2, unsigned int size_dim3 );  //default constructor
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






//*********************************************************************************
//helpers. they do Fourier transform
class FourierTransformer{

private:
    int verbose;

public:
    FourierTransformer();

    //wrapper for the FFTW discrete Fourier Transform
    //the two arrays 'real' and 'imag' are overwritten by
    //resulting arrays for the real and imaginary parts
    int transform( array1D *real, array1D *imag, int direction=1 );
};


#endif
