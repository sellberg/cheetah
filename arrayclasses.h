/*
 *  ArrayClasses.h
 *  xcca
 *
 *  Created by Feldkamp on 2/17/11.
 *  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
 *
 */

#ifndef _arrayclasses_h
#define _arrayclasses_h

#include <string>

class arraydata {
protected:
	double *p_data;												//pointer to the data array
	unsigned int p_size;	
	int p_verbose;												//1 = talk, 0 = silent


	
public:
	arraydata();
	arraydata( unsigned int sizeval );
	~arraydata();
	
	void zero();												//set all elements to zero
	void ones();												//set all elements to 1

	// 'atAbsoluteIndex' functions:
	// the following functions change or return properties
	// assuming a 1D array, no matter what dimension the actual subclass may have
	unsigned int size_absolute();										//total size of the array	
	double get_atAbsoluteIndex( unsigned int index);					// get element value
	void set_atAbsoluteIndex( unsigned int index, double val);			// set element value

	double getMin();
	double getMax();
	
	void readFromRawBinary( std::string filename );
	void writeToRawBinary( std::string filename );

	int verbose();
	void setVerbose( int verbosity );
};



class array1D : public arraydata{
	
private:
	int p_dim1;
	
public:
	array1D();
	array1D( unsigned int size_dim1 );
	~array1D();
	
	unsigned int dim1();
	void setDim1( unsigned int size_dim1 );
	
	double get( unsigned int i );
	void set( unsigned int i, double value );
};



class array2D : public arraydata{
	
private:
	int p_dim1;
	int p_dim2;
	
public:
	array2D();
	array2D( unsigned int size_dim1, unsigned int size_dim2 );
	~array2D();
	
	unsigned int dim1();
	void setDim1( unsigned int size_dim1 );
	unsigned int dim2();
	void setDim2( unsigned int size_dim2 );
	
	double get( unsigned int i, unsigned int j );
	void set( unsigned int i, unsigned int j, double value );
	
	void readFromHDF5( std::string filename );
	
//	int writeToTiff( std::string filename, int scaleFlag = 0 );     //commented out for now (needs libtiff)
	int writeToHDF5( std::string filename );
	int writeToASCII( std::string filename, bool printToConsole = false );
	
	void generateTestPattern( int type );				//for debugging
};




class array3D : public arraydata{
	
private:
	int p_dim1;
	int p_dim2;
	int p_dim3;
	
public:
	array3D();
	array3D( unsigned int size_dim1, unsigned int size_dim2, unsigned int size_dim3 );
	~array3D();
	
	unsigned int dim1();
	void setDim1( unsigned int size_dim1 );
	unsigned int dim2();
	void setDim2( unsigned int size_dim2 );
	unsigned int dim3();
	void setDim3( unsigned int size_dim3 );
	

	double get( unsigned int i, unsigned int j, unsigned int k );
	void set( unsigned int i, unsigned int j, unsigned int k, double value );
};


#endif
