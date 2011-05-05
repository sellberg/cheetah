/*
 *  ArrayClasses.cpp
 *  xcca
 *
 *  Created by Feldkamp on 2/17/11.
 *  Last changed on 04/27/11.
 *  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
 *
 */

#include "arrayclasses.h"

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include <string>
using std::string;

#include <sstream>
using std::ostringstream;

#include <fstream>
using std::ofstream;

#include <cmath>
#include <complex>              // needed to let fftw use regular doubles

//include headers
#include <hdf5.h>
#include <tiffio.h>
#include <fftw3.h>



//*********************************************************************************
//*********************************************************************************
//
// CLASS IMPLEMENTATION OF arraydata
//
//*********************************************************************************
//*********************************************************************************
arraydata::arraydata(){
    init();
}

arraydata::arraydata( unsigned int size_val ){
    init();
	p_size = size_val;
	p_data = new double[size_val];
    if (p_data == 0)
        cerr << "Error in arraydata constructor: could not allocate memory." << endl;

	zero();					// set all elements to zero initially
	setVerbose( 1 );		// set talkativity
}

arraydata::arraydata( int16_t *CArray, unsigned int size_val ){
    init();
    p_size = size_val;
    p_data = new double[size_val];
    if (p_data == 0)
        cerr << "Error in arraydata constructor: could not allocate memory." << endl;
    
    //fill array with data, convert int to double before
    for (int i = 0; i < size_val; i++) {
        p_data[i] = (double)CArray[i];
    }
}

arraydata::arraydata( float *CArray, unsigned int size_val ) {
    init();
    p_size = size_val;
    p_data = new double[size_val];
    if (p_data == 0)
        cerr << "Error in arraydata constructor: could not allocate memory." << endl;
    
    //fill array with data, convert float to double before
    for (int i = 0; i < size_val; i++) {
        p_data[i] = (double)CArray[i];
    }
}

arraydata::arraydata( const arraydata &src ){                       //copy constructor
    init();
    copy( src );
}

arraydata & arraydata::operator=(const arraydata & src){
    if ( this != &src ){
        destroy();
        init();
        copy( src );
    }
    return *this;
}

arraydata::~arraydata(){
	destroy();	
}

//--------------------------------------------------------------------init helper functions
void arraydata::init(){
    p_size = 0;
    p_data = NULL;
}

void arraydata::copy( const double* src_data, int arraysize ){
    p_size = arraysize;
    if (p_size > 0){
        p_data = new double[ arraysize ];
        for (int i = 0; i < p_size; i++) {
            p_data[i] = src_data[i];
        }
    }else{
        p_data = NULL;
    }
}


void arraydata::copy( const arraydata& src ){
    this->copy( src.data(), src.size());
}



void arraydata::destroy(){
    delete []p_data;
}

double *arraydata::data() const{
    return (double *)p_data;
}



//--------------------------------------------------------------------
double arraydata::get_atAbsoluteIndex( unsigned int i) const{
	if (i < size()) {
		return p_data[i];
	} else {
		if (verbose()){
			std::cerr << "Error in arraydata! Index " << i 
			<< " is larger than absolute dimension " << size() << "." << endl;
		}
		
		//fail 'silently', 
		//i.e., return a valid numeric value 0,
		//instead of exiting or throwing an exception
		return 0.;
	}

}


//--------------------------------------------------------------------
void arraydata::set_atAbsoluteIndex( unsigned int i, double val){
	if (i < size()) {
		p_data[i] = val;
	} else {									//fail 'silently', as in: don't do anything
		if (verbose()){
			std::cerr << "Error in arraydata! Index " << i 
			<< " is larger than absolute dimension " << size() << "." << endl;
		}
	}
}


//--------------------------------------------------------------------
unsigned int arraydata::size() const{
	return p_size;
}


//--------------------------------------------------------------------
void arraydata::zero(){				//set all elements to zero
	for (int i = 0; i < size(); i++) {
		set_atAbsoluteIndex(i, 0);
	}	
}

//--------------------------------------------------------------------
double arraydata::getMin() const{
	double tempmin = INFINITY;
	for (int i = 0; i < size(); i++) {
		if (get_atAbsoluteIndex(i) < tempmin) {
			tempmin = get_atAbsoluteIndex(i);
		}
	}
	return tempmin;
}

//--------------------------------------------------------------------
double arraydata::getMax() const{
	double tempmax = -INFINITY;
	for (int i = 0; i < size(); i++) {
		if (get_atAbsoluteIndex(i) > tempmax) {
			tempmax = get_atAbsoluteIndex(i);
		}
	}
	return tempmax;
}


//------------------------------------------------------------- getASCIIdata
string arraydata::getASCIIdataAsRow() const{
    ostringstream osst;
	for (int i = 0; i<size(); i++) {
        osst << get_atAbsoluteIndex(i) << " ";
	}
    osst << endl;
    return osst.str();
}

//------------------------------------------------------------- getASCIIdata
string arraydata::getASCIIdataAsColumn() const{
    ostringstream osst;
	for (int i = 0; i<size(); i++) {
        osst << get_atAbsoluteIndex(i) << endl;
	}
    return osst.str();
}

//--------------------------------------------------------------------arraydata::readFromRawBinary
void arraydata::readFromRawBinary( std::string filename ){
	if ( verbose() ) 
		cout << "reading from raw binary file " << filename << endl;
	
	FILE *filePointerRead;
	long lSize;
	uint16_t *buffer;
	size_t result;
	
	// fiducials for test shots: f909 (black), 13a19 (powder)
	filePointerRead = fopen( filename.c_str(),"r+");
	if (filePointerRead == NULL) 
	{
		fputs ("Error in arraydata::readFromRawBinary. File error\n",stderr); 
		exit (1);
	}
	
	// obtain file size:
	fseek(filePointerRead, 0, SEEK_END);
	lSize = ftell(filePointerRead);
	rewind(filePointerRead);
	lSize /= sizeof(uint16_t);
	
	// allocate memory to contain the whole file:
	buffer = (uint16_t*) malloc(sizeof(uint16_t)*lSize);
	if (buffer == NULL) 
	{
		fputs ("Error in arraydata::readFromRawBinary. Memory error\n",stderr); 
		exit (2);
	}
	
	// copy the file into the buffer:
	result = fread(buffer,sizeof(uint16_t),lSize,filePointerRead);
	if (result != lSize) 
	{
		fputs ("Error in arraydata::readFromRawBinary. Reading error\n",stderr); 
		exit (3);
	}
	
	/* the whole file is now loaded in the memory buffer. */
	printf("data contains %ld unsigned 16-bit integers\n",lSize);
	
	// printRawData(buffer,lSize);
	
	//write read data to array for further processing
	for (int i=0; i<size(); i++) {
		set_atAbsoluteIndex( i, buffer[i] );
	}
	
	
	// terminate reading
	fclose(filePointerRead);
	free(buffer);
}

//--------------------------------------------------------------------array math

//multiply each element by a numerical factor
int arraydata::multiplyByFactor( double factor ){                       
    for (int i = 0; i<this->size(); i++) {
        set_atAbsoluteIndex(i, this->get_atAbsoluteIndex(i)*factor);
    }
    return 0;
}

//multiply each element by an element from a second vector
int arraydata::multiplyByArrayElementwise( arraydata *secondFactor ){
    if (this->size() != secondFactor->size()){
        cerr << "Error in arraydata::multiplyArrayElementwise! Array sizes don't match. ";
        cerr << "(" << this->size() << " != " << secondFactor->size() << "). Operation not performed."<< endl;
        return 1;
    }
    
    for (int i = 0; i<this->size(); i++) {
        set_atAbsoluteIndex(i, this->get_atAbsoluteIndex(i)*secondFactor->get_atAbsoluteIndex(i));
    }
    return 0;
}


//--------------------------------------------------------------------
int arraydata::verbose() const{
	return p_verbose;
}

void arraydata::setVerbose( int verbosity ){
	p_verbose = verbosity;
}











//*********************************************************************************
//*********************************************************************************
//
// CLASS IMPLEMENTATION OF array1D
//
//*********************************************************************************
//*********************************************************************************


//-----------------------------------------------------constructors & destructors
array1D::array1D() 
		: arraydata(){
    setDim1( 0 );
}

array1D::array1D( unsigned int size_dim1 ) 
		: arraydata(size_dim1){
    setDim1( size_dim1 );
}

array1D::array1D( int16_t *CArray, unsigned int size_dim1 )
        : arraydata( CArray, size_dim1 ){
    setDim1( size_dim1 );
}

array1D::array1D( float *CArray, unsigned int size_dim1 )
		: arraydata( CArray, size_dim1 ) {
    setDim1( size_dim1 );
}

//copy constructor
//array1D::array1D( const arraydata &src ){
//}

//constructor to generate a 1D array from a 2D array
array1D::array1D( array2D* dataTwoD ) 
        : arraydata( dataTwoD->size() ){
 	setDim1( dataTwoD->size() );

    //copy contents of dataTwoD to this array1D object
    p_size = dataTwoD->size();
    if (p_size > 0){
        p_data = new double[ dataTwoD->size() ];
        for (int i = 0; i < p_size; i++) {
            p_data[i] = dataTwoD->get_atAbsoluteIndex(i);
        }
    }else{
        p_data = NULL;
    }
}

array1D::~array1D(){
}


//-----------------------------------------------------copy
void array1D::copy( const array1D& src ){
    setDim1( src.size() );
    this->arraydata::copy( src.data(), src.size());
}



//-----------------------------------------------------data accessors
double array1D::get( unsigned int i ) const{
	if ( i < dim1() ) {
		return get_atAbsoluteIndex(i);		
	} else {
		if (verbose()){
			std::cerr << "Error in array1D! Index " << i 
				<< " is larger than dimension " << dim1() << "." << endl;
		}
		return 0.;
	}
}

void array1D::set( unsigned int i, double value ){
	if ( i < dim1() ) {
		set_atAbsoluteIndex(i, value);
	} else {
		if (verbose()){
			std::cerr << "Error in array1D! Index " << i 
				<< " is larger than dimension " << dim1() << "." << endl;
		}
	}


}


//-----------------------------------------------------setters & getters
unsigned int array1D::dim1() const{
	return p_dim1;
}

void array1D::setDim1( unsigned int size_dim1 ){
	p_dim1 = size_dim1;
}


//------------------------------------------------------------- getASCIIdata
string array1D::getASCIIdata() const{
    ostringstream osst;
    osst << "1D data, dim1=" << dim1() << ", size=" << size() << endl;
	osst << " [";
    for (int i = 0; i<dim1(); i++) {
        osst << " " << get(i);
    }
    osst << "]" << endl;
    return osst.str();
}


//------------------------------------------------------------- writeToASCII
int array1D::writeToASCII( std::string filename ) const{
	ofstream fout( filename.c_str() );
	fout << arraydata::getASCIIdataAsColumn();
	fout.close();
	return 0;
}



//*********************************************************************************
//*********************************************************************************
//
// CLASS IMPLEMENTATION OF array2D
//
//*********************************************************************************
//*********************************************************************************

//-----------------------------------------------------constructors & destructors
array2D::array2D()
		: arraydata(){
 	setDim1( 0 );
	setDim2( 0 );
}

array2D::array2D( unsigned int size_dim1, unsigned int size_dim2 )
		: arraydata(size_dim1*size_dim2){
 	setDim1( size_dim1 );
	setDim2( size_dim2 );
}


//constructor to generate a 2D array from a 1D array, given the desired dimensions
array2D::array2D( array1D* dataOneD, unsigned int size_dim1, unsigned int size_dim2) 
        : arraydata( size_dim1*size_dim2 ){
    if (dataOneD->size() != size_dim1*size_dim2) {
        cerr << "Warning in array2D::array2D. Inconsistent array size. ";
        cerr << "size1D=" << dataOneD->size() << ", size2D=" << size_dim1*size_dim2 << endl;
    }
 	setDim1( size_dim1 );
	setDim2( size_dim2 );

    //copy contents of dataOneD to this array2D
    p_size = dataOneD->size();
    if (p_size > 0){
        p_data = new double[ dataOneD->size() ];
        for (int i = 0; i < p_size; i++) {
            p_data[i] = dataOneD->get_atAbsoluteIndex(i);
        }
    }else{
        p_data = NULL;
    }
}


array2D::~array2D(){
}


//-----------------------------------------------------copy
void array2D::copy( const array2D& src ){
    setDim1( src.dim1() );
    setDim2( src.dim2() );
    this->arraydata::copy( src.data(), src.size());
}

//-----------------------------------------------------data accessors
double array2D::get( unsigned int i, unsigned int j ) const{
	if ( i < dim1() && j < dim2()) {
		return get_atAbsoluteIndex( j*dim1() + i );
	} else {
		if (verbose()){
			std::cerr << "Error in array2D::get! Index (" << i << ", " << j 
				<< ") is larger than dimension (" << dim1() << ", " << dim2() << ")." << endl;
		}
		return 0.;
	}
}

void array2D::set( unsigned int i, unsigned int j, double value ){
	if ( i < dim1() && j < dim2()) {
		set_atAbsoluteIndex( j*dim1() + i, value);
	} else {
		if (verbose()){
			std::cerr << "Error in array2D::set! Index (" << i << ", " << j 
				<< ") is larger than dimension (" << dim1() << ", " << dim2() << ")." << endl;
		}
	}
}


//-----------------------------------------------------getRow/getCol
void array2D::getRow( int rownum, array1D *row ) const{
    if (row){
        delete row;
        row = NULL;
    }
    row = new array1D( this->dim1() );
    
    //for a fixed row, i goes through columns (x-values)
    for (int i = 0; i < row->size(); i++){
        row->set( i, this->get(i, rownum) );
    }
}

void array2D::getCol( int colnum, array1D *col ) const{
    if (col){
        delete col;
        col = NULL;
    }
    col = new array1D( this->dim2() );
    
    //for a fixed column number, j goes through the rows (y-values)
    for (int j = 0; j < col->size(); j++){
        col->set( j, this->get(colnum, j) );
    }
}
    
//------------------------------------------------------------- getASCIIdata
std::string array2D::getASCIIdata() const{
    ostringstream osst;
    osst << "2D data, dim1=" << dim1() << ", dim2=" << dim2() << ", size=" << size() << endl;
	for (int j = 0; j<dim2(); j++){
		osst << " [";
	    for (int i = 0; i<dim1(); i++) {
			osst << " " << get(i, j);
		}
		osst << "]" << endl;
	}
    return osst.str();
}


//-------------------------------------------------------------- writeToTiff
// (needs TIFFLIB installation)
//
// scaleFlag 0: direct output of values, 
//				cut off values larger than 2^16-1 = 65535
// scaleFlag 1: scaled output to full tiff scale of 0 to 65535
//
// implementation borrowed following matrixdata::tiff16out 
// of the tomo/matlib package (see http://xray-lens.de )
//--------------------------------------------------------------------------
int array2D::writeToTiff( std::string filename, int scaleFlag ) const{

    if (size() == 0) {
        cerr << "Error in writeToTiff! Array size is zero." << endl;
    }
    
    if (!p_data) {
        cerr << "Error in writeToTiff! No data in array." << endl;
    }

    
	TIFF *out= TIFFOpen(filename.c_str() ,"w");
	if(out)
	{
		uint32 width = (uint32) dim1();
		uint32 height = (uint32) dim2();
		
		//int imagesize = i_width * i_height;
		// float xRes, yRes;
		// uint16 res_unit;
		uint16 spp, bpp, photo;
		
		double MaxValue = getMax();
		double MinValue = getMin();
		double MaxMinRange = MaxValue - MinValue;
		
		uint16 *tifdata = new uint16[width];
		
		
		// define tif-parameters
		spp = 1;      // Samples per pixel
		bpp = 16;     // Bits per sample
		photo = PHOTOMETRIC_MINISBLACK;
		
		// set the width of the image
		TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width/spp);  
		//TIFFSetField(out, TIFFTAG_IMAGEWIDTH, image_width * bpp * 8); 
		
		// set the height of the image
		TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);  
		
		// set the size of the channels   
		TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp); 
		
		// set number of channels per pixel
		TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
		
		// set the origin of the image
		TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    
		
		TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
		TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photo);

		
		
		//write to file
		for (uint32 i = 0; i < height; i++)
		{
			for (unsigned int j = 0; j < width; j++)
			{
				if (scaleFlag) {								//scale image to full range
					tifdata[j] = (uint16) floor(65535.* (get(j,i)-MinValue)/MaxMinRange);
				} 
				else {											//direct output, cut off larger values
					if( get(j,i) < 0 ) {
						tifdata[j] = 0;
					}
					else if ( get(j,i) > 65535 ) {
						tifdata[j] = 65535;
					}
					else{ 
						tifdata[j] = (uint16) floor(get(j,i));
					}
				}
			}
			
			TIFFWriteScanline(out, tifdata, i, 0);
		}
		
		
		cout << "'" << filename << "' written to disc!" << endl;
		
		delete tifdata;
		TIFFClose(out);
		return 0;
	}else{
        cerr << "Error in array2D::writeToTiff! Could not open '" << filename << "' for writing!" << endl;
		return 1;	
    }
}


//-------------------------------------------------------------- writeToHDF5
int array2D::writeToHDF5( std::string filename ) const{
		
	std::string dataset_name = "data";
	int NX = dim1();
	int NY = dim2();
	int image_rank = 2;
	
	cout << "dim1: " << dim1() << ", dim2: " << dim2() << ", sizeof(int)=" << sizeof(int) << endl;
	
    hid_t		fileID, datasetID;          //file and dataset handles
	hid_t		groupID;					//group id
    hid_t		dataspaceID;				// handles 
    hsize_t		dimension[2];				// dataset dimensions
    herr_t		status;
    int			data[NX][NY];				// data to write
    int			i, j;
	
	
	
    // data assignment and output buffer initialization.
	for(i = 0; i < dim1(); i++){			//5
		for(j = 0; j < dim2(); j++){		//10
			data[i][j] = int( get(i, j) );
		}
	}

	//observe original data
	for(i = 0; i < dim1(); i++){
		for(j = 0; j < dim2(); j++){
			cout << "org: " << get(i, j) << " ";
		}
		cout << endl;
	}
	
	//observe new data stream to be written
	for(i = 0; i < dim1(); i++){
		for(j = 0; j < dim2(); j++){
			cout << "new: " << data[i][j] << " ";
		}
		cout << endl;
	}
	
    // create file, ANTON_CHECK
	// returns 'file' handler
    fileID = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
    
	// create group "data", ANTON_CHECK
	std::string groupname = "data";
	groupID = H5Gcreate(fileID, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (groupID < 0){
		std::cerr << "Error in writeToHDF5. Couldn't create group " << groupname << endl;
		H5Fclose(fileID);
		return 1;
	}
	
	
	//describe dataspace dimensions, ANTON_CHECK
    dimension[0] = dim1();
    dimension[1] = dim2();
    dataspaceID = H5Screate_simple(image_rank, dimension, NULL);
	
	
	//create dataset, ANTON_CHECK (anton has H5T_STD_I16LE)
    datasetID = H5Dcreate(groupID, dataset_name.c_str(), H5T_STD_I16LE, dataspaceID,
						 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (datasetID < 0) {
		cerr << "Error in writeToHDF5. Couldn't create dataset" << endl;
		H5Fclose(fileID);
	}
	
	
	// write data to the dataset, ANTON_CHECK (anton has H5T_STD_I16LE)
    status = H5Dwrite(datasetID, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	
	
    //close all open structures 
    H5Dclose(datasetID);
	H5Sclose(dataspaceID);
	
	H5Gclose(groupID);
    H5Fclose(fileID);
	
		
	return 0;
}


//------------------------------------------------------------- writeToASCII
int array2D::writeToASCII( std::string filename ) const{
	
	ofstream fout( filename.c_str() );
	fout << arraydata::getASCIIdataAsColumn();
	fout.close();

	return 0;
}


//--------------------------------------------------------------------arraydata::readFromHDF5
void array2D::readFromHDF5( std::string rfilename ){
	if ( verbose() ) 
		cout << "reading from HDF5 file " << rfilename << endl;
	
	hid_t       file;
//    hid_t       space;
    hid_t       dset;          /* Handles */
    herr_t      status;
//    hsize_t     dims[1] = { size() };
	double		*rdata = new double[size()];
	
	string dataset = "data";
	
	
	
	
	/*
	 * Get datatype and dataspace identifiers,  
	 * then query datatype class, order, and size, and 
	 * then query dataspace rank and dimensions.
	 */
	/*
	 datatype  = H5Dget_type(dataset);     // datatype identifier
	 class     = H5Tget_class(datatype);
	 if (class == H5T_INTEGER) printf("Dataset has INTEGER type \n");
	 order     = H5Tget_order(datatype);
	 if (order == H5T_ORDER_LE) printf("Little endian order \n");
	 
	 size  = H5Tget_size(datatype);
	 printf(" Data size is %d \n", size);
	 
	 dataspace = H5Dget_space(dataset);    // dataspace identifier
	 rank      = H5Sget_simple_extent_ndims(dataspace);
	 status_n  = H5Sget_simple_extent_dims(dataspace, dims_out);
	 printf("rank %d, dimensions %d x %d \n", rank, dims_out[0], dims_out[1]);
	 */
	
	
	
	
	// open file and dataset using the default properties
    file = H5Fopen (rfilename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen (file, dataset.c_str(), H5P_DEFAULT);
	
	
    // read the data using the default properties.
    status = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
					  &rdata[0]);
	
    // close and release resources.
    status = H5Dclose (dset);
    status = H5Fclose (file);
	
	//put data into private variable
	if (p_data){
        delete []p_data;
        p_data = 0;
    }
	p_data = rdata;
	
	//set dimension!!!!!
}


//-----------------------------------------------------setters & getters
unsigned int array2D::dim1() const{
	return p_dim1;
}

void array2D::setDim1( unsigned int size_dim1 ){
	p_dim1 = size_dim1;
}

unsigned int array2D::dim2() const{
	return p_dim2;
}

void array2D::setDim2( unsigned int size_dim2 ){
	p_dim2 = size_dim2;
}




//-----------------------------------------------------test pattern
void array2D::generateTestPattern( int type ){
    
    if (dim1()==0 || dim2()==0)
        cout << "Warning in generateTestPattern. One or more dimensions are zero." << endl;


    cout << "Writing test pattern type " << type << ": ";
	switch (type) {
		case 0:											
			{   
                cout << "2D sinusoidal ";
				double amplitude = 20000;
				double periodX = dim1()/1;
				double periodY = dim2()/2;
				for (int i = 0; i < dim1(); i++) {
					for (int j = 0; j < dim2(); j++) {
						set(i, j, amplitude/2*(sin(2*M_PI*i/periodX)*cos(2*M_PI*j/periodY)+1) );
					}
				}
			}
		case 1:											
			{
                cout << "increment by absolute array index ";
				double val = 0;
				for (int i = 0; i < size(); i++) {
					if (val >= 65535) {
						val = 0; 
					}
					set_atAbsoluteIndex(i, val);
					val += 1;
				}
			}
			break;
		case 2:											
        	{
                cout << "2D centro-symmetric sine ";
				double amplitude = 20000;
				double period = dim1()/5;
                double centerX = dim1()/2;
                double centerY = dim2()/2;
                cout << "amp=" << amplitude << ", period=" << period << ", center=(" << centerX << "," << centerY << ")";
				for (int i = 0; i < dim1(); i++) {
					for (int j = 0; j < dim2(); j++) {
                        double x = i - centerX;
                        double y = j - centerY;
                        double r = sqrt( x*x + y*y );
						set(i, j, amplitude/2*( sin(2*M_PI*r/period)+1 ) );
					}
				}
			}
		case 3:											
        	{
                cout << "2D centro-symmetric sine with circular modulation ";
				double amplitude = 20000;
				double period = dim1()/5;
                double periodMod = dim1()/10;
                double centerX = dim1()/2;
                double centerY = dim2()/2;
                cout << "amp=" << amplitude << ", period=" << period << ", center=(" << centerX << "," << centerY << ")";
				for (int i = 0; i < dim1(); i++) {
					for (int j = 0; j < dim2(); j++) {
                        double x = i - centerX;
                        double y = j - centerY;
                        double r = sqrt( x*x + y*y );
                        double phi = atan(y/x);
						set(i, j, amplitude/2*( sin(2*M_PI*r/period)+1 + 1/10*(tan(2*M_PI*phi/periodMod)+1) ) );
					}
				}
			}
		case 4:											
        	{
                cout << "2D centro-symmetric sine with straight modulation ";
				double amplitude = 20000;
				double period = dim1()/5;
                double periodMod = dim1()/10;
                double centerX = dim1()/2;
                double centerY = dim2()/2;
                cout << "amp=" << amplitude << ", period=" << period << ", center=(" << centerX << "," << centerY << ")";
				for (int i = 0; i < dim1(); i++) {
					for (int j = 0; j < dim2(); j++) {
                        double x = i - centerX;
                        double y = j - centerY;
                        double r = sqrt( x*x + y*y );
						set(i, j, amplitude/2*( sin(2*M_PI*r/period)+1 + cos(2*M_PI*i/periodMod)+1 ) );
					}
				}
			}
            break;
		default:                                        
            cout << "zeros";
			zero();
			break;
	}
    cout << endl;
}




//*********************************************************************************
//*********************************************************************************
//
// CLASS IMPLEMENTATION OF array3D
//
//*********************************************************************************
//*********************************************************************************

//-----------------------------------------------------constructors & destructors
array3D::array3D()
		: arraydata(){
 	setDim1( 0 );
	setDim2( 0 );
	setDim3( 0 );
}

array3D::array3D( unsigned int size_dim1, unsigned int size_dim2, unsigned int size_dim3 )
		: arraydata(size_dim1*size_dim2*size_dim3){
 	setDim1( size_dim1 );
	setDim2( size_dim2 );
	setDim3( size_dim3 );
}

array3D::~array3D(){
}



//-----------------------------------------------------copy
void array3D::copy( const array3D& src ){
    setDim1( src.dim1() );
    setDim2( src.dim2() );
    setDim3( src.dim3() );
    this->arraydata::copy( src.data(), src.size());
}


//-----------------------------------------------------data accessors
double array3D::get( unsigned int i, unsigned int j, unsigned int k ) const{
	if ( i < dim1() && j < dim2() && k < dim3() ) {
		return get_atAbsoluteIndex( k*dim1()*dim2() + j*dim1() + i );
	} else {
		if (verbose()){
			std::cerr << "Error in array3D! Index (" << i << ", " << j << ", " << k 
				<< ") is larger than dimension (" << dim1() << ", " << dim2() << ", " << dim3() << ")." << endl;
		}
		return 0.;
	}
}

void array3D::set( unsigned int i, unsigned int j, unsigned int k, double value ){
	if ( i < dim1() && j < dim2() && k < dim3() ) {
		set_atAbsoluteIndex( k*dim1()*dim2() + j*dim1() + i, value);
	} else {
		if (verbose()){
			std::cerr << "Error in array3D! Index (" << i << ", " << j  << ", " << k 
				<< ") is larger than dimension (" << dim1() << ", " << dim2() << ", " << dim3() << ")." << endl;
		}
	}
}


//-----------------------------------------------------setters & getters
unsigned int array3D::dim1() const{
	return p_dim1;
}

void array3D::setDim1( unsigned int size_dim1 ){
	p_dim1 = size_dim1;
}

unsigned int array3D::dim2() const{
	return p_dim2;
}

void array3D::setDim2( unsigned int size_dim2 ){
	p_dim2 = size_dim2;
}

unsigned int array3D::dim3() const{
	return p_dim3;
}

void array3D::setDim3( unsigned int size_dim3 ){
	p_dim3 = size_dim3;
}


//------------------------------------------------------------- getASCIIdata
string array3D::getASCIIdata() const{
    ostringstream osst;
    osst << "3D data, dim1=" << dim1() << ", dim2=" << dim2() << ", dim3=" << dim3() << ", size=" << size() << endl;
	for (int k = 0; k<dim3(); k++){
		osst << " [[" << endl;
		for (int j = 0; j<dim2(); j++){
            osst << "  [";
            for (int i = 0; i<dim1(); i++) {
                osst << " " << get(i, j, k);
            }//k
            osst << "]" << endl;
        }//j
		osst << "]]" << endl;
	}//i
    return osst.str();
}



//------------------------------------------------------------- writeToASCII
int array3D::writeToASCII( std::string filename ) const{
	
	ofstream fout( filename.c_str() );
	fout << arraydata::getASCIIdataAsColumn();
	fout.close();
    
    return 0;
}







//*********************************************************************************
//lower level FFT routine performs real-to-complex transform
//input: array1D's data, output: both real and imag parts are returned

FourierTransformer::FourierTransformer(){
    verbose = 0;
}


int FourierTransformer::transform( array1D *real, array1D *imag, int direction ){
    int retval = 0;

    if (verbose) {
        cout << "BEFORE TRANSFORM REAL: " << real->getASCIIdata() << endl;
        cout << "BEFORE TRANSFORM IMAG: " << imag->getASCIIdata() << endl;
        real->writeToASCII("/Users/feldkamp/Desktop/before_trafo_real.txt");
        imag->writeToASCII("/Users/feldkamp/Desktop/before_trafo_imag.txt");
        
    }

    //FFTW comments from official documentation:
    //http://www.fftw.org/fftw3_doc/Complex-One_002dDimensional-DFTs.html
    //
    // input data in[i] are purely real numbers,
    // DFT output satisfies the “Hermitian” redundancy: 
    // out[i] is the conjugate of out[n-i]
    // the input is n real numbers, while the output is n/2+1 complex numbers
    //
    // data type fftw_complex is by default a double[2],
    // composed of the real (in[i][0]) and imaginary (in[i][1]) parts of a complex number
    
    /*   ******from FFTW tutorial*******
         fftw_complex *in, *out;
         fftw_plan p;
         ...
         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
         out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
         p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
         ...
         fftw_execute(p); // repeat as needed
         ...
         fftw_destroy_plan(p);
         fftw_free(in); fftw_free(out);
    */
    
    //set up fftw plan input
    int n = real->size();                   //powers of 2 are especially fast!
    //double *in;
    fftw_complex *in;
    fftw_complex *out;
    int sign = 0;
    unsigned int flags = FFTW_ESTIMATE;             //usually FFTW_MEASURE or FFTW_ESTIMATE
                                                // FFTW_MEASURE will be faster, but plans need to be saved for this
                                                // --> implement for final version
    
    if (direction>=0){
        sign = FFTW_FORWARD;
    }else{
        sign = FFTW_BACKWARD;
    }
    
    //allocate fftw-optimized memory for in and out vectors
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    
    //create plan
    //fftw_plan plan = fftw_plan_dft_r2c_1d(n, in, out, flags); 
    fftw_plan plan = fftw_plan_dft_1d(n, in, out, sign, flags);
    //fftw_plan plan = fftw_plan_dft_r2c(int rank, const int *n, double *in, fftw_complex *out, unsigned int flags);

    //fill data into input vector
    for (int i=0; i<n; i++) {
        in[i][0] = real->get(i);
        in[i][1] = 0;        
    }

    //execute fft
    fftw_execute( plan );
    
    //feed 'out' vector back into the argument arrays that were passed
    //note that r2c_1d output is fast, but has changed output size (int)n/2+1;
    int n_out = n;
    delete real;
    real = new array1D(n_out);
    delete imag;
    imag = new array1D(n_out);
    
    for (int i=0; i<n_out; i++) {               
        real->set(i, out[i][0]);                  //real part
        imag->set(i, out[i][1]);                  //imag part
    }
    
    //clean up
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    
    if (verbose) {
        cout << "AFTER TRANSFORM REAL: " << real->getASCIIdata() << endl;
        cout << "AFTER TRANSFORM IMAG: " << imag->getASCIIdata() << endl;
        real->writeToASCII("/Users/feldkamp/Desktop/after_trafo_real.txt");
        imag->writeToASCII("/Users/feldkamp/Desktop/after_trafo_imag.txt");
    }
    
    return retval;
}





