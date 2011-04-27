/*
 *  ArrayClasses.cpp
 *  xcca
 *
 *  Created by Feldkamp on 2/17/11.
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

#include <fstream>
using std::ofstream;

#include <hdf5.h>
//#include <tiffio.h>

#include <cmath>


//*********************************************************************************
//*********************************************************************************
//
// CLASS IMPLEMENTATION OF arraydata
//
//*********************************************************************************
//*********************************************************************************
arraydata::arraydata(){
	arraydata( 0 );
}

arraydata::arraydata( unsigned int size_val ){
	p_size = size_val;
	p_data = new double[size_val];
    if (p_data == 0)
        cerr << "Error in arraydata constructor: could not allocate memory." << endl;

	zero();					// set all elements to zero initially
	setVerbose( 1 );		// set talkativity
}

arraydata::~arraydata(){
	delete []p_data;	
}


//--------------------------------------------------------------------
double arraydata::get_atAbsoluteIndex( unsigned int i){
	if (i < size_absolute()) {
		return p_data[i];
	} else {
		if (verbose()){
			std::cerr << "Error in arraydata! Index " << i 
			<< " is larger than absolute dimension " << size_absolute() << "." << endl;
		}
		
		//fail 'silently', 
		//i.e., return a valid numeric value 0,
		//instead of exiting or throwing an exception
		return 0.;
	}

}


//--------------------------------------------------------------------
void arraydata::set_atAbsoluteIndex( unsigned int i, double val){
	if (i < size_absolute()) {
		p_data[i] = val;
	} else {									//fail 'silently', as in: don't do anything
		if (verbose()){
			std::cerr << "Error in arraydata! Index " << i 
			<< " is larger than absolute dimension " << size_absolute() << "." << endl;
		}
	}
}


//--------------------------------------------------------------------
unsigned int arraydata::size_absolute(){
	return p_size;
}


//--------------------------------------------------------------------
void arraydata::zero(){				//set all elements to zero
	for (int i = 0; i < size_absolute(); i++) {
		set_atAbsoluteIndex(i, 0);
	}	
}

//--------------------------------------------------------------------
double arraydata::getMin(){
	double tempmin = INFINITY;
	for (int i = 0; i < size_absolute(); i++) {
		if (get_atAbsoluteIndex(i) < tempmin) {
			set_atAbsoluteIndex(i, get_atAbsoluteIndex(i));
		}
	}
	return tempmin;
}

//--------------------------------------------------------------------
double arraydata::getMax(){
	double tempmax = -INFINITY;
	for (int i = 0; i < size_absolute(); i++) {
		if (get_atAbsoluteIndex(i) > tempmax) {
			set_atAbsoluteIndex(i, get_atAbsoluteIndex(i));
		}
	}
	return tempmax;
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
	if (filePointerRead==NULL) 
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
		fputs ("Error in arraydata::readFromRawBinary. Memory error",stderr); 
		exit (2);
	}
	
	// copy the file into the buffer:
	result = fread(buffer,sizeof(uint16_t),lSize,filePointerRead);
	if (result != lSize) 
	{
		fputs ("Error in arraydata::readFromRawBinary. Reading error",stderr); 
		exit (3);
	}
	
	/* the whole file is now loaded in the memory buffer. */
	printf("data contains %ld unsigned 16-bit integers\n",lSize);
	
	// printRawData(buffer,lSize);
	
	//write read data to array for further processing
	for (int i=0; i<size_absolute(); i++) {
		set_atAbsoluteIndex( i, buffer[i] );
	}
	
	
	// terminate reading
	fclose(filePointerRead);
	free(buffer);
}


//--------------------------------------------------------------------
int arraydata::verbose(){
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
array1D::array1D() : arraydata(){
	array1D( 0 );
}

array1D::array1D( unsigned int size_dim1 ) 
		: arraydata(size_dim1){
	p_dim1 = size_dim1;
}

array1D::~array1D(){
}




//-----------------------------------------------------data accessors
double array1D::get( unsigned int i ){
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
unsigned int array1D::dim1(){
	return p_dim1;
}

void array1D::setDim1( unsigned int size_dim1 ){
	p_dim1 = size_dim1;
}











//*********************************************************************************
//*********************************************************************************
//
// CLASS IMPLEMENTATION OF array2D
//
//*********************************************************************************
//*********************************************************************************

//-----------------------------------------------------constructors & destructors
array2D::array2D(){
	array2D( 0, 0 );
}

array2D::array2D( unsigned int size_dim1, unsigned int size_dim2 )
		: arraydata(size_dim1*size_dim2){
	p_dim1 = size_dim1;
	p_dim2 = size_dim2;
}

array2D::~array2D(){
}




//-----------------------------------------------------data accessors
double array2D::get( unsigned int i, unsigned int j ){
	if ( i < dim1() && j < dim2()) {
		return get_atAbsoluteIndex( j*dim1() + i );
	} else {
		if (verbose()){
			std::cerr << "Error in array2D! Index (" << i << ", " << j 
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
			std::cerr << "Error in array2D! Index (" << i << ", " << j 
				<< ") is larger than dimension (" << dim1() << ", " << dim2() << ")." << endl;
		}
	}
}



/*
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
int array2D::writeToTiff( std::string filename, int scaleFlag ){
	if (verbose()) 
		cout << "Writing array2D to TIFF file " << filename << endl;

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
		
		if (verbose())
			cout << "\"" << filename << "\" written to disc!" << endl;
		
		delete tifdata;
		TIFFClose(out);
		return 0;
	}
	else 
		return(1);	
}
*/

//-------------------------------------------------------------- writeToHDF5
int array2D::writeToHDF5( std::string filename ){
		
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
int array2D::writeToASCII( std::string filename, bool printToConsole ){
	
	ofstream fout( filename.c_str() );
	for (int i = 0; i<dim1(); i++) {
		fout << " [";
		for (int j = 0; j<dim2(); j++){
			fout << " " << get(i, j);
		}
		fout << "]" << endl;
	}
	fout.close();
	
	
	//option to print to std out as well
	if ( printToConsole ) {
		for (int i = 0; i<dim1(); i++) {
			printf (" [");
			for (int j = 0; j<dim2(); j++){
				printf (" %3f", get(i, j));
			}
		printf ("]\n");
		}
	}
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
//    hsize_t     dims[1] = { size_absolute() };
	double		*rdata = new double[size_absolute()];
	
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
	if (p_data)
        delete []p_data;
	p_data = rdata;
	
	//set dimension!!!!!
}


//-----------------------------------------------------setters & getters
unsigned int array2D::dim1(){
	return p_dim1;
}

void array2D::setDim1( unsigned int size_dim1 ){
	p_dim1 = size_dim1;
}

unsigned int array2D::dim2(){
	return p_dim2;
}

void array2D::setDim2( unsigned int size_dim2 ){
	p_dim2 = size_dim2;
}




//-----------------------------------------------------test pattern
void array2D::generateTestPattern( int type ){
	switch (type) {
		case 0:											// sinusoidal pattern
			{
				double amplitude = 50000;
				double freqX = 100;
				double freqY = 200;
				for (int i = 0; i < dim1(); i++) {
					for (int j = 0; j < dim2(); j++) {
						set(i, j, amplitude/2*(sin(2*M_PI*i/freqX)*cos(2*M_PI*j/freqY)+1) );
					}
				}
			}
		case 1:											// simply increment by index
			{
				double val = 0;
				for (int i = 0; i < size_absolute(); i++) {
					if (val >= 65535) {
						val = 0; 
					}
					set_atAbsoluteIndex(i, val);
					val += 1;
				}
			}
			break;
		default:
			zero();
			break;
	}
}




//*********************************************************************************
//*********************************************************************************
//
// CLASS IMPLEMENTATION OF array3D
//
//*********************************************************************************
//*********************************************************************************

//-----------------------------------------------------constructors & destructors
array3D::array3D(){
	array3D( 0, 0, 0 );
}

array3D::array3D( unsigned int size_dim1, unsigned int size_dim2, unsigned int size_dim3 )
		: arraydata(size_dim1*size_dim2*size_dim3){
	p_dim1 = size_dim1;
	p_dim2 = size_dim2;
	p_dim3 = size_dim3;
}

array3D::~array3D(){
}


//-----------------------------------------------------data accessors
double array3D::get( unsigned int i, unsigned int j, unsigned int k ){
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
unsigned int array3D::dim1(){
	return p_dim1;
}

void array3D::setDim1( unsigned int size_dim1 ){
	p_dim1 = size_dim1;
}

unsigned int array3D::dim2(){
	return p_dim2;
}

void array3D::setDim2( unsigned int size_dim2 ){
	p_dim2 = size_dim2;
}

unsigned int array3D::dim3(){
	return p_dim3;
}

void array3D::setDim3( unsigned int size_dim3 ){
	p_dim3 = size_dim3;
}




