/*
 *  ArrayClasses.cpp
 *  xcca
 *
 *  Created by Feldkamp on 2/17/11.
 *  Last changed on 07/22/11.
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




//=================================================================================
//
// CLASS IMPLEMENTATION OF arraydata
//
//=================================================================================
//arraydata::arraydata(){
//    init();
//	p_data = new double[0];
//}

arraydata::arraydata( unsigned int size_val ){
    init();
	p_size = size_val;
	p_data = new double[size_val];
    if (p_data == 0)
        cerr << "Error in arraydata constructor: could not allocate memory." << endl;

	zero();					// set all elements to zero initially
}

arraydata::arraydata( const int16_t *CArray, const unsigned int size_val ){
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

arraydata::arraydata( const float *CArray, const unsigned int size_val ) {
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
		cout << "DEBUG: ASSIGNMENT OPERATOR FOR arraydata" << endl;
        this->destroy();
        init();
        copy( src );
    }
    return *this;
}

arraydata::~arraydata(){
	this->destroy();	
}

//--------------------------------------------------------------------helper functions
void arraydata::init(){
    p_size = 0;
    p_data = NULL;
}

void arraydata::destroy(){
    delete []p_data;
	p_data = NULL;
}

//-----------------------------------------------------copy
void arraydata::copy( const double* src_data, const unsigned int arraysize ){		//src type: double c-array
    p_size = arraysize;
    if (p_size > 0){
		this->destroy();
        p_data = new double[ arraysize ];
        for (int i = 0; i < p_size; i++) {
            p_data[i] = src_data[i];
        }
    }else{
        p_data = NULL;
    }
}

void arraydata::copy( const float* src_data, const unsigned int arraysize ){		//src type: double c-array
    p_size = arraysize;
    if (p_size > 0){
		this->destroy();
        p_data = new double[ arraysize ];
        for (int i = 0; i < p_size; i++) {
            p_data[i] = (double)src_data[i];
        }
    }else{
        p_data = NULL;
    }
}

void arraydata::copy( const int* src_data, const unsigned int arraysize ){		//src type: int c-array
    p_size = arraysize;
    if (p_size > 0){
		this->destroy();
        p_data = new double[ arraysize ];
        for (int i = 0; i < p_size; i++) {
            p_data[i] = (double) src_data[i];	//convert to double
        }
    }else{
        p_data = NULL;
    }
}

void arraydata::copy( const arraydata& src ){						//src type: arraydata object
    this->copy( src.data(), src.size());
}


double *arraydata::data() const{
    return (double *)p_data;
}



//--------------------------------------------------------------------
double arraydata::get_atAbsoluteIndex( unsigned int i) const{
	//check for errors, but fail 'silently', 
	//i.e., return a valid numeric value 0,
	//instead of exiting or throwing an exception
		
	if (!p_data){
		std::cerr << "Error in arraydata::get_atAbsoluteIndex(" << i << ")! Internal data not allocated." << endl;
		throw;
	} else if (i >= size()) {
		std::cerr << "Error in arraydata::get_atAbsoluteIndex(" << i 
			<< "). Index is larger than absolute dimension " << size() << "." << endl;
		return 0.;
	} else {
		return p_data[i];
	}
}


//--------------------------------------------------------------------
void arraydata::set_atAbsoluteIndex( const unsigned int i, const double val){
	if ( i < size() ) {
		p_data[i] = val;
	} else {									//fail 'silently', as in: don't do anything
		std::cerr << "Error in arraydata::set_atAbsoluteIndex(" << i 
			<< ", " << val <<  "). Index is larger than absolute dimension " << size() << "." << endl;
	}
}


//--------------------------------------------------------------------
unsigned int arraydata::size() const{
	return p_size;
}


//--------------------------------------------------------------------
void arraydata::zero(){					//set all elements to zero
	for (int i = 0; i < size(); i++) {
		set_atAbsoluteIndex(i, 0);
	}
}

void arraydata::zero( unsigned int start, unsigned int stop ){
	if ( start >= stop || stop > size() ){
		cerr << "Error in arraydata::zero("<< start << ", " << stop << "). Check boundaries." << endl;
		return;
	}
	for (int i = start; i < stop; i++) {
		set_atAbsoluteIndex(i, 0);
	}
}

void arraydata::ones(){					//set all elements to one
	for (int i = 0; i < size(); i++) {
		set_atAbsoluteIndex(i, 1);
	}	
}

void arraydata::ones( unsigned int start, unsigned int stop ){
	if ( start >= stop || stop > size() ){
		cerr << "Error in arraydata::ones("<< start << ", " << stop << "). Check boundaries." << endl;
		return;
	}
	for (int i = start; i < stop; i++) {
		set_atAbsoluteIndex(i, 1);
	}
}


void arraydata::range( double neg, double pos ){	//set elements to a range of values, given by the boundaries
	double delta = (pos-neg)/(size()-1);
	for (int i = 0; i < size(); i++) {
		set_atAbsoluteIndex(i, neg+delta*i);
	}	
}

//--------------------------------------------------------------------
double arraydata::calcMin() const{
	double tempmin = +INFINITY;
	for (int i = 0; i < size(); i++) {
		if (get_atAbsoluteIndex(i) < tempmin) {
			tempmin = get_atAbsoluteIndex(i);
		}
	}
	return tempmin;
}

//--------------------------------------------------------------------
double arraydata::calcMax() const{
	double tempmax = -INFINITY;
	for (int i = 0; i < size(); i++) {
		if (get_atAbsoluteIndex(i) > tempmax) {
			tempmax = get_atAbsoluteIndex(i);
		}
	}
	return tempmax;
}

//------------------------------------------------------------- calcSum
double arraydata::calcSum() const{
	double sum = 0.;
	for (int i = 0; i < size(); i++) {
		sum += get_atAbsoluteIndex(i);
	}	
	return sum;
}

//------------------------------------------------------------- calcAvg
double arraydata::calcAvg() const{
	double avg = this->calcSum() / ((double)size());
	return avg;
}

//------------------------------------------------------------- getASCIIdata
string arraydata::getASCIIdataAsRow() const{
    ostringstream osst;
	if (size() == 0) {
  		osst << "Array data has size zero." << endl;
	}else{
		for (int i = 0; i<size(); i++) {
			osst << get_atAbsoluteIndex(i) << " ";
		}
		osst << endl;
	}
    return osst.str();
}

//------------------------------------------------------------- getASCIIdata
string arraydata::getASCIIdataAsColumn() const{
    ostringstream osst;
	if (size() == 0) {
  		osst << "Array data has size zero." << endl;
	}else{
		for (int i = 0; i<size(); i++) {
			osst << get_atAbsoluteIndex(i) << endl;
		}
	}
    return osst.str();
}

//--------------------------------------------------------------------arraydata::readFromRawBinary
void arraydata::readFromRawBinary( std::string filename ){
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
	
	// the whole file is now loaded in the memory buffer.
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
int arraydata::addValue( double val ){
    for (int i = 0; i<this->size(); i++) {
        this->set_atAbsoluteIndex(i, this->get_atAbsoluteIndex(i) + val);
    }
    return 0;
}

//multiply each element by a numerical factor
int arraydata::subtractValue( double val ){
    for (int i = 0; i<this->size(); i++) {
        this->set_atAbsoluteIndex(i, this->get_atAbsoluteIndex(i) - val);
    }
    return 0;
}

//multiply each element by a numerical value
int arraydata::multiplyByValue( double value ){
    for (int i = 0; i<this->size(); i++) {
        this->set_atAbsoluteIndex(i, this->get_atAbsoluteIndex(i) * value);
    }
    return 0;
}

//divide each element by a numerical value (enforcing float division)
int arraydata::divideByValue( double value ){
    for (int i = 0; i<this->size(); i++) {
        this->set_atAbsoluteIndex(i, ((double) this->get_atAbsoluteIndex(i)) / value);
    }
    return 0;
}


//add each element by an element from a second array
int arraydata::addArrayElementwise( const arraydata *secondArray ){
    if (this->size() != secondArray->size()){
        cerr << "Error in arraydata::addArrayElementwise! Array sizes don't match. ";
        cerr << "(" << this->size() << " != " << secondArray->size() << "). Operation not performed."<< endl;
        return 1;
    }
    
    for (int i = 0; i<this->size(); i++) {
        this->set_atAbsoluteIndex(i, this->get_atAbsoluteIndex(i)+secondArray->get_atAbsoluteIndex(i));
    }
    return 0;
}

//subtract each element by an element from a second array
int arraydata::subtractArrayElementwise( const arraydata *secondArray ){
    if (this->size() != secondArray->size()){
        cerr << "Error in arraydata::subtractArrayElementwise! Array sizes don't match. ";
        cerr << "(" << this->size() << " != " << secondArray->size() << "). Operation not performed."<< endl;
        return 1;
    }
    
    for (int i = 0; i<this->size(); i++) {
        this->set_atAbsoluteIndex(i, this->get_atAbsoluteIndex(i)-secondArray->get_atAbsoluteIndex(i));
    }
    return 0;
}

//multiply each element by an element from a second array
int arraydata::multiplyByArrayElementwise( const arraydata *secondArray ){
    if (this->size() != secondArray->size()){
        cerr << "Error in arraydata::multiplyArrayElementwise! Array sizes don't match. ";
        cerr << "(" << this->size() << " != " << secondArray->size() << "). Operation not performed."<< endl;
        return 1;
    }
    
    for (int i = 0; i<this->size(); i++) {
        this->set_atAbsoluteIndex(i, this->get_atAbsoluteIndex(i) * secondArray->get_atAbsoluteIndex(i) );
    }
    return 0;
}

//divide each element by an element from a second array
int arraydata::divideByArrayElementwise( const arraydata *secondArray ){
    if (this->size() != secondArray->size()){
        cerr << "Error in arraydata::divideArrayElementwise! Array sizes don't match. ";
        cerr << "(" << this->size() << " != " << secondArray->size() << "). Operation not performed."<< endl;
        return 1;
    }
    
    for (int i = 0; i<this->size(); i++) {
        this->set_atAbsoluteIndex(i, this->get_atAbsoluteIndex(i)/secondArray->get_atAbsoluteIndex(i));
    }
    return 0;
}


//--------------------------------------------------------------------
//int arraydata::verbose() const{
//	return p_verbose;
//}
//
//void arraydata::setVerbose( int verbosity ){
//	p_verbose = verbosity;
//}











//=================================================================================
//
// CLASS IMPLEMENTATION OF array1D
//
//=================================================================================


//-----------------------------------------------------constructors & destructors
//array1D::array1D() 
//		: arraydata(){
//    setDim1( 0 );
//}

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


//constructor to generate a 1D array from a 2D array
array1D::array1D( array2D* dataTwoD ) 
        : arraydata( dataTwoD->size() ){
 	setDim1( dataTwoD->size() );

    //copy contents of dataTwoD to this array1D object
    p_size = dataTwoD->size();
    if (p_size > 0){
        for (int i = 0; i < p_size; i++) {
            p_data[i] = dataTwoD->get_atAbsoluteIndex(i);
        }
    }else{
        p_data = NULL;
    }
}


//-----------------------------------------------------destructor
array1D::~array1D(){
}


//-----------------------------------------------------copy
void array1D::copy( const array1D& src ){
    setDim1( src.size() );
    this->arraydata::copy( src.data(), src.size() );
}



//-----------------------------------------------------data accessors
double array1D::get( unsigned int i ) const{
	if ( i < dim1() ) {
		return arraydata::get_atAbsoluteIndex(i);		
	} else {
		std::cerr << "Error in array1D! Index " << i 
				<< " is larger than dimension " << dim1() << "." << endl;
		return 0.;
	}
}

void array1D::set( const unsigned int i, const double value ){
	if ( i < dim1() ) {
		arraydata::set_atAbsoluteIndex(i, value);
	} else {
		std::cerr << "Error in array1D! Index " << i 
			<< " is larger than dimension " << dim1() << "." << endl;
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
	if (size() == 0) {
  		osst << "1D data has size zero." << endl;
	}else{
		osst << "1D data, dim1=" << dim1() << ", size=" << size() << endl;
		osst << " [";
		for (int i = 0; i<dim1(); i++) {
			osst << " " << get(i);
		}
		osst << "]" << endl;
	}
    return osst.str();
}


//------------------------------------------------------------- writeToASCII
int array1D::writeToASCII( std::string filename ) const{
	ofstream fout( filename.c_str() );
	fout << arraydata::getASCIIdataAsColumn();
	fout.close();
	return 0;
}



//=================================================================================
//
// CLASS IMPLEMENTATION OF array2D
//
//=================================================================================

//-----------------------------------------------------constructors & destructors
//array2D::array2D()
//		: arraydata(){
// 	setDim1( 0 );
//	setDim2( 0 );
//}

array2D::array2D( unsigned int size_dim1, unsigned int size_dim2 )
		: arraydata(size_dim1*size_dim2){
 	setDim1( size_dim1 );
	setDim2( size_dim2 );
}


//constructor to generate a 2D array from a 1D array, given the desired dimensions
array2D::array2D( array1D* dataOneD, unsigned int size_dim1, unsigned int size_dim2) 
        : arraydata( size_dim1*size_dim2 ){
    
    if (!dataOneD) {
        cerr << "WARNING in array2D::array2D. Input array1D was not allocated! Nothing copied!" << endl;
        return;
    }
    if (dataOneD->size() != size_dim1*size_dim2) {
        cerr << "WARNING in array2D::array2D. Inconsistent array size. ";
        cerr << "size1D=" << dataOneD->size() << ", size2D=" << size_dim1*size_dim2 
			<< "=" << size_dim1 << "*" << size_dim2 << "" << endl;
    }
    
 	setDim1( size_dim1 );
	setDim2( size_dim2 );

    //copy contents of dataOneD to this array2D
    p_size = dataOneD->size();
    if (p_size > 0){
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
    this->arraydata::copy( src.data(), src.size() );
}

//-----------------------------------------------------data accessors
double array2D::get( unsigned int i, unsigned int j ) const{
	if ( i < dim1() && j < dim2()) {
		return arraydata::get_atAbsoluteIndex( j*dim1() + i );
	} else {
		std::cerr << "Error in array2D::get! Index (" << i << ", " << j 
			<< ") is larger than dimension (" << dim1() << ", " << dim2() << ")." << endl;
		return 0.;
	}
}

void array2D::set( const unsigned int i, const unsigned int j, const double value ){
	if ( i < dim1() && j < dim2()) {
		arraydata::set_atAbsoluteIndex( j*dim1() + i, value);
	} else {
		std::cerr << "Error in array2D::set! Index (" << i << ", " << j 
			<< ") is larger than dimension (" << dim1() << ", " << dim2() << ")." << endl;
	}
}



int array2D::getRow( int rownum, array1D *&row ) const{
	if (rownum >= dim2() || rownum < 0){ 
		cerr << "Error in array2D::getRow. row number " << rownum << " too big or below zero." << endl; 
		return 1;
	}
	if (this->dim1()==0){
		cerr << "Error in array2D::getRow. array2D's dimension1 is zero" << endl; 
		return 2;
	}
	
	// create new array if 'row' doesn't have right size
	if (row->size() != this->dim1()) {
		delete row;
		row = new array1D( this->dim1() );
		if (!row){ 
			cerr << "Error in array2D::getRow. Could not allocate row." << endl;
			return 3;
		}
	}

	//for a fixed row, i goes through columns (x-values)
	for (int i = 0; i < row->size(); i++){
		row->set( i, this->get(i, rownum) );
	}
	return 0;
}


int array2D::getCol( int colnum, array1D *&col ) const{
    if (colnum >= dim1() || colnum < 0){ 
		cerr << "Error in array2D::getCol. column number " << colnum << " too big or below zero." << endl; 
		return 1;
	}
	if (this->dim2()==0){
		cerr << "Error in array2D::getCol. array2D's dimension2 is zero" << endl; 
		return 2;
	}

	// create new array if 'col' doesn't have right size
	if (col->size() != this->dim2()) {		
    	delete col;
    	col = new array1D( this->dim2() );
		if (!col){ 
			cerr << "Error in array2D::getCol. Could not allocate column." << endl; 
			return 3;
		}
 	}
	
	//for a fixed column number, j goes through the rows (y-values)
	for (int j = 0; j < col->size(); j++){
		col->set( j, this->get(colnum, j) );
	}
	return 0;
}


//-----------------------------------------------------setRow/setCol
void array2D::setRow( const int rownum, const array1D *row ){
    if (!row){
        cerr << "Error in array2D::setRow. Row not allocated." << endl;
    }else if (row->size() != dim1()){
        cerr << "Error in array2D::setRow. Dimensions don't match. (row size="
			<< row->size() << ", internal dim1=" << dim1() << ")" << endl;
    }else{
        for (int i = 0; i < row->size(); i++){			//for a fixed row, i goes through columns (x-values)
            this->set( i, rownum, row->get(i) );		// copy values of row to this array2D
        }
    }
}

void array2D::setCol( const int colnum, const array1D *col ){
    if (!col){
        cerr << "Error in array2D::setRow. Row not allocated." << endl;
    }else if (col->size() != dim2()){
        cerr << "Error in array2D::setRow. Dimensions don't match. (col size="
			<< col->size() << ", internal dim1=" << dim2() << ")" << endl;
    }else{
        for (int j = 0; j < col->size(); j++){			//for a fixed column number, j goes through the rows (y-values)
            this->set( colnum, j, col->get(j) );
        }
    }
}

//------------------------------------------------------------- xrange
// example of xrange(-2,2)
//	 -2 -1  0  1  2
//	 -2 -1  0  1  2
//	 -2 -1  0  1  2
//	 -2 -1  0  1  2
//	 -2 -1  0  1  2
//
void array2D::xrange( double xneg, double xpos ){	//set elements to a range of values, given by the boundaries
	double xdelta = (xpos-xneg)/(dim1()-1);
	for (int j = 0; j < dim2(); j++) {
		for (int i = 0; i < dim1(); i++) {
			set(i, j, xneg+xdelta*i );
		}
	}
}

//------------------------------------------------------------- range2D
// example of yrange(-2,2)
//	  2  2  2  2  2
//	  1  1  1  1  1
//	  0  0  0  0  0
//	 -1 -1 -1 -1 -1
//	 -2 -2 -2 -2 -2
//
void array2D::yrange( double yneg, double ypos ){	//set elements to a range of values, given by the boundaries
	double ydelta = (ypos-yneg)/(dim2()-1);
	for (int j = 0; j < dim2(); j++) {
		for (int i = 0; i < dim1(); i++) {
			set(i, j, yneg+ydelta*j);
		}
	}
}

    
//------------------------------------------------------------- getASCIIdata
std::string array2D::getASCIIdata() const{
    ostringstream osst;
	if (size() == 0) {
  		osst << "2D data has size zero." << endl;
	}else{
		osst << "2D data, dim1=" << dim1() << ", dim2=" << dim2() << ", size=" << size() << endl;
		for (int j = 0; j<dim2(); j++){
			osst << " [";
			for (int i = 0; i<dim1(); i++) {
				osst << " " << get(i, j);
			}
			osst << "]" << endl;
		}
	}
    return osst.str();
}



//------------------------------------------------------------- writeToASCII
int array2D::writeToASCII( std::string filename ) const{
	ofstream fout( filename.c_str() );
	fout << arraydata::getASCIIdataAsColumn();
	fout.close();
	return 0;
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
        cout << "WARNING in generateTestPattern. One or more dimensions are zero." << endl;


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
            break;
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
            break;
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
            break;
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




//=================================================================================
//
// CLASS IMPLEMENTATION OF array3D
//
//=================================================================================

//-----------------------------------------------------constructors & destructors
//array3D::array3D()
//		: arraydata(){
//	setDim1( 0 );
//	setDim2( 0 );
//	setDim3( 0 );
//}

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
		return arraydata::get_atAbsoluteIndex( k*dim1()*dim2() + j*dim1() + i );
	} else {
		std::cerr << "Error in array3D! Index (" << i << ", " << j << ", " << k 
			<< ") is larger than dimension (" << dim1() << ", " << dim2() << ", " << dim3() << ")." << endl;
		return 0.;
	}
}

void array3D::set( const unsigned int i, const unsigned int j, const unsigned int k, const double value ){
	if ( i < dim1() && j < dim2() && k < dim3() ) {
		arraydata::set_atAbsoluteIndex( k*dim1()*dim2() + j*dim1() + i, value);
	} else {
		std::cerr << "Error in array3D! Index (" << i << ", " << j  << ", " << k 
			<< ") is larger than dimension (" << dim1() << ", " << dim2() << ", " << dim3() << ")." << endl;
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
	if (size() == 0) {
  		osst << "3D data has size zero." << endl;
	}else{
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
	}
    return osst.str();
}



//------------------------------------------------------------- writeToASCII
int array3D::writeToASCII( std::string filename ) const{
	ofstream fout( filename.c_str() );
	fout << arraydata::getASCIIdataAsColumn();
	fout.close();
    return 0;
}




