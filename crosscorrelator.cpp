/*
 *  CrossCorrelator.cpp
 *  xcca
 *
 *  Created by Feldkamp on 2/17/11.
 *  Last changed on 04/27/11.
 *  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
 *
 */

#include "crosscorrelator.h"



//---------------------------------------------------------------------------
/* XCCA analysis code for LCLS
 * The program reads in a binary with raw data,
 * calculates the average SAXS intensity,
 * calculates the angular cross-correlation, and saves the data to binary.
 * Written 2011-01-19 by Jonas Sellberg (version 1)
 */

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <cmath>

#include <string>
using std::string;





//=================================================================================
//
// constructors & destructor
//
//=================================================================================
CrossCorrelator::CrossCorrelator(int arraylength){
    
    initPrivateVariables();
    
    //set basic properties
    setArraySize( arraylength );
    setQmax(1.);

	// create array for data storage
	// data stored COL by COL from (0,0) in agreement with Sections 0,1 for quad 0 in the first column, 
	// x (2nd column) decreases towards control room, 
	// y (3rd column) decreases towards floor, 
	// detector viewed from front (upstream)
	data = new array1D(arraySize());
	
	// create arrays to keep track of pixel locations in x and y directions
	qx = new array1D(arraySize());
	qy = new array1D(arraySize());
    initDefaultQ();
        
    table = new array2D(50, 50);
	
	// create arrays to store polar coordinates for the angular cross-correlation
	q = new array1D(arraySize());
	phi = new array1D(arraySize());
	
	qave = new array1D(samplingLength());
	iave = new array1D(samplingLength());
	phiave = new array1D(samplingAngle());
	
	crossCorrelation = new array3D( samplingLength(), samplingLength(), samplingLag() );
	
}

CrossCorrelator::CrossCorrelator( float *dataCArray, int arraylength ){
    
    initPrivateVariables();

    //set basic properties, just like the default case
    setArraySize( arraylength );
    setQmax(1.);
	
    //special feature: copy data from array over to internal data structure
    data = new array1D( dataCArray, arraylength );
    
    //allocate all other internal objects
    qx = new array1D(arraySize());
	qy = new array1D(arraySize());
    initDefaultQ();

    table = new array2D(50, 50);
    
	q = new array1D(arraySize());
	phi = new array1D(arraySize());
	
	qave = new array1D(samplingLength());
	iave = new array1D(samplingLength());
	phiave = new array1D(samplingAngle());
	
	crossCorrelation = new array3D( samplingLength(), samplingLength(), samplingLag() );
    
}

CrossCorrelator::CrossCorrelator(float *dataCArray, float *qxCArray, float *qyCArray, int arraylength) {
	
    initPrivateVariables();
    
    //set basic properties, just like the default case
    setArraySize(arraylength);
	//jas: calculate qmax of CArray
    setQmax(qmax2CArray(qxCArray, qyCArray, arraySize()));
	
    //special feature: copy data from array over to internal data structure
    data = new array1D(dataCArray, arraySize());
    
    //allocate all other internal objects
    qx = new array1D(qxCArray, arraySize());
	qy = new array1D(qyCArray, arraySize());
    table = new array2D(50, 50);
	
	q = new array1D(arraySize());
	phi = new array1D(arraySize());
	
	qave = new array1D(samplingLength());
	iave = new array1D(samplingLength());
	phiave = new array1D(samplingAngle());
	
	crossCorrelation = new array3D( samplingLength(), samplingLength(), samplingLag() );

	p_dataCArray = dataCArray;	
}

CrossCorrelator::CrossCorrelator(float *dataCArray, float *qxCArray, float *qyCArray, int arraylength, double qMax, double qMin) {
	
    initPrivateVariables();
    
    //set basic properties, just like the default case
    setArraySize(arraylength);
    setQmaxmin(qMax, qMin);
	
    //special feature: copy data from array over to internal data structure
    data = new array1D(dataCArray, arraySize());
    
    //allocate all other internal objects
    qx = new array1D(qxCArray, arraySize());
	qy = new array1D(qyCArray, arraySize());
    table = new array2D(50, 50);
	
	q = new array1D(arraySize());
	phi = new array1D(arraySize());
	
	qave = new array1D(samplingLength());
	iave = new array1D(samplingLength());
	phiave = new array1D(samplingAngle());
	
	crossCorrelation = new array3D( samplingLength(), samplingLength(), samplingLag() );
	
	p_dataCArray = dataCArray;
	
}

CrossCorrelator::CrossCorrelator(float *dataCArray, float *qxCArray, float *qyCArray, int arraylength, int nq, int nphi) {
	
    initPrivateVariables();
    
    //set basic properties for the size of the arrays
    setArraySize(arraylength);
	p_samplingLength = nq;
	p_samplingAngle = nphi;
	p_samplingLag = (int) ceil(p_samplingAngle/2.0+1);
	
    //special feature: copy data from array over to internal data structure
    data = new array1D(dataCArray, arraySize());
    
    //allocate all other internal objects
    qx = new array1D(qxCArray, arraySize());
	qy = new array1D(qyCArray, arraySize());
    table = new array2D(50, 50);
	
	q = new array1D(arraySize());
	phi = new array1D(arraySize());
	
	qave = new array1D(samplingLength());
	iave = new array1D(samplingLength());
	phiave = new array1D(samplingAngle());
	
	crossCorrelation = new array3D( samplingLength(), samplingLength(), samplingLag() );
	
	p_dataCArray = dataCArray;
	
}

CrossCorrelator::~CrossCorrelator(){
	//free memory for objects
	delete data;
	
	delete qx;
	delete qy;
    delete table;
	
	delete q;
	delete phi;
	
	delete qave;
	delete iave;
	delete phiave;
	
	delete crossCorrelation;
}


//=================================================================================
//
// initialize the class internals and set some defaults
//
//=================================================================================

//----------------------------------------------------------------------------initPrivateVariables
//make sure all private variables are initialized
//so that they don't contain or point to random memory
//----------------------------------------------------------------------------
void CrossCorrelator::initPrivateVariables(){
	p_arraySize = 1;
	p_qmin = 0;
	p_qmax = 0;
	p_deltaq = 0;
	p_phimin = 0;
	p_phimax = 0;
	p_deltaphi = 0;
	p_samplingLength = 0;
	p_samplingAngle = 0;
	p_samplingLag = 0;
    p_outputdir = "";

	data = NULL;
	qx = NULL;
	qy = NULL;	
	q = NULL;
	phi = NULL;	
	qave = NULL;
	iave = NULL;
	phiave = NULL;
	crossCorrelation = NULL;	    
    table = NULL;
	
    p_debug = 0;      
//    check1D = NULL;
	
}


//----------------------------------------------------------------------------initDefaultQ
//if a q-calibration was not given to this class from the outside
//create a default one here
//----------------------------------------------------------------------------
void CrossCorrelator::initDefaultQ(){

    cout << "Initializing qx and qy vectors with default values." << endl;

    //make sure that qx and qy have been allocated previously, otherwise, do so now
    delete qx;
    qx = new array1D(arraySize());
    delete qy;
    qy = new array1D(arraySize());
    
    //set new values for deltaq and qmax
    p_deltaq = 1;                            
    p_qmax = arraySize()/2.*deltaq();

    for (int i=0; i<arraySize(); i++){
        qx->set(i, -qmax()+deltaq()*i );
        qy->set(i, -qmax()+deltaq()*i );
    }
    
    int onedim = (int) floor(sqrt(arraySize()));
    for (int i = 0; i < onedim; i++){
        for (int j = 0; j < onedim; j++){
            qx->set( i*onedim+j,    (i - (onedim-1)/2.)*deltaq() );
            qy->set( i*onedim+j,    (j - (onedim-1)/2.)*deltaq() );
        }
    } 

}



//----------------------------------------------------------------------------initWithTestPattern
void CrossCorrelator::initWithTestPattern( int sizex, int sizey, int type ){

    //create a nice test pattern
    array2D *test = new array2D(sizex, sizey);
    test->generateTestPattern(type);

    //convert it to 1D object and feed it to 'data'
    delete data;
    data = new array1D( test );
    this->setArraySize( sizex*sizey );
    
    //create test q-calibration (should have been done already, but just to make sure)
    initDefaultQ();
    
    //write tiff image to check what the pattern looks like
    test->writeToTiff( outputdir() + "testpattern.tif" );
    
	if (debug()>1) {
		cout << "CrossCorrelator::initWithTestPattern done." << endl;
		cout << "qx = " << qx->getASCIIdata() << endl;
		cout << "qy = " << qy->getASCIIdata() << endl;
		cout << "data = " << data->getASCIIdata() << endl;
    }
    delete test;
}
    

//----------------------------------------------------------------------------initFromFile
// load data from file
//		type 0: HDF5 file (output from Anton's hit finder code
//		type 1: raw 16-bit binary file
void CrossCorrelator::initFromFile( std::string filename, int type ){
	
	//read from file
	switch (type) {
		case 0:
			cout << "FILE READ FOR HDF5 UNDER CONSTRUCTION!!!!!!!" << endl;
			//data->readFromHDF5( filename );
			break;
		case 1:
			data->readFromRawBinary( filename );
			break;
		default:
			std::cerr << "Error in CrossCorrelator::loadFromFile. Type '" << type << "' is no valid file type." << endl;
			break;
	}
		

	printf("maximum Q: %f\n",qmax());
	printf("deltaQ: %f\n",deltaq());
	
	
	//set up q-vectors
	int counter = 0;
	for (int i=0; i<matrixSize(); i++) {
		double qx_val, qy_val;
		qx_val = i*deltaq() - 1;
		for (int j=0; j<matrixSize(); j++) {
			qy_val = j*deltaq() - 1;
			qx->set( counter, qx_val );
			qy->set( counter, qy_val );
			counter++;
		}
	}
}

//=================================================================================
//
// CROSS-CORRELATION ALGORITHM 0
//
//=================================================================================

//----------------------------------------------------------------------------calculatePolarCoordinates
// calculate polar coordinates from cartesian coordinate system
void CrossCorrelator::calculatePolarCoordinates()
{	
	// calculate phi for each pixel and bin angles with correct deltaphi
	for (int i=0; i<arraySize(); i++) {
		
		double phii = phi->get(i);
		double qxi = qx->get(i);
		double qyi = qy->get(i);
		
		// setup UHP
		if (qxi == 0) { // make sure that the column with qx = 0 has angle 0 (q = 0 is assumed to have phi = 0)
			phii = 0;
		} else {
			phii = atan(qxi/qyi); // If qy = 0 and qx != 0, atan gives the correct result, but only for the UHP! Need to add PI for all LHP!
		}

		// correct LHP by adding PI
		if (qyi < 0) {
			phii += M_PI;
		}
		
		if (phii < -deltaphi()/2) { // make sure the binned angle is between 0 and 2PI-deltaphi()
			phii += 2*M_PI;
		}
		
		if (phii < 0) {
			if (p_debug >= 2) cout << "phii: " << phii << ", sampleAngle(phii): " << round(phii/deltaphi()) << endl;
		} else if (phii > (samplingAngle()-1)*deltaphi()) {
			if (p_debug >= 2) cout << "phii: " << phii << ", sampleAngle(phii): " << round(phii/deltaphi()) << endl;
		}
		
		phi->set( i, round(phii/deltaphi()) * deltaphi() );
		
	}
}

void CrossCorrelator::calculatePolarCoordinates(double start_q, double stop_q)
{
	// calculate dependent variables for cheetah.ini control
	p_qmin = start_q;
	p_qmax = stop_q;
	p_deltaq = (qmax()-qmin())/(samplingLength()-1);	// make sure deltaq samples start and stop
	p_deltaphi = (double) 2.0*M_PI/(p_samplingAngle);	// make sure deltaphi samples exactly an interval of 2PI
	
	if (p_debug >= 1) cout << "qmin: " << qmin() << ", qmax: " << qmax() << ", p_deltaq: " << p_deltaq << ", p_samplingLength: " << p_samplingLength << ", p_deltaphi: " << p_deltaphi << ", p_samplingAngle: " << p_samplingAngle << ", p_samplingLag: " << p_samplingLag << endl;
	
	// calculate phi for each pixel and bin angles with correct deltaphi
	for (int i=0; i<arraySize(); i++) {
				
		double phii;
		double qxi = qx->get(i);
		double qyi = qy->get(i);
		
		// setup UHP
		if (qxi == 0) { // make sure that the column with qx = 0 has angle 0 (q = 0 is assumed to have phi = 0)
			phii = 0;
		} else {
			phii = atan(qxi/qyi); // If qy = 0 and qx != 0, atan gives the correct result, but only for the UHP! Need to add PI for all LHP!
		}
		
		// correct LHP by adding PI
		if (qyi < 0) {
			phii += M_PI;
		}
		
		if (phii < -deltaphi()/2) { // make sure the binned angle is between 0 and 2PI-deltaphi()
			phii += 2*M_PI;
		}
		
		if (phii < 0) {
			if (p_debug >= 2) cout << "phii: " << phii << ", sampleAngle(phii): " << round(phii/deltaphi()) << endl;
		} else if (phii > (samplingAngle()-1)*deltaphi()) {
			if (p_debug >= 2) cout << "phii: " << phii << ", sampleAngle(phii): " << round(phii/deltaphi()) << endl;
		}
		
		phi->set( i, round(phii/deltaphi()) * deltaphi() );
		
	}
}	



//----------------------------------------------------------------------------calculateSAXS
void CrossCorrelator::calculateSAXS() {
	
	// using SAXS average for all shots to calculate cross-correlation 
    // or just the SAXS from the specific shots will give different results. 
    // The second choice is probably preferable and is performed here.
	if (p_debug >= 1) printf("calculating average SAXS intensity...\n");
	
	// calculate |q| for each pixel and bin lengths with correct resolution
	for (int i=0; i<arraySize(); i++) {
		q->set(i, round(sqrt( (qx->get(i)*qx->get(i))+(qy->get(i)*qy->get(i)) ) / deltaq()) * deltaq() );
		if (p_debug >= 3) printf("q[%d]: %f\n",i,q->get(i));
	}
	
	// angular average for each |q|
	if (p_debug >= 1) printf("# of steps: %d\n",samplingLength());
	if (p_debug >= 2) printf("average SAXS intensity:\n");
	
	for (int i=0; i<samplingLength(); i++) {
		qave->set( i, qmin()+i*deltaq() );
		double itot = 0; // reset summed intensity
		int counter = 0; // reset counter
		for (int j=0; j<arraySize(); j++) {
			if ( q->get(j) == qave->get(i) ) {
				itot += data->get(j);
				counter++;
			}
		}
		iave->set( i, itot/counter );
		if (p_debug >= 2) cout << "Q: " << qave->get(i) << ",   \t# pixels: " << counter << ",\tI: " << iave->get(i) << endl;
	}
}



//----------------------------------------------------------------------------calculateXCCA
void CrossCorrelator::calculateXCCA(){

	if (p_debug >= 1) cout << "deltaPhi: " << deltaphi() << endl;
	if (p_debug >= 1) cout << "# of angles: " << samplingAngle() << endl;
	for (int i=0; i<samplingAngle(); i++) {
		phiave->set( i, phimin()+i*deltaphi() );
	}
	
	// create array of the speckle pattern with the correct binning
	array2D *speckle = new array2D( samplingLength(), samplingAngle() );
	
	// create array over pixel counts for each sampled q and phi
	array2D *pixelCount = new array2D( samplingLength(), samplingAngle() );
	array2D *pixelBool = new array2D( samplingLength(), samplingAngle() );

	if (p_debug >= 1) printf("calculating speckle arrays...\n");

	for (int i=0; i<arraySize(); i++) {
		int qIndex = (int) round((q->get(i)-qmin())/deltaq()); // the index in qave[] that corresponds to q[i]
		int phiIndex = (int) round((phi->get(i)-phimin())/deltaphi()); // the index in phiave[] that corresponds to phi[i]
		if (p_debug >= 3) printf("qIndex: %d, phiIndex: %d\n", qIndex, phiIndex);
		if (qIndex >= 0 && qIndex < samplingLength() && phiIndex >= 0 && phiIndex < samplingAngle()) { // make sure qIndex and phiIndex is not out of array bounds
			speckle->set(qIndex, phiIndex, speckle->get(qIndex, phiIndex) + data->get(i) );
			pixelCount->set(qIndex, phiIndex, pixelCount->get(qIndex,phiIndex)+1);
			if (pixelBool->get(qIndex, phiIndex) != 1) {
				pixelBool->set(qIndex, phiIndex, 1);
			}
		} else if (p_debug >= 2) printf("POINT EXCLUDED! qIndex: %d, phiIndex: %d\n", qIndex, phiIndex);
	}
		
	// subtract the average SAXS intensity from the speckle array
	array2D *speckleNorm = new array2D( samplingLength(), samplingAngle() );
	
	for (int i=0; i<samplingLength(); i++) { // make sure each element is initially zero
		for (int j=0; j<samplingAngle(); j++) {
			if (pixelBool->get(i,j) != 0) {
				speckle->set(i, j, speckle->get(i,j) / pixelCount->get(i,j) );
				speckleNorm->set(i, j, speckle->get(i,j) - iave->get(i) );
			}
			if (p_debug >= 2) printf("q: %f, phi: %f --> bool: %f, count: %f\n", qave->get(i), phiave->get(j), pixelBool->get(i, j), pixelCount->get(i, j));
		}
	}
	
	// calculate cross-correlation array and normalization array for cross-correlation
	if (p_debug >= 1) cout << "# of angular lags: " << samplingLag() << endl;

	array3D *normalization = new array3D( samplingLength(), samplingLength(), samplingLag() );
	
	if (p_debug >= 1) printf("starting main loop to calculate cross-correlation...\n");
	
	for (int i=0; i<samplingLength(); i++) { // q1 index
		for (int j=0; j<samplingLength(); j++) { // q2 index 
			for (int k=0; k<samplingLag(); k++) { // phi lag => phi2 index = (l+k)%samplingAngle
				for (int l=0; l<samplingAngle(); l++) { // phi1 index
					crossCorrelation->set(i,j,k, crossCorrelation->get(i,j,k) + speckleNorm->get(i,l)*speckleNorm->get(j, (l+k)%samplingAngle()) );
					normalization->set(i, j, k, normalization->get(i, j, k) + pixelBool->get(i,l)*pixelBool->get(j, (l+k)%samplingAngle()) );
					if (p_debug >= 3) printf("phi2: %d\n",(l+k)%samplingAngle());
				}
			}
		}
	}
	
	// normalize the cross-correlation array with the average SAXS intensity and the calculated normalization constant
	for (int i=0; i<samplingLength(); i++) { // make sure each element is initially zero
		for (int j=0; j<samplingLength(); j++) {
			for (int k=0; k<samplingLag(); k++) {
				if (normalization->get(i,j,k) != 0) {
					crossCorrelation->set(i, j, k, crossCorrelation->get(i,j,k) / (normalization->get(i,j,k)*iave->get(i)*iave->get(j)) );
				}
			}
		}
	}
	
	delete pixelBool;
	delete pixelCount;
	
	delete normalization;
	
	delete speckle;
	delete speckleNorm;
	
	if (p_debug >= 1) printf("done calculating cross-correlation...\n");
}


///////////////OLD STUFF////////////////
/*
//----------------------------------------------------------------------------writeSAXS
void CrossCorrelator::writeSAXS()
{	
	// jas: writeSAXS is currently just used to check that the algorithm works
	// this was done for r0003 and it has NOT been generalized yet!
	
	// write cross-correlation and average SAXS intensity to binary
	printf("writing data to file...\n");
	FILE *filePointerWrite;
	double samplingLengthD = (double) samplingLength(); // save everything as doubles
	double samplingLagD = (double) samplingLag();
	double samplingAngleD = (double) samplingAngle();
	double *buffer;
	buffer = (double*) calloc(samplingLength()*samplingLength()*samplingLag(), sizeof(double));
	
    string filename = outputdir()+"r0003-xcca.bin";
	filePointerWrite = fopen(filename.c_str(),"w+"); // jas: TEST FILE, need to change this to a general string name later
	
	// angular averages
	for (int i=0; i<samplingLength(); i++) {
		buffer[i] = iave->get(i);
	}
	fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite); // saving dimensions of array before the actual data
	fwrite(&buffer[0],sizeof(double),samplingLength(),filePointerWrite);
	
	// q binning
	for (int i=0; i<samplingLength(); i++) {
		buffer[i] = qave->get(i);
	}
	fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
	fwrite(&buffer[0],sizeof(double),samplingLength(),filePointerWrite);
	
	// angle binning
	for (int i=0; i<samplingAngle(); i++) {
		buffer[i] = phiave->get(i);
	}
	fwrite(&samplingAngleD,sizeof(double),1,filePointerWrite);
	fwrite(&buffer[0],sizeof(double),samplingAngle(),filePointerWrite);
	
	// cross-correlation - full version
	for (int i=0; i<samplingLength(); i++) {
		for (int j=0; j<samplingLength(); j++) {
			for (int k=0; k<samplingLag(); k++) {
				buffer[i*samplingLength()*samplingLag()+j*samplingLag()+k] = crossCorrelation->get(i,j,k);
			}
		}
	}
	fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
	fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
	fwrite(&samplingLagD,sizeof(double),1,filePointerWrite);
	fwrite(&buffer[0],sizeof(double),samplingLength()*samplingLength()*samplingLag(),filePointerWrite);
	
	// cross-correlation - autocorrelation only (q1=q2)
//	for (int i=0; i<samplingLength(); i++) {
//		for (int k=0; k<samplingLag(); k++) {
//			buffer[i*samplingLag()+k] = crossCorrelation->get(i,i,k);
//		}
//	}
//	fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
//	fwrite(&samplingLagD,sizeof(double),1,filePointerWrite);
//	fwrite(&buffer[0],sizeof(double),samplingLength()*samplingLag(),filePointerWrite);
	
	fclose(filePointerWrite);
	free(buffer);
	
	cout << "writeSAXS done" << endl;
}



//----------------------------------------------------------------------------writeXCCA
void CrossCorrelator::writeXCCA(){
	printf("writing data to file...\n");
	
	// All saving is currently handled by writeSAXS()
	
	//jas: saving &array1D->get(0) does NOT work with fwrite, need to loop through array1D and save into Carray before saving to file...
	//jas: have to call the pointer to the data array in array1D through &array1D->data() instead.
	
	FILE *filePointerWrite;
	
	filePointerWrite = fopen("f909-q0-xcca.bin","w+");
	
	double samplingLengthD = (double) samplingLength(); // save everything as doubles
	double samplingLagD = (double) samplingLag();
	double samplingAngleD = (double) samplingAngle();
	
	fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite); // saving dimensions of array before the actual data
	fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
	fwrite(&samplingLagD,sizeof(double),1,filePointerWrite);
	fwrite(&crossCorrelation->get(0,0,0),sizeof(double),samplingLength()*samplingLength()*samplingLag(),filePointerWrite); // saving data as arrays of LAG in the following order [0][0][LAG], [0][1][LAG], ... , [0][LENGTH][LAG], [1][0][LAG], [1][1][LAG], and so on
	fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
	fwrite(&iave->get(0),sizeof(double),samplingLength(),filePointerWrite);
	fwrite(&samplingLengthD,sizeof(double),1,filePointerWrite);
	fwrite(&qave->get(0),sizeof(double),samplingLength(),filePointerWrite);
	fwrite(&samplingAngleD,sizeof(double),1,filePointerWrite);
	fwrite(&phiave->get(0),sizeof(double),samplingAngle,filePointerWrite);
	
	fclose(filePointerWrite);
	
	cout << "writeXCCA done" << endl;
}
*/


//----------------------------------------------------------------------------dumpResults
void CrossCorrelator::dumpResults( std::string filename ){
	cout << "Writing results to TIFF file '" << filename << "'." << endl;
    
    //convert data to 2D array class (using known dimensions)
    int dim1 = matrixSize();                // not exact.... will have to change this soon
    int dim2 = matrixSize();
    array2D *dataTwoD = new array2D( data, dim1, dim2 );
    
    //write to tiff image
    dataTwoD->writeToTiff( outputdir()+filename );
    delete dataTwoD;
}



//----------------------------------------------------------------------------printRawData
// print raw data after having read from file
void CrossCorrelator::printRawData(uint16_t *buffer,long lSize) {
	for (int i=0; i<lSize; i++) {
		printf("%u ",buffer[i]);
	}
	printf("\n");
}


//=================================================================================
//
// setters and getters for private variables
//
//=================================================================================
int CrossCorrelator::arraySize() const{
	return p_arraySize;
}

void CrossCorrelator::setArraySize( int arraySize_val ){
	p_arraySize = arraySize_val;
}

int CrossCorrelator::matrixSize() const{
	return (int)floor(sqrt( (double)arraySize() ));            //assuming a square image
}

void CrossCorrelator::setMatrixSize( int matrixSize_val ){
	setArraySize( matrixSize_val*matrixSize_val );
	updateDependentVariables();
}

void CrossCorrelator::setQmaxmin( double qmax_val, double qmin_val ){
	p_qmax = qmax_val;
	p_qmin = qmin_val;
	updateDependentVariables();
}

double CrossCorrelator::qmin() const{
	return p_qmin;
}

void CrossCorrelator::setQmin( double qmin_val ){
	p_qmin = qmin_val;
}

double CrossCorrelator::qmax() const{
	return p_qmax;
}

void CrossCorrelator::setQmax( double qmax_val ){
	p_qmax = qmax_val;
	updateDependentVariables();
}

double CrossCorrelator::phimin() const{
	return p_phimin;
}

void CrossCorrelator::setPhimin( double phimin_val ){
	p_phimin = phimin_val;
}

double CrossCorrelator::phimax() const{
	return p_phimax;
}

void CrossCorrelator::setPhimax( double phimax_val ){
	p_phimax = phimax_val;
}

double CrossCorrelator::qmax2CArray( float *qxCArray, float *qyCArray, int arraylength ) {
	double qmax = 0;
	for (int i=0; i<arraylength; i++) {
		double qtemp = (double) sqrt(qxCArray[i]*qxCArray[i] + qyCArray[i]*qyCArray[i]);
		if (qtemp > qmax) qmax = qtemp;
	}
	return qmax;
}

double CrossCorrelator::qmax1CArray( float *qCArray, int arraylength ) {
	double qmax = 0;
	for (int i=0; i<arraylength; i++) {
		double qtemp; 
		if (qCArray[i] > 0) qtemp = (double) qCArray[i];
		else qtemp = (double) -qCArray[i];
		if (qtemp > qmax) qmax = qtemp;
	}
	return qmax;
}

void CrossCorrelator::setOutputdir( std::string dir ){
    p_outputdir = dir;
}

string CrossCorrelator::outputdir(){
    return p_outputdir;
}

int CrossCorrelator::debug(){
    return p_debug;
}

void CrossCorrelator::setDebug( int debuglevel ){
    p_debug = debuglevel;
}

    
double CrossCorrelator::deltaq() const{						//getter only, dependent variable
	return p_deltaq;
}

double CrossCorrelator::deltaphi() const{						//getter only, dependent variable
	return p_deltaphi;
}

int CrossCorrelator::samplingLength() const{						//getter only, dependent variable
	return p_samplingLength;
}

int CrossCorrelator::samplingAngle() const{						//getter only, dependent variable
	return p_samplingAngle;
}

int CrossCorrelator::samplingLag() const{					//getter only, dependent variable
	return p_samplingLag;
}

void CrossCorrelator::updateDependentVariables(){		//update the values that depend on qmax and matrixSize
	// FINE BINNING
//	p_deltaq = 2*qmax()/(matrixSize()-1);
//	p_samplingLength = int(qmax()/p_deltaq+1+0.001);
//	p_deltaphi = 2*atan(1/(2*(p_samplingLength-1.0)));
//	p_samplingAngle = int(2*round(M_PI/p_deltaphi)); // make sure p_samplingAngle is even (exclude 2PI)
//	p_deltaphi = (double) 2.0*M_PI/(p_samplingAngle); // make sure deltaphi samples exactly an interval of 2PI
//	p_samplingLag = (int) round(p_samplingAngle/2.0+1);
	
	// COARSE BINNING
	p_deltaq = 20*qmax()/(matrixSize()-1);
	p_samplingLength = int(qmax()/p_deltaq+1+0.001);
	p_deltaphi = 2*atan(1/(2*(p_samplingLength-1.0)));
	p_samplingAngle = int(2*round(M_PI/p_deltaphi)); // make sure p_samplingAngle is even (exclude 2PI)
	p_deltaphi = (double) 2.0*M_PI/(p_samplingAngle); // make sure deltaphi samples exactly an interval of 2PI
	p_samplingLag = (int) round(p_samplingAngle/2.0+1);
	
	if (p_debug >= 1) {
		cout << "CrossCorrelator::updateDependentVariables done. qmax: " << qmax() << ", p_deltaq: " << p_deltaq << ", p_samplingLength: " << p_samplingLength << ", p_deltaphi: " << p_deltaphi << ", p_samplingAngle: " << p_samplingAngle << ", p_samplingLag: " << p_samplingLag << endl;
	}
}

double CrossCorrelator::getQave(unsigned index) const {
	return qave->get(index);
}

double CrossCorrelator::getPhiave(unsigned index) const {
	return phiave->get(index);
}

double CrossCorrelator::getIave(unsigned index) const {
	return iave->get(index);
}

double CrossCorrelator::getCrossCorrelation(unsigned index1, unsigned index2, unsigned index3) const {
	return crossCorrelation->get(index1,index2,index3);
}



//=================================================================================
//
// CROSS-CORRELATION ALGORITHM 1 (_FAST)
//
//=================================================================================
int CrossCorrelator::calculatePolarCoordinates_FAST(array2D *&polar, int number_q, int number_phi){

    //call the more specific function with some reasonable default values
    double start_q = 3*deltaq();    //default: start at a close-to but non-zero value
    double stop_q = qmax();         //default: go out to the maximum q
    double start_phi = 0;           //default: full circle
    double stop_phi = 360;    
    return calculatePolarCoordinates_FAST(polar, number_q, start_q, stop_q, number_phi, start_phi, stop_phi);
}


//----------------------------------------------------------------------------transform to polar coordinates
int CrossCorrelator::calculatePolarCoordinates_FAST(array2D *&polar, 
                                                        int number_q, double start_q, double stop_q,
                                                        int number_phi, double start_phi, double stop_phi ){
	if( debug()>1 ){ 
		cout << "CrossCorrelator::calculatePolarCoordinates_FAST" << endl; 
		cout << "varying scattering vector q from " << start_q << " to " <<  stop_q << " in " << number_q << " steps, "
			<< "and angle phi from " << start_phi << " to " <<  stop_phi << " in " << number_phi << " steps." << endl;
	}
    int retval = 0;

    double step_q = (stop_q - start_q)/number_q;
    double step_phi = (stop_phi - start_phi)/number_phi;
    
    if (step_q < 0)
        cerr << "Error in CrossCorrelator::calculatePolarCoordinates_FAST -- step_r value " 
            << step_q << " is smaller than zero." << endl;
    if (step_phi < 0)
        cerr << "Error in CrossCorrelator::calculatePolarCoordinates_FAST -- step_phi value " 
            << step_phi << " is smaller than zero." << endl;
    
	delete polar;
    polar = new array2D( (int) floor(number_phi+0.5), (int) floor(number_q+0.5) );
	if (!polar){
		cerr << "Error in CrossCorrelator::calculatePolarCoordinates_FAST. polar couldn't be allocated." << endl;
		return 1;
	}

    double xcoord = 0.;
    double ycoord = 0.;
    double value = 0.;
    double q = 0.;
    double p = 0.;
    int qcounter = 0;
    int pcounter = 0;
    
	for(q = start_q, qcounter=0; qcounter < number_q; q+=step_q, qcounter++){                        // q: for all rings/q-values
        if(debug()){cout << "#" << qcounter << ",q=" << q << "  " << std::flush;}
		for(p = start_phi, pcounter=0; pcounter < number_phi; p+=step_phi, pcounter++){				// phi: go through all angles

            //find lookup coordinates
			xcoord = q * cos(p*M_PI/180);
			ycoord = q * sin(p*M_PI/180);

            //lookup that value in original scattering data
            value = lookup( xcoord, ycoord );
            
			//assign the new values (note the functional determinant q to account for bins growing as q grows)
			polar->set(pcounter, qcounter, value * q);
		}
	}
    
    if (debug()) {
		cout << "polarCoordinates done. dimensions=(" << polar->dim1() << " x " << polar->dim2() << ")" << endl;
	    //write output of the intermediate files? (all zero by default, turn on for debugging or whatever)
		int output_data_ASCII = 0;
		int output_polar_ASCII = 0;
		int output_polar_TIF = 0;
		if (output_data_ASCII){ cout << "data: " << data->getASCIIdata() << endl; }
		if (output_polar_ASCII){ cout << "polar: " << polar->getASCIIdata() << endl; }
		if (output_polar_TIF){ polar->writeToTiff( outputdir()+"polar.tif" ); }            
	}
    return retval;
}



//----------------------------------------------------------------------------calculate XCCA
int CrossCorrelator::calculateXCCA_FAST( array2D *&polar, array2D *&corr, int writeToInternalDataStructure ){
	if(debug()>1){ 
		cout << "CrossCorrelator::calculateXCCA_FAST ("<< polar->dim1() << ", " << polar->dim2() << ")" << endl; 
	}
	
	//do some sanity checks first
	if (!polar) {
  		cerr << "Error in CrossCorrelator::calculateXCCA_FAST. Polar coordinate matrix not allocated. Aborting." << endl;
		return 1;
	}
	if (polar->dim1()==0 || polar->dim2()==0) {
		cerr << "Error in CrossCorrelator::calculateXCCA_FAST. Polar coordinate matrix dimension is zero. Aborting." << endl;
		return 2;
	}
	
	//create a new output array 'corr'
	delete corr;
	corr = new array2D( polar->dim1(), polar->dim2() );
	if (!corr) {
		cerr << "Error in CrossCorrelator::calculateXCCA_FAST. Output 'corr' could not be allocated. Aborting." << endl;
		return 3;
	}
	
	//calculate the auto-correlation for all rings
	for(int q_ct=0; q_ct < polar->dim2(); q_ct++){

		array1D *f = new array1D( polar->dim1() );
		if (!f){ 
			cerr << "Error in CrossCorrelator::calculateXCCA_FAST. Single row 'f' could not be allocated. Aborting." << endl;
			return 4;
		}

		//get one row out of the polar coordinates (fixed q)
		polar->getRow( q_ct, f );		
		if (f->size() == 0) {
			cerr << "Error in CrossCorrelator::calculateXCCA_FAST. Could not get single row. Aborting." << endl;
			return 5;
		}
		double avg = f->calcAvg();
		
		
		if (debug()>1){ cout << "   #" << q_ct << ", f before FFT: " << f->getASCIIdata() << endl; }
		
		//perform autocorrelation --> compute via FFT
		if (debug()) { cout << q_ct << " " << std::flush; }
		
		int fail = autocorrelateFFT( f );          // should yield the same result as correlateFFT( f, f );
		if (fail){
			cerr << "Error in CrossCorrelator::calculateXCCA_FAST. Could calculate correlation. Aborting." << endl;
			return 6;
		}
		
		//normalize
		f->multiplyByFactor(1/avg);
		
		//feed result into corr
		corr->setRow( q_ct, f );
		
		if (debug()>1){ cout << "   #" << q_ct << ", f after FFT: " << f->getASCIIdata() << endl; }
		delete f;
		f = NULL;
	}//for

	if(debug()){
		cout << "Translating to internal data structure";
		cout << " (" << samplingLength() << ", " << samplingLength() << ", " << samplingLag() << ")" << endl;
	}
	if (writeToInternalDataStructure){
		for (int i=0; i<samplingLength(); i++) { // q1=q2 index (only doing autocorrelation at this point)
			for (int k=0; k<samplingLag(); k++) { // phi lag => phi2 index = (l+k)%samplingAngle
				crossCorrelation->set(i,i,k, corr->get(k, i) );
			}
		}
		cout << endl;
	}
	
	if(debug()){ cout << "CrossCorrelator::calculateXCCA_FAST done." << endl; }			
    return 0;
}




//----------------------------------------------------------------------------correlate
// compute 1D correlation corr(f,g) using FFT, result is written to f
int CrossCorrelator::correlateFFT( array1D *f, array1D *g ){
    int retval = 0;

	if (!f || !g) {
  		cerr << "CrossCorrelator::correlateFFT. Input not allocated." << endl;
		return 1;
	}
	if (f->size()==0 || g->size()==0) {
  		cerr << "CrossCorrelator::correlateFFT. Input has size zero." << endl;
		return 2;
	}
	    
    //-------------------------------------------------------------------------
    //   Correlation Theorem:
    //   multiplying the FT of one function by the complex conjugate 
    //   of the FT of the other gives the FT of their correlation
    //
    //   http://mathworld.wolfram.com/Cross-CorrelationTheorem.html
    //-------------------------------------------------------------------------
    
    FourierTransformer *ft = NULL;
        
    // transform f -> F
    array1D *f_real = new array1D( *f );
    array1D *f_imag = new array1D(f->size());				// initialized to zeros --> f is real
    ft = new FourierTransformer( f_real, f_imag );
    int f_fail = ft->transform( 1 );
	ft->getData( f_real, f_imag );
	delete ft;
    if (f_fail){
        cerr << "Error in CrossCorrelator::correlateFFT. Transform (f->F) failed." << endl;
        retval++;
    }
        
    // transform g -> G
    array1D *g_real = new array1D( *g );
    array1D *g_imag = new array1D(g->size());				// initialized to zeros --> g is real
	ft = new FourierTransformer( g_real, g_imag );
    int g_fail = ft->transform( 1 );
	ft->getData( g_real, g_imag );
	delete ft;
    if (g_fail){
        cerr << "Error in CrossCorrelator::correlateFFT. Transform (g->G) failed." << endl;
        retval++;
    }

    // compute F * G_cc (complex conjugate)
    // if F = a+ib, G = c+id, then FG_cc = ac + bd + ibc - iad
    array1D *FG_real = new array1D( g_real->size() );
    array1D *FG_imag = new array1D( g_real->size() );
    for (int i=0; i<f_real->size(); i++) {
        FG_real->set( i,   f_real->get(i)*g_real->get(i) + f_imag->get(i)*g_imag->get(i)   );   // ac + bd
        FG_imag->set( i,   f_imag->get(i)*g_real->get(i) - f_real->get(i)*g_imag->get(i)   );   // i(bc - ad)
    }
    
    // transform the result back to find the correlation
    // transform FG -> corr(f,g)
	ft = new FourierTransformer( FG_real, FG_imag );
    int FG_fail = ft->transform( -1 );
	ft->getData( FG_real, FG_imag );
	delete ft;
    if (FG_fail){
        cerr << "Error in CrossCorrelator::correlateFFT. Transform (FG->corr) failed." << endl;
        retval++;
    }

    
    // return result in original argument arrays
    f->copy( *FG_real );
    g->copy( *FG_imag );
    
    //normalize
    f->multiplyByFactor( 1/((double)f->size()) );
    g->multiplyByFactor( 1/((double)f->size()) );

    
    delete f_real;
    delete f_imag;    
    delete g_real;
    delete g_imag;
    delete FG_real;
    delete FG_imag;
    
    
    return retval;
}


//----------------------------------------------------------------------------autocorrelate
// compute 1D autocorrelation corr(f,f) using FFT, result is written to f
int CrossCorrelator::autocorrelateFFT( array1D *f ){         
    int retval = 0;
    
	if (!f || !f->data()) {
  		cerr << "CrossCorrelator::autocorrelateFFT. Input not allocated." << endl;
		return 1;
	}
	if (f->size()==0) {
  		cerr << "CrossCorrelator::autocorrelateFFT. Input has size zero." << endl;
		return 2;
	}
	
    //-------------------------------------------------------------------------
    //   Wiener-Khinchin Theorem:
    //   the autocorrelation of f is simply given by the Fourier transform 
    //   of the absolute square of F
    //   http://mathworld.wolfram.com/Wiener-KhinchinTheorem.html
    //-------------------------------------------------------------------------
    
    FourierTransformer *ft = NULL;
	
    //transform forward
    array1D *f_imag = new array1D(f->size());							// initialized to zeros --> f is real
	
	ft = new FourierTransformer(f, f_imag);
    int fail = ft->transform( 1 );
	ft->getData( f, f_imag );
	delete ft;
    if (fail){
        cerr << "Error in CrossCorrelator::autocorrelateFFT. Transform (forward) failed." << endl;
        retval++;
    }
    
    //calculate the magnitude squared
    // if F = a+ib, then |F|^2 = a^2 + b^2
    for (int i=0; i < f->size(); i++) {
        f->set( i,   f->get(i)*f->get(i) + f_imag->get(i)*f_imag->get(i)   );
    }    
    f_imag->zero();                 //set to zero for back transform

    //transform back
    // after inverse transform, result is stored in original argument array f
    ft = new FourierTransformer( f, f_imag );
	int fail_inv = ft->transform( -1 );    
	ft->getData( f, f_imag );
	delete ft;
    if (fail_inv){
        cerr << "Error in CrossCorrelator::autocorrelateFFT. Transform (backward) failed." << endl;
        retval++;
    }
    
	//normalize
    f->multiplyByFactor( 1/((double)f->size()) ); 
    
    delete f_imag;
    return retval;
}






//----------------------------------------------------------------------------lookup
double CrossCorrelator::lookup( double xcoord, double ycoord ) {

    double val = 0.;    //return data value at the given coordinates
    int index = 0;      //to do that, the index in the data is determined first
    
    //create lookup index from original coordinates
    //we assume the data to be centered!
    //(the add-one-half->floor trick is to achieve reasonable rounded integers)
	double deltaLUT = deltaq();
	//double deltaLUT2 = 1;
	double qmin = -qmax();
    double xc = (xcoord-qmin) / deltaLUT;
    double yc = (ycoord-qmin) / deltaLUT;
    int ix = (int) floor( xc + 0.5 );		// round to nearest integer
    int iy = (int) floor( yc + 0.5 );
	    
    if ( !table ){
        cerr << "Error in lookup! No lookup table was allocated." << endl;
    } else if ( (ix < 0) || (ix >= table->dim1()) ){
        cerr << "Error in lookup! xcoord=" << xcoord << " is too large or too small.";
        cerr << " (ix=" << ix << ", table dimx=" << table->dim1() << ")" << endl;
    } else if ( (iy < 0) || (iy >= table->dim2()) ){
        cerr << "Error in lookup! ycoord=" << ycoord << " is too large or too small.";
        cerr << " (iy=" << iy << ", table dimy=" << table->dim2() << ")" << endl;
    } else {
        index = (int) floor( table->get(ix, iy) + 0.5 );
        val = data->get( index );


		//-----------------------!!!DEBUG!!!
		//make a hot pixel out of the one that was just looked up
		//careful! this actually changes the data in the input array!!
		//---> therefore, uncomment the following line only for debugging
		//p_dataCArray[index] = 100000.;
		//-----------------------!!!DEBUG!!!
    }

    //keep track of where the value was read in a separate 'check1D' array
    if( debug()>2 ){
        //check1D->set(index, 65535);         
        cout << "lookup (" << xcoord << ", " << ycoord 
            << ") --> LUT: (xc,yc)=(" << xc << ", " << yc 
            << ") ==> (" << ix << ", " << iy << ") "
            << "--> index=" << index << ", --> val=" << val << endl;
    }

    return val;
}



//----------------------------------------------------------------------------setlookupTable
void CrossCorrelator::setLookupTable( array2D *LUT ){
	double qx_min, qx_stepsize, qy_min, qy_stepsize = 0.;
	calcLUTvariables( LUT->dim1(), LUT->dim2(), qx_min, qx_stepsize, qy_min, qy_stepsize );

	table->copy(*LUT);										//store a copy locally
}

void CrossCorrelator::setLookupTable( const int *cLUT, unsigned int LUT_dim1, unsigned int LUT_dim2 ){
	double qx_min, qx_stepsize, qy_min, qy_stepsize = 0.;
	calcLUTvariables( LUT_dim1, LUT_dim2, qx_min, qx_stepsize, qy_min, qy_stepsize );
	
	unsigned int tablelength = LUT_dim1*LUT_dim2;
	table->setDim1(LUT_dim1);						//set dimensions
	table->setDim2(LUT_dim2);
	table->arraydata::copy(cLUT, tablelength);		//store a copy in 'table' locally
}

//----------------------------------------------------------------------------calcLUTvariables
// calculate some variables needed to fill the lookup table
// and update the private variables p_deltaq and p_qmax
// which are important parameters for the lookup() function
//----------------------------------------------------------------------------
void CrossCorrelator::calcLUTvariables( int lutNx, int lutNy, double &qx_min, double &qx_stepsize, double &qy_min, double &qy_stepsize ){
    double qx_max, qx_range, qy_max, qy_range = 0.;
	qx_min = qx->calcMin();
    qx_max = qx->calcMax();
    qx_range = fabs(qx_max - qx_min);
    qx_stepsize = qx_range/(double)(lutNx-1);
    
    qy_min = qy->calcMin();
    qy_max = qy->calcMax();
    qy_range = fabs(qy_max - qy_min);    
    qy_stepsize = qy_range/(double)(lutNy-1);    
    
    //update class variables to reflect these changes
    //(if values should differ, use smaller one to not lose accuracy)
    p_deltaq = (qx_stepsize<=qy_stepsize) ? qx_stepsize : qy_stepsize;         
    p_qmax = (qx_max<=qy_max) ? qx_max : qy_max;
	
	if (debug()>1) {
		cout << "qx: min=" << qx_min << ", max=" << qx_max << ", range=" << qx_range << ", step size=" << qx_stepsize << endl;  
		cout << "qy: min=" << qy_min << ", max=" << qy_max << ", range=" << qy_range << ", step size=" << qy_stepsize << endl;                    
		cout << "deltaq=" << deltaq() << ", qmax=" << qmax() << endl;
	}
}

//----------------------------------------------------------------------------createLookupTable
// rearrange data into a fast lookup table to get values fast using 
// val=lookup(x,y)
// dimensions of the argument 'table' determines the accuracy of the lookup
//----------------------------------------------------------------------------
int CrossCorrelator::createLookupTable( int lutNx, int lutNy ){
    int retval = 0;
    
	if (debug()){ cout << "CrossCorrelator::createLookupTable() begin." << endl; }

	double qx_min, qx_stepsize, qy_min, qy_stepsize = 0.;
	calcLUTvariables( lutNx, lutNy, qx_min, qx_stepsize, qy_min, qy_stepsize );
	
	array2D *myTable = new array2D( lutNx, lutNy );
    myTable->zero();
    
    if ( qx->size()!=qy->size() ) {
        cerr << "Error in createLookupTable! Array sizes don't match: " 
            << "qx=" << qx->size() << ", qy=" << qy->size() << endl;
        retval++;
    } else {
        double ix = 0;
        double iy = 0;
        for (int i = 0; i < qx->size(); i++){           //go through all the data
            //get q-values from qx and qy arrays
            //and determine at what index (ix, iy) to put them in the lookup table
            ix = (qx->get(i)-qx_min) / qx_stepsize;
            iy = (qy->get(i)-qy_min) / qy_stepsize;
            
            //and fill table at the found coordinates with the data index
            //overwriting whatever value it had before
		    //(the add-one-half->floor trick is to achieve reasonable rounded integers)
            myTable->set( (int) floor(ix+0.5), (int) floor(iy+0.5), i );
            
            /////////////////////////////////////////////////////////////////////////////////////////
            //ATTENTION: THIS METHOD WILL LEAD TO A LOSS OF DATA,
            //ESPECIALLY FOR SMALL TABLE SIZES,
            //BUT IT WILL BUY A LOT OF SPEED IN THE LOOKUP PROCESS
            //--> this should be improved to a more precise version, 
            //    maybe even one that allows the lookup(x,y) to interpolate
            //    for that to work, we need to find the four closest points in the data or so
            //    (for instance, instead of one index, the table could contain 
            //    a vector of all applicable indices)
            /////////////////////////////////////////////////////////////////////////////////////////
            
        }//for
		
		//after all is done, set the zero element to zero 
		//to make sure a failure of lookup() doesn't result in an acutal value
		myTable->set(0, 0, 0);		
    }//if

		
	if (debug()>2){ cout << "table = " << myTable->getASCIIdata() << endl; }
    if (debug()){ cout << "CrossCorrelator::createLookupTable() done." << endl; }
	
	//replace table
	setLookupTable(myTable);
	
	delete myTable;	
    return retval;
}



