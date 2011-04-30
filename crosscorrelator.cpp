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





// ***********************************************************************************
// constructors & destructor
// ***********************************************************************************
CrossCorrelator::CrossCorrelator(int arraylength){
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
	
	// create arrays to store polar coordinates for the angular cross-correlation
	q = new array1D(arraySize());
	phi = new array1D(arraySize());
	
	qave = new array1D(arraySize());
	iave = new array1D(arraySize());
	phiave = new array1D(samplingAngle());

}

CrossCorrelator::CrossCorrelator( int16_t *dataCArray, int arraylength ){

    //set basic properties, just like the default case
    setArraySize( arraylength );
    setQmax(1.);
	
    //special feature: copy data from array over to internal data structure
    data = new array1D( dataCArray, arraylength );
    
    //allocate all other internal objects
    qx = new array1D(arraySize());
	qy = new array1D(arraySize());
    
	q = new array1D(arraySize());
	phi = new array1D(arraySize());
	
	qave = new array1D(arraySize());
	iave = new array1D(arraySize());
	phiave = new array1D(samplingAngle());    
    
}

CrossCorrelator::~CrossCorrelator(){
	//free memory for objects
	delete data;
	
	delete qx;
	delete qy;
	
	delete q;
	delete phi;
	
	delete qave;
	delete iave;
	delete phiave;
}



// ***********************************************************************************
// load data from file
//		type 0: HDF5 file (output from Anton's hit finder code
//		type 1: raw 16-bit binary file
// ***********************************************************************************
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
//			data[counter] = double(buffer[counter]);
//			qx[counter] = qx_val;
//			qy[counter] = qy_val;
			counter++;
		}
	}
}

//----------------------------------------------------------------------------initWithTestPattern
void CrossCorrelator::initWithTestPattern( int type ){

    //create a nice test pattern
    int sizex = 100;
    int sizey = 100;
    array2D *test = new array2D(sizex, sizey);
    test->generateTestPattern(type);

    //convert it to 1D object and feed it to 'data'
    delete data;
    data = new array1D( test );
    
    //set the right center
    setCenterX( (double)sizex/2 );
    setCenterY( (double)sizey/2 );
    
    //write tiff image to check what the pattern looks like
    std::string filename = "/Users/feldkamp/Desktop/testpattern.tif";
    test->writeToTiff( filename );
    
    delete test;
}
    



// ***********************************************************************************
// CROSS-CORRELATION ALGORITHM 0
// ***********************************************************************************

//----------------------------------------------------------------------------calculatePolarCoordinates
// calculate polar coordinates from cartesian coordinate system
void CrossCorrelator::calculatePolarCoordinates()
{	
	// calculate phi for each pixel and bin angles with correct resolution
	for (int i=0; i<arraySize(); i++) {
		
		double phii = phi->get(i);
		double qxi = qx->get(i);
		double qyi = qy->get(i);
		
		// setup UHP
		if (qxi == 0) { // make sure that the column with qx = 0 has angle 0 (q = 0 is assumed to have phi = 0)
			phi->set(i, 0);
		} else {
			phi->set(i, atan(qxi/qyi) ); // pixel 5050-5099 have -0 and 5101-5150 have 0 in angle (5100 has NaN), WHY??? ANSWER: That column has qx = 0, so the angle should be 0 or PI. If qy = 0 and qx != 0, atan gives the correct result, but only for the UHP! Need to add PI for all LHP!
		}

		// correct LHP by adding PI
		if (qyi < 0) {
			phi->set(i, phii + M_PI);
		}
		
		if (phii<0) { // make sure the angle is between 0 and 2PI
			phi->set(i, phii+2*M_PI);
			//phi[i] += 2*M_PI;
		}
		
		phi->set( i, round(phii/deltaphi()) * deltaphi() );
		// printf("phi[%d]: %f\n",i,phi[i]);
	}
}	


//----------------------------------------------------------------------------calculateSAXS
void CrossCorrelator::calculateSAXS(){
	
	
	// using SAXS average for all shots to calculate cross-correlation 
    //or just the SAXS from the specific shots will give different results. 
    //The first choice is probably preferable. Here, the second one is performed.
	printf("calculating average SAXS intensity...\n");
	
	// calculate |q| for each pixel and bin lengths with correct resolution
	for (int i=0; i<arraySize(); i++) {
		q->set(i, round(sqrt( (qx->get(i)*qx->get(i))+(qy->get(i)*qy->get(i)) ) / deltaq()) * deltaq() );
		// printf("q[%d]: %f\n",i,q[i]);
	}
	
	// angular average for each |q|
	printf("# of steps: %d\n",samplingLength());
	printf("average SAXS intensity:\n");
	
	int counter = 0;
	for (unsigned int i=0; i<samplingLength(); i++) {
		qave->set(i, i*deltaq());
		double itot = 0;
		counter = 0; // reset counter
		for (int j=0; j<arraySize(); j++) {
			if ( q->get(j) == qave->get(i) ) {
				itot += data->get(j);
				counter++;
			}
		}
		iave->set( i, itot/counter );
		cout << "Q: " << qave->get(i) << ",   \t# pixels: " << counter << ",\tI: " << iave->get(i) << endl;
	}
}



//----------------------------------------------------------------------------calculateXCCA
void CrossCorrelator::calculateXCCA(){

	cout << "deltaPhi: " << deltaphi() << endl;
	cout << "# of angles: " << samplingAngle() << endl;
	for (int i=0; i<samplingAngle(); i++) {
		phiave->set( i, i*deltaphi() );
	}
	
	// create array over pixel counts for each sampled q and phi
	//unsigned pixelCount[(int)samplingLength()][samplingAngle];
	//unsigned pixelBool[samplingLength()][samplingAngle];
	array2D *pixelCount = new array2D( samplingLength(), samplingAngle() );
	array2D *pixelBool = new array2D( samplingLength(), samplingAngle() );
	
	for (int i=0; i<arraySize(); i++) {
		int qIndex = int(q->get(i)/deltaq()+0.001); // the index in qave[] that corresponds to q[i]
		int phiIndex = int(phi->get(i)/deltaphi()+0.001); // the index in phiave[] that corresponds to phi[i]
		// printf("qIndex: %d, phiIndex: %d\n",qIndex,phiIndex);
		if (qIndex < samplingLength() && phiIndex < samplingAngle()) { // make sure qIndex is not larger 
                                                                        //than the samplingLength 
                                                                        //(corners where q > 1 are excluded)
			pixelCount->set(qIndex, phiIndex, pixelCount->get(qIndex,phiIndex)+1);
			if (pixelBool->get(qIndex, phiIndex) != 1) {
				pixelBool->set(qIndex, phiIndex, 1);
			}
		} // else printf("POINT EXCLUDED! qIndex: %d, phiIndex: %d\n",qIndex,phiIndex);
	}
	//  for (int i=0; i<samplingLength(); i++) {
	//    for (int j=0; j<samplingAngle; j++) { 
	//      printf("q: %f, phi: %f --> bool: %u, count: %u\n",qave[i],phiave[j],pixelBool[i][j],pixelCount[i][j]);
	//    }
	//  }
	
	// calculate normalization constant for cross-correlation
	cout << "# of angular lags: " << samplingLag() << endl;
	printf("calculating normalization array...\n");
	array3D *normalization = new array3D( samplingLength(), samplingLength(), samplingLag() );
	
	for (int i=0; i<samplingLength(); i++) { // q1 index
		for (int j=0; j<samplingLength(); j++) { // q2 index 
			for (int k=0; k<samplingLag(); k++) { // phi lag => phi2 index = (l+k)%samplingAngle
				for (int l=0; l<samplingAngle(); l++) { // phi1 index
					// printf("phi2: %d\n",(l+k)%samplingAngle);
					//normalization[i][j][k] += pixelBool[i][l]*pixelBool[j][(l+k)%samplingAngle()];
					normalization->set(i, j, k, normalization->get(i, j, k) + pixelBool->get(i,l)*pixelBool->get(j, (l+k)%samplingAngle()) );
				}
			}
		}
	}

	
	// create array with average SAXS intensity
	// NO NEED since we already have average intensity stored in iave[]
	
	// create array of the speckle pattern with the correct binning
	array2D *speckle = new array2D( samplingLength(), samplingAngle() );
	
	// *** LOOPS COPIED FROM pixelCount/pixelBool ***
	for (int i=0; i<arraySize(); i++) {
		int qIndex = int(q->get(i)/deltaq()+0.001); // the index in qave[] that corresponds to q[i]
		int phiIndex = int(phi->get(i)/deltaphi()+0.001); // the index in phiave[] that corresponds to phi[i]
		// printf("qIndex: %d, phiIndex: %d\n",qIndex,phiIndex);
		if (qIndex < samplingLength() && phiIndex < samplingAngle()) { // make sure qIndex is not larger than the samplingLength 
                                                                        //(corners where q > 1 are excluded)
			speckle->set(qIndex, phiIndex, speckle->get(qIndex, phiIndex) + data->get(i) );
		} // else printf("POINT EXCLUDED! qIndex: %d, phiIndex: %d\n",qIndex,phiIndex);
	}
	// *** END OF LOOP COPIED FROM pixelCount/pixelBool ***
	
	// subtract the average SAXS intensity from the speckle array
	array2D *speckleNorm = new array2D( samplingLength(), samplingAngle() );
	
	for (int i=0; i<samplingLength(); i++) { // make sure each element is initially zero
		for (int j=0; j<samplingAngle(); j++) {
			if (pixelBool->get(i,j) != 0) {
				speckle->set(i, j, speckle->get(i,j) / pixelCount->get(i,j) );
				speckleNorm->set(i, j, speckle->get(i,j) - iave->get(i) );
			}
		}
	}
	
	// create cross-correlation array
	printf("starting main loop to calculate cross-correlation...\n");
	array3D *crossCorrelation = new array3D( samplingLength(), samplingLength(), samplingLag() );
	
	// calculate cross-correlation
	for (int i=0; i<samplingLength(); i++) { // q1 index
		for (int j=0; j<samplingLength(); j++) { // q2 index 
			for (int k=0; k<samplingLag(); k++) { // phi lag => phi2 index = (l+k)%samplingAngle
				for (int l=0; l<samplingAngle(); l++) { // phi1 index
					crossCorrelation->set(i,j,k, crossCorrelation->get(i,j,k) + speckleNorm->get(i,l)*speckleNorm->get(j, (l+k)%samplingAngle()) );
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
	
	delete speckle;
	delete speckleNorm;
	delete crossCorrelation;
	
	printf("done calculating cross-correlation...\n");
}



//----------------------------------------------------------------------------writeSAXS
void CrossCorrelator::writeSAXS(){
	// write cross-correlation and average SAXS intensity to binary
	printf("writing data to file...\n");
	FILE *filePointerWrite;
	
	filePointerWrite = fopen("f909-q0-xcca.bin","w+");
	
	/*
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
	*/
	
	fclose(filePointerWrite);
	
	cout << "writeSAXS done" << endl;
}



//----------------------------------------------------------------------------dumpResults
void CrossCorrelator::dumpResults( std::string filename ){
	cout << "Writing results to TIFF file '" << filename << "'." << endl;
    
    //convert data to 2D array class (using known dimensions)
    int dim1 = matrixSize();                // not exact.... will have to change this soon
    int dim2 = matrixSize();
    array2D *dataTwoD = new array2D( data, dim1, dim2 );
    
    //write to tiff image
    dataTwoD->writeToTiff( filename );
    delete dataTwoD;
}



//----------------------------------------------------------------------------writeXCCA
void CrossCorrelator::writeXCCA(){
	//not implemented yet. 
	//right now, all output is handled by writeSAXS()
}




//----------------------------------------------------------------------------printRawData
// print raw data after having read from file
void CrossCorrelator::printRawData(uint16_t *buffer,long lSize) {
	for (int i=0; i<lSize; i++) {
		printf("%u ",buffer[i]);
	}
	printf("\n");
}



// ***********************************************************************************
// setters and getters for private variables
// ***********************************************************************************
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

double CrossCorrelator::qmax() const{
	return p_qmax;
}

void CrossCorrelator::setQmax( double qmax_val ){
	p_qmax = qmax_val;
	updateDependentVariables();
}


double CrossCorrelator::centerX() const{
    return p_centerX;
}

void CrossCorrelator::setCenterX( double cen_x ){
    p_centerX = cen_x;
}

double CrossCorrelator::centerY() const{
    return p_centerY;
}

void CrossCorrelator::setCenterY( double cen_y ){
    p_centerY = cen_y;
}


double CrossCorrelator::deltaq() const{						//getter only, dependent variable
	return p_deltaq;
}

double CrossCorrelator::deltaphi() const{						//getter only, dependent variable
	return p_deltaphi;
}

int CrossCorrelator::samplingLength() const{					//getter only, dependent variable
	return p_samplingLength;
}

int CrossCorrelator::samplingAngle() const{					//getter only, dependent variable
	return p_samplingAngle;
}

int CrossCorrelator::samplingLag() const{					//getter only, dependent variable
	return p_samplingLag;
}

void CrossCorrelator::updateDependentVariables(){		//update the values that depend on qmax and matrixSize
	p_deltaq = 2*qmax()/(matrixSize()-1);
    
	p_samplingLength = int(1/p_deltaq+1+0.001);
	p_deltaphi = 2*atan(1/(2*(p_samplingLength-1.0)));
	p_samplingAngle = (int) floor(2*M_PI/p_deltaphi);
	p_samplingLag = (int) ceil(p_samplingAngle/2.0)+2;
}




// ***********************************************************************************
// CROSS-CORRELATION ALGORITHM 1 (_FAST)
// ***********************************************************************************

//----------------------------------------------------------------------------transform to polar coordinates
int CrossCorrelator::calculatePolarCoordinates_FAST(array2D* polar){

    cout << "calculatePolarCoordinates_FAST" << endl;

    int retval = 0;
    
    
    //write output of the intermediate files?
    int output_data_flag = 0;               
    int output_polar_flag = 1;

    //some of the following variables should be set specifically for each dataset
    //this is currently still in testing mode
    //needs to be changed, soon -----> put into ini file (or so)
    double start_r = 10;
    double stop_r = 45;
    int number_r = 35;
    double step_r = (stop_r - start_r)/number_r;

    double start_phi = 0;
    double stop_phi = 360;
    int number_phi =  90;
    double step_phi = (stop_phi - start_phi)/number_phi;
    
    if (step_r < 0)
        cerr << "Error in CrossCorrelator::calculatePolarCoordinates_Jan -- step_r value " 
            << step_r << " is smaller than zero." << endl;
    if (step_phi < 0)
        cerr << "Error in CrossCorrelator::calculatePolarCoordinates_Jan -- step_phi value " 
            << step_phi << " is smaller than zero." << endl;
    
    if (polar){
        delete polar;
        polar = NULL;
    }
    polar = new array2D(number_phi, number_r);

    //DEBUG!!!!!!!!!!!!
    check1D = new array1D(*data);
    array2D *polarSampling = new array2D(2*(stop_r+1), 2*(stop_r+1));

    double xcoord = 0.;
    double ycoord = 0.;
    double value = 0.;
    double r = 0.;
    double p = 0.;
    int rcounter = 0;
    int pcounter = 0;
    
	for(r = start_r, rcounter=0; rcounter < number_r; r+=step_r, rcounter++){                        // r: for all rings/radii
        cout << "ring r=" << r << ", #" << rcounter << endl;
		for(p = start_phi, pcounter=0; pcounter < number_phi; p+=step_phi, pcounter++){				// phi: go through all angles

            //find lookup coordinates
			xcoord = r * cos(p*M_PI/180);
			ycoord = r * sin(p*M_PI/180);
            polarSampling->set(xcoord+stop_r, ycoord+stop_r, 65535);
			
            //lookup that value in original scattering data
            value = lookup( xcoord + centerX(), ycoord + centerY() );                 //COMMENT THIS BACK IN FOR REAL CASE
            
			//assign the new values (note the functional determinant r)
			polar->set(pcounter, rcounter, value * r);
		}
	}

    if (output_data_flag) {
//        data->writeToTiff("/Users/feldkamp/Desktop/polar.tif");
        cout << "data: " << data->getASCIIdata() << endl;
    }
    
    if (output_polar_flag) {
        polar->writeToTiff("/Users/feldkamp/Desktop/polar.tif");
//        cout << "polar: " << polar->getASCIIdata() << endl;
    }
    
    array2D *check2D = new array2D( check1D, 100, 100 );
    check2D->writeToTiff("/Users/feldkamp/Desktop/check2D.tif");
    polarSampling->writeToTiff("/Users/feldkamp/Desktop/polarSampling.tif");
    delete check2D;
    delete check1D;
    delete polarSampling;
    
    return retval;
}


//---------------------------------------------------------------------------- lookup
double CrossCorrelator::lookup( double xcoord, double ycoord ) const{
    //YET TO IMPLEMENT ACCORDING TO CSPAD GEOMETRY
    
    
    //test case: feed everything back into 2D and be done with it
    int xdim = 100;
    xcoord = round(xcoord);
    ycoord = round(ycoord);
    int index = (int) xcoord + (xdim*ycoord);
    
    double val = data->get( index );
    
    //DEBUG!!!!!!!!!!
    check1D->set(index, 65535);         //keep track of where the value was read

    return val;
}


//----------------------------------------------------------------------------calculate XCCA
int CrossCorrelator::calculateXCCA_FAST( array2D *polar, array2D *corr ){
    int retval = 0;

    cout << "calculateXCCA_FAST" << endl;
    
    if (corr)
        delete corr;
    corr = new array2D( polar->dim1(), polar->dim2() );

    for(int r_ct=0; r_ct < polar->dim2(); r_ct++){							// r_ct: for all rings

        //perform autocorrelation
        array1D *f = new array1D;
        polar->getRow( r_ct, f);
        autocorrelateFFT( f );
        
        cout << "AUTOCORRELATION -- " << f->getASCIIdata() << endl;
	}
    
    corr->writeToTiff("/Users/feldkamp/Desktop/corr_scaled.tif", 1);
    
    return retval;
}




//----------------------------------------------------------------------------correlate
// compute 1D correlation corr(f,g) using FFT, result is written to f
int CrossCorrelator::correlateFFT( array1D *f, array1D *g ){
    int retval = 0;
    
    //Correlation Theorem:
    //multiplying the FT of one function by the complex conjugate 
    //of the FT of the other gives the FT of their correlation
    
    
    FourierTransformer *ft = new FourierTransformer;
        
    // transform f -> F
    array1D *f_real = new array1D( *f );
    array1D *f_imag = new array1D;
    int f_fail = ft->transform( f_real, f_imag );
    delete ft;
    
    if (f_fail)
        cerr << "Error in CrossCorrelator::autocorrelateFFT. Transform (f) failed." << endl;
    
    // transform g -> G
    array1D *g_real = new array1D( *g );
    array1D *g_imag = new array1D;
    int g_fail = ft->transform( g_real, g_imag );
    if (g_fail)
        cerr << "Error in CrossCorrelator::autocorrelateFFT. Transform (g) failed." << endl;
        
    delete ft;
            
    
    // compute F * G_cc (complex conjugate)
    //WARNING::::::::NOT SURE THE MATH IS ALL CORRECT AT THIS POINT! (by JF, 2011-04-29)
    for (int i=0; i<f_real->size(); i++) {
        double val = f_real->get(i)*g_real->get(i) + f_imag->get(i)*(-g_imag->get(i));
        f->set(i, val);
    }

    delete f_real;
    delete f_imag;    
    delete g_real;
    delete g_imag;
    
    
    return retval;
}


//----------------------------------------------------------------------------autocorrelate
// compute 1D autocorrelation corr(f,f) using FFT, result is written to f
int CrossCorrelator::autocorrelateFFT( array1D *f ){         
    int retval = 0;
    
    //-------------------------------------------------------------------------
    //   Wiener-Khinchin Theorem:
    //   the autocorrelation of a function f with itself 
    //   is found by computing the magnitude squared of its Fourier transform
    //-------------------------------------------------------------------------
    
    //transform f -> F
    array1D *real = new array1D( *f );
    array1D *imag = new array1D;
    FourierTransformer *ft = new FourierTransformer;
    int fail = ft->transform( real, imag );
    delete ft;
    
    if (fail)
        cerr << "Error in CrossCorrelator::autocorrelateFFT. Transform failed." << endl;
    
    //compute absolute square |F|^2
    for (int i=0; i<real->size(); i++) {
        double m2 = real->get(i)*real->get(i) + imag->get(i)*imag->get(i);
        f->set(i, m2);
    }
    
    delete real;
    delete imag;
    return retval;   
}


