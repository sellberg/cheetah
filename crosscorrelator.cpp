/*
 *  crosscorrelator.cpp
 *
 *  XCCA analysis code 
 *  
 *  calculates the average SAXS intensity
 *  calculates the angular auto-correlation or cross-correlation
 *
 *  algorithm 1 (direct calculation) written by Jonas Sellberg, 2011-08-20 
 *  algorithm 2 (Fourier method) written by Jan Feldkamp, 2011-08-22
 *
 *  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
 */

#include "crosscorrelator.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <cmath>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "fouriertransformer.h"
#include "arraydataIO.h"				// can be used, isn't be mandatory here


//=================================================================================
//
// constructors & destructor
//
//=================================================================================

//----------------------------------------------------------------------------constructor with C-style arrays
CrossCorrelator::CrossCorrelator( float *dataCArray, float *qxCArray, float *qyCArray, 
									int arraylength, int nphi, int nq1, 
									int nq2, int16_t *maskCArray){
    initPrivateVariables();
    
    //set basic properties for the size of the arrays
    setArraySize(arraylength);
	p_nPhi = nphi;
	p_nLag = (int) ceil(p_nPhi/2.0+1);
	
	//check, user wants to run in 2D xaca or full 3D xcca mode
	if (nq2 == 0){
		p_xcca_enable = false;
		p_nQ = nq1;
		p_autoCorrelation = new array2D( nQ(), nLag() );
	} else {
		p_xcca_enable = true;
		p_nQ = (nq1 > nq2) ? nq1 : nq2;
		p_crossCorrelation = new array3D( nQ(), nQ(), nLag() );
	}

	//check, if a mask was given
	if (maskCArray == NULL){
		p_mask_enable = false;
	}else{
		p_mask_enable = true;
		setMask(maskCArray, arraylength);	
	}
	
    //copy data from array over to internal data structure
	setData(dataCArray, arraylength);

	//set q-calibration arrays (or generate a default)
	if ( qxCArray && qyCArray ){
		setQx(qxCArray, arraylength);
		setQy(qyCArray, arraylength);	
	}else{
		initDefaultQ();
	}

	initInternalArrays();
}


//----------------------------------------------------------------------------constructor with arraydata objects
CrossCorrelator::CrossCorrelator( array2D *dataArray, array2D *qxArray, array2D *qyArray, 
									int nphi, int nq1,
									int nq2, int16_t *maskCArray ){
	initPrivateVariables();
	
	if ( (dataArray->size() != qyArray->size()) || (dataArray->size() != qyArray->size()) ){
		cerr << "Warning in CrossCorrelator constructor! Array sizes don't match" << endl;
		cerr << "qxArray size = " << qxArray->size() << endl;
		cerr << "qyArray size = " << qyArray->size() << endl;
		cerr << "dataArray size = " << dataArray->size() << endl;
	}
	
    //set basic properties for the size of the arrays
    setArraySize(dataArray->size());
	p_nPhi = nphi;
	p_nLag = (int) ceil(p_nPhi/2.0+1);
	
	//check, user wants to run in 2D xaca or full 3D xcca mode
	if (nq2 == 0){
		p_xcca_enable = false;
		p_nQ = nq1;
		p_autoCorrelation = new array2D( nQ(), nLag() );
	} else {
		p_xcca_enable = true;
		p_nQ = (nq1 > nq2) ? nq1 : nq2;
		p_crossCorrelation = new array3D( nQ(), nQ(), nLag() );
	}

	//check, if a mask was given
	if (maskCArray == NULL){
		p_mask_enable = false;
	}else{
		p_mask_enable = true;
		setMask(maskCArray, dataArray->size());	
	}
	
    //copy data from array over to internal data structure
    setData( dataArray );

	//set q-calibration arrays (or generate a default)
	if ( qxArray && qyArray ){
		setQx( qxArray );
		setQy( qyArray );	  
	}else{
		initDefaultQ();
	}
	
	initInternalArrays();
}


//----------------------------------------------------------------------------destructor
CrossCorrelator::~CrossCorrelator(){
	//free memory for objects
	delete p_data;
	delete p_qx;
	delete p_qy;
	delete p_mask;
	
    delete p_polar;
	delete p_corr;
	delete p_mask_polar;
	delete p_mask_corr;
	delete p_table;
	
	delete p_q;
	delete p_phi;	
	delete p_qave;
	delete p_iave;
	delete p_phiave;
	
	delete p_crossCorrelation;
	delete p_autoCorrelation;
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
	p_nQ = 0;
	p_nPhi = 0;
	p_nLag = 0;
	p_mask_enable = false;
	p_xcca_enable = true;
    p_outputdir = "";

	p_data = NULL;
	p_qx = NULL;
	p_qy = NULL;
	p_mask = NULL;

	p_polar = NULL;
	p_corr = NULL;
	p_mask_polar = NULL;
	p_mask_corr = NULL;
    p_table = NULL;

	p_q = NULL;
	p_phi = NULL;	
	p_qave = NULL;
	p_iave = NULL;
	p_phiave = NULL;
	
	p_crossCorrelation = NULL;
	p_autoCorrelation = NULL;
	
	p_qxmin = 0;
	p_qymin = 0;
	p_qxdelta = 0;
	p_qydelta = 0;
	
    p_debug = 0;      
	//check1D = NULL;
}

void CrossCorrelator::initInternalArrays(){
	//allocate all other internal objects
	p_table = new array2D(50, 50);
	
	p_q = new array1D(arraySize());
	p_phi = new array1D(arraySize());
	
	p_qave = new array1D(nQ());
	p_iave = new array1D(nQ());
	p_phiave = new array1D(nPhi());
}


//----------------------------------------------------------------------------initDefaultQ
//if a q-calibration was not given to this class from the outside
//create a default one here
//----------------------------------------------------------------------------
void CrossCorrelator::initDefaultQ(){
    cout << "Initializing p_qx and p_qx vectors with default values." << endl;

    array1D *default_qx = new array1D(arraySize());
    array1D *default_qy = new array1D(arraySize());
    
    //set new values for deltaq and qmax
    p_deltaq = 1;                            
    p_qmax = arraySize()/2.*deltaq();

    for (int i=0; i<arraySize(); i++){
        default_qx->set(i, -qmax()+deltaq()*i );
        default_qy->set(i, -qmax()+deltaq()*i );
    }
    
    int onedim = (int) floor(sqrt(arraySize()));
    for (int i = 0; i < onedim; i++){
        for (int j = 0; j < onedim; j++){
            default_qx->set( i*onedim+j,    (i - (onedim-1)/2.)*deltaq() );
            default_qy->set( i*onedim+j,    (j - (onedim-1)/2.)*deltaq() );
        }
    }
	
	setQx(default_qx);
	setQy(default_qy);
		
	delete default_qx;
    delete default_qy;
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
			//data()->readFromHDF5( filename );
			break;
		case 1:
			data()->readFromRawBinary( filename );
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
			qx()->set( counter, qx_val );
			qy()->set( counter, qy_val );
			counter++;
		}
	}
}

//=================================================================================
//
// CROSS-CORRELATION ALGORITHM 1
//
//=================================================================================

//----------------------------------------------------------------------------calculatePolarCoordinates
// calculate polar coordinates from cartesian coordinate system
void CrossCorrelator::calculatePolarCoordinates()
{	
	// calculate phi for each pixel and bin angles with correct deltaphi
	for (int i=0; i<arraySize(); i++) {
		
		double phii = p_phi->get(i);
		double qxi = qx()->get(i);
		double qyi = qy()->get(i);
		
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
		} else if (phii > (nPhi()-1)*deltaphi()) {
			if (p_debug >= 2) cout << "phii: " << phii << ", sampleAngle(phii): " << round(phii/deltaphi()) << endl;
		}
		
		p_phi->set( i, round(phii/deltaphi()) * deltaphi() );
		
	}
}

//----------------------------------------------------------------------------calculatePolarCoordinates
void CrossCorrelator::calculatePolarCoordinates(double start_q, double stop_q)
{
	// calculate dependent variables for cheetah.ini control
	p_qmin = start_q;
	p_qmax = stop_q;
	p_deltaq = (qmax()-qmin())/(nQ()-1);	// make sure deltaq samples start and stop
	p_deltaphi = (double) 2.0*M_PI/(p_nPhi);	// make sure deltaphi samples exactly an interval of 2PI
	
	if (p_debug >= 1) cout << "qmin: " << qmin() << ", qmax: " << qmax() << ", p_deltaq: " << p_deltaq << ", p_nQ: " << p_nQ << ", p_deltaphi: " << p_deltaphi << ", p_nPhi: " << p_nPhi << ", p_nLag: " << p_nLag << endl;
	
	// calculate phi for each pixel and bin angles with correct deltaphi
	for (int i=0; i<arraySize(); i++) {
		double phii = 0;
		double qxi = qx()->get(i);
		double qyi = qy()->get(i);
		
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
		} else if (phii > (nPhi()-1)*deltaphi()) {
			if (p_debug >= 2) cout << "phii: " << phii << ", sampleAngle(phii): " << round(phii/deltaphi()) << endl;
		}
		
		p_phi->set( i, round(phii/deltaphi()) * deltaphi() );
		
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
		p_q->set(i, round(sqrt( (qx()->get(i)*qx()->get(i))+(qy()->get(i)*qy()->get(i)) ) / deltaq()) * deltaq() );
		if (p_debug >= 3) printf("q[%d]: %f\n",i,p_q->get(i));
	}
	
	// angular average for each |q|
	if (p_debug >= 1) printf("# of steps: %d\n",nQ());
	if (p_debug >= 2) printf("average SAXS intensity:\n");
	
	for (int i=0; i<nQ(); i++) {
		p_qave->set( i, qmin()+i*deltaq() );
		double itot = 0; // reset summed intensity
		int counter = 0; // reset counter
		for (int j=0; j<arraySize(); j++) {
			if ( p_q->get(j) == p_qave->get(i) && (!p_mask_enable || mask()->get(j))) {
				itot += data()->get(j);
				counter++;
			}
		}
		p_iave->set( i, itot/counter );
		if (p_debug >= 2) cout << "Q: " << p_qave->get(i) << ",   \t# pixels: " << counter << ",\tI: " << p_iave->get(i) << endl;
	}
}



//----------------------------------------------------------------------------calculateXCCA
void CrossCorrelator::calculateXCCA(){
	
	if (p_debug >= 1) cout << "deltaPhi: " << deltaphi() << endl;
	if (p_debug >= 1) cout << "# of angles: " << nPhi() << endl;
	for (int i=0; i<nPhi(); i++) {
		p_phiave->set( i, phimin()+i*deltaphi() );
	}
	
	// create array of the speckle pattern with the correct binning
	array2D *speckle = new array2D( nQ(), nPhi() );
	
	// create array over pixel counts for each sampled q and phi
	array2D *pixelCount = new array2D( nQ(), nPhi() );
	array2D *pixelBool = new array2D( nQ(), nPhi() );

	if (p_debug >= 1) printf("calculating speckle arrays...\n");

	for (int i=0; i<arraySize(); i++) {
		if (!p_mask_enable || mask()->get(i)) {
			int qIndex = (int) round((p_q->get(i)-qmin())/deltaq()); // the index in qave[] that corresponds to q[i]
			int phiIndex = (int) round((p_phi->get(i)-phimin())/deltaphi()); // the index in phiave[] that corresponds to phi[i]
			if (p_debug >= 3) printf("qIndex: %d, phiIndex: %d\n", qIndex, phiIndex);
			if (qIndex >= 0 && qIndex < nQ() && phiIndex >= 0 && phiIndex < nPhi()) { // make sure qIndex and phiIndex is not out of array bounds
				speckle->set(qIndex, phiIndex, speckle->get(qIndex, phiIndex) + data()->get(i) );
				pixelCount->set(qIndex, phiIndex, pixelCount->get(qIndex,phiIndex)+1);
				if (pixelBool->get(qIndex, phiIndex) != 1) {
					pixelBool->set(qIndex, phiIndex, 1);
				}
			} else if (p_debug >= 2) printf("POINT EXCLUDED! qIndex: %d, phiIndex: %d\n", qIndex, phiIndex);
		}
	}
		
	// subtract the average SAXS intensity from the speckle array
	array2D *speckleNorm = new array2D( nQ(), nPhi() );
	
	for (int i=0; i<nQ(); i++) { // make sure each element is initially zero
		for (int j=0; j<nPhi(); j++) {
			if (pixelBool->get(i,j) != 0) {
				speckle->set(i, j, speckle->get(i,j) / pixelCount->get(i,j) );
				speckleNorm->set(i, j, speckle->get(i,j) - p_iave->get(i) );
			}
			if (p_debug >= 2) printf("q: %f, phi: %f --> bool: %f, count: %f\n", p_qave->get(i), p_phiave->get(j), pixelBool->get(i, j), pixelCount->get(i, j));
		}
	}
	
	// calculate cross-correlation array and normalization array for cross-correlation
	if (p_debug >= 1) cout << "# of angular lags: " << nLag() << endl;
	
	if (p_xcca_enable) {
		
		if (p_debug >= 1) printf("starting main loop to calculate cross-correlation...\n");
		
		for (int i=0; i<nQ(); i++) { // q1 index
			for (int j=0; j<nQ(); j++) { // q2 index 
				for (int k=0; k<nLag(); k++) { // phi lag => phi2 index = (l+k)%nPhi
					double norm = 0;
					for (int l=0; l<nPhi(); l++) { // phi1 index
						int phi2_index = (l+k)%nPhi();
						p_crossCorrelation->set(i,j,k, p_crossCorrelation->get(i,j,k) + speckleNorm->get(i,l)*speckleNorm->get(j, phi2_index) );
						norm += pixelBool->get(i,l)*pixelBool->get(j, phi2_index);
					}
			
					if (norm) {
						// normalize the cross-correlation array with the average SAXS intensity and the calculated normalization constant					
						// p_crossCorrelation->set(i, j, k, p_crossCorrelation->get(i,j,k) / ( norm*p_iave->get(i)*p_iave->get(j)) );
					
						//normalize by standard deviation (or the zeroth element of the correlation)
						p_crossCorrelation->set(i, j, k, p_crossCorrelation->get(i,j,k) / ( norm*p_crossCorrelation->get(i, j, 0)) );
					}
				}
			}
		}
		
	} else {
		if (p_debug >= 1) printf("starting main loop to calculate cross-correlation...\n");
		
		for (int i=0; i<nQ(); i++) { // q index
			double stdev = 0;
			for (int k=0; k<nLag(); k++) { // phi lag => phi2 index = (l+k)%nPhi
				double norm = 0;
				for (int l=0; l<nPhi(); l++) { // phi1 index
					int phi2_index = (l+k)%nPhi();
					p_autoCorrelation->set(i,k, p_autoCorrelation->get(i,k) + speckleNorm->get(i,l)*speckleNorm->get(i, phi2_index) );
					norm += pixelBool->get(i,l)*pixelBool->get(i, phi2_index);
				}
				if (norm != 0) {
					if (k == 0) {
						stdev = p_autoCorrelation->get(i, 0)/norm;
					}
					if (stdev != 0) {
						//p_autoCorrelation->set(i, k, p_autoCorrelation->get(i,k) / (norm*p_iave->get(i)*p_iave->get(i)) );
						p_autoCorrelation->set(i, k, p_autoCorrelation->get(i,k) / (norm*stdev) );
					}
				} else {
					p_autoCorrelation->set(i, k, -2);
				}
			}		
		}
	}
	
	delete pixelBool;
	delete pixelCount;
	
	delete speckle;
	delete speckleNorm;
	
	if (p_debug >= 1) printf("done calculating cross-correlation...\n");
}



//----------------------------------------------------------------------------calculateXACA
void CrossCorrelator::calculateXACA() {
	
	if (p_debug >= 1) cout << "deltaPhi: " << deltaphi() << endl;
	if (p_debug >= 1) cout << "# of angles: " << nPhi() << endl;
	for (int i=0; i<nPhi(); i++) {
		p_phiave->set( i, phimin()+i*deltaphi() );
	}
	
	// create array of the speckle pattern with the correct binning
	array2D *speckle = new array2D( nQ(), nPhi() );
	
	// create array over pixel counts for each sampled q and phi
	array2D *pixelCount = new array2D( nQ(), nPhi() );
	array2D *pixelBool = new array2D( nQ(), nPhi() );
	
	if (p_debug >= 1) printf("calculating speckle arrays...\n");
	
	for (int i=0; i<arraySize(); i++) {
		if (!p_mask_enable || mask()->get(i)) {
			int qIndex = (int) round((p_q->get(i)-qmin())/deltaq()); // the index in qave[] that corresponds to q[i]
			int phiIndex = (int) round((p_phi->get(i)-phimin())/deltaphi()); // the index in phiave[] that corresponds to phi[i]
			if (p_debug >= 3) printf("qIndex: %d, phiIndex: %d\n", qIndex, phiIndex);
			if (qIndex >= 0 && qIndex < nQ() && phiIndex >= 0 && phiIndex < nPhi()) { // make sure qIndex and phiIndex is not out of array bounds
				speckle->set(qIndex, phiIndex, speckle->get(qIndex, phiIndex) + data()->get(i) );
				pixelCount->set(qIndex, phiIndex, pixelCount->get(qIndex,phiIndex)+1);
				if (pixelBool->get(qIndex, phiIndex) != 1) {
					pixelBool->set(qIndex, phiIndex, 1);
				}
			} else if (p_debug >= 2) printf("POINT EXCLUDED! qIndex: %d, phiIndex: %d\n", qIndex, phiIndex);
		}
	}
	
	// subtract the average SAXS intensity from the speckle array
	array2D *speckleNorm = new array2D( nQ(), nPhi() );
	
	for (int i=0; i<nQ(); i++) { // make sure each element is initially zero
		for (int j=0; j<nPhi(); j++) {
			if (pixelBool->get(i,j) != 0) {
				speckle->set(i, j, speckle->get(i,j) / pixelCount->get(i,j) );
				speckleNorm->set(i, j, speckle->get(i,j) - p_iave->get(i) );
			}
			if (p_debug >= 2) printf("q: %f, phi: %f --> bool: %f, count: %f\n", p_qave->get(i), p_phiave->get(j), pixelBool->get(i, j), pixelCount->get(i, j));
		}
	}
	
	// calculate cross-correlation array and normalization array for cross-correlation
	if (p_debug >= 1) cout << "# of angular lags: " << nLag() << endl;
	
	array2D *normalization = new array2D( nQ(), nLag() );
	
	if (p_debug >= 1) printf("starting main loop to calculate cross-correlation...\n");
	
	for (int i=0; i<nQ(); i++) { // q index
		for (int k=0; k<nLag(); k++) { // phi lag => phi2 index = (l+k)%nPhi
			for (int l=0; l<nPhi(); l++) { // phi1 index
				p_crossCorrelation->set(i,i,k, p_crossCorrelation->get(i,i,k) + speckleNorm->get(i,l)*speckleNorm->get(i, (l+k)%nPhi()) );
				normalization->set(i, k, normalization->get(i, k) + pixelBool->get(i,l)*pixelBool->get(i, (l+k)%nPhi()) );
			}
		}		
	}
	
	// normalize the cross-correlation array with the average SAXS intensity and the calculated normalization constant
	for (int i=0; i<nQ(); i++) {
		for (int k=0; k<nLag(); k++) {
			if (normalization->get(i,k) != 0) {
				p_crossCorrelation->set(i, i, k, p_crossCorrelation->get(i,i,k) / (normalization->get(i,k)*p_iave->get(i)*p_iave->get(i)) );
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

//----------------------------------------------------------------------------data
array1D *CrossCorrelator::data() const {
	return p_data;
}

void CrossCorrelator::setData( array1D *data ) {
	if (p_data) {
		delete p_data;
	}
	p_data = new array1D( *data );
	setArraySize( data->size() );
}

void CrossCorrelator::setData( array2D *data2D ) {
	if (p_data) {
		delete p_data;
	}
	p_data = new array1D( data2D );
	setArraySize( data2D->dim1() * data2D->dim2() );
}

void CrossCorrelator::setData( float *dataCArray, unsigned int size ){
	if (p_data) {
		delete p_data;
	}
	p_data = new array1D( dataCArray, size );
	setArraySize( size );
}

//----------------------------------------------------------------------------qx
array1D *CrossCorrelator::qx() const {
	return p_qx;
}

void CrossCorrelator::setQx( array1D *qx ){
	if (p_qx) {
		delete p_qx;
	}
	p_qx = new array1D( *qx );
}

void CrossCorrelator::setQx( array2D *qx ){
	if (p_qx) {
		delete p_qx;
	}
	p_qx = new array1D( qx );
}

void CrossCorrelator::setQx( float *qxArray, unsigned int size ){
	if (p_qx) {
		delete p_qx;
	}
	p_qx = new array1D( qxArray, size );
}


//----------------------------------------------------------------------------qy
array1D *CrossCorrelator::qy() const {
	return p_qy;
}

void CrossCorrelator::setQy( array1D *qy ) {
	if (p_qy) {
		delete p_qy;
	}
	p_qy = new array1D( *qy );
}

void CrossCorrelator::setQy( array2D *qy ) {
	if (p_qy) {
		delete p_qy;
	}
	p_qy = new array1D( qy );
}

void CrossCorrelator::setQy( float *qyArray, unsigned int size ){
	if (p_qy) {
		delete p_qy;
	}
	p_qy = new array1D( qyArray, size );
}

//----------------------------------------------------------------------------mask
array1D *CrossCorrelator::mask() const {
	return p_mask;
}

void CrossCorrelator::setMask( array1D *mask ) {
	if (p_mask) {
		delete p_mask;
	}
	p_mask = new array1D( *mask );
	p_mask_enable = true;
	normalizeMask();
}

void CrossCorrelator::setMask( array2D *mask ) {
	if (p_mask) {
		delete p_mask;
	}
	p_mask = new array1D( mask );
	p_mask_enable = true;
	normalizeMask();
}

void CrossCorrelator::setMask(int16_t *maskCArray, unsigned int size){
	if (p_mask) {
		delete p_mask;
	}
	p_mask = new array1D( maskCArray, size );
	p_mask_enable = true;
	normalizeMask();
}

void CrossCorrelator::normalizeMask(){
	//if the mask has a value other than zero (any value, pos. or neg.), set that to exactly 1
	for (int i = 0; i < mask()->size(); i++){
		double tolerance = 0.000001;	// float comparison, so allow for some error
		if ( mask()->get_atAbsoluteIndex(i) < tolerance && mask()->get_atAbsoluteIndex(i) > -tolerance ){
			mask()->set_atAbsoluteIndex(i, 0.);
		}else{
			mask()->set_atAbsoluteIndex(i, 1.);
		}
	}
}

void CrossCorrelator::setLookupTable( array2D *LUT ){
	calcLUTvariables( LUT->dim1(), LUT->dim2() );
	p_table->copy(*LUT);								//store a copy locally
}

void CrossCorrelator::setLookupTable( const int *cLUT, unsigned int LUT_dim1, unsigned int LUT_dim2 ){
	calcLUTvariables( LUT_dim1, LUT_dim2 );	
	unsigned int tablelength = LUT_dim1*LUT_dim2;
	p_table->setDim1(LUT_dim1);						//set dimensions
	p_table->setDim2(LUT_dim2);
	p_table->arraydata::copy(cLUT, tablelength);		//store a copy in 'table' locally
}


array2D *CrossCorrelator::polar() const {
	return p_polar;
}

array2D *CrossCorrelator::corr() const {
	return p_corr;
}

array2D *CrossCorrelator::mask_polar() const {
	return p_mask_polar;
}

array2D *CrossCorrelator::mask_corr() const {
	return p_mask_corr;
}

array2D *CrossCorrelator::lookupTable() const{
	return p_table;
}

//----------------------------------------------------------------------------arraySize
int CrossCorrelator::arraySize() const{
	return p_arraySize;
}

void CrossCorrelator::setArraySize( int arraySize_val ){
	p_arraySize = arraySize_val;
}

//----------------------------------------------------------------------------matrixSize
int CrossCorrelator::matrixSize() const{
	return (int)floor(sqrt( (double)arraySize() ));            //assuming a square image
}

void CrossCorrelator::setMatrixSize( int matrixSize_val ){
	setArraySize( matrixSize_val*matrixSize_val );
	updateDependentVariables();
}

//----------------------------------------------------------------------------qmin/qmax
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

//----------------------------------------------------------------------------phimin/phimax
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

//----------------------------------------------------------------------------outputdir
string CrossCorrelator::outputdir(){
    return p_outputdir;
}

void CrossCorrelator::setOutputdir( std::string dir ){
	const char lastchar = dir.at( dir.size()-1 );
	if( lastchar != '/' ){		//if last character is not a slash, append one
		dir += '/';
	}
    p_outputdir = dir;
}

//----------------------------------------------------------------------------debug
int CrossCorrelator::debug() const {
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

int CrossCorrelator::nQ() const{						//getter only, dependent variable
	return p_nQ;
}

int CrossCorrelator::nPhi() const{						//getter only, dependent variable
	return p_nPhi;
}

int CrossCorrelator::nLag() const{					//getter only, dependent variable
	return p_nLag;
}

void CrossCorrelator::updateDependentVariables(){		//update the values that depend on qmax and matrixSize
	// FINE BINNING
//	p_deltaq = 2*qmax()/(matrixSize()-1);
//	p_nQ = int(qmax()/p_deltaq+1+0.001);
//	p_deltaphi = 2*atan(1/(2*(p_nQ-1.0)));
//	p_nPhi = int(2*round(M_PI/p_deltaphi)); // make sure p_nPhi is even (exclude 2PI)
//	p_deltaphi = (double) 2.0*M_PI/(p_nPhi); // make sure deltaphi samples exactly an interval of 2PI
//	p_nLag = (int) round(p_nPhi/2.0+1);
	
	// COARSE BINNING
	p_deltaq = 20*qmax()/(matrixSize()-1);
	p_nQ = int(qmax()/p_deltaq+1+0.001);
	p_deltaphi = 2*atan(1/(2*(p_nQ-1.0)));
	p_nPhi = int(2*round(M_PI/p_deltaphi)); // make sure p_nPhi is even (exclude 2PI)
	p_deltaphi = (double) 2.0*M_PI/(p_nPhi); // make sure deltaphi samples exactly an interval of 2PI
	p_nLag = (int) round(p_nPhi/2.0+1);
	
	if (p_debug >= 1) {
		cout << "CrossCorrelator::updateDependentVariables done. qmax: " << qmax() << ", p_deltaq: " << p_deltaq << ", p_nQ: " << p_nQ << ", p_deltaphi: " << p_deltaphi << ", p_nPhi: " << p_nPhi << ", p_nLag: " << p_nLag << endl;
	}
}

double CrossCorrelator::getQave(unsigned index) const {
	return p_qave->get(index);
}

double CrossCorrelator::getPhiave(unsigned index) const {
	return p_phiave->get(index);
}

double CrossCorrelator::getIave(unsigned index) const {
	return p_iave->get(index);
}

double CrossCorrelator::getCrossCorrelation(unsigned index1, unsigned index2) const {
	return p_autoCorrelation->get(index1,index2);
}

double CrossCorrelator::getCrossCorrelation(unsigned index1, unsigned index2, unsigned index3) const {
	return p_crossCorrelation->get(index1,index2,index3);
}










//=================================================================================
//
// CROSS-CORRELATION ALGORITHM 2 (_FAST)
//
//=================================================================================



//----------------------------------------------------------------------------createLookupTable
// re-arrange data into a fast lookup table to get values fast using 
// val=lookup(x,y)
// dimensions of the argument 'table' determines the accuracy of the lookup
//----------------------------------------------------------------------------
int CrossCorrelator::createLookupTable( int lutNx, int lutNy ){
	if (debug()){ cout << "CrossCorrelator::createLookupTable() begin." << endl; }

    int retval = 0;
	
	calcLUTvariables( lutNx, lutNy );
	
	array2D *myTable = new array2D( lutNx, lutNy );
    myTable->zero();
	
	//initialize table with -1, which can never be a real index
	const double initval = -1;
	myTable->addValue(initval);			
    
    if ( qx()->size()!=qy()->size() ) {
        cerr << "Error in createLookupTable! Array sizes don't match: " 
            << "qx=" << qx()->size() << ", qy=" << qy()->size() << endl;
        retval++;
    } else {
        double ix = 0;
        double iy = 0;
        for (int i = 0; i < qx()->size(); i++){           //go through all the data
            //get q-values from qx and qy arrays
            //and determine at what index (ix, iy) to put them in the lookup table
			double qxi = qx()->get(i);
			double qyi = qy()->get(i);
            ix = (qxi-p_qxmin) / p_qxdelta;
            iy = (qyi-p_qymin) / p_qydelta;
            
            //and fill table at the found coordinates with the data index
            //overwriting whatever value it had before
		    //(the add-one-half->floor trick is to achieve reasonable rounded integers)
            myTable->set( (int) floor(ix+0.5), (int) floor(iy+0.5), double(i) );
            
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

		
	if (debug()>1){ cout << "table = " << myTable->getASCIIdata() << endl; }
    if (debug()){ cout << "CrossCorrelator::createLookupTable() done." << endl; }
	
	//sanity check: how many, if any places are still at the default value?
	int initcount = 0;
	for (int i = 0; i < myTable->size(); i++){
		if ( myTable->get_atAbsoluteIndex(i) == initval ){
			initcount++;
		}
	}
	if (initcount > 0) {
		cout << "==============================================================================================" << endl;
		cout << "Warning in createLookupTable! There are still " << initcount << " unassigned values (" 
			<< ((double)initcount)/myTable->size()*100 << "%) in the lookup table." << endl;
		cout << "==============================================================================================" << endl;
	}
	
	//replace current table
	setLookupTable(myTable);
	
	delete myTable;	
    return retval;
}



//----------------------------------------------------------------------------calcLUTvariables
// calculate some variables needed to fill the lookup table
// and update the private variables
// which are important parameters for the lookup() function
//----------------------------------------------------------------------------
void CrossCorrelator::calcLUTvariables( int lutNx, int lutNy ){
	p_qxmin = qx()->calcMin();
    p_qxmax = qx()->calcMax();
    double qx_range = fabs(p_qxmax - p_qxmin);
    p_qxdelta = qx_range/(double)(lutNx-1);
    
    p_qymin = qy()->calcMin();
    p_qymax = qy()->calcMax();
    double qy_range = fabs(p_qymax - p_qymin);    
    p_qydelta = qy_range/(double)(lutNy-1);    
    
	if (debug()>1) {
		cout << "qx: min=" << p_qxmin << ", max=" << p_qxmax << ", range=" << qx_range << ", p_qxdelta=" << p_qxdelta << endl;  
		cout << "qy: min=" << p_qymin << ", max=" << p_qymax << ", range=" << qy_range << ", p_qxdelta=" << p_qydelta << endl;                    
	}
}



//---------------------------------------------------------------------------- normalizeToSAXS
int CrossCorrelator::normalizeToSAXS(){
	for(int ct=0; ct < polar()->dim2(); ct++){
		
		//get one row out of the polar coordinates
		array1D *row = new array1D( polar()->dim1() );
		polar()->getRow( ct, row );		
		
		if (row->size() == 0) { 
			cerr << "ERROR in CrossCorrelator::normalizeToSAXS. Row has size zero." << endl; 
			return 1; 
		}

		//calculate average intensity in that row
		double avg = 0;
		if (!p_mask_enable){			//no bad pixels --> just take the average of the row
			avg = row->calcAvg();
		}else{							//mask --> leave out bad pixels
			double sum = 0.;
			double valid = 0;
			for (int i = 0; i < row->size(); i++) {
				sum += row->get_atAbsoluteIndex(i);
				valid += mask_polar()->get(i, ct);		// if mask has a 1 here, the point is valid
			}
			avg = (valid > 0) ? (sum/((double)valid)) : 0;
		}
		
					
		row->subtractValue(avg);				//subtract
		
//		if (avg != 0){
//			//row->multiplyByValue(1/avg/avg);		
//			row->divideByValue(avg);				//divide
//		} else {
//			row->zero();							//set to zero (harsh, but avoids catastrophe)
//		}
		polar()->setRow(ct, row);
		
		delete row;
	}
	return 0;
}

//----------------------------------------------------------------------------transform to polar coordinates
int CrossCorrelator::calculatePolarCoordinates_FAST(){
	return calculatePolarCoordinates_FAST( 0, nQ() );
}

int CrossCorrelator::calculatePolarCoordinates_FAST( double start_q, double stop_q ) {
	
	int number_phi = nPhi();
	int number_q = nQ();
	
	//find smallest common q in all directions, use that as upper limit, if necessary
	double abs_q_max = HUGE_VALF;
	if ( abs_q_max > fabs(p_qxmax)) { abs_q_max = fabs(p_qxmax); }
	if ( abs_q_max > fabs(p_qymax)) { abs_q_max = fabs(p_qymax); }
	if ( abs_q_max > fabs(p_qxmin)) { abs_q_max = fabs(p_qxmin); }
	if ( abs_q_max > fabs(p_qymin)) { abs_q_max = fabs(p_qymin); }
    stop_q = abs_q_max<stop_q ? abs_q_max : stop_q;  	//take the smallest value of the two

	if( debug()>1 ){ 
		cout << "CrossCorrelator::calculatePolarCoordinates_FAST" << endl; 
		cout << "varying scattering vector q from " << start_q << " to " <<  stop_q << " in " << number_q << " steps, "
			<< "and angle phi full circle in " << number_phi << " steps." << endl;
	}
    int retval = 0;
	
	//apply mask, if there is one
	//mask should be a map of ones (good pixels) and zeros (bad pixels)
	//multiplication with the real data then represents masking
	if ( p_mask_enable ){
		data()->multiplyByArrayElementwise( mask() );
	}


	//create new array2D to store polar coordinate representation
	delete p_polar;
    p_polar = new array2D( number_phi, number_q );
	if (!p_polar){
		cerr << "Error in CrossCorrelator::calculatePolarCoordinates_FAST. polar couldn't be allocated." << endl;
		return 1;
	}
	
	if (p_mask_enable) {
		//create new array2D to store polar coordinate representation
		delete p_mask_polar;
		p_mask_polar = new array2D( number_phi, number_q );
		if (!p_mask_polar){
			cerr << "Error in CrossCorrelator::calculatePolarCoordinates_FAST. p_mask_polar couldn't be allocated." << endl;
			return 1;
		}
	}
	
	int novalue_count = calculatePolarCoordinates_FAST( data(), polar(), number_phi, number_q, start_q, stop_q );
	if ( novalue_count > 0 ){
		cout << "Couldn't assign a true value in " << novalue_count << " cases (" 
			<< ((double)novalue_count)/number_q/number_phi*100 << "%) in the polar coordinate image of " 
			<< number_phi << " x " << number_q << " pixels." << endl;
    }
	
	if (p_mask_enable) {
		calculatePolarCoordinates_FAST( mask(), mask_polar(), number_phi, number_q, start_q, stop_q  );
	}

	//normalize data
	normalizeToSAXS();

	return retval;
}


//----------------------------------------------------------------------------calculatePolarCoordinates_FAST
// returns the number of times the lookup has failed
int CrossCorrelator::calculatePolarCoordinates_FAST( array1D* cartesian1D, array2D* polar2D, 
													int number_phi, int number_q, double start_q, double stop_q) const {
    
	double xcoord = 0.;
    double ycoord = 0.;
    double data_value = 0.;
    double q = 0.;
    double p = 0.;
    int qcounter = 0;
    int pcounter = 0;
	int novalue_count = 0;
 	
	//in principle, these could be arbitrary, but the FFT approach currently assumes wrap-around data
	const int start_phi = 0;
	const int stop_phi = 360;
	
	//leave out the upper boundary, e.g. 200 <= q < 400, 0 <= phi < 360
	double step_q = (stop_q - start_q)/((double)number_q-1);
    double step_phi = (stop_phi - start_phi)/((double)number_phi-1);
	
	if (step_q < 0)
        cerr << "Error in CrossCorrelator::calculatePolarCoordinates_FAST -- step_r value " 
            << step_q << " is smaller than zero." << endl;
    if (step_phi < 0)
        cerr << "Error in CrossCorrelator::calculatePolarCoordinates_FAST -- step_phi value " 
            << step_phi << " is smaller than zero." << endl;
    
	
	for(q = start_q, qcounter=0; qcounter < number_q; q+=step_q, qcounter++){                        // q: for all rings/q-values
        if(debug()){cout << "#" << qcounter << ",q=" << q << "  " << std::flush;}
		for(p = start_phi, pcounter=0; pcounter < number_phi; p+=step_phi, pcounter++){				// phi: go through all angles

            //find lookup coordinates
			xcoord = q * cos(p*M_PI/180.);
			ycoord = q * sin(p*M_PI/180.);

            //lookup that value in original scattering data
			try {
            	data_value = lookup( xcoord, ycoord, cartesian1D );
            } catch (int e) {
				//cout << "Exception: " << e << endl;
				data_value = 0;									//setting the 0 here will introduce false correlations!!!
				novalue_count++;
				if ( e == -1){
					cout << "Exception " << e << ": interpolation in lookup failed." << endl;
				}else if ( e == 0){
					cout << "Exception " << e << ": table not allocated." << endl;
				}else{					
					//cases e = 1-4: asked for (x,y) values that are larger or smaller than the actual image
					//fail silently... data value remains zero. 
					
					if (p_mask_enable){
						//additionally, include this point in the mask to ignore it later in the analysis
						mask_polar()->set(pcounter, qcounter, 1);
					}
				}
			}
			
			//assign the new values 
			polar2D->set(pcounter, qcounter, data_value);
			//polar->set(pcounter, qcounter, data_value * q);	//(functional determinant q to account for bins growing as q grows)
		}
	}
	
    if (debug()) {
		cout << "polarCoordinates done. dimensions=(" << polar()->dim1() << " x " << polar()->dim2() << ")" << endl;
	    //write output of the intermediate files? (all zero by default, turn on for debugging or whatever)
		if (false){ cout << "data: " << cartesian1D->getASCIIdata() << endl; }
		if (false){ cout << "polar: " << polar2D->getASCIIdata() << endl; }
		//if (false){ io->writeToTiff( outputdir()+"polar.tif", polar ); }      //(arraydataIO needed for this)      
	}
	
	return novalue_count;		
}




//----------------------------------------------------------------------------lookup
double CrossCorrelator::lookup( double xcoord, double ycoord, array1D *dataArray ) const {

    double value = 0.;    //return data value at the given coordinates
    int index = 0;      //to do that, the index in the data is determined first

    double xc = (xcoord-p_qxmin) / ((double)p_qxdelta);
    double yc = (ycoord-p_qymin) / ((double)p_qydelta);

    int ix = (int) floor( xc + 0.5 );		// round to nearest integer
    int iy = (int) floor( yc + 0.5 );
	
	bool enable_warnings = false;
    if ( !lookupTable() ){
        cerr << "Error in lookup! No lookup table was allocated." << endl;
		throw 0;
    } else if ( ix < 0 ){
        if (enable_warnings){
			cerr << "Error in lookup! xcoord=" << xcoord << " is too small.";
        	cerr << "   ix=" << ix << " < 0)" << endl;
		}
		throw 1;
    } else if ( ix >= lookupTable()->dim1() ){
        if (enable_warnings){
	        cerr << "Error in lookup! xcoord=" << xcoord << " is too large.";
	        cerr << "   ix=" << ix << " >= " << lookupTable()->dim1() << " (table dimx)" << endl;
		}
		throw 2;
    } else if ( iy < 0 ){
        if (enable_warnings){
	        cerr << "Error in lookup! ycoord=" << ycoord << " is too small.";
    	    cerr << "   iy=" << iy << " < 0)" << endl;
		}
		throw 3;
    } else if ( iy >= lookupTable()->dim2() ){
        if (enable_warnings){
			cerr << "Error in lookup! ycoord=" << ycoord << " is too large.";
			cerr << "   iy=" << iy << " >= " << lookupTable()->dim2() << " (table dimy)" << endl;
    	}
		throw 4;
	} else {
	    //create lookup index from original coordinates (assuming the data is properly centered)
        index = (int) lookupTable()->get(ix, iy);
		if ( index > 0 ){
			value = dataArray->get( index );
		} else {
			// if the value returned was negative, that means there is no lookup value assigned to this pair of (ix, iy)
			// this could be because the table is too large
			// solution: create a binned value from the surrounding pixels
			int valid = 0;			//number of valid indices
			double sum = 0.;		//sum of retrieved values (at valid indices)
			vector<int> indices;
			indices.push_back( (int) lookupTable()->get(ix+1, iy) );
			indices.push_back( (int) lookupTable()->get(ix-1, iy) );
			indices.push_back( (int) lookupTable()->get(ix, iy+1) );
			indices.push_back( (int) lookupTable()->get(ix, iy-1) );
			for (int i = 0; i < indices.size(); i++){
				int tmpindx = indices.at(i);
				if ( tmpindx > 0 ){
					valid++;
					sum += dataArray->get( tmpindx );
				}
			}
			if ( valid > 0 ){
				value = sum / ((double)valid);
			} else {
				cout << "WARNING in lookup()! Couldn't find a valid value in LUT for coordinates (" << xcoord << ", " << ycoord << "). " << endl;
				// search region would have to be expanded in a complete implementation....... 
				// the way it is should be enough for the moment, though, if the table isn't way too big......
				throw -1;
			}
			
		}
		
		//-----------------------!!!DEBUG!!!
		//make a hot pixel out of the one that was just looked up
		//careful! this actually changes the data in the input array!!
		//---> therefore, uncomment one of the following line only for debugging
		//p_dataCArray[index] = 100000.;
		//p_data->set( index, 1000);
		//-----------------------!!!DEBUG!!!
    }

    if( debug()>2 ){
        cout << "lookup (" << xcoord << ", " << ycoord 
            << ") --> LUT: (xc,yc)=(" << xc << ", " << yc 
            << ") ==> (" << ix << ", " << iy << ") "
            << "--> index=" << index << ", --> val=" << value << endl;
    }

    return value;
}






//----------------------------------------------------------------------------calculateXCCA_FAST
//
//normalize, according to eq. (3) in Altarelli, et al., PRB, 82, 104207 (2010)
// general, eq. (3):  C = ( <I(q1, phi) * I(q2, phi+deltaphi)> - <I(q1)>*<I(q2)> ) / ( <I(q1)>*<I(q2)> )
// autocorr, eq. (1): C =    <I(q, phi) * I(q, phi+deltaphi)>  /  ( <I(q)>^2 )  -  1
//
int CrossCorrelator::calculateXCCA_FAST(){
	if(debug()>1){ 
		cout << "CrossCorrelator::calculateXCCA_FAST ("<< polar()->dim1() << ", " << polar()->dim2() << ")" << endl; 
	}
	
	//do some sanity checks first
	if (!polar()) { 
		cerr << "No polar coordinate representation found. Call calculatePolarCoordinates_FAST() first." << endl; 
		throw 1; 
	}
	if (polar()->dim1()==0 || polar()->dim2()==0) { throw 2; }
	
	//create a new output array 'corr'
	delete p_corr;
	p_corr = new array2D( polar()->dim1(), polar()->dim2() );
	if (!p_corr) { throw 3; }
	
	//create a new output array 'mask_corr', if needed
	if ( p_mask_enable ){
		delete p_mask_corr;
		p_mask_corr = new array2D( polar()->dim1(), polar()->dim2() );
		if (!p_mask_corr) { throw 4; }
	}
	
	//correlate data in polar, write result to corr
	try {
		if (!p_xcca_enable){		// autocorrelation only
			autocorrelateFFT( polar(), corr() );
			
			//write to 2D dataset: only auto-correlation output
			for (int i=0; i<nQ(); i++) { 			// q1 index 
				for (int k=0; k<nLag(); k++) { 			// phi lag => phi2 index = (l+k)%nPhi
					p_autoCorrelation->set(i, k, corr()->get(k, i) );
				}
			}
		}else{						// full-blown cross-correlation
			array3D *corr3D = new array3D;
			crosscorrelateFFT( polar(), corr3D );
			
			//write to 3D dataset: full cross-correlation output (WARNING!!! NOT TESTED YET!!!)
			for (int i=0; i<nQ(); i++) { 			// q1 index
				for (int j=0; i<nQ(); i++) { 		// q2 index
					for (int k=0; k<nLag(); k++) { 		// phi lag => phi2 index = (l+k)%nPhi
						p_crossCorrelation->set(i, j, k, corr3D->get(k, j, i) );
					}
				}
			}
			delete corr3D;
		}
	} catch(int e) {
		cerr << "Exception " << e << " thrown by autocorrelateFFT_usingArrays: ";
		switch (e) {
			case 1:
				cerr << "Single row 'f' could not be allocated.";
				break;
			case 2:
				cerr << "Could not get single row.";
				break;
			case 3:
				cerr << "Could not calculate correlation.";
				break;
			default:
				cerr << "Unkown exception.";
		}//switch
		cerr << " Aborting. " << endl;
		throw e;	//re-throw
	}
	
	//---------------------------------debug: this if-loop is unneccessary and is just for academic purposes right now
	if ( p_mask_enable ){
		//correlate mask, write result to mask_corr
		try {
			autocorrelateFFT( mask_polar(), mask_corr() );
		} catch(int e) {
			cerr << "Exception caught in autocorrelateFFT_usingArrays (mask enabled)" << endl;
			throw e;	//re-throw
		}
	}
	
	if(debug()){ cout << endl << "CrossCorrelator::calculateXCCA_FAST done." << endl; }			
    return 0;
}




//----------------------------------------------------------------------------autocorrelateFFT_byRow
int CrossCorrelator::autocorrelateFFT(array2D *polar2D, array2D *corr2D) const {
	int retval, fail = 0;
	
	//calculate the auto-correlation for all rings
	for(int q_ct=0; q_ct < polar2D->dim2(); q_ct++){

		//get one row out of the polar coordinates (fixed q)
		array1D *f = new array1D( polar2D->dim1() );
		if (!f){ throw 1; } 
		polar2D->getRow( q_ct, f );		
		if (f->size() == 0) { throw 2; }
				
		if (debug()>1){ cout << "   #" << q_ct << ", f before FFT: " << f->getASCIIdata() << endl; }
		if (debug()){ cout << q_ct << " " << std::flush; }
		
		//perform autocorrelation --> compute via FFT
		if ( !p_mask_enable ){												// no mask

			// initialize imaginary part to zero --> f is real
		    array1D *f_imag = new array1D(f->size());						
	
			FourierTransformer *ft = new FourierTransformer();
			fail = ft->autocorrelation( f, f_imag );
			
			delete ft;
			delete f_imag;
			if (fail){ cerr << "Error in CrossCorrelator::autocorrelateFFT. Transform (forward) failed." << endl; throw 3; }

		}else{
			//get one mask row, corresponding to the data one
			
			if (!f || !f->data() ) {
				cerr << "CrossCorrelator::autocorrelateFFT. Input not allocated." << endl;
				return 1;
			}
			if (f->size()==0 ) {
				cerr << "CrossCorrelator::autocorrelateFFT. Input has size zero." << endl;
				return 2;
			}
			
			FourierTransformer *ft = new FourierTransformer();
			
			//-----calculate 
			array1D *f_imag = new array1D(f->size());							// initialized to zeros --> f is real
			int f_fail = ft->transformForward( f, f_imag );
			if (f_fail){ cerr << "Error in CrossCorrelator::autocorrelateFFT. autocorr failed." << endl; throw 4; }
			
			
			//-----prepare for back transform
			array1D *filter = new array1D( f->size() );
			int maxorder = 700;
			int winstart = int( f->size()/2. - maxorder - 1 );
			int winstop = int( f->size()/2. + maxorder + 1 );
			filter->zero();
			filter->ones( winstart, winstop );
			//cout << "filter: " << filter->getASCIIdata();
			f->multiplyByArrayElementwise( filter );	//apply filter in Fourier space
			f_imag->multiplyByArrayElementwise( filter );
			
			//mag squared
			for (int i=0; i < f->size(); i++) {
				f->set( i,   ( f->get(i)*f->get(i) + f_imag->get(i)*f_imag->get(i) )  );
				f_imag->set( i, 0 );		    
			}
			
			//-----transform back
			// after inverse transform, result is stored in original argument array mf
			int inv_fail = ft->transformInverse( f, f_imag );
			if (inv_fail){ cerr << "Error in CrossCorrelator::autocorrelateFFT(masked). transformInverse failed." << endl; throw 4; }

			delete ft;    
			delete f_imag;
			delete filter;
		}
		
		//normalize to zero-correlation
		f->divideByValue( f->get(0) );
		
		//feed result into corr
		corr2D->setRow( q_ct, f );
		
		if (debug()>1){ cout << "   #" << q_ct << ", f after FFT: " << f->getASCIIdata() << endl; }
		delete f;
		f = NULL;
	}//for
	
	return retval;
}



//----------------------------------------------------------------------------autocorrelateFFT_byRow
int CrossCorrelator::crosscorrelateFFT(array2D *polar2D, array3D *corr3D) const {
	int retval, fail = 0;

	//calculate the cross-correlation for all combination of rings
	for(int fq_ct=0; fq_ct < polar2D->dim2(); fq_ct++){
		
		//get one row out of the polar coordinates
		array1D *f = new array1D( polar2D->dim1() );
		if (!f){ throw 1; } 
		polar2D->getRow( fq_ct, f );		
		if (f->size() == 0) { throw 2; }
		
		//inner loop!
		for(int gq_ct=0; gq_ct < polar2D->dim2(); gq_ct++){

			//get second row out of the polar coordinates
			array1D *g = new array1D( polar2D->dim1() );
			if (!f){ throw 1; } 
			polar2D->getRow( gq_ct, g );		
			if (g->size() == 0) { throw 2; }
					
			if (debug()>1){ 
				cout << "   #" << fq_ct << ", f before FFT: " << f->getASCIIdata() << endl; 
			}
			if (debug()){ cout << fq_ct << " " << std::flush; }
			
			//perform CROSS-correlation --> compute via FFT
			if (!f || !g) {
				cerr << "CrossCorrelator::crosscorrelateFFT. Input not allocated." << endl;
				return 1;
			}
			if (f->size()==0 || g->size()==0) {
				cerr << "CrossCorrelator::crosscorrelateFFT. Input has size zero." << endl;
				return 2;
			}
			
			array1D *f_real = new array1D( *f );
			array1D *f_imag = new array1D(f->size());				// initialized to zeros --> f is real
			array1D *g_real = new array1D( *g );
			array1D *g_imag = new array1D(g->size());				// initialized to zeros --> g is real   
			
			// calculate the cross correlation, result is returned in f_real, f_imag
			FourierTransformer *ft = new FourierTransformer();
			fail = ft->crosscorrelation( f_real, f_imag, g_real, g_imag );
			delete ft;
			
			if (fail){ cerr << "Error in CrossCorrelator::crosscorrelateFFT. Transform (forward) failed." << endl; fail = 3; }
			
			// return result in original argument arrays (not really interested in imaginary part right now)
			f->copy( *f_real );
			//g->copy( *f_imag );
			
			delete f_real;
			delete f_imag;    
			delete g_real;
			delete g_imag;

			//normalize to zero-correlation
			f->divideByValue( f->get(0) );
			//g->divideByValue( f->get(0) );
	
			//feed result into corr
			//corr3D->setRow( fq_ct, gq_ct, f );		//(no such function exists, yet --> revisit later)
			for (int i = 0; i < f->size(); i++){
				corr3D->set( fq_ct, gq_ct, i, f->get(i) );
			}
			
			if (debug()>1){
				cout << "   #" << fq_ct << ", f before FFT: " << f->getASCIIdata() << endl; 
			}
			
			delete g;
			g = NULL;
		}//for: inner loop
		
		delete f;
		f = NULL;
	}//for: outer loop
	
	return retval;
}

