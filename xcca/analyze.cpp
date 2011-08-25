//
//  test.cpp
//  xcca_commandline
//
//  Created by Feldkamp on 7/21/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//

#include "analyze.h"

#include <iostream>
using std::cout;
using std::endl;

#include <sstream>

#include "crosscorrelator.h"
#include "arraydataIO.h"

Analyzer::Analyzer(){
	p_back = new array2D();
	p_mask = new array2D();
	
	p_back_weight = 1;
	
	flag_subtract_background = false;
	flag_single_correlation_output = false;
	flag_use_mask = false;
}

Analyzer::~Analyzer(){
	delete p_back;
}

// =============================================================
// assuming all files are similar
// i.e. same dimensions, same beam center
//
//	int number_phi = 256;
//	int number_q = 200;
//	int cenX = 260;			//in detector pixels
//	int cenY = 308;	
// =============================================================
int Analyzer::processFiles( vector<string> files, int shiftX, int shiftY, int num_phi, int num_q){
	double start_q = 0;
	double stop_q = num_q;
	double LUTx = 200;
	double LUTy = 200;
	return processFiles(files, shiftX, shiftY, num_phi, num_q, start_q, stop_q, LUTx, LUTy);
}

int Analyzer::processFiles( std::vector<string> files, int shiftX, int shiftY, int num_phi, int num_q, 
					double start_q, double stop_q, int LUTx, int LUTy ){
					
	//prepare array2D's to hold overall averages
	arraydataIO *io = new arraydataIO;
	array2D *detavg = new array2D;
	io->readFromEDF( files.at(0), detavg );
	detavg->zero();							// this is now an image of correct dimensions, all zeros
	int imgdim1 = detavg->dim1();
	int imgdim2 = detavg->dim2();
	
	array2D *detavg_copy = new array2D( imgdim1, imgdim2 );
	array2D *backavg = new array2D( imgdim1, imgdim2 );
	
	//make a centered q-range
	int detX = detavg->dim1();		//detector size	(for example, pilatus = 487 * 619)
	int detY = detavg->dim2();
	
	array2D *qx = new array2D( imgdim1, imgdim2 );
	array2D *qy = new array2D( imgdim1, imgdim2 );
	qx->xrange(-detX/2+shiftX, +detX/2+shiftX);
	qy->yrange(-detY/2+shiftY, +detY/2+shiftY);
	//cout << "qx: " << qx->getASCIIdata();
	//cout << "qy: " << qy->getASCIIdata();	


	//prepare lookup table once, so it doesn't have to be done every time
	CrossCorrelator *lutcc = new CrossCorrelator(detavg, qx, qy, num_phi, num_q);
	lutcc->createLookupTable(LUTx, LUTy);
	array2D *LUT = new array2D( *(lutcc->lookupTable()) );
	io->writeToEDF( outputDirectory()+"LUT.edf", LUT );
	
	array2D *polaravg = new array2D( lutcc->nPhi(), lutcc->nQ() );
	array2D *corravg = new array2D( lutcc->nQ(), lutcc->nLag() );
	
	delete lutcc;
	
	//prepare background file
	background()->multiplyByValue( backgroundWeight() );
	
	//process all files
	unsigned int num_files = (unsigned int)files.size();
	for (int k = 0; k < num_files; k++){
		std::ostringstream osst_num;
		osst_num << k;
		string single_desc = osst_num.str();
		string fn = files.at(k);
		cout << "#" << k << ": ";// << fn << endl;
		
		array2D *image = new array2D;
		io->readFromEDF( fn, image );
		array2D *image_copy = new array2D( *image );
		
		
		CrossCorrelator *cc = new CrossCorrelator(image, qx, qy, num_phi, num_q);
		cc->setOutputdir( outputDirectory() );
		cc->setDebug(0);
		
		if (flag_subtract_background){
			image->subtractArrayElementwise( background() );
			backavg->addArrayElementwise( background() );
		}

		if ( flag_use_mask ){
			cc->setMask( this->mask() );
		}
		
		int alg = 1;
		switch (alg) {
			case 0:
				cout << "XCCA regular" << endl;
				cc->calculatePolarCoordinates();
				cc->calculateSAXS();
				cc->calculateXCCA();	
				
			break;
			case 1:
				cout << "XCCA FAST" << endl;

				cc->setLookupTable( LUT );
				cc->calculatePolarCoordinates_FAST( start_q, stop_q );
				cc->calculateXCCA_FAST();

				if ( flag_single_correlation_output ){
					io->writeToEDF( cc->outputdir()+"polar"+single_desc+".edf", cc->polar() );           		
					io->writeToEDF( cc->outputdir()+"corr"+single_desc+".edf", cc->autoCorr() );
				}
				
				polaravg->addArrayElementwise( cc->polar() );
				corravg->addArrayElementwise( cc->autoCorr() );
			break;
		}
		
		detavg->addArrayElementwise( image );			//sum up
		detavg_copy->addArrayElementwise( image_copy );
		
		delete image;
		delete image_copy;
		delete cc;
	}
	
	detavg->divideByValue( num_files );			//normalize
	detavg_copy->divideByValue( num_files );	//normalize
	backavg->divideByValue( num_files );		//normalize
	polaravg->divideByValue( num_files );		//normalize	
	corravg->divideByValue( num_files );		//normalize

	corravg->transpose();
	
	io->writeToEDF( outputDirectory()+"det_avg.edf", detavg);			// average background-subtracted detector image
	io->writeToEDF( outputDirectory()+"polar_avg.edf", polaravg);		// average image in polar coordinates
	io->writeToEDF( outputDirectory()+"corr_avg.edf", corravg);			// average autocorrelation
	if ( flag_subtract_background ){
//		io->writeToEDF( outputDirectory()+"det_avg_original.edf", detavg_copy);		// no background subtraction
		io->writeToEDF( outputDirectory()+"det_background_avg.edf", backavg);					// just the background
	}


	delete qx;
	delete qy;
	delete detavg;
	delete detavg_copy;
	delete backavg;
	delete polaravg;
	delete corravg;		
	delete LUT;
	delete io;
	
	return 0;
}

void Analyzer::setBackground( array2D *back ){
	if (p_back) {
  		delete p_back;
	}
	p_back = new array2D(*back);
}

array2D* Analyzer::background(){
	return p_back;
}

void Analyzer::setBackgroundWeight( double weight ){
	p_back_weight = weight;
}

double Analyzer::backgroundWeight(){
	return p_back_weight;
}

void Analyzer::setMask( array2D *newmask ){
	if (p_mask) {
  		delete p_mask;
	}
	p_mask = new array2D(*newmask);
	flag_use_mask = true;
}

array2D* Analyzer::mask(){
	return p_mask;
}

void Analyzer::setOutputDirectory( string outdir ){
	const char lastchar = outdir.at( outdir.size()-1 );
	if( lastchar != '/' ){		//if last character is not a slash, append one
		outdir += '/';
	}
	p_out_dir = outdir;
}

string Analyzer::outputDirectory(){
	return p_out_dir;
}



