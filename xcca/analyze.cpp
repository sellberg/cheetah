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

#include <cmath>


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
	
	array2D *polaravg = new array2D( num_phi, num_q );
	array2D *corravg = new array2D( num_phi, num_q );
	
	array2D *mask_polar, *mask_corr = NULL;
	
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
	delete lutcc;
	
	//prepare background file
	background()->multiplyByValue( backgroundWeight() );
	
	//process all files
	unsigned int num_files = files.size();
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

		if (flag_subtract_background){
			image->subtractArrayElementwise( background() );
			backavg->addArrayElementwise( background() );
		}

		if ( flag_use_mask ){
			cc->setMask( this->mask() );
		}
		
		cc->setDebug(0);        
		
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
//					io->writeToTiff( cc->outputdir()+"polar"+single_desc+".tif", polar, 1 );		//dump scaled output from correlation			
//					io->writeToTiff( cc->outputdir()+"corr"+single_desc+".tif", corr, 1 );			//dump scaled output from correlation
					io->writeToEDF( cc->outputdir()+"polar"+single_desc+".edf", cc->polar() );           		
					io->writeToEDF( cc->outputdir()+"corr"+single_desc+".edf", cc->corr() );
				}
				
				polaravg->addArrayElementwise( cc->polar() );
				corravg->addArrayElementwise( cc->corr() );
			break;
		}
		
		detavg->addArrayElementwise( image );			//sum up
		detavg_copy->addArrayElementwise( image_copy );
		
		if ( flag_use_mask ){
			//debug: at the end, write to disk the last mask that was used 
			mask_polar = new array2D( *(cc->mask_polar()) );
			mask_corr = new array2D( *(cc->mask_corr()) );
		}
		
		delete image;
		delete image_copy;
		delete cc;
	}
	
	detavg->divideByValue( num_files );			//normalize
	detavg_copy->divideByValue( num_files );	//normalize
	backavg->divideByValue( num_files );		//normalize
	polaravg->divideByValue( num_files );		//normalize	
	corravg->divideByValue( num_files );		//normalize
	
//	io->writeToTiff("outputDirectory()+"det_avg.edf", detavg, 0);
//	io->writeToTiff("outputDirectory()+"det_avg_original.edf", detavg_copy, 0);
//	io->writeToTiff("outputDirectory()+"back_avg.edf", detavg, 0);
//	io->writeToTiff("outputDirectory()+"polar"+fileinfo+"_avg.edf", polaravg, 1);
//	io->writeToTiff("outputDirectory()+"corr"+fileinfo+"_avg.edf", corravg, 1);

	io->writeToEDF( outputDirectory()+"det_avg.edf", detavg);						// average background-subtracted detector image
	io->writeToEDF( outputDirectory()+"polar_avg.edf", polaravg);		// average image in polar coordinates
	io->writeToEDF( outputDirectory()+"corr_avg.edf", corravg);			// average autocorrelation
	if ( flag_subtract_background ){
//		io->writeToEDF( outputDirectory()+"det_avg_original.edf", detavg_copy);		// no background subtraction
		io->writeToEDF( outputDirectory()+"det_background_avg.edf", backavg);					// just the background
	}
	
	if ( flag_use_mask ){
		io->writeToEDF( outputDirectory()+"mask_polar.edf", mask_polar );           		
		io->writeToEDF( outputDirectory()+"mask_corr.edf", mask_corr );
	}

	delete qx;
	delete qy;
	delete detavg;
	delete detavg_copy;
	delete backavg;
	delete polaravg;
	delete corravg;		
	delete mask_corr;
	delete LUT;
	delete io;
	if (flag_use_mask){delete mask_polar;}
	
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


//=================================================================================================
//=================================================================================================
//=================================================================================================
//=================================================================================================
Test::Test(){
	const char * homepath = getenv("HOME");
	cout << "Initializing Test environment with path '" << homepath << "'" << endl;
	setBase(homepath);
}

Test::Test( string base ){
	setBase(base);
}

Test::~Test(){
}

void Test::setBase( string base ){
	const char lastchar = base.at( base.size()-1 );
	if( lastchar != '/' ){		//if last character is not a slash, append one
		base += '/';
	}
	p_base = base;
}

string Test::base(){
	return p_base;
}

//-------------------------------------------------------- -t1
int Test::testCrossCorrelator( int alg ){						
    cout << "testCrossCorrelator(" << alg << ")" << endl;
    

    CrossCorrelator *cc = new CrossCorrelator();
	int number_q = 200;
	int number_phi = 256;

    cc->setOutputdir( base() );
    cc->initWithTestPattern(500, 500, 4 );		//cases 0 - 4 available
//	cc->initFromFile(base()+"/LCLS_2011_Feb27_r0079_054306_2e71_cspad.h5");
    
    cc->setDebug(1);        //TURN ON DEBUGGING FOR NOW --> a lot of output
    
    switch (alg) {
        case 0:
            cout << "XCCA regular" << endl;
            cc->calculatePolarCoordinates();
            cc->calculateSAXS();
            cc->calculateXCCA();	
            
        break;
        case 1:
            cout << "XCCA FAST" << endl;
			
            cc->createLookupTable(500, 500);
            double start_q = 5*cc->deltaq();
            double stop_q = cc->qmax();

			cc->calculatePolarCoordinates_FAST(start_q, stop_q);
			cc->calculateXCCA_FAST();

			arraydataIO *io = new arraydataIO;
			io->writeToTiff( cc->outputdir()+"polar.tif", cc->polar(), 1 );            //dump scaled output from correlation			
			io->writeToTiff( cc->outputdir()+"corr.tif", cc->corr(), 1 );            //dump scaled output from correlation
			delete io;
			
        break;
    }

    
	delete cc;
    return 0;
}


//-------------------------------------------------------- -t2
int Test::testArrayClasses(){
	cout << "TESTING THE ARRAY CLASSES" << endl;

	cout << "----------array1D------------" << endl;
	array1D *my1Darray = new array1D(5);
	my1Darray->set(0, 2);
	my1Darray->set(2, 5);
	my1Darray->set(6, 5);
	cout << my1Darray->getASCIIdata();
	cout << my1Darray->get(2) << endl;
	cout << my1Darray->get(3) << endl;	
    delete my1Darray;
	
	cout << "----------array2D------------" << endl;    
	array2D *my2Darray = new array2D(5, 10);
	my2Darray->set(0, 0, 2);
	my2Darray->set(2, 8, 5.5);
	my2Darray->set(6, 11, -2);
	cout << my2Darray->getASCIIdata();
	cout << my2Darray->get(2, 8) << endl;
    cout << "filling my2Darray" << endl;	
	for (int j = 0; j < my2Darray->dim2(); j++) {
		for (int i = 0; i < my2Darray->dim1(); i++) {
			my2Darray->set(i, j, j*my2Darray->dim1() + i);
		}
	}	
	cout << my2Darray->getASCIIdata();
	my2Darray->addValue(-1.5);
	cout << "after adding -1.5 --- " << my2Darray->getASCIIdata();
	my2Darray->multiplyByValue(1/2.);
	cout << "after multiplying by 1/2 --- " << my2Darray->getASCIIdata();
	my2Darray->xrange(-12, 8);
	cout << "after xrange(-12, 8) --- " << my2Darray->getASCIIdata();
	my2Darray->yrange(-12, 8);
	cout << "after yrange(-12, 8) --- " << my2Darray->getASCIIdata();
	
	array2D *mask = new array2D( my2Darray->dim1(), my2Darray->dim2() );
	mask->set(1,1,1);
	mask->set(4,4,1);
	my2Darray->multiplyByArrayElementwise( mask );
	cout << "after multiplication (1,1)=1 and (4,4)=1, else=0 --- " << my2Darray->getASCIIdata();
	delete mask;
	delete my2Darray;
    
    
	cout << "----------array2D------------" << endl;
	array3D *my3Darray = new array3D(20, 10, 3);
	my3Darray->set(0, 0, 0, 2);
    my3Darray->set(1, 1, 1, -2);
	my3Darray->set(2, 8, 2, 5.5);
	my3Darray->set(6, 11, 33, -2);
	cout << my3Darray->getASCIIdata();
	cout << my3Darray->get(2, 8, 2) << endl;
	cout << my3Darray->get(0, 1, 1) << endl;
	delete my3Darray;
    
/*    
    cout << "----------forging array1D ---> array2D----------" << endl;
	array1D *one = new array1D(25);
	one->set(0, 2);
	one->set(2, 2);
	one->set(6, 6);
    one->set(24, 24);
    array2D *two = new array2D(one, 5, 5);
    cout << "one: " << one->getASCIIdata() << endl;
    cout << "two: " << two->getASCIIdata() << endl;
    cout << "-->deleting one" << endl;
    delete one;
//    cout << "one: " << one->getASCIIdata() << endl;            //THIS IS SUPPOSED TO CRASH!!!
    cout << "two: " << two->getASCIIdata() << endl;
    cout << "-->deleting two" << endl;
    delete two;                                                
 //   cout << "one: " << one->getASCIIdata() << endl;            //THIS IS SUPPOSED TO CRASH!!!
 //   cout << "two: " << two->getASCIIdata() << endl;            //THIS IS SUPPOSED TO CRASH!!!
 
 
 
    cout << "----------forging array2D ---> array1D----------" << endl;
    two = new array2D(5, 5);
    two->set(0,0,1);
    two->set(1,1,7);
    two->set(2,0,-2);
    two->set(4,4,16);
    one = new array1D(two);
    cout << "two: " << two->getASCIIdata() << endl;  
    cout << "one: " << one->getASCIIdata() << endl; 
    delete two;
//    cout << "two: " << two->getASCIIdata() << endl; 			//THIS IS SUPPOSED TO CRASH!!!
//    cout << "one: " << one->getASCIIdata() << endl; 			
    delete one;  
*/     
    return 0;
}



//-------------------------------------------------------- -t3
int Test::testFourierTrafo(){

    int size = 50;
    array1D *f = new array1D(size);
    array1D *g = new array1D(size);
    array1D *model = new array1D(size);
	
	
    //create typical 2D model
	int delta = 5;
	std::vector<int> seeds;
	seeds.push_back(6);
	seeds.push_back(24);
	seeds.push_back(32);
	seeds.push_back(3);
	seeds.push_back(22);
	seeds.push_back(1);
	for (int i = 0; i < seeds.size(); i++){
	    model->set(seeds.at(i), 1);
    	model->set(seeds.at(i)+delta, 1);
	}
	//f->addValue(1);	//put it on a pedestal
    
	//another model
	for (int i = 0; i < size; i++){
		int amp = 10.;
		int steps = 20.;	//degrees per 'i'
	    model->set(i , amp*sin( i*steps * (M_PI/180) ) );
	}
	
	
	
    cout << "--- BEFORE ALL ---" << endl;
    cout << "model: " << model->getASCIIdata() << endl;
    cout << "g: " << g->getASCIIdata() << endl;
    //f->writeToASCII(base()+"/f_model.txt");
    
	double sum_model = f->calcSum();	
	double avg_model = f->calcAvg();
	
	//subtract avg SAXS intensity
	//	f->subtractValue(avg_model);
	
	
	FourierTransformer *trafo = new FourierTransformer();

    cout << "--- TRANSFORM FORWARD ---" << endl;
	f->copy( *model ); g->zero();			//refresh arrays
	trafo->transformForward( f, g );
    cout << "f: " << f->getASCIIdata() << endl;
    cout << "g: " << g->getASCIIdata() << endl; 

    cout << "--- TRANSFORM INVERSE ---" << endl;
	f->copy( *model ); g->zero();			//refresh arrays
	trafo->transformInverse( f, g );
    cout << "f: " << f->getASCIIdata() << endl;
    cout << "g: " << g->getASCIIdata() << endl; 
	
    cout << "--- MAGNITUDE SQUARED ---" << endl;
	f->copy( *model ); g->zero();			//refresh arrays
	trafo->magnitudeSquared( f, g );
    cout << "f: " << f->getASCIIdata() << endl;
    cout << "g: " << g->getASCIIdata() << endl;       
    //f->writeToASCII(base()+"/f_power.txt");
    
    cout << "--- CORRELATION ---" << endl;
	f->copy( *model ); g->zero();			//refresh arrays
    trafo->autocorrelation(f, g);
    cout << "f: " << f->getASCIIdata() << endl;
    cout << "g: " << g->getASCIIdata() << endl;      
    //f->writeToASCII(base()+"/f_corr.txt");
	
	
	//normalize according to Altarelli
	//f->multiplyByValue( 1/avg_model/avg_model);
	
  	double sum_corr = f->calcSum();
	double avg_corr = f->calcAvg();
	cout << "sum_model=" << sum_model << ", avg_model=" << avg_model << endl;
	cout << "sum_corr=" << sum_corr << ", avg_corr=" << avg_corr << endl;  
	
	delete trafo;
    delete f;
    delete g;
    
	return 0;
}


//-------------------------------------------------------- -t4
int Test::testIO(){
	cout << "TESTING THE ARRAYDATA INPUT/OUTPUT" << endl;
	
	string outfn = "test";
	string infn = "test";
	
	arraydataIO *io = new arraydataIO;

	array2D *wmat = new array2D(10, 6);
	wmat->generateTestPattern( 4 );		//cases 0 - 4 available
	wmat->writeToASCII(base()+outfn+".txt");
	
	cout << "writing to EDF" << endl;
	io->writeToEDF(base()+outfn+".edf", wmat);
	
	cout << "writing to TIFF scaled" << endl;
	io->writeToTiff(base()+outfn+"_scaled.tif", wmat, 1);

	cout << "writing to TIFF unscaled" << endl;
	io->writeToTiff(base()+outfn+"_unscaled.tif", wmat, 0);

	
	cout << "writing to HDF5 (doubles)" << endl;
	io->writeToHDF5(base()+outfn+"_d.h5", wmat, 0);

	cout << "writing to HDF5 (float)" << endl;
	io->writeToHDF5(base()+outfn+"_f.h5", wmat, 1);
	
	cout << "writing to HDF5 (int)" << endl;
	io->writeToHDF5(base()+outfn+"_i.h5", wmat, 2);
	
	cout << "writing to HDF5 (default)" << endl;
	io->writeToHDF5(base()+outfn+".h5", wmat, 3);
	
	delete wmat;
	
	array2D *rmat = new array2D;

	cout << "reading from EDF" << endl;
	io->readFromEDF(base()+infn+".edf", rmat);
	cout << rmat->getASCIIdata();
	
	cout << "reading from TIFF scaled" << endl;
	io->readFromTiff(base()+infn+"_scaled.tif", rmat);
	cout << rmat->getASCIIdata();

	cout << "reading from TIFF unscaled" << endl;
	io->readFromTiff(base()+infn+"_unscaled.tif", rmat);
	cout << rmat->getASCIIdata();

	cout << "reading from HDF5" << endl;
	io->readFromHDF5(base()+infn+".h5", rmat);
	cout << rmat->getASCIIdata();
		
	delete rmat;
	
	delete io;
	
	return 0;
}


//-------------------------------------------------------- -t5
int Test::testDataTypes(){
	cout << "--- TESTING DATA TYPES ON THIS SYSTEM ---" << endl;
	cout << "sizeof(int) =                " << sizeof(int) << endl;
	cout << "sizeof(short int) =          " << sizeof(short int) << endl;
	cout << "sizeof(long int) =           " << sizeof(long int) << endl;
	cout << "sizeof(unsigned int) =       " << sizeof(unsigned int) << endl;	
	cout << "sizeof(unsigned short int) = " << sizeof(unsigned short int) << endl;	
	cout << "sizeof(unsigned long int)  = " << sizeof(unsigned long int) << endl;		
	cout << "sizeof(float) =              " << sizeof(float) << endl;
	cout << "sizeof(double) =             " << sizeof(double) << endl;
	cout << "-----------------------------------------" << endl;
	return 0;
}
