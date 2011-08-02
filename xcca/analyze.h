#ifndef ANALYZE_H
#define ANALYZE_H
//
//  analyze.h
//  xcca_commandline
//
//  Created by Feldkamp on 7/21/11.
//  Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
//

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "arrayclasses.h"

class Analyzer {
public:	
	Analyzer();
	~Analyzer();

	int processFiles( vector<string> files, int cenX, int cenY, int num_phi, int num_q);
	int processFiles( vector<string> files, int cenX, int cenY, int num_phi, int num_q, 
						double start_q, double stop_q, int LUTx, int LUTy );
						
	void setBackground( array2D *back );
	array2D *background();
	void setBackgroundWeight( double weight );
	double backgroundWeight();
	void setMask( array2D *newmask );
	array2D *mask();
	void setOutputDirectory( string outdir );
	string outputDirectory();
	
	bool flag_subtract_background;				// subtract background?
	bool flag_single_correlation_output;		// write every single correlation?
	bool flag_use_mask;
	
private:
	array2D *p_back;
	array2D *p_mask;
	double p_back_weight;
	string p_out_dir;
};

class Test {
	public:
		Test();
		Test( string base );
		~Test();
		
		void setBase( string base );
		string base();

		int testArrayClasses();
		int testCrossCorrelator( int alg = 1 );
		int testFourierTrafo();
		int testIO();
		int testDataTypes();
		
	private:
		string p_base;
};

#endif
