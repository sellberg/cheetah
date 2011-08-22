
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include <string>
using std::string;

#include <sstream>		//strings
#include <fstream>		//files
#include <vector>

#include "crosscorrelator.h"
#include "arrayclasses.h"
#include "arraydataIO.h"

#include "analyze.h"

void usage(){
	cout << "====================usage============================" << endl;
	cout << " -l<filelist>   : input from a file list " << endl;
	cout << " -o<directory>  : specify output directory" << endl;
	cout << " -s             : output of separate correlation for each individual image" << endl;
	cout << " -b<filename>   : subtract background file " << endl;
	cout << " -m<filename>   : use mask file " << endl;
	cout << " -B<factor>     : weighting of background" << endl;
	cout << "-----------------------------------------------------" << endl;
	cout << " -t<num>        : test different features" << endl;
	cout << "                      num=1 : testCrossCorrelator" << endl;	
	cout << "                      num=2 : testArrayClasses" << endl;	
	cout << "                      num=3 : testFourierTrafo" << endl;	
	cout << "                      num=4 : testIO" << endl;	
	cout << "                      num=5 : testDataTypes" << endl;	
	cout << "                      num=0 : test all" << endl;	
	cout << "=====================================================" << endl;
}

string argToString( char * const argument ){
	//fill argument after the -<option letter> into a string
	std::ostringstream osst;
	int j = 2;
	while ( argument[j] != '\0' ){
		osst << argument[j];
		j++;
	}
	string retstring = osst.str();
	return retstring;
}

int main (int argc, char * const argv[]) {

    cout << "Hello, World!\n";

//	cout << "argc = " << argc << endl;
//	for (int i = 0; i < argc; i++){
//		cout << "argv[" << i << "] = '" << argv[i] << "'" << endl;
//	}

	Analyzer *ana = new Analyzer;
	std::vector<string> files;								// vector to be filled with individual files to process
				
	if (argc < 2){
		usage();
		return 1;
	} else {
		for(int i = 1; i < argc; i++){						// check if all options are valid first
			if (argv[i][0] != '-')
			{
				cout << "'" << argv[i] << "' does not start with a dash: not a valid option. Exiting." << endl;
				usage();
				return 1;
			}
		}

		for(int i = 1; i < argc; i++){						//if all options are valid, proceed and evaluate
			
				if (argv[i][1] == 'l'){
				string list_fn = argToString(argv[i]);
				
				//read from file list
				std::ifstream fin;
				fin.open(list_fn.c_str());					
				if( fin.is_open() ){
					string line;
					while (fin.good()) {
						getline(fin, line);
						if (line != ""){
							files.push_back(line);
						}
					}
				}else{
					cerr << "Could not open file '" << list_fn << "'" << endl;
					exit( 1 );
				}
				fin.close();
				

				cout << "--> using file list in " << list_fn << endl;

			}else if (argv[i][1] == 'o'){
				string outdir = argToString(argv[i]);
				ana->setOutputDirectory( outdir );
				cout << "--> using output directory " << outdir << endl;	
				
			}else if (argv[i][1] == 's'){
				ana->flag_single_correlation_output = true;
				cout << "--> using single image output " << endl;
				
			} else if (argv[i][1] == 'b') {				
				string back_fn = argToString(argv[i]);
				
				array2D *back = new array2D;
				arraydataIO *io = new arraydataIO;
				io->readFromEDF( back_fn, back);
				ana->setBackground( back );
				ana->flag_subtract_background = true;				
				delete io;
				delete back;
				cout << "--> using background file " << back_fn << endl;
				
			} else if (argv[i][1] == 'B') {
				double weight = 1;
				std::stringstream sst;
				int j = 2;
				while ( argv[i][j] != '\0' ){				
					sst << argv[i][j];
					j++;
				}
				sst >> weight;
				ana->setBackgroundWeight( weight );
				cout << "--> using background weighting factor " << weight << endl;	
			
			}else if (argv[i][1] == 'm') {				
				string mask_fn = argToString(argv[i]);
				
				array2D *mask = new array2D;
				arraydataIO *io = new arraydataIO;
				io->readFromEDF( mask_fn, mask);
				ana->setMask( mask );			
				delete io;
				delete mask;
				cout << "--> using mask file " << mask_fn << endl;
				
			} else if (argv[i][1] == 't') {
				cout << "---- testing option '" << argv[i][2] << "' ----" << endl;
				//---------------------------------------------------
				//run various tests
				//---------------------------------------------------
				string base = "/Users/feldkamp/Desktop/test/";
				cout << "output directory '" << base << "'" << endl;
				Test *t = new Test();
				t->setBase(base);
				
				switch(argv[i][2]){
					case '1':
						t->testCrossCorrelator( 1 );					
						break;
					case '2':
						t->testArrayClasses();
						break;
					case '3':
						t->testFourierTrafo();
						break;
					case '4':
						t->testIO();
						break;
					case '5':
						t->testDataTypes();					
						break;
					case '0':							// fall through to default
					default:
						t->testCrossCorrelator( 1 );
						t->testArrayClasses();
						t->testFourierTrafo();
						t->testIO();
						t->testDataTypes();
						break;
				}//end switch
				cout << "---- testing done ----" << endl;
				delete t;
				return 0;
			} else {
				cout << "-" << argv[i][1] << " is not a valid option." << endl;
				usage();
				return 2;
			}
		}//end for i
	}//end if
	
	
	//define set of variables to pass to the processFiles function, should probabaly go into an ini file or so at some point
	double shiftX = -7;
	double shiftY = 1;
	int num_phi = 2048;
	int num_q = 250;
	double start_q = 0;
	double stop_q = num_q;
	int LUTx = 487;				//pilatus detector x and y values
	int LUTy = 619;
	
	//prototype: processFiles( files, shiftX, shiftY, num_phi, num_q);
	//prototype: processFiles( files, shiftX, shiftY, num_phi, num_q, start_phi, stop_phi, start_q, stop_q, LUTx, LUTy);
	ana->processFiles( files, shiftX, shiftY, num_phi, num_q, start_q, stop_q, LUTx, LUTy);
//	ana->processFiles( files, -7, 1, 2048, 250, 0, 360, 0, 250, 400, 400);
//	ana->processFiles( files, -7, 1, 2048, 250 );

    delete ana;
	
    cout << "Goodbye, World" << endl;
    return 0;
}



