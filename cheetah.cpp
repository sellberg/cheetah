/*
 *  cheetah.cpp
 *  cheetah
 *
 *  Created by Anton Barty on 7/2/11.
 *  Copyright 2011 CFEL. All rights reserved.
 *	
 *	This file based on myana_cspad.cc,v 1.6 2010/10/21 22:09:43 by Matt Weaver, SLAC 
 *
 *  The contents of the directories 
 *		myana,
 *		cspad-gjw, and 
 *		build 
 *	can be obtained from SLAC, are copyright by the SLAC instrument development team, 
 *	and are covered by whatever license that software is distributed under. 
 *
 *	The rest of this code is free software: you can modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "myana/myana.hh"
#include "myana/main.hh"
//#include "myana/XtcRun.hh"
#include "release/pdsdata/cspad/ConfigV1.hh"
#include "release/pdsdata/cspad/ConfigV2.hh"
#include "release/pdsdata/cspad/ConfigV3.hh"
#include "release/pdsdata/cspad/ConfigV4.hh"
#include "release/pdsdata/cspad/ElementHeader.hh"
#include "release/pdsdata/cspad/ElementIterator.hh"
//#include "cspad-gjw/CspadTemp.hh"
//#include "cspad-gjw/CspadCorrector.hh"
//#include "cspad-gjw/CspadGeometry.hh"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::string;

#include <string.h>
#include <sys/time.h>
#include <hdf5.h>
#include <pthread.h>
#include <stdlib.h>
#include <limits>
#include <cmath>

#include "setup.h"
#include "worker.h"
#include "attenuation.h"


static cGlobal		global;
static long			frameNumber;


// Quad class definition
/*
class MyQuad {
	public:
	MyQuad(unsigned q) : _quad(q) {
		char buff[64];
		for(unsigned i=0; i<16; i++) {
	  		sprintf(buff,"Q%d_ASIC%d_values",q,i);
	  		sprintf(buff,"Q%d_ASIC%d_map",q,i);
		}
	}
	
  	void Fill(Pds::CsPad::ElementIterator& iter) {
    	unsigned section_id;
    	const Pds::CsPad::Section* s;
    	while((s=iter.next(section_id))) {
      		for(unsigned col=0; col<COLS; col++)
				for(unsigned row=0; row<ROWS; row++) {	
				}
    	}
  	}    
  
  	void write() {
    	for(unsigned i=0; i<16; i++) {
    	}
  	}

	private:
  		unsigned _quad;
  		CspadSection _s;
};

static MyQuad* quads[4];
*/

using namespace Pds;


/*
 *	Beam on or off??
 */
static bool beamOn()
{
	int nfifo = getEvrDataNumber();
	for(int i=0; i<nfifo; i++) {
		unsigned eventCode, fiducial, timestamp;
		if (getEvrData(i,eventCode,fiducial,timestamp)) 
			printf("Failed to fetch evr fifo data\n");
		else if (eventCode==140)
			return true;
	}
	return false;;
}


// This function is called once at the beginning of the analysis job,
// You can ask for detector "configuration" information here.
void beginjob() {
	printf("beginjob()\n");
	
	/*
	 * Get time information
	 */
//	int fail = 0;
//	const char* time;
//  fail = getLocalTime( time );
    
	int seconds, nanoSeconds;
	getTime( seconds, nanoSeconds );	

	
	/*
	 *	New csPad corrector
	 */
	//corrector = new CspadCorrector(Pds::DetInfo::CxiDs1,0,CspadCorrector::DarkFrameOffset);

	/*
	for(unsigned i=0; i<4; i++)
		quads[i] = new MyQuad(i);
	*/
	
	/*
	 *	Stuff for worker thread management
	 */
	global.defaultConfiguration();
	global.parseConfigFile(global.configFile);
	global.readDetectorGeometry(global.geometryFile);
	global.setup();
	global.readDarkcal(global.darkcalFile);
	global.readBadpixelMask(global.badpixelFile);
	global.readGaincal(global.gaincalFile);
	global.readPeakmask(global.hitfinder.peaksearchFile);
	global.readIcemask(global.icefinder.peaksearchFile);
	global.readWatermask(global.waterfinder.peaksearchFile);
	global.readBackgroundmask(global.backgroundfinder.peaksearchFile);
	if (global.useAttenuationCorrection >= 0) global.readAttenuations(global.attenuationFile);
	if (global.usePixelStatistics) global.readPixels(global.pixelFile);
	if (global.useCorrelation) global.createLookupTable();	// <-- important that this is done after detector geometry is determined
}




void fetchConfig()
{
	if (getCspadConfig( Pds::DetInfo::CxiDs1, configV1 )==0) {
		configVsn= 1;
		quadMask = configV1.quadMask();
		asicMask = configV1.asicMask();
		printf("CSPAD configuration: quadMask %x  asicMask %x  runDelay %d\n", quadMask, asicMask, configV1.runDelay());
		printf("\tintTime %d/%d/%d/%d\n", configV1.quads()[0].intTime(), configV1.quads()[1].intTime(), configV1.quads()[2].intTime(), configV1.quads()[3].intTime());
	}
	else if (getCspadConfig( Pds::DetInfo::CxiDs1, configV2 )==0) {
		configVsn= 2;
		quadMask = configV2.quadMask();
		asicMask = configV2.asicMask();
		printf("CSPAD configuration: quadMask %x  asicMask %x  runDelay %d\n", quadMask, asicMask, configV2.runDelay());
		printf("\tintTime %d/%d/%d/%d\n", configV2.quads()[0].intTime(), configV2.quads()[1].intTime(), configV2.quads()[2].intTime(), configV2.quads()[3].intTime());
	}
	else if (getCspadConfig( Pds::DetInfo::CxiDs1, configV3 )==0) {
		configVsn= 3;
		quadMask = configV3.quadMask();
		asicMask = configV3.asicMask();
		printf("CSPAD configuration: quadMask %x  asicMask %x  runDelay %d\n", quadMask, asicMask, configV3.runDelay());
		printf("\tintTime %d/%d/%d/%d\n", configV3.quads()[0].intTime(), configV3.quads()[1].intTime(), configV3.quads()[2].intTime(), configV3.quads()[3].intTime());
	}
	else if (getCspadConfig( Pds::DetInfo::CxiDs1, configV4 )==0) {
		configVsn= 4;
		quadMask = configV4.quadMask();
		asicMask = configV4.asicMask();
		printf("CSPAD configuration: quadMask %x  asicMask %x  runDelay %d\n", quadMask, asicMask, configV4.runDelay());
		printf("\tintTime %d/%d/%d/%d\n", configV4.quads()[0].intTime(), configV4.quads()[1].intTime(), configV4.quads()[2].intTime(), configV4.quads()[3].intTime());
		}
	else {
		configVsn= 0;
		printf("Failed to get CspadConfig\n");
		exit(1);
	}
}


/*
 * 	This function is called once for each run.  You should check to see
 *  	if detector configuration information has changed.
 */
void beginrun() 
{
	//printf("beginrun\n");
	printf("Processing r%04u\n",getRunNumber());
	fetchConfig();
	//corrector->loadConstants(getRunNumber());
	global.runNumber = getRunNumber();
	global.writeInitialLog();	
	frameNumber = 0;

}

/*
 *	Calibration
 */
void begincalib()
{
	printf("begincalib\n");
	fetchConfig();
}



/*
 *	This function is called once per shot
 */
void event() {
	
	// Variables
	frameNumber++;
	//printf("Processing event %i\n", frameNumber);
//	FILE* fp;
//	char filename[1024];
	int fail = 0;

	
	
	
	/*
	 *	Get run number
	 */
	unsigned runNumber;
	runNumber = getRunNumber();
	
	
	/*
	 *	Get event information
	 */	 
	int 			numEvrData;
	unsigned		fiducial;
	unsigned int 	eventCode;
	unsigned int 	timeStamp;	
	numEvrData = getEvrDataNumber();
	for (int i=0; i<numEvrData; i++) {
		fail = getEvrData( i, eventCode, fiducial, timeStamp );
	}
	// EventCode==140 = Beam On
	//printf("fiducial: %x\n", fiducial);
	
	
	/*
	 *	How quickly are we processing the data? (average over last 10 events)
	 */	
	timeval	now;
	float dt, dt_us, datarate;
	gettimeofday(&now, NULL);
	dt_us = (float) (now.tv_usec - global.lasttime.tv_usec);
	dt = clock() - global.lastclock;
	if(dt_us != 0) {
		datarate = ((float)CLOCKS_PER_SEC)/dt;
		//datarate = 1/(1e6*dt_us);
		gettimeofday(&global.lasttime, NULL);
		global.lastclock = clock();
		global.datarate = (datarate+9*global.datarate)/10.;
	}
	
	
	

	/*
	 * Get time information
	 */
	int seconds, nanoSeconds;
	getTime( seconds, nanoSeconds );
	
    //	const char* timestring;
	//fail = getLocalTime( timestring );
	//printf("%s\n",timestring);
	
	/*
	 *	Get fiducials
	 */
	fail = getFiducials(fiducial);

	
	/* 
	 *	Is the beam on?
	 */
	bool beam = beamOn();
	//printf("Beam %s : fiducial %x\n", beam ? "On":"Off", fiducial);
	if(!beam && !global.generateDarkcal)
		return;
	
	
	/*
	 * Get electron beam parameters from beamline data
	 */     
	double fEbeamCharge;    // in nC
	double fEbeamL3Energy;  // in MeV
	double fEbeamLTUPosX;   // in mm
	double fEbeamLTUPosY;   // in mm
	double fEbeamLTUAngX;   // in mrad
	double fEbeamLTUAngY;   // in mrad
	double fEbeamPkCurrBC2; // in Amps
	double photonEnergyeV;
	double wavelengthA;
	
	if ( getEBeam(fEbeamCharge, fEbeamL3Energy, fEbeamLTUPosX, fEbeamLTUPosY,
	              fEbeamLTUAngX, fEbeamLTUAngY, fEbeamPkCurrBC2) ) {
		
		wavelengthA = std::numeric_limits<double>::quiet_NaN();
		photonEnergyeV = std::numeric_limits<double>::quiet_NaN();
		
	} else {
		
		/* Calculate the resonant photon energy (ie: photon wavelength) */
		// Get the present peak current in Amps
		double peakCurrent = fEbeamPkCurrBC2;
		// Get present beam energy [GeV]
		double DL2energyGeV = 0.001*fEbeamL3Energy;
		// wakeloss prior to undulators
		double LTUwakeLoss = 0.0016293*peakCurrent;
		// Spontaneous radiation loss per segment
		double SRlossPerSegment = 0.63*DL2energyGeV;
		// wakeloss in an undulator segment
		double wakeLossPerSegment = 0.0003*peakCurrent;
		// energy loss per segment
		double energyLossPerSegment = SRlossPerSegment + wakeLossPerSegment;
		// energy in first active undulator segment [GeV]
		double energyProfile = DL2energyGeV - 0.001*LTUwakeLoss
		- 0.0005*energyLossPerSegment;
		// Calculate the resonant photon energy of the first active segment
		photonEnergyeV = 44.42*energyProfile*energyProfile;
		// Calculate wavelength in Angstrom
		wavelengthA = 12398.42/photonEnergyeV;
	}
	
	
	/*
	 * 	FEE gas detectors (pulse energy in mJ)
	 */     
	double gasdet[4];
	double gmd1;
	double gmd2;
	if ( getFeeGasDet(gasdet) ) {
		gmd1 = std::numeric_limits<double>::quiet_NaN();
		gmd2 = std::numeric_limits<double>::quiet_NaN();
	} else {
		gmd1 = (gasdet[0]+gasdet[1])/2;
		gmd2 = (gasdet[2]+gasdet[3])/2;
	}
	
	
	/*
	 * Phase cavity data
	 *	(we probably won't need this info)
	 */     
	double 	phaseCavityTime1;
	double	phaseCavityTime2;
	double	phaseCavityCharge1;
	double	phaseCavityCharge2;
	fail = getPhaseCavity(phaseCavityTime1, phaseCavityTime2, phaseCavityCharge1, phaseCavityCharge2);
	
	
	/*
	 *	CXI detector position (Z)
	 */
	float detposnew;
	if ( getPvFloat(global.detectorZpvname, detposnew) == 0 ) {
		// sanity check of detector readout, corrects artificial 'jumps' in the detector position if detector change is larger than 0.83 mm between events, corresponding to a speed larger than 100 mm/s
		if (fabs(detposnew - global.detposold) > 102.5 && nevents != 0) { // 102.5 corresponds to 100 mm/s, which is far above the upper limit of the detector speed
			printf("Old detector position = %2.2f mm\n", global.detposold);
			printf("New detector position = %2.2f mm\n", detposnew);
			cout << "Detector change of " << fabs(detposnew - global.detposold)*120/123 << " mm/s" << endl;
			cout << "\tArtificial detector jump detected! Keeping the old detector position." << endl;
			detposnew = global.detposold;
		} else {
			global.detposold = detposnew;
		}
		/* When encoder reads -500mm, detector is at its closest possible
		 * position to the specimen, and is 79mm from the centre of the 
		 * 8" flange where the injector is mounted.  The injector itself is
		 * about 4mm further away from the detector than this. */
		global.detectorZ = 500.0 + detposnew + 79.0 + global.detectorOffset;
	} else if ((nevents-global.attenuationOffset) % 123 == 0) {
		cout << "Failed to retrieve detector position for EVENT #" << nevents+1 << endl;
		if (nevents == 0) {
			global.detectorZ = 500.0 + global.detectorZpos + 79.0 + global.detectorOffset;
			cerr << "WARNING: using pre-set detector position = " << global.detectorZpos << " mm" << endl;
		}
	}
	
	
	/*
	 *	Attenuation information is only saved for every 1/123 event
	 */
	if ((nevents-global.attenuationOffset) % 123 == 0 && global.useAttenuationCorrection >= 0) {
		
		/*
		 *	Get total thickness for Si filters in XRT
		 */
		unsigned totalThickness;
		fail = getSiThickness(totalThickness, global.nFilters, global.filterThicknesses, global.filterPositionIn, global.filterPositionOut);
		// cout << "Total thickness: " << totalThickness << endl;
		if (fail) {
			cout << "Failed to retrieve attenuation for EVENT #" << nevents+1 << " [failcode " << fail << "]" << endl;
			if (fail == 30) global.attenuationOffset++;
		} else {
		
			/*
			 *	Calculate attenuation from Si filters
			 */
			double attenuation;
			attenuation = getAttenuation(totalThickness, global.nThicknesses, global.possibleThicknesses, global.possibleAttenuations);
			// cout << "Calculated attenuation is: " << attenuation << endl;
		
			/*
			 *	Add attenuation to dynamic array
			 */
			if (global.nAttenuations >= global.attenuationCapacity) global.expandAttenuationCapacity();
			if (global.nAttenuations == 0 || global.attenuations[global.nAttenuations-1] != attenuation) {
				global.attenuations[global.nAttenuations] = attenuation;
				global.changedAttenuationEvents[global.nAttenuations] = nevents+1;
				global.totalThicknesses[global.nAttenuations++] = totalThickness;
			}
		}
	}
	
	
	/*
	 *	Energy calibration
	 */
	if (global.useEnergyCalibration && photonEnergyeV == photonEnergyeV) { // Check if photonEnergyeV is NAN
		if (global.nEnergies >= global.energyCapacity) global.expandEnergyCapacity();
		global.energies[global.nEnergies] = photonEnergyeV;
		global.wavelengths[global.nEnergies++] = wavelengthA;
		if (photonEnergyeV > global.Emax) global.Emax = photonEnergyeV;
		if (photonEnergyeV < global.Emin) global.Emin = photonEnergyeV;
		if (wavelengthA > global.Lmax) global.Lmax = wavelengthA;
		if (wavelengthA < global.Lmin) global.Lmin = wavelengthA;
	}
	
	
	/*
	 *	Create a new threadInfo structure in which to place all information
	 */
	tThreadInfo	*threadInfo;
	threadInfo = (tThreadInfo*) malloc(sizeof(tThreadInfo));
		
	
	/*
	 *	Copy all interesting information into worker thread structure if we got this far.
	 *	(ie: presume that myana itself is NOT thread safe and any event info may get overwritten)
	 */
	threadInfo->seconds = seconds;
	threadInfo->nanoSeconds = nanoSeconds;
	threadInfo->fiducial = fiducial;
	threadInfo->runNumber = getRunNumber();
	threadInfo->beamOn = beam;
	threadInfo->nPeaks = 0;
	
	threadInfo->gmd11 = gasdet[0];
	threadInfo->gmd12 = gasdet[1];
	threadInfo->gmd21 = gasdet[2];
	threadInfo->gmd22 = gasdet[3];
	
	threadInfo->fEbeamCharge = fEbeamCharge;		// in nC
	threadInfo->fEbeamL3Energy = fEbeamL3Energy;	// in MeV
	threadInfo->fEbeamLTUPosX = fEbeamLTUPosX;		// in mm
	threadInfo->fEbeamLTUPosY = fEbeamLTUPosY;		// in mm
	threadInfo->fEbeamLTUAngX = fEbeamLTUAngX;		// in mrad
	threadInfo->fEbeamLTUAngY = fEbeamLTUAngY;		// in mrad
	threadInfo->fEbeamPkCurrBC2 = fEbeamPkCurrBC2;	// in Amps
	threadInfo->photonEnergyeV = photonEnergyeV;	// in eV
	threadInfo->wavelengthA = wavelengthA;			// in Angstrom
	
	threadInfo->phaseCavityTime1 = phaseCavityTime1;
	threadInfo->phaseCavityTime2 = phaseCavityTime2;
	threadInfo->phaseCavityCharge1 = phaseCavityCharge1;
	threadInfo->phaseCavityCharge1 = phaseCavityCharge2;
	
	threadInfo->detectorPosition = global.detectorZ;	// in mm
	
	threadInfo->pixelCenterX = global.pixelCenterX;
	threadInfo->pixelCenterY = global.pixelCenterY;
	
	if (global.useAttenuationCorrection >= 2 && global.nAttenuations >= 1) {
		threadInfo->attenuation = global.attenuations[global.nAttenuations-1]/((gmd1+gmd2)/2);	// 1/(transmission*gmd) taken from last successful readouts of the Si filters and FEE gas detectors
	}
	else if (global.useAttenuationCorrection >= 0 && global.nAttenuations >= 1) {
			threadInfo->attenuation = global.attenuations[global.nAttenuations-1];	// 1/transmission taken from last successful readout of the Si filters
	} else {
		threadInfo->attenuation = 0;
	}
	
	threadInfo->pGlobal = &global;
	
	for(int quadrant=0; quadrant<4; quadrant++) {
		threadInfo->quad_data[quadrant] = (uint16_t*) calloc(ROWS*COLS*16, sizeof(uint16_t));
		memset(threadInfo->quad_data[quadrant], 0, ROWS*COLS*16*sizeof(uint16_t));
	}
	
	
	/*
	 *	Check if the list contains the event if the listfinder is enabled
	 */
	if (global.listfinder.use) {
		nameEvent(threadInfo, &global);
		global.eventIsHit = containsEvent((string) threadInfo->eventname, &global);
		if (global.eventIsHit) global.nhits++;
	}
	
	/*
	 *	Copy raw cspad image data into worker thread structure for processing
	 */
	Pds::CsPad::ElementIterator iter;
	fail=getCspadData(DetInfo::CxiDs1, iter);


	if (fail || (global.listfinder.use && !global.eventIsHit)) { 
		if (fail) {
			printf("getCspadData fail=%d, skipping this event (%x).\n",fail,fiducial);
			threadInfo->cspad_fail = fail;
		}
		
		// if the event is not a hit, update counters, save log files, and flush memory
		if (!fail) {
			nevents++;
			global.nprocessedframes++;
			
			/*
			 *	Save periodic powder patterns and update log file
			 */
			if(global.saveInterval!=0 && (global.nprocessedframes%global.saveInterval)==0 && (global.nprocessedframes > global.startFrames+50) ){
				saveRunningSums(&global);
				global.updateLogfile();
			}
			
			/*
			 *	Periodically flush all HDF5 data and memory
			 */
			if(global.flushInterval != 0 && (global.nprocessedframes%global.flushInterval) == 0) {
				cout << "Preparing for HDF5 memory flush..." << endl;
				// Wait for threads to finish before flushing
				while(global.nActiveThreads > 0) {
					printf("\twaiting for %i worker threads to terminate\n", (int)global.nActiveThreads);
					usleep(100000);
				}
				flushHDF5();
			}
		}
		
		// cleanup allocated arrays
		for(int quadrant=0; quadrant<4; quadrant++) 
			free(threadInfo->quad_data[quadrant]);
		free(threadInfo);
		
		return;
	} else {
		nevents++;
		const Pds::CsPad::ElementHeader* element;

		// loop over elements (quadrants)
		while(( element=iter.next() )) {  
			if(element->quad() < 4) {
				// Which quadrant is this?
				int quadrant = element->quad();
				
				// Have we jumped to a new fiducial (event??)
				//if (fiducial != element->fiducials())
				//	printf("Fiducial jump: %x/%d:%x\n",fiducial,element->quad(),element->fiducials());
				
				// Get temperature on strong back 
				//float	temperature = CspadTemp::instance().getTemp(element->sb_temp(2));
				//float	temperature = CspadTemp::instance().getTemp(element->sb_temp((element->quad()%2==0)?3:0));
				//printf("Temperature on quadrant %i: %3.1fC\n",quadrant, temperature);
				//printf("Temperature: %3.1fC\n",CspadTemp::instance().getTemp(element->sb_temp((element->quad()%2==0)?3:0)));
				//threadInfo->quad_temperature[quadrant] = temperature;
				
				
				// Read 2x1 "sections" into data array in DAQ format, i.e., 2x8 array of asics (two bytes / pixel)
				const Pds::CsPad::Section* s;
				unsigned section_id;
				while(( s=iter.next(section_id) )) {  
					//printf("\tQuadrant %d, Section %d  { %04x %04x %04x %04x }\n", quadrant, section_id, s->pixel[0][0], s->pixel[0][1], s->pixel[0][2], s->pixel[0][3]);
					memcpy(&threadInfo->quad_data[quadrant][section_id*2*ROWS*COLS],s->pixel[0],2*ROWS*COLS*sizeof(uint16_t));
				}
			}
		}
	}//if (fail)
	
	
	
	/*
	 *	Spawn worker thread to process this frame
	 *	Threads are created detached so we don't have to wait for anything to happen before returning
	 *		(each thread is responsible for cleaning up its own threadInfo structure when done)
	 */
	pthread_t		thread;
	pthread_attr_t	threadAttribute;
	int				returnStatus;
	

	// Avoid fork-bombing the system: wait until we have a spare thread in the thread pool
	while(global.nActiveThreads >= global.nThreads) {
		usleep(100);
	}
	
	
	// Increment threadpool counter
	pthread_mutex_lock(&global.nActiveThreads_mutex);
	global.nActiveThreads += 1;
	threadInfo->threadNum = ++global.threadCounter;
	pthread_mutex_unlock(&global.nActiveThreads_mutex);
	
	// Set detached state
	pthread_attr_init(&threadAttribute);
	//pthread_attr_setdetachstate(&threadAttribute, PTHREAD_CREATE_JOINABLE);
	pthread_attr_setdetachstate(&threadAttribute, PTHREAD_CREATE_DETACHED);
	
	// Create a new worker thread for this data frame
	returnStatus = pthread_create(&thread, &threadAttribute, worker, (void *)threadInfo); 
	pthread_attr_destroy(&threadAttribute);
	//pthread_detach(thread);
	global.nprocessedframes += 1;
	
	
	/*
	 *	Save periodic powder patterns and update log file
	 */
	if(global.saveInterval!=0 && (global.nprocessedframes%global.saveInterval)==0 && (global.nprocessedframes > global.startFrames+50) ){
		saveRunningSums(&global);
		global.updateLogfile();
	}
	
	
	/*
	 *	Periodically flush all HDF5 data and memory
	 */
	if(global.flushInterval != 0 && (global.nprocessedframes%global.flushInterval) == 0) {
		cout << "Preparing for HDF5 memory flush..." << endl;
		// Wait for threads to finish before flushing
		while(global.nActiveThreads > 0) {
			printf("\twaiting for %i worker threads to terminate\n", (int)global.nActiveThreads);
			usleep(100000);
		}
		flushHDF5();
	}
	
}
// End of event data processing block




/*
 *	Stuff that happens at the end
 */
void endcalib() {
	printf("endcalib()\n");
}

void endrun() 
{
	printf("User analysis endrun() routine called.\n");
}

void endjob()
{
	printf("User analysis endjob() routine called.\n");


	// Wait for threads to finish
	while(global.nActiveThreads > 0) {
		printf("Waiting for %i worker threads to terminate\n", (int)global.nActiveThreads);
		usleep(100000);
	}
	
	
	// Calculate center correction from regular powder pattern
	if (global.hitfinder.use && global.powdersum && global.calculateCenterCorrectionPowder) {
		calculateCenterCorrection(&global, global.powderRaw, global.npowder);
		shiftPixelCenter(&global);
		updatePixelArrays(&global);
		updateImageArrays(&global);
		updateAngularAvgArrays(&global);
	}
	
	
	// Calculate center correction from regular powder pattern for each quad
	if (global.hitfinder.use && global.powdersum && global.calculateCenterCorrectionQuad) {
		calculateCenterCorrectionQuad(&global, global.powderRaw, global.npowder); // function also shifts Quads
		updatePixelArrays(&global);
		updateImageArrays(&global);
		updateAngularAvgArrays(&global);
	}
	
	
	// Refine metrology by translating/rotating quads w.r.t each other from regular powder pattern
	if (global.hitfinder.use && global.powdersum && global.refineMetrology) {
		cout << "Refining metrology..." << endl;
		translateQuads(&global);
		rotateQuads(&global);
		optimizeGap(&global);
		updatePixelArrays(&global);
		updateImageArrays(&global);
		updateAngularAvgArrays(&global);
	}
	
	
	// Save intensity statistics
	if (global.useIntensityStatistics) {
		saveIntensities(&global);
		makeIntensityHistograms(&global);
	}
	
	
	// Save energy calibration
	if (global.useEnergyCalibration) {
		saveEnergies(&global);
		makeEnergyHistograms(&global);
		if (global.hitAngularAvg || (global.powdersum && global.powderAngularAvg)) makePowderQcalibration(&global);
	}
	
	
	// Calculate angular average of powder pattern
	if (global.powdersum && global.powderAngularAvg) {
		calculatePowderAngularAvg(&global);
		savePowderAngularAvg(&global);
	}
	
	
	// Save powder patterns
	saveRunningSums(&global);
	
	
	global.writeFinalLog();
	
	
	// Attenuation?
	if (global.useAttenuationCorrection >= 0) {
		printf("%i attenuations calculated:\n", global.nAttenuations);
		for (int i=0; i<global.nAttenuations; i++) {
			cout << "\t" << global.attenuations[i] << " [first recorded for EVENT #" << global.changedAttenuationEvents[i] << "] corresponding to " << global.totalThicknesses[i] << " um Si" << endl;
		}
	}
	
	
	// Hitrate?
	printf("%i files processed, %i hits (%2.2f%%)\n",(int)global.nprocessedframes, (int)global.nhits, 100.*( global.nhits / (float) global.nprocessedframes));

	
	// Cleanup
	printf("Cheetah ");
	
	free(global.darkcal);
	free(global.powderAssembled);
	free(global.powderRaw);
	free(global.powderAverage);
	free(global.powderCorrelation);
	free(global.powderVariance);
	free(global.iceAssembled);
	free(global.iceRaw);
	free(global.iceAverage);
	free(global.iceCorrelation);
	free(global.waterAssembled);
	free(global.waterRaw);
	free(global.waterAverage);
	free(global.waterCorrelation);
	free(global.pix_x);
	free(global.pix_y);
	free(global.pix_z);
	free(global.pix_r);
	free(global.quad_dx);
	free(global.quad_dy);
	free(global.hotpixelmask);
	free(global.selfdark);
	free(global.gaincal);
	free(global.badpixelmask);
	free(global.hitfinder.peakmask);
	free(global.icefinder.peakmask);
	free(global.waterfinder.peakmask);
	free(global.backgroundfinder.peakmask);
	
	delete[] global.filterThicknesses;
	delete[] global.possibleThicknesses;
	delete[] global.possibleAttenuations;
	delete[] global.attenuations;
	delete[] global.changedAttenuationEvents;
	delete[] global.totalThicknesses;
	delete[] global.energies;
	delete[] global.intensities;
	delete[] global.hits;
	delete[] global.pixels;
	delete[] global.pixelXYList;	
	delete[] global.wavelengths;
	delete[] global.phi;
	delete[] global.angularAvg_i;
	delete[] global.angularAvgQ;
	delete[] global.angularAvgQcal;
	delete[] global.correlationLUT;
	
	pthread_mutex_destroy(&global.nActiveThreads_mutex);
	pthread_mutex_destroy(&global.powdersumraw_mutex);
	pthread_mutex_destroy(&global.powdersumassembled_mutex);
	pthread_mutex_destroy(&global.powdersumcorrelation_mutex);
	pthread_mutex_destroy(&global.powdersumvariance_mutex);
	pthread_mutex_destroy(&global.icesumraw_mutex);
	pthread_mutex_destroy(&global.icesumassembled_mutex);
	pthread_mutex_destroy(&global.icesumcorrelation_mutex);
	pthread_mutex_destroy(&global.watersumraw_mutex);
	pthread_mutex_destroy(&global.watersumassembled_mutex);
	pthread_mutex_destroy(&global.watersumcorrelation_mutex);
	pthread_mutex_destroy(&global.correlation_mutex);
	pthread_mutex_destroy(&global.correlationFFT_mutex);
	pthread_mutex_destroy(&global.pixelcenter_mutex);
	pthread_mutex_destroy(&global.image_mutex);	
	pthread_mutex_destroy(&global.intensities_mutex);
	pthread_mutex_destroy(&global.selfdark_mutex);
	pthread_mutex_destroy(&global.hotpixel_mutex);
	pthread_mutex_destroy(&global.nhits_mutex);
	pthread_mutex_destroy(&global.framefp_mutex);
	
	printf("done!\n");
}
