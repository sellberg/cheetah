/*
 *  worker.h
 *  cheetah
 *
 *  Created by Anton Barty on 6/2/11.
 *  Copyright 2011 CFEL. All rights reserved.
 *
 *	You can modify this software under the terms of the GNU General Public License 
 *	as published by the Free Software Foundation, either version 3 of the License, or
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






/*
 *	Structure used for passing information to worker threads
 */
typedef struct {
	
	// Reference to common global structure
	cGlobal		*pGlobal;
	int			busy;
	long		threadNum;
	
	// cspad data
	int			cspad_fail;
	float		quad_temperature[4];
	uint16_t	*quad_data[4];
	uint16_t	*raw_data;
	int16_t		*corrected_data;
	int16_t		*image;
	int			nPeaks;
	int			nHot;
	
	
	// Beamline data, etc
	int			seconds;
	int			nanoSeconds;
	unsigned	fiducial;
	char		timeString[1024];
	char		eventname[1024];
	bool		beamOn;
	unsigned	runNumber;
	
	double		gmd11;
	double		gmd12;
	double		gmd21;
	double		gmd22;
	
	double		fEbeamCharge;		// in nC
	double		fEbeamL3Energy;		// in MeV
	double		fEbeamLTUPosX;		// in mm
	double		fEbeamLTUPosY;		// in mm
	double		fEbeamLTUAngX;		// in mrad
	double		fEbeamLTUAngY;		// in mrad
	double		fEbeamPkCurrBC2;	// in Amps
	
	double		photonEnergyeV;		// in eV
	double		wavelengthA;		// in Angstrom
	
	double		detectorPosition; 
	
	double		phaseCavityTime1;
	double		phaseCavityTime2;
	double		phaseCavityCharge1;
	double		phaseCavityCharge2;
	
	// Thread management
	int		threadID;
	
	
	
} tThreadInfo;



/*
 *	Stuff from Garth's original code
 */

// Static variables
using namespace std;
static CspadCorrector*      corrector;
static Pds::CsPad::ConfigV1 configV1;
static Pds::CsPad::ConfigV2 configV2;
static unsigned             configVsn;
static unsigned             quadMask;
static unsigned             asicMask;

static const unsigned  ROWS = 194;
static const unsigned  COLS = 185;

static const unsigned int cbufsize = 1024;


static uint32_t nevents = 0;

#define ERROR(...) fprintf(stderr, __VA_ARGS__)
#define STATUS(...) fprintf(stderr, __VA_ARGS__)

#define DEBUGL1_ONLY if(global->debugLevel >= 1)
#define DEBUGL2_ONLY if(global->debugLevel >= 2)

/*
 *	Function prototypes
 */
void *worker(void *);
void cmModuleSubtract(tThreadInfo*, cGlobal*);
void cmSubModuleSubtract(tThreadInfo*, cGlobal*);
void subtractDarkcal(tThreadInfo*, cGlobal*);
void applyGainCorrection(tThreadInfo*, cGlobal*);
void applyBadPixelMask(tThreadInfo*, cGlobal*);
void subtractPersistentBackground(tThreadInfo*, cGlobal*);
void killHotpixels(tThreadInfo*, cGlobal*);
int  hitfinder(tThreadInfo*, cGlobal*);
void addToPowder(tThreadInfo*, cGlobal*);
void assemble2Dimage(tThreadInfo*, cGlobal*);
void nameEvent(tThreadInfo*, cGlobal*);
void writeHDF5(tThreadInfo*, cGlobal*);
void writeSimpleHDF5(const char*, const void*, int, int, int);
void saveRunningSums(cGlobal*);




