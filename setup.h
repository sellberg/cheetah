/*
 *  setup.cpp
 *  cspad_cryst
 *
 *  Created by Anton Barty on 7/2/11.
 *  Copyright 2011 CFEL. All rights reserved.
 *
 */



/*
 *	Global variables
 */
class cGlobal {
	
public:
	// Switches for processing options
	int			startFrames;
	int			cmModule;
	int			cmColumn;
	int			subtractBg;
	int			subtractDarkcal;
	int			selfDarkcal;
	int			hitfinder;
	int			savehits;
	int			powdersum;
	int			saveRaw;
	int			debugLevel;
	int			hdf5dump;
	int			autohotpixel;
	
	
	// Power user settings
	float		cmFloor;
	int			saveInterval;
	int			powderthresh;
	int			hitfinderADC;
	int			hitfinderNAT;
	int			scaleDarkcal;
	float		hotpixFreq;
	float		hotpixADC;
	float		hotpixMemory;
	float		selfDarkMemory;
	
	
	// Configuration files
	char		configFile[1024];
	char		geometryFile[1024];
	char		darkcalFile[1024];
	
	// Run information
	unsigned	runNumber;
	
	
	// Thread management
	long			nThreads;
	long			nActiveThreads;
	long			threadCounter;
	pthread_t		*threadID;
	pthread_mutex_t	nActiveThreads_mutex;
	pthread_mutex_t	hotpixel_mutex;
	pthread_mutex_t	selfdark_mutex;
	pthread_mutex_t	powdersum1_mutex;
	pthread_mutex_t	powdersum2_mutex;
	
	
	// Detector geometry
	long			pix_nx;
	long			pix_ny;
	long			pix_nn;
	float			*pix_x;
	float			*pix_y;
	float			*pix_z;
	float			pix_dx;
	unsigned		module_rows;
	unsigned		module_cols;
	long			image_nx;
	long			image_nn;
	
	
	// Common variables
	int16_t			*darkcal;
	int32_t			*powderRaw;
	int32_t			*powderAssembled;
	float			*hotpixelmask;
	float			*selfdark;
	long			npowder;
	long			nprocessedframes;
	clock_t			lastclock;
	float			datarate;

	
	
	
	
public:
	void defaultConfiguration(void);
	void parseConfigFile(char *);
	void parseCommandLineArguments(int, char**);
	void setupThreads(void);
	void readDetectorGeometry(char *);
	void readDarkcal(char *);

	
private:
	void parseConfigTag(char*, char*);

	
};

