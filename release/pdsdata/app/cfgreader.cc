
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>

#include "pdsdata/xtc/DetInfo.hh"
#include "pdsdata/xtc/BldInfo.hh"
#include "pdsdata/xtc/ProcInfo.hh"
#include "pdsdata/xtc/XtcIterator.hh"
#include "pdsdata/xtc/XtcFileIterator.hh"
#include "pdsdata/acqiris/ConfigV1.hh"
#include "pdsdata/ipimb/ConfigV1.hh"
#include "pdsdata/ipimb/ConfigV2.hh"
#include "pdsdata/encoder/ConfigV1.hh"
#include "pdsdata/camera/FrameFexConfigV1.hh"
#include "pdsdata/camera/FrameFccdConfigV1.hh"
#include "pdsdata/fccd/FccdConfigV1.hh"
#include "pdsdata/opal1k/ConfigV1.hh"
#include "pdsdata/pulnix/TM6740ConfigV1.hh"
#include "pdsdata/pnCCD/ConfigV1.hh"
#include "pdsdata/pnCCD/ConfigV2.hh"
#include "pdsdata/evr/IOConfigV1.hh"
#include "pdsdata/evr/ConfigV1.hh"
#include "pdsdata/evr/ConfigV2.hh"
#include "pdsdata/evr/ConfigV3.hh"
#include "pdsdata/evr/ConfigV4.hh"
#include "pdsdata/evr/ConfigV7.hh"
#include "pdsdata/control/ConfigV1.hh"
#include "pdsdata/control/PVControl.hh"
#include "pdsdata/control/PVMonitor.hh"
#include "pdsdata/epics/EpicsPvData.hh"
#include "pdsdata/epics/EpicsXtcSettings.hh"
#include "pdsdata/bld/bldData.hh"
#include "pdsdata/princeton/ConfigV1.hh"
#include "pdsdata/cspad/ConfigV4.hh"

using namespace Pds;

static bool noEpics;

class myLevelIter : public XtcIterator {
public:
  enum {Stop, Continue};
  myLevelIter(Xtc* xtc, unsigned depth) : 
    XtcIterator(xtc), _depth(depth) {}

  void process(const Src&, const Acqiris::ConfigV1&) {
    printf("*** Processing Acqiris config object\n");
  }
  void process(const Src&, const Ipimb::ConfigV1& o) {
    printf("*** Processing Ipimb config object\n");
    o.dump();
  }
  void process(const Src&, const Ipimb::ConfigV2& o) {
    printf("*** Processing Ipimb config object\n");
    o.dump();
  }
  void process(const Src&, const Encoder::ConfigV1&) {
    printf("*** Processing Encoder config object\n");
  }
  void process(const Src&, const Opal1k::ConfigV1&) {
    printf("*** Processing Opal1000 config object\n");
  }
  void process(const Src&, const Pulnix::TM6740ConfigV1&) {
    printf("*** Processing TM6740 config object\n");
  }
  void process(const Src&, const Camera::FrameFexConfigV1& c) {
    printf("*** Processing frame feature extraction config object\n");
    printf("roiBegin (%d,%d)  roiEnd(%d,%d)\n",
     c.roiBegin().column, c.roiBegin().row,
     c.roiEnd().column, c.roiEnd().row);
  }
  void process(const Src&, const Camera::FrameFccdConfigV1&) {
    printf("*** Processing FCCD Frame ConfigV1 object\n");
  }
  void process(const Src&, const FCCD::FccdConfigV1&) {
    printf("*** Processing FCCD ConfigV1 object\n");
  }
  void process(const Src& info, const PNCCD::ConfigV1& config) {
    const DetInfo& det = static_cast<const DetInfo&>(info);
    if ( det.detId() != 0 )
    {
      printf( "myLevelIter::process(...,PNCCD::ConfigV1&): pnCCD detector Id (%d) is not 0\n", det.detId() );
      return;
    }
    if ( det.devId() < 0 || det.devId() > 1)
    {
      printf( "myLevelIter::process(...,PNCCD::ConfigV1&): pnCCD device Id (%d) is out of range (0..1)\n", det.devId() );
      return;
    }
    
    _pnCcdCfgListV1[det.devId()] = config;
    printf("*** Processing pnCCD config.  Number of Links: %d, PayloadSize per Link: %d\n",
           config.numLinks(),config.payloadSizePerLink());
  }  
  void process(const Src& info, const PNCCD::ConfigV2& config) {
    const DetInfo& det = static_cast<const DetInfo&>(info);
    if ( det.detId() != 0 )
    {
      printf( "myLevelIter::process(...,PNCCD::ConfigV2&): pnCCD detector Id (%d) is not 0\n", det.detId() );
      return;
    }
    if ( det.devId() < 0 || det.devId() > 1)
    {
      printf( "myLevelIter::process(...,PNCCD::ConfigV2&): pnCCD device Id (%d) is out of range (0..1)\n", det.devId() );
      return;
    }

    _pnCcdCfgListV2[det.devId()] = config;
    printf("*** Processing pnCCD config.  Number of Links: %u, PayloadSize per Link: %u\n",
           config.numLinks(),config.payloadSizePerLink());
    printf("\tNumber of Channels %u, Number of Rows %u, Number of SubModule Channels %u\n\tNumber of SubModule Rows %u, Number of SubModules, %u\n",
        config.numChannels(),config.numRows(), config.numSubmoduleChannels(),config.numSubmoduleRows(),config.numSubmodules());
    printf("\tCamex Magic 0x%x, info %s, Timing File Name %s\n", config.camexMagic(),config.info(),config.timingFName());
  }
  void process(const Src&, const ControlData::ConfigV1& config) {
    printf("*** Processing Control config object\n");    
    
    printf( "Control PV Number = %d, Monitor PV Number = %d\n", config.npvControls(), config.npvMonitors() );
    for(unsigned int iPvControl=0; iPvControl < config.npvControls(); iPvControl++) {      
      const Pds::ControlData::PVControl& pvControlCur = config.pvControl(iPvControl);
      if (pvControlCur.array())
        printf( "%s[%d] = ", pvControlCur.name(), pvControlCur.index() );
      else
        printf( "%s = ", pvControlCur.name() );
      printf( "%lf\n", pvControlCur.value() );
    }
    
    for(unsigned int iPvMonitor=0; iPvMonitor < config.npvMonitors(); iPvMonitor++) {      
      const Pds::ControlData::PVMonitor& pvMonitorCur = config.pvMonitor(iPvMonitor);
      if (pvMonitorCur.array())
        printf( "%s[%d]  ", pvMonitorCur.name(), pvMonitorCur.index() );
      else
        printf( "%s  ", pvMonitorCur.name() );
      printf( "Low %lf  High %lf\n", pvMonitorCur.loValue(), pvMonitorCur.hiValue() );
    }
          
  }  
  void process(const Src&, const EpicsPvHeader& epicsPv)
  {    
    printf("*** Processing Epics object\n");
    epicsPv.printPv();
    printf( "\n" );
  }
  void process(const Src&, const BldDataFEEGasDetEnergy& bldData) {
    printf("*** Processing FEEGasDetEnergy object\n");
    bldData.print();
    printf( "\n" );    
  }  
  void process(const Src&, const BldDataEBeamV0& bldData) {
    printf("*** Processing EBeamV0 object\n");
    bldData.print();
    printf( "\n" );    
  }  
  void process(const Src&, const BldDataEBeamV1& bldData) {
    printf("*** Processing EBeamV1 object\n");
    bldData.print();
    printf( "\n" );    
  }  
  void process(const Src&, const BldDataEBeam& bldData) {
    printf("*** Processing EBeam object\n");
    bldData.print();
    printf( "\n" );    
  }  
  void process(const Src&, const BldDataPhaseCavity& bldData) {
    printf("*** Processing PhaseCavity object\n");
    bldData.print();
    printf( "\n" );    
  }
  void process(const Src&, const BldDataIpimbV0& bldData) {
    printf("*** Processing Bld-Ipimb V0 object\n");
    bldData.print();
    printf( "\n" );    
  } 

  void process(const Src&, const BldDataIpimb& bldData) {
    printf("*** Processing Bld-Ipimb V1 object\n");
    bldData.print();
    printf( "\n" );    
  } 

  void process(const Src&, const BldDataGMDV0& bldData) {
    printf("*** Processing Bld-GMD V0 object\n");
    bldData.print();
    printf( "\n" );    
  }   

  void process(const Src&, const BldDataGMDV1& bldData) {
    printf("*** Processing Bld-GMD V1 object\n");
    bldData.print();
    printf( "\n" );
  }
  
  void process(const Src&, const EvrData::IOConfigV1&) {
    printf("*** Processing EVR IOconfig V1 object\n");
  }
  void process(const Src&, const EvrData::ConfigV1&) {
    printf("*** Processing EVR config V1 object\n");
  }
  void process(const Src&, const EvrData::ConfigV2&) {
    printf("*** Processing EVR config V2 object\n");
  }
  void process(const Src&, const EvrData::ConfigV3&) {
    printf("*** Processing EVR config V3 object\n");
  }
  void process(const Src&, const EvrData::ConfigV4&) {
    printf("*** Processing EVR config V4 object\n");
  }
  void process(const Src&, const EvrData::ConfigV7& c) {
    printf("*** Processing EVR config V4 object\n");
    c.print();
    const EvrData::ConfigV7::SeqConfigType& s = c.seq_config();
    printf(" seq src %d/%d : len %d : cycles %d\n",
           s.sync_source(), s.beam_source(), s.length(), s.cycles());
  }
  void process(const Src&, const Princeton::ConfigV1&) {
    printf("*** Processing Princeton ConfigV1 object\n");
  }
  void process(const Src&, const CsPad::ConfigV4& c) {
    printf("*** Processing Cspad ConfigV4 object\n");
    printf("  runDelay %x  intTime %x\n",
           c.runDelay(), c.quads()[0].intTime());
  }
  int process(Xtc* xtc) {
    unsigned      i         =_depth; while (i--) printf("  ");
    Level::Type   level     = xtc->src.level();
    printf("%s level, payload size %d contains: %s: ",
     Level::name(level), xtc->sizeofPayload(), TypeId::name(xtc->contains.id()));
     
    const Src& info = xtc->src;
    switch(xtc->src.level()) {
    case Level::Source:
      printf("%s\n", DetInfo::name(static_cast<const DetInfo&>(xtc->src)));
      break;
    case Level::Reporter:
      printf("%s\n", BldInfo::name(static_cast<const BldInfo&>(xtc->src)));
      break;
    default: {
      const ProcInfo& pinfo = static_cast<const ProcInfo&>(xtc->src);
      printf("IpAddress 0x%x ProcessId 0x%x\n",pinfo.ipAddr(),pinfo.processId());
      break;
    }
    }
    if (level < 0 || level >= Level::NumberOfLevels )
    {
        printf("Unsupported Level %d\n", (int) level);
        return Continue;
    }    
    switch (xtc->contains.id()) {
    case (TypeId::Id_Xtc) : {
      myLevelIter iter(xtc,_depth+1);
      iter.iterate();
      break;
    }
    case (TypeId::Id_AcqConfig) :
    {      
      unsigned version = xtc->contains.version();
      switch (version) {
      case 1:
        process(info,*(const Acqiris::ConfigV1*)(xtc->payload()));
        break;
      default:
        printf("Unsupported acqiris configuration version %d\n",version);
        break;
      }
      break;      
    }
    case (TypeId::Id_IpimbConfig) :
    {      
      unsigned version = xtc->contains.version();
      switch (version) {
      case 1:
        process(info,*(const Ipimb::ConfigV1*)(xtc->payload()));
        break;
      case 2:
        process(info,*(const Ipimb::ConfigV2*)(xtc->payload()));
        break;
      default:
        printf("Unsupported ipimb configuration version %d\n",version);
        break;
      }
      break;      
    }
    case (TypeId::Id_EncoderConfig) :
    {      
      unsigned version = xtc->contains.version();
      switch (version) {
      case 1:
        process(info,*(const Encoder::ConfigV1*)(xtc->payload()));
        break;
      default:
        printf("Unsupported encoder configuration version %d\n",version);
        break;
      }
      break;      
    }
    case (TypeId::Id_Opal1kConfig) :
      process(info, *(const Opal1k::ConfigV1*)(xtc->payload()));
      break;
    case (TypeId::Id_FrameFexConfig) :
      process(info, *(const Camera::FrameFexConfigV1*)(xtc->payload()));
      break;
    case (TypeId::Id_pnCCDconfig) :
      {
      unsigned version = xtc->contains.version();
      switch (version) {
        case 1:
          process(info, *(const PNCCD::ConfigV1*)(xtc->payload()));
          break;
        case 2:
          process(info, *(const PNCCD::ConfigV2*)(xtc->payload()));
          break;
        default:
          printf("Unsupported pnCCD configuration version %d\n",version);
      }
      break;
      }
    case (TypeId::Id_EvrIOConfig) :
      {      
      process(info, *(const EvrData::IOConfigV1*)(xtc->payload()));
      break;
      }
    case (TypeId::Id_EvrConfig) :
    {      
      unsigned version = xtc->contains.version();
      switch (version) {
      case 1:
        process(info, *(const EvrData::ConfigV1*)(xtc->payload()));
        break;
      case 2:
        process(info, *(const EvrData::ConfigV2*)(xtc->payload()));
        break;
      case 3:
        process(info, *(const EvrData::ConfigV3*)(xtc->payload()));
        break;
      case 4:
        process(info, *(const EvrData::ConfigV4*)(xtc->payload()));
        break;
      case 7:
        process(info, *(const EvrData::ConfigV7*)(xtc->payload()));
        break;
      default:
        printf("Unsupported evr configuration version %d\n",version);
        break;
      }
      break;      
    }
    case (TypeId::Id_ControlConfig) :
      process(info, *(const ControlData::ConfigV1*)(xtc->payload()));
      break;
    case (TypeId::Id_Epics) :      
    {
      if (!noEpics) {
        int iVersion = xtc->contains.version();
        if ( iVersion != EpicsXtcSettings::iXtcVersion ) 
          {
            printf( "Xtc Epics version (%d) is not compatible with reader supported version (%d)", iVersion, EpicsXtcSettings::iXtcVersion );
            break;
          }
        process(info, *(const EpicsPvHeader*)(xtc->payload()));
      }
      break;
    }
    /*
     * BLD data
     */
    case (TypeId::Id_FEEGasDetEnergy) :
    {
      process(info, *(const BldDataFEEGasDetEnergy*) xtc->payload() );
      break;        
    }
    case (TypeId::Id_EBeam) :
    {
      switch(xtc->contains.version()) {
      case 0:
        process(info, *(const BldDataEBeamV0*) xtc->payload() );
        break; 
      case 1:
        process(info, *(const BldDataEBeamV1*) xtc->payload() );
        break; 
      case 2:
        process(info, *(const BldDataEBeam*) xtc->payload() );
        break; 
      default:
        break;
      }       
    }    
    case (TypeId::Id_PhaseCavity) :
    {
      process(info, *(const BldDataPhaseCavity*) xtc->payload() );
      break;        
    }
    case (TypeId::Id_GMD) :
    {
      switch(xtc->contains.version()) {
        case 0:
          process(info, *(const BldDataGMDV0*) xtc->payload() );
          break;
        case 1:
          process(info, *(const BldDataGMDV1*) xtc->payload() );
          break;
        default:
          break;
      }
      break;
      break;
    }
    case (TypeId::Id_SharedIpimb) :
    {

     switch(xtc->contains.version()) {
      case 0:
        process(info, *(const BldDataIpimbV0*) xtc->payload() );
        break; 
      case 1:
        process(info, *(const BldDataIpimb*) xtc->payload() );
        break; 
      default:
        break;
      }   
      break; 
       
    }
    case (TypeId::Id_PrincetonConfig) :
    {
      process(info, *(const Princeton::ConfigV1*)(xtc->payload()));
      break;
    }
    case (TypeId::Id_CspadConfig) :
    {
      process(info, *(const CsPad::ConfigV4*)(xtc->payload()));
      break;
    }
    default :
      break;
    }
    return Continue;
  }
private:
  unsigned       _depth;
  const char*    _hdr;

  /* static private data */
  static PNCCD::ConfigV1 _pnCcdCfgListV1[2];
  static PNCCD::ConfigV2 _pnCcdCfgListV2[2];
};

PNCCD::ConfigV1 myLevelIter::_pnCcdCfgListV1[2] = { PNCCD::ConfigV1(), PNCCD::ConfigV1() };
PNCCD::ConfigV2 myLevelIter::_pnCcdCfgListV2[2] = { PNCCD::ConfigV2(), PNCCD::ConfigV2() };

void usage(char* progname) {
  fprintf(stderr,"Usage: %s -f <filename> [-h]\n", progname);
}

int main(int argc, char* argv[]) {
  int c;
  char* xtcname=0;
  int parseErr = 0;
  noEpics=false;
  unsigned nevents=0;

  while ((c = getopt(argc, argv, "hef:n:")) != -1) {
    switch (c) {
    case 'e':
      noEpics=true;
      break;
    case 'h':
      usage(argv[0]);
      exit(0);
    case 'f':
      xtcname = optarg;
      break;
    case 'n':
      nevents = strtoul(optarg,NULL,0);
      break;
    default:
      parseErr++;
    }
  }
  
  if (!xtcname) {
    usage(argv[0]);
    exit(2);
  }

  int fd = open(xtcname, O_RDONLY | O_LARGEFILE);
  if (fd < 0) {
    perror("Unable to open file %s\n");
    exit(2);
  }

  unsigned count(0);
  XtcFileIterator iter(fd,0x900000);
  Dgram* dg;
  while( dg = iter.next() ) {
    if (dg->seq.service() == TransitionId::L1Accept)
      if (++count>nevents) break;

    printf("%s transition: time 0x%x/0x%x, payloadSize %d\n",TransitionId::name(dg->seq.service()),
            dg->seq.stamp().fiducials(),dg->seq.stamp().ticks(), dg->xtc.sizeofPayload());
    myLevelIter liter(&(dg->xtc),0);
    liter.iterate();
  }
  ::close(fd);
  return 0;
}
