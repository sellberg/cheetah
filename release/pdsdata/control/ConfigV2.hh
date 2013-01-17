//
//  Class for DAQ Control
//
#ifndef Pds_ConfigV2_hh
#define Pds_ConfigV2_hh

#include <list>
#include <stdint.h>

#include "pdsdata/xtc/ClockTime.hh"
using Pds::ClockTime;

namespace Pds {

  namespace ControlData {

    class PVControl;
    class PVMonitor;
    class PVLabel;

    class ConfigV2 {
    public:
      enum { Version=2 };
    public:
      enum Initialize { Default };
      ConfigV2();
      ConfigV2(Initialize);
      ConfigV2(const std::list<PVControl>&, const std::list<PVMonitor>&, const std::list<PVLabel>&);
      ConfigV2(const std::list<PVControl>&, const std::list<PVMonitor>&, const std::list<PVLabel>&, const ClockTime&);
      ConfigV2(const std::list<PVControl>&, const std::list<PVMonitor>&, const std::list<PVLabel>&, unsigned events );
      ConfigV2(const ConfigV2&);
    private:
      ~ConfigV2();  // Should not be built on the stack (placement new only)
    public:
      bool             uses_duration()       const;
      bool             uses_events  ()       const;
      const ClockTime& duration   ()         const;
      unsigned         events     ()         const;
      unsigned         npvControls()         const;
      const PVControl& pvControl  (unsigned) const;
      unsigned         npvMonitors()         const;
      const PVMonitor& pvMonitor  (unsigned) const;
      unsigned         npvLabels  ()         const;
      const PVLabel&   pvLabel    (unsigned) const;
      unsigned         size       ()         const;
    private:
      uint32_t  _control;
      uint32_t  _reserved;
      ClockTime _duration;
      uint32_t  _npvControls;
      uint32_t  _npvMonitors;
      uint32_t  _npvLabels;
    };

  };

};

#endif
