#include "pdsdata/control/ConfigV2.hh"

#include "pdsdata/control/PVControl.hh"
#include "pdsdata/control/PVMonitor.hh"
#include "pdsdata/control/PVLabel.hh"

using namespace Pds::ControlData;

enum { EventsMask   = 0x3fffffff,
       UsesDuration = 0x40000000,
       UsesEvents   = 0x80000000 };

static void appendChannels(const std::list<PVControl>& pvcs, 
         const std::list<PVMonitor>& pvms, 
         const std::list<PVLabel > & pvls, 
         ConfigV2* s)
{
  std::list<PVControl> sorted_pvcs (pvcs); sorted_pvcs.sort();
  std::list<PVMonitor> sorted_pvms (pvms); sorted_pvms.sort();
  std::list<PVLabel > sorted_pvls (pvls); sorted_pvls.sort();

  PVControl* c = reinterpret_cast<PVControl*>(s+1);
  for(std::list<PVControl>::const_iterator iter = sorted_pvcs.begin();
      iter != sorted_pvcs.end(); ++iter)
    *c++ = *iter;
  PVMonitor* m = reinterpret_cast<PVMonitor*>(c);
  for(std::list<PVMonitor>::const_iterator iter = sorted_pvms.begin();
      iter != sorted_pvms.end(); ++iter)
    *m++ = *iter;
  PVLabel*  p = reinterpret_cast<PVLabel*>(m);
  for(std::list<PVLabel >::const_iterator iter = sorted_pvls.begin();
      iter != sorted_pvls.end(); ++iter)
    *p++ = *iter;
}

ConfigV2::ConfigV2() {}

ConfigV2::ConfigV2(Initialize) :
  _control    (0),
  _reserved   (0),
  _npvControls(0),
  _npvMonitors(0),
  _npvLabels (0)
{
}

ConfigV2::ConfigV2(const std::list<PVControl>& pvcs,const std::list<PVMonitor>& pvms,const std::list<PVLabel>& pvls) :
  _control    (0),
  _reserved   (0),
  _duration   (0,0),
  _npvControls(pvcs.size()),
  _npvMonitors(pvms.size()),
  _npvLabels  (pvls.size())
{
  appendChannels(pvcs,pvms,pvls,this);
}

ConfigV2::ConfigV2(const std::list<PVControl>& pvcs,const std::list<PVMonitor>& pvms,const std::list<PVLabel>& pvls,
       const ClockTime& t) :
  _control    (UsesDuration),
  _reserved   (0),
  _duration   (t),
  _npvControls(pvcs.size()),
  _npvMonitors(pvms.size()),
  _npvLabels (pvls.size())
{
  appendChannels(pvcs,pvms,pvls,this);
}

ConfigV2::ConfigV2(const std::list<PVControl>& pvcs,const std::list<PVMonitor>& pvms,const std::list<PVLabel>& pvls,
       unsigned events ) :
  _control    (UsesEvents | (events&EventsMask)),
  _reserved   (0),
  _duration   (0,0),
  _npvControls(pvcs.size()),
  _npvMonitors(pvms.size()),
  _npvLabels  (pvls.size())
{
  appendChannels(pvcs,pvms,pvls,this);
}

ConfigV2::ConfigV2(const ConfigV2& s) :
  _control    (s._control),
  _reserved   (0),
  _npvControls(s._npvControls),
  _npvMonitors(s._npvMonitors),
  _npvLabels  (s._npvLabels )
{
  PVControl* c = reinterpret_cast<PVControl*>(this+1);
  for(unsigned k=0; k<_npvControls; k++)
    *c++ = s.pvControl(k);
  PVMonitor* m = reinterpret_cast<PVMonitor*>(c);
  for(unsigned k=0; k<_npvMonitors; k++)
    *m++ = s.pvMonitor(k);
  PVLabel*   p = reinterpret_cast<PVLabel*  >(m);
  for(unsigned k=0; k<_npvLabels; k++)
    *p++ = s.pvLabel  (k);
}

ConfigV2::~ConfigV2() {}

bool             ConfigV2::uses_duration()       const
{ return _control & UsesDuration; }

bool             ConfigV2::uses_events  ()       const
{ return _control & UsesEvents; }

const ClockTime& ConfigV2::duration   ()         const 
{ return _duration; }

unsigned         ConfigV2::events     ()         const 
{ return _control & EventsMask; }

unsigned         ConfigV2::npvControls()         const
{ return _npvControls; }

const PVControl& ConfigV2::pvControl  (unsigned i) const
{ return reinterpret_cast<const PVControl*>(this+1)[i]; }

unsigned         ConfigV2::npvMonitors()         const
{ return _npvMonitors; }

const PVMonitor& ConfigV2::pvMonitor  (unsigned i) const
{ return reinterpret_cast<const PVMonitor*>(&pvControl(_npvControls))[i]; }

unsigned         ConfigV2::npvLabels()         const
{ return _npvLabels; }

const PVLabel& ConfigV2::pvLabel      (unsigned i) const
{ return reinterpret_cast<const PVLabel*>(&pvMonitor(_npvMonitors))[i]; }

unsigned         ConfigV2::size       ()         const
{ return sizeof(ConfigV2) + _npvControls*sizeof(PVControl) + _npvMonitors*sizeof(PVMonitor) + _npvLabels*sizeof(PVLabel); }
