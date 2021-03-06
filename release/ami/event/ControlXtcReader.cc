#include "ControlXtcReader.hh"
#include "ami/data/FeatureCache.hh"

#include "pdsdata/xtc/Dgram.hh"
#include "pdsdata/xtc/Level.hh"
#include "pdsdata/xtc/Xtc.hh"
#include "pdsdata/xtc/ProcInfo.hh"
#include "pdsdata/control/ConfigV1.hh"
#include "pdsdata/control/ConfigV2.hh"
#include "pdsdata/control/PVControl.hh"

#include <stdio.h>

using namespace Ami;

ControlXtcReader::ControlXtcReader(FeatureCache& f)  : 
  EventHandler(Pds::ProcInfo(Pds::Level::Control,0,-1UL),
	       Pds::TypeId::NumberOf,
	       Pds::TypeId::Id_ControlConfig),
  _cache(f),
  _index(-1),
  _values(0)
{
}

ControlXtcReader::~ControlXtcReader()
{
}

void   ControlXtcReader::_calibrate(Pds::TypeId id, const void* payload, const Pds::ClockTime& t)
{ _configure(id, payload,t); }

void   ControlXtcReader::_configure(Pds::TypeId id, const void* payload, const Pds::ClockTime& t)
{
  _index = -1;
  _values.resize(0);

  if (id.version()==Pds::ControlData::ConfigV1::Version) {
    const Pds::ControlData::ConfigV1& c =
      *reinterpret_cast<const Pds::ControlData::ConfigV1*>(payload);

    //  Add the PVControl variables
    char nbuf[FeatureCache::FEATURE_NAMELEN];
    for(unsigned k=0; k<c.npvControls(); k++) {
      const Pds::ControlData::PVControl& pv = c.pvControl(k);
      int index;
      if (pv.array()) {
        snprintf(nbuf, FeatureCache::FEATURE_NAMELEN, "%s[%d](SCAN)", pv.name(), pv.index());
        index = _cache.add(nbuf);
      }
      else {
        snprintf(nbuf, FeatureCache::FEATURE_NAMELEN, "%s(SCAN)", pv.name());
        index = _cache.add(nbuf);
      }
      _cache.cache(index,pv.value());
      if (_index<0) _index = index;
      _values.push_back(pv.value());
    }
  }
  else if (id.version()==Pds::ControlData::ConfigV2::Version) {
    const Pds::ControlData::ConfigV2& c =
      *reinterpret_cast<const Pds::ControlData::ConfigV2*>(payload);

    //  Add the PVControl variables
    char nbuf[FeatureCache::FEATURE_NAMELEN];
    for(unsigned k=0; k<c.npvControls(); k++) {
      const Pds::ControlData::PVControl& pv = c.pvControl(k);
      int index;
      if (pv.array()) {
        snprintf(nbuf, FeatureCache::FEATURE_NAMELEN, "%s[%d](SCAN)", pv.name(), pv.index());
        index = _cache.add(nbuf);
      }
      else {
        snprintf(nbuf, FeatureCache::FEATURE_NAMELEN, "%s(SCAN)", pv.name());
        index = _cache.add(nbuf);
      }
      _cache.cache(index,pv.value());
      if (_index<0) _index = index;
      _values.push_back(pv.value());
    }
  }
}

//  no L1 data will appear from Control
void   ControlXtcReader::_event    (Pds::TypeId, const void* payload, const Pds::ClockTime& t) 
{
  for(unsigned i=0, index=_index; i<_values.size(); i++,index++)
    _cache.cache(index,_values[i]);
}

void   ControlXtcReader::_damaged  () {}

//  No Entry data
unsigned     ControlXtcReader::nentries() const { return 0; }
const Entry* ControlXtcReader::entry   (unsigned) const { return 0; }
void         ControlXtcReader::reset   () { _index=-1; }
