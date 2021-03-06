#include "IpmFexHandler.hh"

#include "ami/data/FeatureCache.hh"

#include "pdsdata/lusi/IpmFexV1.hh"
#include "pdsdata/xtc/BldInfo.hh"
#include "pdsdata/xtc/DetInfo.hh"

#include <stdio.h>
#include <string.h>

using namespace Ami;

IpmFexHandler::IpmFexHandler(const Pds::DetInfo& info, FeatureCache& f) :
  EventHandler(info,
	       Pds::TypeId::Id_IpmFex,
	       Pds::TypeId::Id_IpmFexConfig),
  _cache(f)
{
}

IpmFexHandler::~IpmFexHandler()
{
}

void   IpmFexHandler::_calibrate(const void* payload, const Pds::ClockTime& t) {}
void   IpmFexHandler::_configure(const void* payload, const Pds::ClockTime& t)
{
  char buffer[64];
  char* iptr;
  unsigned i=0;

  switch(info().level()) {
  case Level::Reporter:
    strncpy(buffer,Pds::BldInfo::name(static_cast<const Pds::BldInfo&>(info())),60);
    iptr = buffer+strlen(buffer);
    sprintf(iptr,":FEX:CH0");  _index[i] = _cache.add(buffer);  i++;
    sprintf(iptr,":FEX:CH1");  _index[i] = _cache.add(buffer);  i++;
    sprintf(iptr,":FEX:CH2");  _index[i] = _cache.add(buffer);  i++;
    sprintf(iptr,":FEX:CH3");  _index[i] = _cache.add(buffer);  i++;
    sprintf(iptr,":FEX:SUM");  _index[i] = _cache.add(buffer);  i++;
    sprintf(iptr,":FEX:XPOS"); _index[i] = _cache.add(buffer);  i++;
    sprintf(iptr,":FEX:YPOS"); _index[i] = _cache.add(buffer);  i++;
    break;
  case Level::Source:
  default:
    strncpy(buffer,Pds::DetInfo::name(static_cast<const Pds::DetInfo&>(info())),60);
    iptr = buffer+strlen(buffer);
    sprintf(iptr,":CH0"); _index[i++] = _cache.add(buffer);
    sprintf(iptr,":CH1"); _index[i++] = _cache.add(buffer);
    sprintf(iptr,":CH2"); _index[i++] = _cache.add(buffer);
    sprintf(iptr,":CH3"); _index[i++] = _cache.add(buffer);
    sprintf(iptr,":SUM");
    _index[i++] = _cache.add(buffer);
    sprintf(iptr,":XPOS");
    _index[i++] = _cache.add(buffer);
    sprintf(iptr,":YPOS");
    _index[i++] = _cache.add(buffer);
    break;
  }
}

void   IpmFexHandler::_event    (const void* payload, const Pds::ClockTime& t)
{
  const Pds::Lusi::IpmFexV1& d = *reinterpret_cast<const Pds::Lusi::IpmFexV1*>(payload);

  unsigned i=0;
  _cache.cache(_index[i], d.channel[i]); i++;
  _cache.cache(_index[i], d.channel[i]); i++;
  _cache.cache(_index[i], d.channel[i]); i++;
  _cache.cache(_index[i], d.channel[i]); i++;
  _cache.cache(_index[i], d.sum       ); i++;
  _cache.cache(_index[i], d.xpos      ); i++;
  _cache.cache(_index[i], d.ypos      ); i++;
}

void   IpmFexHandler::_damaged  ()
{
  for(unsigned i=0; i<NChannels; i++)
    _cache.cache(_index[i], 0, true);
}

//  No Entry data
unsigned     IpmFexHandler::nentries() const { return 0; }
const Entry* IpmFexHandler::entry   (unsigned) const { return 0; }
void         IpmFexHandler::reset   () 
{
}
