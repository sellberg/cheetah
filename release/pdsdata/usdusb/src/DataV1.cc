#include "pdsdata/usdusb/DataV1.hh"

using namespace Pds::UsdUsb;

//
//  Convert to 24-bit signed integer
//
int      DataV1::encoder_count(unsigned i) const
{ 
  int32_t v(_count[i]<<8); 
  return v>>8;
}

unsigned DataV1::analog_in    (unsigned i) const
{ return _ain[i]; }

unsigned DataV1::digital_in   () const
{ return _din; }

unsigned DataV1::timestamp    () const
{ return _timestamp; }
