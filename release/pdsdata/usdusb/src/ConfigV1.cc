#include "pdsdata/usdusb/ConfigV1.hh"

#include <stdlib.h>

using namespace Pds::UsdUsb;

ConfigV1::ConfigV1( Count_Mode cm[],
		    Quad_Mode  qm[] ) :
  _count_mode( (uint32_t[4]){cm[0], cm[1], cm[2], cm[3]} ),  
  _quad_mode ( (uint32_t[4]){qm[0], qm[1], qm[2], qm[3]} )
{
}

ConfigV1::Count_Mode ConfigV1::counting_mode   (unsigned channel) const
{ 
  return Count_Mode(_count_mode[channel]);
}

ConfigV1::Quad_Mode  ConfigV1::quadrature_mode (unsigned channel) const
{
  return Quad_Mode(_quad_mode[channel]);
}

static const char* _count_mode_labels[] = { "WRAP_FULL",
					    "LIMIT",
					    "HALT",
					    "WRAP_PRESET",
					    NULL };

static const char* _quad_mode_labels[] = { "CLOCK_DIR",
					   "X1",
					   "X2",
					   "X4",
					   NULL };

const char** ConfigV1::count_mode_labels() 
{
  return _count_mode_labels;
}

const char** ConfigV1::quad_mode_labels()
{
  return _quad_mode_labels;
}
