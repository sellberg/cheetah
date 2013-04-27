#include "pdsdata/orca/ConfigV1.hh"
#include "pdsdata/camera/FrameCoord.hh"
#include "pdsdata/xtc/DetInfo.hh"

#include <string.h>

using namespace Pds;
using namespace Orca;

enum { Mode_Offset               = 0,
       Cooling_Offset            = 2,
       Pixel_Corr_Enable_Offset  = 4 };

ConfigV1::ConfigV1() {}

ConfigV1::ConfigV1(ReadoutMode rmode,
		   unsigned    rows,
		   Cooling     cmode,
		   bool        enable_pixel_correction ) :
  _options         ( ( rmode                << Mode_Offset) |
		     ( cmode                << Cooling_Offset) |
                     ( enable_pixel_correction ? 1 << Pixel_Corr_Enable_Offset : 0) ),
  _rows            ( rows )
{
}

ConfigV1::ReadoutMode     ConfigV1::mode() const
{ return ReadoutMode((_options >> Mode_Offset)&0x3); }

unsigned                  ConfigV1::rows   () const
{ return _rows; }

ConfigV1::Cooling         ConfigV1::cooling() const
{ return Cooling((_options >> Cooling_Offset)&0x3); }

bool            ConfigV1::defect_pixel_correction_enabled() const
{ return _options & (1<<Pixel_Corr_Enable_Offset); }

unsigned        ConfigV1::size() const
{ 
  return sizeof(*this);
}

