#ifndef Evr_EventCodeV4_hh
#define Evr_EventCodeV4_hh

#include <stdint.h>

#pragma pack(4)

namespace Pds
{
namespace EvrData
{

class EventCodeV4
{
public:
  EventCodeV4(
    uint16_t u16Code,
    bool     bReadout,
    bool     bTerminator,    
    uint32_t u32ReportDelay = 0,
    uint32_t u32ReportWidth = 1,
    uint32_t u32MaskTrigger = 0,
    uint32_t u32MaskSet     = 0, 
    uint32_t u32MaskClear   = 0);
    
  EventCodeV4() {} // For array initialization

  uint16_t  code        () const { return _u16Code; }
  bool      isReadout   () const { return ( _u16MaskEventAttr & (1<<EventAttrBitReadout) )     != 0; }
  bool      isTerminator() const { return ( _u16MaskEventAttr & (1<<EventAttrBitTerminator) )  != 0; }
  uint32_t  reportDelay () const { return _u32ReportDelay; }
  uint32_t  reportWidth () const { return _u32ReportWidth; }
  uint32_t  maskTrigger () const { return _u32MaskTrigger; }
  uint32_t  maskSet     () const { return _u32MaskSet; }
  uint32_t  maskClear   () const { return _u32MaskClear; }

private:
  enum EventAttrBitEnum { EventAttrBitReadout = 0, EventAttrBitTerminator = 1 };
    
  uint16_t _u16Code;
  uint16_t _u16MaskEventAttr;
  uint32_t _u32ReportDelay;
  uint32_t _u32ReportWidth;
  uint32_t _u32MaskTrigger;
  uint32_t _u32MaskSet;
  uint32_t _u32MaskClear;
};

} // namespace EvrData
} // namespace Pds

#pragma pack()

#endif
