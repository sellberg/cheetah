#ifndef ANDOR_CONFIG_V1_HH
#define ANDOR_CONFIG_V1_HH

#include <stdint.h>
#include "ConfigV1.hh"

#pragma pack(4)

namespace Pds 
{

namespace Andor 
{

class ConfigV1
{
public:
  enum { Version = 1 };
  enum EnumFanMode
  {
    ENUM_FAN_FULL     = 0,
    ENUM_FAN_LOW      = 1,
    ENUM_FAN_OFF      = 2,
    ENUM_FAN_ACQOFF   = 3,
    ENUM_FAN_NUM      = 4
  };

  ConfigV1()  {}  
  ConfigV1(
   uint32_t         uWidth, 
   uint32_t         uHeight, 
   uint32_t         uOrgX, 
   uint32_t         uOrgY, 
   uint32_t         uBinX, 
   uint32_t         uBinY,
   float            f32ExposureTime, 
   float            f32CoolingTemp, 
   uint8_t          u8FanMode,
   uint8_t          u8BaselineClamp,
   uint8_t          u8HighCapacity,
   uint8_t          u8GainIndex,
   uint16_t         u16ReadoutSpeedIndex,
   uint16_t         u16ExposureEventCode, 
   uint32_t         u32NumDelayShots
   );  

  uint32_t          width ()            const         { return _uWidth; }
  uint32_t          height()            const         { return _uHeight; }
  uint32_t          orgX  ()            const         { return _uOrgX; }
  uint32_t          orgY  ()            const         { return _uOrgY; }    
  uint32_t          binX  ()            const         { return _uBinX; }
  uint32_t          binY  ()            const         { return _uBinY; }    
  float             exposureTime  ()    const         { return _f32ExposureTime; }
  float             coolingTemp   ()    const         { return _f32CoolingTemp; }
  uint8_t           fanMode       ()    const         { return _u8FanMode; }
  uint8_t           baselineClamp ()    const         { return _u8BaselineClamp; }
  uint8_t           highCapacity  ()    const         { return _u8HighCapacity; }
  uint8_t           gainIndex     ()    const         { return _u8GainIndex; }
  uint16_t          readoutSpeedIndex() const         { return _u16ReadoutSpeedIndex; }
  
  uint16_t          exposureEventCode()  const         { return _u16ExposureEventCode; }
  uint32_t          numDelayShots()         const         { return _u32NumDelayShots; }

  uint32_t          setWidth    (uint32_t uWidth)     { return _uWidth = uWidth; }
  uint32_t          setHeight   (uint32_t uHeight)    { return _uHeight = uHeight; }
  uint16_t          setReadoutSpeedIndex
                            (uint16_t uSpeedIndex)    { return _u16ReadoutSpeedIndex = uSpeedIndex; }
  uint32_t          setNumDelayShots
                            (uint32_t uNumDelayShots) { return _u32NumDelayShots = uNumDelayShots;  }
  
  int               size      ()        const         { return sizeof(*this); }
  int               frameSize ()        const; // calculate the frame size based on the current ROI and binning settings
    
private:
  uint32_t          _uWidth, _uHeight;
  uint32_t          _uOrgX,  _uOrgY;
  uint32_t          _uBinX,  _uBinY;
  float             _f32ExposureTime;
  float             _f32CoolingTemp;
  uint8_t           _u8FanMode;
  uint8_t           _u8BaselineClamp;
  uint8_t           _u8HighCapacity;
  uint8_t           _u8GainIndex;
  uint16_t          _u16ReadoutSpeedIndex;
  uint16_t          _u16ExposureEventCode;
  uint32_t          _u32NumDelayShots;
};

} // namespace Andor

} // namespace Pds 

#pragma pack()

#endif //#ifndef ANDOR_CONFIG_V1_HH
