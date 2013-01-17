#include "pdsdata/andor/ConfigV1.hh"
#include "pdsdata/andor/FrameV1.hh"

#include <string.h>

using namespace Pds;
using namespace Andor;

ConfigV1::ConfigV1(
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
 uint32_t         u32NumDelayShots) :
 _uWidth  (uWidth), 
 _uHeight (uHeight), 
 _uOrgX   (uOrgX), 
 _uOrgY   (uOrgY), 
 _uBinX   (uBinX), 
 _uBinY   (uBinY),
 _f32ExposureTime       (f32ExposureTime),
 _f32CoolingTemp        (f32CoolingTemp), 
 _u8FanMode             (u8FanMode), 
 _u8BaselineClamp       (u8BaselineClamp), 
 _u8HighCapacity        (u8HighCapacity), 
 _u8GainIndex           (u8GainIndex),  
 _u16ReadoutSpeedIndex  (u16ReadoutSpeedIndex), 
 _u16ExposureEventCode  (u16ExposureEventCode),
 _u32NumDelayShots      (u32NumDelayShots)
 {} 
 
int ConfigV1::frameSize() const
{
  // Note: This formula is different from Princeton image
  return sizeof(FrameV1) + 
    (int) (_uWidth / _uBinX ) * 
    (int) (_uHeight/ _uBinY ) * 2; // 2 -> 16 bit color depth  
  //return sizeof(FrameV1) + 4*1024*1024*2; // 2 -> 16 bit color depth // !! debug
}
