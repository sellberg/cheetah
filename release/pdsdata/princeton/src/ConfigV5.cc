#include "pdsdata/princeton/ConfigV5.hh"
#include "pdsdata/princeton/FrameV2.hh"

#include <string.h>

using namespace Pds;
using namespace Princeton;

ConfigV5::ConfigV5(
 uint32_t         u32Width,
 uint32_t         u32Height,
 uint32_t         u32OrgX,
 uint32_t         u32OrgY,
 uint32_t         u32BinX,
 uint32_t         u32BinY,
 float            f32ExposureTime,
 float            f32CoolingTemp,
 uint16_t         u16GainIndex,
 uint16_t         u16ReadoutSpeedIndex,
 uint32_t         u32MaskedHeight,
 uint32_t         u32KineticHeight,
 float            f32VsSpeed,
 int16_t          i16InfoReportInterval,
 uint16_t         u16ExposureEventCode,
 uint32_t         u32NumDelayShots) :
 _u32Width  (u32Width),
 _u32Height (u32Height),
 _uOrgX     (u32OrgX),
 _uOrgY     (u32OrgY),
 _uBinX     (u32BinX),
 _uBinY     (u32BinY),
 _f32ExposureTime       (f32ExposureTime),
 _f32CoolingTemp        (f32CoolingTemp),
 _u16GainIndex          (u16GainIndex),
 _u16ReadoutSpeedIndex  (u16ReadoutSpeedIndex),
 _u32MaskedHeight       (u32MaskedHeight),
 _u32KineticHeight      (u32KineticHeight),
 _f32VsSpeed            (f32VsSpeed),
 _i16InfoReportInterval (i16InfoReportInterval),
 _u16ExposureEventCode  (u16ExposureEventCode),
 _u32NumDelayShots      (u32NumDelayShots)
 {}

ConfigV5::ConfigV5(const ConfigV1& config) :
 _u32Width  (config._uWidth),
 _u32Height (config._uHeight),
 _uOrgX   (config._uOrgX),
 _uOrgY   (config._uOrgY),
 _uBinX   (config._uBinX),
 _uBinY   (config._uBinY),
 _f32ExposureTime       (config._f32ExposureTime),
 _f32CoolingTemp        (config._f32CoolingTemp),
 _u16GainIndex          (0),
 _u16ReadoutSpeedIndex  (config._u32ReadoutSpeedIndex),
 _u32MaskedHeight       (0),
 _u32KineticHeight      (0),
 _f32VsSpeed            (0),
 _i16InfoReportInterval (0),
 _u16ExposureEventCode  (config._u16ReadoutEventCode),
 _u32NumDelayShots      (config._u16DelayMode)
{
}

ConfigV5::ConfigV5(const ConfigV2& config) :
 _u32Width  (config._uWidth),
 _u32Height (config._uHeight),
 _uOrgX   (config._uOrgX),
 _uOrgY   (config._uOrgY),
 _uBinX   (config._uBinX),
 _uBinY   (config._uBinY),
 _f32ExposureTime       (config._f32ExposureTime),
 _f32CoolingTemp        (config._f32CoolingTemp),
 _u16GainIndex          (config._u16GainIndex),
 _u16ReadoutSpeedIndex  (config._u16ReadoutSpeedIndex),
 _u32MaskedHeight       (0),
 _u32KineticHeight      (0),
 _f32VsSpeed            (0),
 _i16InfoReportInterval (0),
 _u16ExposureEventCode  (config._u16ReadoutEventCode),
 _u32NumDelayShots      (config._u16DelayMode)
{
}

ConfigV5::ConfigV5(const ConfigV3& config) :
 _u32Width  (config._uWidth),
 _u32Height (config._uHeight),
 _uOrgX   (config._uOrgX),
 _uOrgY   (config._uOrgY),
 _uBinX   (config._uBinX),
 _uBinY   (config._uBinY),
 _f32ExposureTime       (config._f32ExposureTime),
 _f32CoolingTemp        (config._f32CoolingTemp),
 _u16GainIndex          (config._u8GainIndex),
 _u16ReadoutSpeedIndex  (config._u8ReadoutSpeedIndex),
 _u32MaskedHeight       (0),
 _u32KineticHeight      (0),
 _f32VsSpeed            (0),
 _i16InfoReportInterval (0),
 _u16ExposureEventCode  (config._u16ExposureEventCode),
 _u32NumDelayShots      (config._u32NumDelayShots)
{
}

ConfigV5::ConfigV5(const ConfigV4& config) :
  _u32Width   (config._u32Width),
  _u32Height  (config._u32Width),
  _uOrgX      (config._uOrgX),
  _uOrgY      (config._uOrgY),
  _uBinX      (config._uBinX),
  _uBinY      (config._uBinY),
  _f32ExposureTime       (config._f32ExposureTime),
  _f32CoolingTemp        (config._f32CoolingTemp),
  _u16GainIndex          (config._u8GainIndex),
  _u16ReadoutSpeedIndex  (config._u8ReadoutSpeedIndex),
  _u32MaskedHeight       (config._u32MaskedHeight),
  _u32KineticHeight      (config._u32KineticHeight),
  _f32VsSpeed            (config._f32VsSpeed),
  _i16InfoReportInterval (0),
  _u16ExposureEventCode  (config._u16ExposureEventCode),
  _u32NumDelayShots      (config._u32NumDelayShots)
{
}

int ConfigV5::frameSize() const
{
  return sizeof(FrameV2) +
    (int) ((_u32Width + _uBinX-1)/ _uBinX ) *
    (int) ((_u32Height+ _uBinY-1)/ _uBinY ) * 2; // 2 -> 16 bit color depth
}
