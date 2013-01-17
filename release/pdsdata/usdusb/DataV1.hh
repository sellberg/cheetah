#ifndef PdsUsdUsb_DataV1_hh
#define PdsUsdUsb_DataV1_hh

#include <stdint.h>

namespace Pds {
  namespace UsdUsb {

    class DataV1 {
    public:					
      enum { Version = 1 };
      enum { Encoder_Inputs = 4 };
      enum { Analog_Inputs  = 4 };
      enum { Digital_Inputs = 8 };
    public:
      DataV1() {}
    public:
      int      encoder_count(unsigned) const;
      unsigned analog_in    (unsigned) const;
      unsigned digital_in   () const;
      unsigned timestamp    () const;
    private:
      uint8_t  _header[6];
      uint8_t  _din;
      uint8_t  _estop;
      uint32_t _timestamp;
      uint32_t _count[4];
      uint8_t  _status[4];
      uint16_t _ain[4];
    };
  };
};

#endif
