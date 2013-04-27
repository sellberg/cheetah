//
//  Class for configuration of Adimec Opal-1000 monochrome camera
//
#ifndef Orca_ConfigV1_hh
#define Orca_ConfigV1_hh

#include <stdint.h>

namespace Pds {

  class DetInfo;

  namespace Camera {
    class FrameCoord;
  };

  namespace Orca {

    class ConfigV1 {
    public:
      enum { Version=1 };
      enum ReadoutMode { x1, x2, x4, Subarray };
      enum Cooling { Off, On, Max };
      enum { Row_Pixels=2048 };
      enum { Column_Pixels=2048 };

      ConfigV1();
      ConfigV1(ReadoutMode     rmode,
	       unsigned        rows,
	       Cooling         cmode,
	       bool            enable_pixel_correction );

      //
      //  Accessors
      //
      
      ReadoutMode      mode   () const;
      unsigned         rows   () const;
      Cooling          cooling() const;
      bool             defect_pixel_correction_enabled() const;

      //  total size of this structure 
      //  (including defective pixel coords and output lookup table)
      unsigned         size() const;

    private:
      uint32_t _options;        // bit mask of enumerations
      uint32_t _rows;           // sub array dimensions
    };

  };
};

#endif
