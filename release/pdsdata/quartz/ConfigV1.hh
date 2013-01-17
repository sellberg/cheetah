//
//  Class for configuration of Adimec Opal-1000 monochrome camera
//
#ifndef Quartz_ConfigV1_hh
#define Quartz_ConfigV1_hh

#include <stdint.h>

namespace Pds {

  class DetInfo;

  namespace Camera {
    class FrameCoord;
  };

  namespace Quartz {

    class ConfigV1 {
    public:
      enum { Version=1 };
      enum Depth     { Eight_bit, Ten_bit };
      enum Binning   { x1, x2, x4 };
      enum Mirroring { None, HFlip, VFlip, HVFlip };
      enum { LUT_Size=4096 };
      enum { Row_Pixels=2048 };
      enum { Column_Pixels=2048 };

      ConfigV1();
      ConfigV1(unsigned short  black,
	       unsigned short  gain,   // 100-800
	       Depth           depth,
	       Binning         vbinning,
	       Binning         hbinning,
	       Mirroring       mirroring,
	       bool            enable_pixel_correction );

      ConfigV1(unsigned short   black,
	       unsigned short   gain,
	       Depth            depth,
	       Binning          vbinning,
	       Binning          hbinning,
	       Mirroring        mirroring,
	       bool             enable_pixel_correction,
	       const uint16_t*  lut);

      //
      //  Accessors
      //

      //  offset/pedestal setting for camera (before gain)
      unsigned short   black_level() const;

      //  camera gain setting in percentile [100-800] = [1x-8x]
      unsigned short   gain_percent() const;

      //  offset/pedestal value in pixel counts
      unsigned short   output_offset() const;

      //  bit-depth of pixel counts
      Depth            output_resolution() const;
      unsigned         output_resolution_bits() const;

      //  horizontal re-binning of output (consecutive columns summed)
      Binning          horizontal_binning() const;

      //  vertical re-binning of output (consecutive rows summed)
      Binning          vertical_binning() const;

      //  geometric transformation of the image
      Mirroring        output_mirroring() const;

      //  correct defective pixels internally
      bool             defect_pixel_correction_enabled() const;

      //  apply output lookup table corrections
      bool             output_lookup_table_enabled() const;
      //  location of output lookup table: output_value[input_value]
      //  (appended to the end of this structure)
      const
      uint16_t*        output_lookup_table() const;

      //  defective pixel count
      unsigned         number_of_defect_pixels() const;
      void             set_number_of_defect_pixels(unsigned);

      //  location of defective pixel coordinates appended to
      //  the end of this structure
      const 
      Camera::FrameCoord* defect_pixel_coordinates() const;
      Camera::FrameCoord* defect_pixel_coordinates();

      //  total size of this structure 
      //  (including defective pixel coords and output lookup table)
      unsigned         size() const;

      static unsigned max_row_pixels   (const DetInfo&);
      static unsigned max_column_pixels(const DetInfo&);

    private:
      uint32_t _offsetAndGain; // offset and gain
      uint32_t _outputOptions; // bit mask of output formatting options
      uint32_t _defectPixelCount;
    };

  };
};

#endif
