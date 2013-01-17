#ifndef PdsUsdUsb_ConfigV1_hh
#define PdsUsdUsb_ConfigV1_hh

#include <stdint.h>

namespace Pds {
  namespace UsdUsb {
    class ConfigV1 {
    public:
      enum { Version = 1 };
      enum { NCHANNELS = 4 };
      enum Count_Mode {
	WRAP_FULL,
	LIMIT,
	HALT,
	WRAP_PRESET,
      };
      enum Quad_Mode {
	CLOCK_DIR,
	X1,
	X2,
	X4,
      };

      ConfigV1() {}
      ConfigV1( Count_Mode  cm[NCHANNELS],
		Quad_Mode   qm[NCHANNELS] );
      ~ConfigV1() {}

      void dump() const;

      Count_Mode counting_mode   (unsigned channel) const;
      Quad_Mode  quadrature_mode (unsigned channel) const;

    public:
      static const char** count_mode_labels();
      static const char** quad_mode_labels ();

    private:
      uint32_t _count_mode[NCHANNELS];
      uint32_t _quad_mode [NCHANNELS];
    };
  };
};

#endif
