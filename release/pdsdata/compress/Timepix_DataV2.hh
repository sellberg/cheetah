#ifndef Pds_TimePix_CompressedDataV2_hh
#define Pds_TimePix_CompressedDataV2_hh

#include "pdsdata/timepix/DataV2.hh"
#include "pdsdata/compress/CompressedPayload.hh"

#include <stdint.h>

namespace boost { template<class T> class shared_ptr; };

namespace Pds {
  namespace Timepix {
    class CompressedDataV2 {
    public:
      CompressedDataV2() {}
      CompressedDataV2(const DataV2&);      
    public:
      unsigned                width       () const;
      unsigned                height      () const;
      unsigned                timestamp   () const;
      unsigned                frameCounter() const;
      unsigned                lostRows    () const;
    public:
      const CompressedPayload& pd() const;
      boost::shared_ptr<DataV2> uncompressed() const;
    private:
      DataV2            _data;
      CompressedPayload _pd;
    };
  };
};

#endif
