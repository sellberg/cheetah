#ifndef Pds_TimePix_CompressedDataV1_hh
#define Pds_TimePix_CompressedDataV1_hh

#include "pdsdata/timepix/DataV1.hh"
#include "pdsdata/compress/CompressedPayload.hh"

#include <stdint.h>
#include "pdsdata/xtc/TypeId.hh"

namespace boost { template<class T> class shared_ptr; };

namespace Pds
{
   namespace Timepix
   {
     class CompressedDataV1
     {
     public:
       enum { Version = 1 };

       CompressedDataV1() {}

       CompressedDataV1(const DataV1&);

       unsigned width() const { return _data.width(); }

       unsigned height() const { return _data.height(); }

       unsigned depth() const { return _data.depth(); }

       // number of bytes per pixel
       unsigned depth_bytes() const { return _data.depth_bytes(); }

       uint32_t timestamp(void) const;
       uint16_t frameCounter(void) const;
       uint16_t lostRows(void) const;
     public:
       const CompressedPayload& pd() const;
       boost::shared_ptr<DataV1> uncompressed() const;
    private:
       DataV1            _data;
       CompressedPayload _pd;
     };
   };
};

#endif
