#ifndef Pds_CompressedFrameV1_hh
#define Pds_CompressedFrameV1_hh

#include "pdsdata/camera/FrameV1.hh"
#include "pdsdata/compress/CompressedPayload.hh"

#include <stdint.h>
#include <memory>

namespace boost { template<class T> class shared_ptr; };

namespace Pds {
  namespace Camera {
    class CompressedFrameV1 {
    public:
      enum { Version = 1 };
      CompressedFrameV1() {}
      CompressedFrameV1(const FrameV1&);
    public:
      unsigned                width () const;
      unsigned                height() const;
      unsigned                depth () const;
      unsigned                offset() const;
    public:
      const CompressedPayload& pd() const;
      boost::shared_ptr<FrameV1> uncompressed() const;
    private:
      FrameV1           _frame;
      CompressedPayload _pd;
    };
  };
};

#endif
