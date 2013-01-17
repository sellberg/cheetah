#ifndef Pds_CompressedMiniElementV1_hh
#define Pds_CompressedMiniElementV1_hh

#include "pdsdata/cspad/ElementHeader.hh"

#include "pdsdata/compress/CompressedPayload.hh"

namespace boost { template<class T> class shared_ptr; };

namespace Pds {
  namespace CsPad {
    class MiniElementV1;
    class CompressedMiniElementV1 : public ElementHeader {
    public:
      CompressedMiniElementV1();
      CompressedMiniElementV1(const MiniElementV1&);
    public:
      const CompressedPayload& pd() const;
      boost::shared_ptr<MiniElementV1> uncompressed() const;
    private:
      CompressedPayload _pd;
    };
  };
};
#endif
