#ifndef Pds_CompressedElementV1_hh
#define Pds_CompressedElementV1_hh

#include "pdsdata/cspad/ElementHeader.hh"

#include "pdsdata/compress/CompressedPayload.hh"

namespace boost { template<class T> class shared_ptr; };

namespace Pds {
  namespace CsPad {
    class ElementV1;
    class CompressedElementV1 : public ElementHeader {
    public:
      CompressedElementV1();
      CompressedElementV1(const ElementV1&);
    public:
      const CompressedPayload& pd() const;
      boost::shared_ptr<ElementV1> uncompressed() const;
    private:
      CompressedPayload _pd;
    };
  };
};
#endif
