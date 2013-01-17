#ifndef Pds_CompressedElementV2_hh
#define Pds_CompressedElementV2_hh

#include "pdsdata/cspad/ElementHeader.hh"

#include "pdsdata/compress/CompressedPayload.hh"

namespace boost { template<class T> class shared_ptr; };

namespace Pds {
  namespace CsPad {
    class ElementV2;
    class CompressedElementV2 : public ElementHeader {
    public:
      CompressedElementV2();
      CompressedElementV2(const ElementV2&);
    public:
      const CompressedPayload& pd() const;
      boost::shared_ptr<ElementV2> uncompressed() const;
   private:
      CompressedPayload _pd;
    };
  };
};
#endif
