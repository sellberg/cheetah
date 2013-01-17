#ifndef Pds_CompressedElementV1_2x2_hh
#define Pds_CompressedElementV1_2x2_hh

#include "pdsdata/cspad2x2/ElementHeader.hh"

#include "pdsdata/compress/CompressedPayload.hh"

namespace boost { template<class T> class shared_ptr; };

namespace Pds {
  namespace CsPad2x2 {
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
