#include "pdsdata/compress/Cspad_ElementV1.hh"
#include "pdsdata/cspad/ElementV1.hh"

#include <boost/shared_ptr.hpp>
#include <new>

using namespace Pds::CsPad;

CompressedElementV1::CompressedElementV1() 
{
}

CompressedElementV1::CompressedElementV1(const ElementV1& o) :
  ElementHeader(o)
{
}

const Pds::CompressedPayload& CompressedElementV1::pd() const
{ return _pd; }

static void Destroy(ElementV1* p) { delete[](char*)p; }

//
//  This method only uncompresses one element(quad)
//
boost::shared_ptr<ElementV1> CompressedElementV1::uncompressed() const
{
  char* p = new char[sizeof(ElementV1)+_pd.dsize()+sizeof(uint32_t)];
  ElementV1* v = static_cast<ElementV1*>(new(p) ElementHeader(*this));
  if (!_pd.uncompress(p+sizeof(ElementV1))) {
    delete[] p;
    v = 0;
  }
  else {
    //  remember the trailer word
    *reinterpret_cast<uint32_t*>(p+sizeof(ElementV1)+_pd.dsize()) = 
      *reinterpret_cast<uint32_t*>((char*)(this+1)+_pd.csize());
  }
  boost::shared_ptr<ElementV1> q(v,Destroy);
  return q;
}
