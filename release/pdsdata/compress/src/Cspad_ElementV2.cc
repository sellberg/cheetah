#include "pdsdata/compress/Cspad_ElementV2.hh"
#include "pdsdata/cspad/ElementV2.hh"

#include <boost/shared_ptr.hpp>
#include <new>

using namespace Pds::CsPad;

CompressedElementV2::CompressedElementV2() 
{
}

CompressedElementV2::CompressedElementV2(const ElementV2& o) :
  ElementHeader(o)
{
}

const Pds::CompressedPayload& CompressedElementV2::pd() const
{ return _pd; }

static void Destroy(ElementV2* p) { delete[](char*)p; }

//
//  This method only uncompresses one element(quad)
//
boost::shared_ptr<ElementV2> CompressedElementV2::uncompressed() const
{
  char* p = new char[sizeof(ElementV2)+_pd.dsize()+sizeof(uint32_t)];
  ElementV2* v = static_cast<ElementV2*>(new(p) ElementHeader(*this));
  if (!_pd.uncompress(p+sizeof(ElementV2))) {
    delete[] p;
    v = 0;
  }
  else {
    //  remember the trailer word
    *reinterpret_cast<uint32_t*>(p+sizeof(ElementV2)+_pd.dsize()) = 
      *reinterpret_cast<uint32_t*>((char*)(this+1)+_pd.csize());
  }
  boost::shared_ptr<ElementV2> q(v,Destroy);
  return q;
}
