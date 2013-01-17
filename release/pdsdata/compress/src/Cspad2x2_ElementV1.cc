#include "pdsdata/compress/Cspad2x2_ElementV1.hh"
#include "pdsdata/cspad2x2/ElementV1.hh"

#include <boost/shared_ptr.hpp>
#include <new>

using namespace Pds::CsPad2x2;

CompressedElementV1::CompressedElementV1() 
{
}

CompressedElementV1::CompressedElementV1(const ElementV1& o) :
  ElementHeader(o)
{
}

const Pds::CompressedPayload& CompressedElementV1::pd() const
{ return _pd; }

boost::shared_ptr<ElementV1> CompressedElementV1::uncompressed() const
{
  ElementV1* v = new ElementV1;
  char* p = reinterpret_cast<char*>(v);
  new(p) ElementHeader(*this);
  if (!_pd.uncompress(p+sizeof(ElementHeader))) {
    delete v;
    v = 0;
  }
  boost::shared_ptr<ElementV1> q(v);
  return q;
}
