#include "pdsdata/compress/Cspad_MiniElementV1.hh"
#include "pdsdata/cspad/MiniElementV1.hh"

#include <boost/shared_ptr.hpp>
#include <new>

using namespace Pds::CsPad;

CompressedMiniElementV1::CompressedMiniElementV1() 
{
}

CompressedMiniElementV1::CompressedMiniElementV1(const MiniElementV1& o) :
  ElementHeader(o)
{
}

const Pds::CompressedPayload& CompressedMiniElementV1::pd() const
{ return _pd; }

boost::shared_ptr<MiniElementV1> CompressedMiniElementV1::uncompressed() const
{
  MiniElementV1* v = new MiniElementV1;
  char* p = reinterpret_cast<char*>(v);
  new(p) ElementHeader(*this);
  if (!_pd.uncompress(p+sizeof(ElementHeader))) {
    delete v;
    v = 0;
  }
  boost::shared_ptr<MiniElementV1> q(v);
  return q;
}
