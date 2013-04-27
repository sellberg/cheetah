#include "pdsdata/compress/Camera_FrameV1.hh"

#include <boost/shared_ptr.hpp>
#include <new>

using namespace Pds::Camera;

CompressedFrameV1::CompressedFrameV1(const FrameV1& f) :
  _frame(f.width(),
         f.height(),
         f.depth(),
         f.offset())
{
}

unsigned                 CompressedFrameV1::width () const { return _frame.width (); }
unsigned                 CompressedFrameV1::height() const { return _frame.height(); }
unsigned                 CompressedFrameV1::depth () const { return _frame.depth (); }
unsigned                 CompressedFrameV1::offset() const { return _frame.offset(); }

const Pds::CompressedPayload& CompressedFrameV1::pd    () const { return _pd; }

static void Destroy(FrameV1* p) { delete[](char*)p; }

boost::shared_ptr<FrameV1> CompressedFrameV1::uncompressed() const
{
  char* p = new char[sizeof(_frame)+_pd.dsize()];
  FrameV1* v = new(p) FrameV1(_frame);
  if (!_pd.uncompress(p + sizeof(_frame))) {
    delete[] p;
    v = 0;
  }
  boost::shared_ptr<FrameV1> q(v,Destroy);
  return q;
}
