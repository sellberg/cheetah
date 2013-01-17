#include "pdsdata/compress/Timepix_DataV2.hh"
#include "pdsdata/timepix/DataV2.hh"

#include <boost/shared_ptr.hpp>
#include <new>

using namespace Pds::Timepix;

CompressedDataV2::CompressedDataV2(const DataV2& d) :
  _data(d)
{
}

unsigned                CompressedDataV2::width       () const { return _data.width(); }
unsigned                CompressedDataV2::height      () const { return _data.height(); }
unsigned                CompressedDataV2::timestamp   () const { return _data.timestamp(); }
unsigned                CompressedDataV2::frameCounter() const { return _data.frameCounter(); }
unsigned                CompressedDataV2::lostRows    () const { return _data.lostRows(); }

const Pds::CompressedPayload& CompressedDataV2::pd() const { return _pd; }

static void Destroy(Pds::Timepix::DataV2* p) { delete[](char*)p; }

boost::shared_ptr<DataV2> CompressedDataV2::uncompressed() const
{
  char* p = new char[sizeof(_data)+_pd.dsize()];
  DataV2* v = new(p) DataV2(_data);
  if (!_pd.uncompress(p+sizeof(_data))) {
    delete[] p;
    v = 0;
  }
  boost::shared_ptr<DataV2> q(v,Destroy);
  return q;
}
