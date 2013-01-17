#include "pdsdata/compress/Timepix_DataV1.hh"
#include "pdsdata/timepix/DataV1.hh"

#include <boost/shared_ptr.hpp>
#include <new>

Pds::Timepix::CompressedDataV1::CompressedDataV1(const Pds::Timepix::DataV1& c) :
  _data(c)
{
}

uint32_t Pds::Timepix::CompressedDataV1::timestamp(void) const
{
  return (_data.timestamp());
}

uint16_t Pds::Timepix::CompressedDataV1::frameCounter(void) const
{
  return (_data.frameCounter());
}

uint16_t Pds::Timepix::CompressedDataV1::lostRows(void) const
{
  return (_data.lostRows());
}

const Pds::CompressedPayload& Pds::Timepix::CompressedDataV1::pd() const
{
  return _pd;
}

static void Destroy(Pds::Timepix::DataV1* p) { delete[](char*)p; }

boost::shared_ptr<Pds::Timepix::DataV1> Pds::Timepix::CompressedDataV1::uncompressed() const
{
  char* p = new char[sizeof(_data)+_pd.dsize()];
  Pds::Timepix::DataV1* v = new(p) Pds::Timepix::DataV1(_data);
  if (!_pd.uncompress(p+sizeof(_data))) {
    delete[] p;
    v = 0;
  }
  boost::shared_ptr<Pds::Timepix::DataV1> q(v,Destroy);
  return q;
}
