#include "pdsdata/compress/CompressedPayload.hh"

using namespace Pds;

CompressedPayload::CompressedPayload() {}

CompressedPayload::CompressedPayload(Engine   e,
                                     unsigned d,
                                     unsigned c) :
  _engine  (e),
  _reserved(0),
  _dsize   (d),
  _csize   (c)
{
}

CompressedPayload::Engine CompressedPayload::compressor() const { return Engine(_engine); }

unsigned CompressedPayload::dsize() const { return _dsize; }

unsigned CompressedPayload::csize() const { return _csize; }

const void* CompressedPayload::cdata() const { return (const void*)(this+1); }

bool CompressedPayload::uncompress(void* outbuf) const
{
  return false;
}
