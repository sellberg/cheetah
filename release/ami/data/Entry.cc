#include <string.h>
#include <sys/uio.h>
#include <math.h>
#include <new>

#include "ami/data/Entry.hh"
#include "ami/data/DescEntry.hh"
#include "pdsdata/xtc/ClockTime.hh"

using namespace Ami;

Entry::Entry() : 
  _payloadsize(0),
  _payload(0)
{}

Entry::~Entry()
{
  if (_payload)
    delete [] _payload;  
}

void Entry::payload(iovec& iov)
{
  iov.iov_base = _payload;
  iov.iov_len = _payloadsize;
}

void Entry::payload(iovec& iov) const
{
  iov.iov_base = static_cast<void*>(_payload);
  iov.iov_len = _payloadsize;
}

void Entry::reset()
{
  memset(_payload, 0, _payloadsize);
}

void* Entry::allocate(unsigned size)
{
  if (_payload)
    delete [] _payload;  

  _payloadsize = sizeof(unsigned long long)+size;
  _payload = new unsigned long long[(_payloadsize>>3)+1];

  reset();

  return (_payload+1);
}

double Entry::last() const 
{
  const Pds::ClockTime& t = time();
  return t.seconds() + 1.e-9 * t.nanoseconds();
}

const Pds::ClockTime& Entry::time() const 
{
  return *reinterpret_cast<const Pds::ClockTime*>(_payload);
}

static const Pds::ClockTime _invalid(0,1);

void Entry::valid(const Pds::ClockTime& t) 
{
  *_payload = *(reinterpret_cast<const unsigned long long*>(&t)) & ~1ULL;
}

void Entry::valid(double v)
{
  new ((char*)_payload) Pds::ClockTime(unsigned(v),unsigned(drem(v,1)*1.e9)&~1);
}

void Entry::invalid() 
{ 
  *_payload |= 1ULL;
}

bool Entry::valid() const { return ((*_payload)&1)==0; }

void Entry::merge(char* p) const
{
  uint64_t* u = reinterpret_cast<uint64_t*>(p);
  if (*u > *_payload) {
    if (desc().aggregate()) {
      *u = (*u) | (*_payload & 1ULL);
      _merge((char*)(u+1));
    }
  }
  else {
    if (desc().aggregate()) {
      *u = (*_payload) | (*u & 1ULL);
      _merge((char*)(u+1));
    }
    else {
      memcpy(p,_payload,_payloadsize);
    }
  }
}
