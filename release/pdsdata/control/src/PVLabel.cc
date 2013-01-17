#include "pdsdata/control/PVLabel.hh"

using namespace Pds::ControlData;

#include <string.h>

PVLabel::PVLabel() {}

PVLabel::PVLabel(const char* pvname, const char* value)
{
  strncpy(_name , pvname, NameSize);
  strncpy(_value, value , ValueSize);
}

PVLabel::PVLabel(const PVLabel& c)
{
  strncpy(_name , c._name , NameSize );
  strncpy(_value, c._value, ValueSize);
}

PVLabel::~PVLabel() {}

bool PVLabel::operator<(const PVLabel& m) const
{
  int nt = strncmp(_name, m._name, NameSize);
  return (nt<0);
}

const char* PVLabel::name() const { return _name; }

const char* PVLabel::value() const { return _value; }
