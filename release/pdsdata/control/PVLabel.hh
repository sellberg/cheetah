//
//  Class for Process Variable String
//
#ifndef PdsData_PVLabel_hh
#define PdsData_PVLabel_hh

#include <stdint.h>

namespace Pds {

#pragma pack(4)

  namespace ControlData {

    class PVLabel {
    public:
      enum { NameSize=32 };
      enum { ValueSize=64 };
    public:
      PVLabel();
      PVLabel(const char* pvname, const char* value);
      PVLabel(const PVLabel&);
      ~PVLabel();
    public:
      bool operator<(const PVLabel&) const;
    public:
      const char* name () const;
      const char* value() const;
    private:
      char     _name [NameSize];
      char     _value[ValueSize];
    };

  };

#pragma pack()
};

#endif
