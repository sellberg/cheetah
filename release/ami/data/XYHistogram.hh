#ifndef Ami_XYHistogram_hh
#define Ami_XYHistogram_hh

//
//  class XYHistogram : an operator that generates a histogram of
//    an EntryImage pixel values.
//

#include "ami/data/AbsOperator.hh"
#include "ami/data/DescScalarRange.hh"

namespace Ami {

  class DescEntry;
  class Entry;

  class XYHistogram : public AbsOperator {
  public:
    XYHistogram(const DescEntry& output);
    XYHistogram(const char*&, const DescEntry&);
    XYHistogram(const char*&);
    ~XYHistogram();
  private:
    DescEntry& _routput   () const;
    Entry&     _operate  (const Entry&) const;
    void*      _serialize(void*) const;
    bool       _valid    () const { return true; }
  private:
    enum { DESC_LEN = sizeof(DescScalarRange) };
    char             _desc_buffer[DESC_LEN];
    Entry*           _output;
  };

};

#endif
