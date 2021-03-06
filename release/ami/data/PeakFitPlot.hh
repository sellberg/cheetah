#ifndef Ami_PeakFitPlot_hh
#define Ami_PeakFitPlot_hh

//
//  class PeakFitPlot : an operator that extracts peak parameters from a profile
//    and generates a mean value (Scalar), a distribution (TH1F), or a profile (Prof)
//    versus a BLD or PV quantity.
//

#include "ami/data/AbsOperator.hh"

namespace Ami {

  class DescEntry;
  class Entry;
  class Term;
  class FeatureCache;

  class PeakFitPlot : public AbsOperator {
  public:
    enum Parameter { Position, Height, FWHM, RMS, NumberOf };
    static const char* name(Parameter);
  public:
    //  Defined by the output entry's description,
    //    the BLD/PV, and any BLD/PV dependence for profiles.
    PeakFitPlot(const DescEntry& output,
		double    baseline,
		Parameter prm);
    PeakFitPlot(const DescEntry& output,
		int       nbins,
                int      *bins,
		Parameter prm);
    //  Reconstituted from the input serial stream, the BldState and PvState
    //    accessors, and the Cds input entry accessor.
    PeakFitPlot(const char*&, FeatureCache&);
    PeakFitPlot(const char*&);
    ~PeakFitPlot();
  public:
    unsigned   input    () const;
    Parameter  prm      () const;
    const char* feature() const;
  private:
    DescEntry& _routput   () const;
    Entry&     _operate  (const Entry&) const;
    void*      _serialize(void*) const;
    bool       _valid    () const { return _v; }
  private:
    enum { DESC_LEN = 1024, MAX_BINS = 10 };
    char             _desc_buffer[DESC_LEN];
    int              _nbins;
    int              _bins[MAX_BINS];
    double           _baseline;
    Parameter        _prm;

    FeatureCache* _cache;
    Term*         _term;
    Entry*        _entry;
    bool          _v;
  };

};

#endif
