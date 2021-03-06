 #ifndef AmiQt_QtTH1F_hh
#define AmiQt_QtTH1F_hh

#include "ami/qt/QtBase.hh"
#include "ami/qt/QtPlotCurve.hh"

class QwtPlot;
class QColor;

namespace Ami {
  class EntryTH1F;
  class AbsTransform;
  namespace Qt {
    class QtTH1F : public QtBase {
    public:
      QtTH1F(const QString&   title,
	     const Ami::EntryTH1F&,
	     const AbsTransform& x,
	     const AbsTransform& y,
	     const QColor&);
      ~QtTH1F();
    public:
      void        dump  (FILE*   ) const;
      void        attach(QwtPlot*);
      void        update()        ;
      void        xscale_update() ;
      void        yscale_update() ;
      void        set_color(unsigned);
      const AxisInfo* xinfo() const;
      double      normalization() const;
    private:
      const AbsTransform&     _xscale;
      const AbsTransform&     _yscale;
      QtPlotCurve      _curve;
      double*          _x;
      double*          _y;
      AxisInfo*     _xinfo;
    };
  };
};

#endif
