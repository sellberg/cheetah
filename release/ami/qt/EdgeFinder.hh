#ifndef AmiQt_EdgeFinder_hh
#define AmiQt_EdgeFinder_hh

#include "ami/qt/QtPWidget.hh"
#include "ami/qt/OverlayParent.hh"

#include <QtCore/QString>
#include <QtCore/QStringList>

class QLineEdit;
class QButtonGroup;
class QVBoxLayout;
class QComboBox;
class QCheckBox;

#include <list>

namespace Ami {

  class AbsOperator;
  class Cds;
  class DescEntry;
  class Entry;

  namespace Qt {
    class AxisInfo;
    class ChannelDefinition;
    class EdgeCursor;
    class EdgePlot;
    class EdgeOverlay;
    class DescTH1F;
    class DescTH2T;
    class WaveformDisplay;

    class EdgeFinder : public QtPWidget,
                       public OverlayParent {
      Q_OBJECT
    public:
      EdgeFinder(QWidget* parent, ChannelDefinition* channels[], unsigned nchannels, 
                 WaveformDisplay&, QtPWidget* =0);
      ~EdgeFinder();
    public:
      void save(char*& p) const;
      void load(const char*& p);
      void save_plots(const QString&) const;
    public:
      void configure(char*& p, unsigned input, unsigned& output,
		     ChannelDefinition* ch[], int* signatures, unsigned nchannels);
      void setup_payload(Cds&);
      void update();
      void initialize(const Ami::DescEntry&);
    public slots:
      void set_channel   (int); // set the source
      void plot          ();   // configure the plot
      void overlay       ();   // configure the plot
      void remove_plot   (QObject*);
    signals:
      void changed();
    private:
      ChannelDefinition** _channels;
      unsigned _nchannels;
      unsigned _channel;

      WaveformDisplay&  _frame;
      EdgeCursor* _baseline;
      EdgeCursor* _threshold;
 
      QLineEdit* _title;
      QButtonGroup* _vtbutton;
      DescTH1F*  _hist;
      DescTH2T*  _hist2d;
      QCheckBox *_leading, *_trailing;
      QLineEdit* _dead;

      std::list<EdgePlot*> _plots;

    public:
      void add_overlay   (DescEntry*, QtPlot*, SharedData*);
      void remove_overlay(QtOverlay*);
    private:
      std::list<EdgeOverlay*> _ovls;
    };
  };
};

#endif
