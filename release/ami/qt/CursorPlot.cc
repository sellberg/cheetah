#include "CursorPlot.hh"

#include "ami/qt/AxisInfo.hh"
#include "ami/qt/ChannelDefinition.hh"
#include "ami/qt/Filter.hh"
#include "ami/qt/PlotFactory.hh"
#include "ami/qt/QtTH1F.hh"
#include "ami/qt/QtTH2F.hh"
#include "ami/qt/QtChart.hh"
#include "ami/qt/QtProf.hh"
#include "ami/qt/QtScan.hh"
#include "ami/qt/QtEmpty.hh"

#include "ami/data/BinMath.hh"
#include "ami/data/Cds.hh"
#include "ami/data/ConfigureRequest.hh"
#include "ami/data/AbsTransform.hh"
#include "ami/data/DescEntry.hh"
#include "ami/data/RawFilter.hh"

#include "ami/data/EntryTH1F.hh"
#include "ami/data/EntryTH2F.hh"
#include "ami/data/EntryProf.hh"
#include "ami/data/EntryScan.hh"
#include "ami/data/EntryScalar.hh"
#include "ami/data/EntryScalarRange.hh"
#include "ami/data/EntryScalarDRange.hh"

#include <QtGui/QLabel>
#include "qwt_plot.h"

namespace Ami {
  namespace Qt {
    class NullTransform : public Ami::AbsTransform {
    public:
      ~NullTransform() {}
      double operator()(double x) const { return x; }
    };
  };
};

using namespace Ami::Qt;

static NullTransform noTransform;


CursorPlot::CursorPlot(QWidget* parent,
		       const QString&   name,
		       unsigned         channel,
		       BinMath*         input) :
  QtPlot   (parent, name),
  _channel (channel),
  _input   (input),
  _output_signature  (0),
  _plot    (new QtEmpty),
  _auto_range(0)
{
  _plot->attach(_frame);
  setPlotType(_input->output().type());
}

CursorPlot::CursorPlot(QWidget* parent,
		       const char*& p) :
  QtPlot   (parent),
  _plot    (new QtEmpty),
  _auto_range(0)
{
  load(p);
}

CursorPlot::~CursorPlot()
{
  delete _input;
  if (_plot    ) delete _plot;
}

void CursorPlot::save(char*& p) const
{
  char* buff = new char[8*1024];
  XML_insert(p, "QtPlot", "self", QtPlot::save(p) );
  XML_insert(p, "int", "_channel", QtPersistent::insert(p,(int)_channel) );
  XML_insert(p, "BinMath", "_input", QtPersistent::insert(p,buff,(char*)_input->serialize(buff)-buff) );
  delete[] buff;
}

void CursorPlot::load(const char*& p)
{
  XML_iterate_open(p,tag)
    if (tag.element == "QtPlot")
      QtPlot::load(p);
    else if (tag.name == "_channel")
      _channel = QtPersistent::extract_i(p);
    else if (tag.name == "_input") {
      const char* b = (const char*)QtPersistent::extract_op(p);
      b += 2*sizeof(uint32_t);
      _input = new BinMath(b);
    }
  XML_iterate_close(CursorPlot,tag);

  _output_signature=0;
  _plot->attach(_frame);
  setPlotType(_input->output().type());
}

void CursorPlot::dump(FILE* f) const { _plot->dump(f); }

#include "ami/data/Entry.hh"
#include "ami/data/DescEntry.hh"

void CursorPlot::setup_payload(Cds& cds)
{
  Ami::Entry* entry = cds.entry(_output_signature);
  if (entry) {

    if (_plot && !_req.changed() && !_auto_range) {
      _plot->entry(*entry);
    }
    else {
      if (_plot)
        delete _plot;
    
      _auto_range = 0;

      edit_xrange(true);
      edit_yrange(true);

      switch(entry->desc().type()) {
      case Ami::DescEntry::TH1F: 
        _plot = new QtTH1F(_name,*static_cast<const Ami::EntryTH1F*>(entry),
                           noTransform,noTransform,QColor(0,0,0));
        break;
      case Ami::DescEntry::TH2F: 
        _plot = new QtTH2F(_name,*static_cast<const Ami::EntryTH2F*>(entry),
                           noTransform,noTransform,QColor(0,0,0));
        edit_xrange(false);
        edit_yrange(false);
        break;
      case Ami::DescEntry::Scalar:  // create a chart from a scalar
        _plot = new QtChart(_name,*static_cast<const Ami::EntryScalar*>(entry),
                            QColor(0,0,0));
        edit_xrange(false);
        break;
      case Ami::DescEntry::ScalarRange:
        _auto_range = static_cast<const Ami::EntryScalarRange*>(entry);
        _plot = new QtEmpty;
        return;
      case Ami::DescEntry::ScalarDRange:
        _auto_range = static_cast<const Ami::EntryScalarDRange*>(entry);
        _plot = new QtEmpty;
        return;
      case Ami::DescEntry::Prof: 
        _plot = new QtProf(_name,*static_cast<const Ami::EntryProf*>(entry),
                           noTransform,noTransform,QColor(0,0,0));
        break;
      case Ami::DescEntry::Scan: 
        _plot = new QtScan(_name,*static_cast<const Ami::EntryScan*>(entry),
                           noTransform,noTransform,QColor(0,0,0),
                           _style.symbol_size());
        break;
      default:
        printf("CursorPlot type %d not implemented yet\n",entry->desc().type()); 
        _plot = new QtEmpty;
        return;
      }
      _plot->attach(_frame);
    }
  }
  else {
    if (_output_signature>=0)
      printf("%s output_signature %d not found\n",qPrintable(_name),_output_signature);
    if (_plot) {
      delete _plot;
      _plot = new QtEmpty;
    }
  }
}

void CursorPlot::configure(char*& p, unsigned input, unsigned& output,
			   ChannelDefinition* channels[], int* signatures, unsigned nchannels,
			   const AxisInfo& xinfo, ConfigureRequest::Source source)
{
  unsigned channel = _channel;
  unsigned input_signature = signatures[channel];
  configure(p, input_signature, output, xinfo, source);
}

void CursorPlot::configure(char*& p, unsigned input, unsigned& output,
			   const AxisInfo& xinfo, ConfigureRequest::Source source)
{
  // replace cursor values with bin indices
  QString expr(_input->expression());
  QString new_expr;
  { QRegExp match("\\[[^\\]]*\\]");
    int last=0;
    int pos=0;
    while( (pos=match.indexIn(expr,pos)) != -1) {
      QString use = expr.mid(pos+1,match.matchedLength()-2);
      bool ok;
      double v = use.toDouble(&ok);
      unsigned bin=0;
      if (!ok)
	printf("error parsing double %s\n",qPrintable(use));
      else {
	bin = xinfo.tick(v);
      }
      new_expr.append(expr.mid(last,pos-last));
      new_expr.append(QString("[%1]").arg(bin));
      pos += match.matchedLength();
      last = pos;
    }
    new_expr.append(expr.mid(last));
    new_expr.replace(QString("]%1[").arg(BinMath::integrate()),QString(BinMath::integrate()));
    new_expr.replace(QString("]%1[").arg(BinMath::moment1  ()),QString(BinMath::moment1  ()));
    new_expr.replace(QString("]%1[").arg(BinMath::moment2  ()),QString(BinMath::moment2  ()));
    new_expr.replace(QString("]%1[").arg(BinMath::range    ()),QString(BinMath::range    ()));
    new_expr.replace(QString("]%1[").arg(BinMath::contrast ()),QString(BinMath::contrast ()));
    new_expr.replace(QString("]%1[").arg(BinMath::xmoment  ()),QString(BinMath::xmoment  ()));
    new_expr.replace(QString("]%1[").arg(BinMath::ymoment  ()),QString(BinMath::ymoment  ()));
    new_expr.replace(QString("]%1[").arg(BinMath::mean     ()),QString(BinMath::mean     ()));
    new_expr.replace(QString("]%1[").arg(BinMath::variance ()),QString(BinMath::variance ()));
  }
  QString end_expr;
  { int last=0, next=0, pos=0;
    while( (pos=new_expr.indexOf(BinMath::range(),pos)) != -1) {
      if ( (next=new_expr.lastIndexOf("[",pos))==-1 )
	printf("error parsing range in %s\n",qPrintable(expr));
      else {
	end_expr.append(new_expr.mid(last,next-last));
	last  = new_expr.indexOf("]",pos);
	int a = new_expr.mid(next+1,pos -next-1).toInt();
	int b = new_expr.mid(pos +1,last-pos -1).toInt();
	printf("%s/%d %s/%d\n",
	       qPrintable(new_expr.mid(next+1,pos -next-1)),a,
	       qPrintable(new_expr.mid(pos +1,last-pos -1)),b);
	end_expr.append(QString("(%1)").arg(QString::number(abs(a-b)+1)));
	pos  = ++last;
      }
    }
    end_expr.append(new_expr.mid(last));
  }

  Ami::BinMath op(_input->output(), qPrintable(end_expr));
  
  ConfigureRequest& r = *new (p) ConfigureRequest(ConfigureRequest::Create,
						  source,
						  input,
						  -1,
                                                  RawFilter(),
						  op);
  p += r.size();
  _req.request(r,output);
  _output_signature = r.output();
}

void CursorPlot::update()
{
  if (_plot) {
    _plot->update();
    emit counts_changed(_plot->normalization());
    emit redraw();
  }
  if (_auto_range) {
    double v = _auto_range->entries();
    emit counts_changed(v);
    if (v >= 0) {
      _auto_range->result(&_input->output());
      _auto_range = 0;
      emit changed();
    }
  }
}
