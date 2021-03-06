#include "SummaryClient.hh"

#include "ami/qt/Control.hh"
#include "ami/qt/Status.hh"
#include "ami/qt/QtTH1F.hh"
#include "ami/qt/QtTH2F.hh"
#include "ami/qt/QtProf.hh"
#include "ami/qt/QtChart.hh"
#include "ami/qt/QtScan.hh"
#include "ami/qt/QtImage.hh"
#include "ami/qt/QtPlot.hh"
#include "ami/qt/ImageDisplay.hh"
#include "ami/qt/QtPWidget.hh"
#include "ami/qt/PWidgetManager.hh"

#include "ami/client/ClientManager.hh"

#include "ami/data/AbsTransform.hh"
#include "ami/data/ConfigureRequest.hh"
#include "ami/data/Discovery.hh"
#include "ami/data/DescEntry.hh"
#include "ami/data/EntryTH1F.hh"
#include "ami/data/EntryTH2F.hh"
#include "ami/data/EntryScalar.hh"
#include "ami/data/EntryProf.hh"
#include "ami/data/EntryScan.hh"
#include "ami/data/EntryImage.hh"
#include "ami/data/EntryFactory.hh"
#include "ami/data/RawFilter.hh"

#include "ami/service/Socket.hh"
#include "ami/service/Semaphore.hh"

#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QTabWidget>

#include <sys/types.h>
#include <sys/socket.h>

#include <list>

static const QChar PageIndex('#');

namespace Ami {
  namespace Qt {
    class QtBasePlot : public QtPlot {
    public:
      QtBasePlot(QtBase* b) : QtPlot(NULL,b->title()) 
      {
        PWidgetManager::remove(this);
        _base.push_back(b); 
        b->attach(_frame); 
      }
      virtual ~QtBasePlot() 
      {
        for(std::list<QtBase*>::iterator it=_base.begin(); it!=_base.end(); it++) {
          QtBase* b = *it;
          delete b;
        }
      }
    public:
      void add   (QtBase* b)
      {
        _base.push_back(b);
        b->attach(_frame);
      }
      void update() 
      { 
        for(std::list<QtBase*>::iterator it=_base.begin(); it!=_base.end(); it++)
          (*it)->update(); 
        emit redraw(); 
      }
      void dump(FILE* f) const 
      {
        for(std::list<QtBase*>::const_iterator it=_base.begin(); it!=_base.end(); it++)
          (*it)->dump(f); 
      }
    private:
      std::list<QtBase*> _base; 
    };

    class PagePlot : public QtPWidget {
    public:
      PagePlot(const QString& title) : QtPWidget(0), _title (title),_layout(new QGridLayout) {}
      ~PagePlot() {}
    public:
      const QString& title() const { return _title; }
      void add   (QtBase* plot, int row, int col)
      {
        if (row<0 || col<0) {
          row = _layout->rowCount();
          col = 0;
        }
        { int trow, tcol, trows, tcols;
          const unsigned count(_layout->count());
          for(unsigned i=0; i<count; i++) {
            _layout->getItemPosition(i, &trow, &tcol, &trows, &tcols);
            if (trow==row &&
                tcol==col) {
              QtBasePlot* frame = reinterpret_cast<QtBasePlot*>(_layout->itemAt(i)->widget());
              frame->add(plot);
              return;
            }
          }
        }
	QtBasePlot* frame = new QtBasePlot(plot);
	_layout->addWidget(frame, row, col);
	_plots.push_back(frame);
      }
      void add   (QtImage* plot, int row, int col)
      {
        ImageDisplay* frame = new ImageDisplay;
        frame->add(plot,true);
        if (row<0 || col<0) {
          row = _layout->rowCount();
          col = 0;
        }
        _layout->addWidget(frame, row, col);
        _images.push_back(frame);
      }
      void layout() { setLayout(_layout); }
      void update() 
      {
	for(std::list<QtBasePlot*>::iterator it = _plots.begin(); it!=_plots.end(); it++)
	  (*it)->update(); 
	for(std::list<ImageDisplay*>::iterator it = _images.begin(); it!=_images.end(); it++)
	  (*it)->update(); 
      }
      void save(char*& p) const
      {
	for(std::list<QtBasePlot*>::const_iterator it = _plots.begin(); it!=_plots.end(); it++)
	  XML_insert(p, "QtBasePlot", "_plots", (*it)->save(p));
	for(std::list<ImageDisplay*>::const_iterator it = _images.begin(); it!=_images.end(); it++)
	  XML_insert(p, "ImageDisplay", "_images", (*it)->save(p));
      }
      void load(const char*& p) 
      {
	std::list<QtBasePlot*  >::iterator it1 = _plots.begin();
	std::list<ImageDisplay*>::iterator it2 = _images.begin();
	XML_iterate_open(p,tag)
	  if (tag.element == "QtBasePlot")
	    (*it1++)->load(p);
	  else if (tag.element == "ImageDisplay")
	    (*it2++)->load(p);
	XML_iterate_close(PagePlot,tag);
      }
    private:
      QString                  _title;
      QGridLayout*             _layout;
      std::list<QtBasePlot*>   _plots;
      std::list<ImageDisplay*> _images;
    };
  };
};

using namespace Ami::Qt;

static const int BufferSize = 0x8000;
static Ami::AbsTransform& noTransform = Ami::AbsTransform::null();

SummaryClient::SummaryClient(QWidget* parent, const Pds::DetInfo& info, unsigned channel,
			     const QString& title, ConfigureRequest::Source source) :
  AbsClient        (parent, info, channel),
  _title           (title),
  _source          (source),
  _request         (new char[BufferSize]),
  _description     (new char[BufferSize]),
  _cds             ("Client"),
  _manager         (0),
  _niovload        (5),
  _niovread        (5),
  _iovload         (new iovec[_niovload]),
  _sem             (new Semaphore(Semaphore::EMPTY))
{
  PWidgetManager::add(this, _title);

  setWindowTitle(_title);
  setAttribute(::Qt::WA_DeleteOnClose, false);

  _control = new Control(*this,1);
  _status  = new Status;

  QVBoxLayout* layout = new QVBoxLayout;
  { QHBoxLayout* layout1 = new QHBoxLayout;
    layout1->addWidget(_control);
    layout1->addStretch();
    layout1->addWidget(_status);
    layout->addLayout(layout1); }
  layout->addWidget(_tab = new QTabWidget);
  setLayout(layout);

  connect(this, SIGNAL(description_changed(int)), this, SLOT(_read_description(int)));
}

SummaryClient::~SummaryClient() { PWidgetManager::remove(this); }

void SummaryClient::save(char*& p) const
{
  XML_insert(p, "QtPWidget", "self", QtPWidget::save(p) );
  XML_insert(p, "Control", "_control", _control->save(p) );
#if 0
  for(int i=0; i<_tab->count(); i++) {
    XML_insert( p, "PagePlot", "_tab",
		static_cast<PagePlot*>(_tab->widget(i))->save(p) );
  }
#endif
}

void SummaryClient::load(const char*& p)
{
  int i=0;
  XML_iterate_open(p,tag)
    if (tag.element == "QtPWidget")
      QtPWidget::load(p);
    else if (tag.name == "_control")
      _control->load(p);
#if 0
    else if (tag.element == "PagePlot") {
      static_cast<PagePlot*>(_tab->widget(i++))->load(p);
    }
#endif
  XML_iterate_close(SummaryClient,tag);
}

void SummaryClient::connected()
{
  _status->set_state(Status::Connected);
  _manager->discover();
}

void SummaryClient::discovered(const DiscoveryRx& rx)
{
  _status->set_state(Status::Discovered);
  printf("%s Discovered\n",qPrintable(_title));

  _manager->configure();
}

int  SummaryClient::configure       (iovec* iov) 
{
  _status->set_state(Status::Configured);
  printf("%s Configure\n",qPrintable(_title));

  char* p = _request;
  ConfigureRequest& r = *new (p) ConfigureRequest(ConfigureRequest::Create, _source, info.phy());
  p += r.size();

  iov[0].iov_base = _request;
  iov[0].iov_len  = p - _request;
  return 1;
}

int  SummaryClient::configured      () 
{
  printf("%s Configured\n",qPrintable(_title));
  return 0; 
}

int  SummaryClient::read_description(Socket& socket,int len)
{
  printf("%s Described so\n",qPrintable(_title));
  int size = socket.read(_description,len);

  if (size<0) {
    printf("Read error in Ami::Qt::Client::read_description.\n");
    return 0;
  }

  if (size==BufferSize) {
    printf("Buffer overflow in Ami::Qt::Client::read_description.  Dying...\n");
    abort();
  }

  //  printf("emit description\n");
  emit description_changed(size);

  _sem->take();

  return size;
}

void SummaryClient::_read_description(int size)
{
  printf("%s Described si\n",qPrintable(_title));

#if 0
  char* cache = new char[32*1024];
  { char* cache_p = cache;
    save(cache_p);
    sprintf(cache_p,"</SummaryClient>");
    printf("Saved %d bytes\n", cache_p-cache); 
    printf("%s\n",cache); }
#endif
    
  while(_tab->count()) {
    QWidget* w = _tab->widget(0);
    _tab->removeTab(0);
    delete w;
  }

  _cds.reset();

  std::list<PagePlot*> pages;

  const char* payload = _description;
  const char* const end = payload + size;

  payload += sizeof(Desc);

  while( payload < end ) {
    const DescEntry* desc = reinterpret_cast<const DescEntry*>(payload);
    if (desc->size()==0) {
      printf("read_description size==0\n");
      break;
    }
    Entry* entry = EntryFactory::entry(*desc);
    printf("%s found entry %s of type %d\n",
           qPrintable(_title),
	   desc->name(), desc->type());

    QString page_title, plot_title;
    int row=-1, col=-1;
    int rgb=0;
    { QString title(desc->name());
      int index = title.indexOf( PageIndex );
      if (index >= 0) {
	plot_title = title.left(index);
        int rIndex = title.indexOf( PageIndex, index+1 );
        if (rIndex >= 0) {
          page_title = title.mid(index+1,rIndex-index-1);
          int cIndex = title.indexOf( PageIndex, rIndex+1 );
          if (cIndex >= 0) {
            row = title.mid( rIndex+1, cIndex-rIndex-1).toInt();
            int rgbIndex = title.indexOf( PageIndex, cIndex+1 );
            if (rgbIndex >= 0) {
              col = title.mid( cIndex+1, rgbIndex-cIndex-1).toInt();
              rgb = title.mid( rgbIndex+1 ).toInt(0,16);
            }
            else
              col = title.mid( cIndex+1 ).toInt();
          }
          else
            row = title.mid( rIndex+1 ).toInt();
        }
        else
          page_title = title.mid(index+1);
      }
      else {
	plot_title = page_title = title;
      }
    }

    QtBase* plot = 0;
    QtImage* img = 0;
    switch(desc->type()) {
    case Ami::DescEntry::TH1F: 
      plot = new QtTH1F(plot_title,*static_cast<const Ami::EntryTH1F*>(entry),
			noTransform,noTransform,QColor((rgb>>16)&0xff,(rgb>>8)&0xff,(rgb>>0)&0xff));
      break;
    case Ami::DescEntry::Scalar:  // create a chart from a scalar
      plot = new QtChart(plot_title,*static_cast<const Ami::EntryScalar*>(entry),
                         QColor((rgb>>16)&0xff,(rgb>>8)&0xff,(rgb>>0)&0xff));
      break;
    case Ami::DescEntry::Prof: 
      plot = new QtProf(plot_title,*static_cast<const Ami::EntryProf*>(entry),
			noTransform,noTransform,QColor((rgb>>16)&0xff,(rgb>>8)&0xff,(rgb>>0)&0xff));
      break;
    case Ami::DescEntry::Scan: 
      plot = new QtScan(plot_title,*static_cast<const Ami::EntryScan*>(entry),
			noTransform,noTransform,QColor((rgb>>16)&0xff,(rgb>>8)&0xff,(rgb>>0)&0xff));
      break;
    case Ami::DescEntry::Image:
      img  = new QtImage(plot_title,*static_cast<const Ami::EntryImage*>(entry),
                         noTransform,noTransform,QColor(0,0,0));
      break;
    case Ami::DescEntry::TH2F:
      plot = new QtTH2F(plot_title,*static_cast<const Ami::EntryTH2F*>(entry),
                        noTransform,noTransform,QColor(0,0,0));
      break;
    default:
      printf("%s type %d not implemented yet\n",qPrintable(_title),desc->type()); 
      return;
    }

    bool lFound=false;
    for(std::list<PagePlot*>::iterator it = pages.begin(); it != pages.end(); it++) {
      if ((*it)->title() == page_title) {
        if (plot)
          (*it)->add(plot,row,col);
        else
          (*it)->add(img ,row,col);
        lFound=true;
        break;
      }
    }
    if (!lFound) {
      PagePlot* page = new PagePlot(page_title);
      if (plot)
        page->add(plot,row,col);
      else
        page->add(img ,row,col);
      pages.push_back(page);
    }

    _cds.add(entry, desc->signature());
    payload += desc->size();
  }
  if (_cds.totalentries()>_niovload) {
    delete[] _iovload;
    _iovload = new iovec[_niovload=_cds.totalentries()];
  }
  _niovread = _cds.payload(_iovload, _cds.request());

  for(std::list<PagePlot*>::iterator it = pages.begin(); it != pages.end(); it++) {
    (*it)->layout();
    _tab->addTab(*it, (*it)->title());
  }

#if 0
  { const char* cache_p = cache;
    load(cache_p); 
    printf("Loaded %d bytes\n",cache_p-cache); }
#endif

  _status->set_state(Status::Described);

  _sem->give();
}

int  SummaryClient::read_payload     (Socket& socket,int)
{
  int nbytes = 0;
  if (_niovread==0)
    ;
  else if (_status->state() == Status::Requested) {
    nbytes = socket.readv(_iovload,_niovread);
  }
  else {
    printf("Ami::Qt::Client::read_payload unrequested payload\n");
  }
  _status->set_state(Status::Received);
  return nbytes;
}

bool SummaryClient::svc             () const { return false; }

void SummaryClient::process         () 
{
//   QWidget* w = _tab->currentWidget();
//   if (w)
//     static_cast<PagePlot*>(w)->update();

  for(int i=0; i<_tab->count(); i++)
    static_cast<PagePlot*>(_tab->widget(i))->update();

  _status->set_state(Status::Processed);
}

void SummaryClient::managed(ClientManager& mgr)
{
  _manager = &mgr;
  show();
  printf("%s connecting\n",qPrintable(_title));
  _manager->connect();
}

void SummaryClient::request_payload()
{
  if (_status->state() >= Status::Described) {
    _cds.request(isVisible() ? Cds::All : Cds::None);
    _niovread = _cds.payload(_iovload, _cds.request());
    _manager->request_payload(_cds.request());
    _status->set_state(Status::Requested);
  }
}

const QString& SummaryClient::title() const { return _title; }

void SummaryClient::save_plots(const QString&) const {}
void SummaryClient::reset_plots() { if (_manager) _manager->configure(); }
