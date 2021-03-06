#include "CspadMiniHandler.hh"
#include "CspadAlignment.hh"
#include "CspadAlignment_Commissioning.hh"

#include "ami/event/CspadTemp.hh"
#include "ami/event/CspadCalib.hh"
#include "ami/data/DescImage.hh"
#include "ami/data/EntryImage.hh"
#include "ami/data/ChannelID.hh"
#include "ami/data/FeatureCache.hh"
#include "ami/data/PeakFinder.hh"
#include "ami/data/PeakFinderFn.hh"

#include "pdsdata/cspad/MiniElementV1.hh"
#include "pdsdata/xtc/Xtc.hh"

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define UNBINNED
#define DO_PED_CORR
#define POST_INTEGRAL

using Ami::CspadCalib;

typedef Pds::CsPad::MiniElementV1 CspadElement;

static const unsigned Offset = 0x4000;
static const double pixel_size = 110e-6;

enum Rotation { D0, D90, D180, D270, NPHI=4 };

static void _transform(double& x,double& y,double dx,double dy,Rotation r)
{
  switch(r) {
  case D0  :    x += dx; y += dy; break;
  case D90 :    x += dy; y -= dx; break;
  case D180:    x -= dx; y -= dy; break;
  case D270:    x -= dy; y += dx; break;
  default:                        break;
  }
}

static inline unsigned sum1(const uint16_t*& data,
                            const uint16_t*& off,
                            const uint16_t* const*& psp,
                            const double fn,
                            const float*& gn)
{ 
  unsigned v;
  if (off==*psp) { psp++; v = Offset; }
  else { 
    double d = (*gn)*(double(*data + *off - Offset) - fn);
    d += Offset; 
    v = unsigned(d+0.5);
  }
  off++;
  gn++;
  return v;
}

static const unsigned no_threshold = 0x00ffffff;

static inline unsigned thr1(double v0, double v1,
                            const uint16_t*& off,
                            const uint16_t* const*& psp,
                            const float*& rms)
{ 
  unsigned v;
  if (off==*psp) { psp++; v = no_threshold; }
  else {  v = unsigned(v0 + *rms*v1 + 0.5); }
  off++;
  rms++;
  return v;
}


static double frameNoise(const uint16_t*  data,
                         const uint16_t*  off,
                         const uint16_t* const* sta)
{
  const unsigned ColBins = CsPad::ColumnsPerASIC;
  const unsigned RowBins = CsPad::MaxRowsPerASIC<<1;
  const int fnPixelMin = -100 + Offset;
  const int fnPixelMax =  100 + Offset;
  const int fnPixelBins = fnPixelMax - fnPixelMin;
  const int peakSpace   = 5;
  
  //  histogram the pixel values
  unsigned hist[fnPixelBins];
  { memset(hist, 0, fnPixelBins*sizeof(unsigned));
    const uint16_t* d(data);
    const uint16_t* o(off );
    for(unsigned i=0; i<ColBins; i++) {
      for(unsigned j=0; j<RowBins; j++, d+=2, o++) {
        if (*sta == o)
          sta++;
        else {
          int v = *d + *o - fnPixelMin;
          if (v >= 0 && v < int(fnPixelBins))
            hist[v]++;
        }
      }
    }
  }

  double v = 0;
  // the first peak from the left above this is the pedestal
  { const int fnPeakBins = 5;
    const int fnPixelRange = fnPixelBins-fnPeakBins-1;
    const unsigned fnPedestalThreshold = 1000;
    
    unsigned i=fnPeakBins;
    while( int(i)<fnPixelRange ) {
      if (hist[i]>fnPedestalThreshold) break;
      i++;
    }

    unsigned thresholdPeakBin=i;
    unsigned thresholdPeakBinContent=hist[i];
    while( int(++i)<fnPixelRange ) {
      if (hist[i]<thresholdPeakBinContent) {
        if (i > thresholdPeakBin+peakSpace)
          break;
      }
      else {
        thresholdPeakBin = i;
        thresholdPeakBinContent = hist[i];
      }
    }

    i = thresholdPeakBin;
    if ( int(i)+fnPeakBins<=fnPixelRange ) {
      unsigned s0 = 0;
      unsigned s1 = 0;
      for(unsigned j=i-fnPeakBins-1; j<i+fnPeakBins; j++) {
        s0 += hist[j];
        s1 += hist[j]*j;
      }
      
      double binMean = double(s1)/double(s0);
      v =  binMean + fnPixelMin - Offset;
      
      s0 = 0;
      unsigned s2 = 0;
      for(unsigned j=i-10; j<i+fnPeakBins; j++) {
        s0 += hist[j];
	s2 += hist[j]*(j-int(binMean))*(j-int(binMean));
      }
//      const double allowedPedestalWidthSquared = 2.5*2.5;
      //      printf("frameNoise finds mean %f, variance %f\n", v, double(s2)/double(s0));
//      if (double(s2)/double(s0)>allowedPedestalWidthSquared) v = 0;
      // this isn't the standard rms around the mean, but should be similar if rms_real < 3
      //      printf("frameNoise finds mean %f, variance %f\n", v, double(s2)/double(s0));

    }
    else {
      static unsigned nPrint=0;
      nPrint++;
      if ((nPrint<10) || (nPrint%100)==0)
        printf("frameNoise : peak not found [%d]\n",nPrint);
    }
//    printf("CspadMiniHandler::frameNoise v=%lf\n", v);
  }

  return v;
}

static FILE *fopen_dual(char *path1,char * path2, char *description)
{
  const int ErrMsgSize=200;
  char errmsg[ErrMsgSize];
  FILE *f = fopen(path1, "r");

  if (f) {
    printf("Loaded %s from %s\n", description, path1);
  } else {
    snprintf(errmsg, ErrMsgSize, "fopen: Failed to load %s from %s", description, path1);
    perror(errmsg);
    f = fopen(path2, "r");
    if (f) {
      printf("Loaded %s from %s\n", description, path2);
    } else {
      snprintf(errmsg, ErrMsgSize, "fopen: Failed to load %s from %s", description, path2);
      perror(errmsg);
    }
  }
  return (f);
}

namespace CspadMiniGeometry {

  //
  //  When filling the image, compensate data which
  //    only partially fills a pixel (at the edges)
  //
#define FRAME_BOUNDS 							\
  const unsigned ColLen   =   CsPad::ColumnsPerASIC/ppb-1;              \
    const unsigned RowLen = 2*CsPad::MaxRowsPerASIC/ppb-1;		\
    unsigned x0 = CALC_X(column,0,0);					\
    unsigned x1 = CALC_X(column,ColLen,RowLen);                         \
    unsigned y0 = CALC_Y(row,0,0);					\
    unsigned y1 = CALC_Y(row,ColLen,RowLen);				\
    if (x0 > x1) { unsigned t=x0; x0=x1; x1=t; }			\
    if (y0 > y1) { unsigned t=y0; y0=y1; y1=t; }			


#define BIN_ITER1(F1) {							\
    const unsigned ColBins = CsPad::ColumnsPerASIC;			\
    const unsigned RowBins = CsPad::MaxRowsPerASIC<<1;			\
    /*  fill the target region  */					\
    for(unsigned i=0; i<ColBins; i++) {					\
      for(unsigned j=0; j<RowBins; j++, data+=2) {                      \
	const unsigned x = CALC_X(column,i,j);				\
	const unsigned y = CALC_Y(row   ,i,j);				\
	image.content(F1,x,y);                                          \
      }									\
    }									\
  }

  //
  //  This class locates the ASIC data to the binned image grid
  //
  class Asic {
  public:
    Asic(double x, double y, unsigned ppbin) :
      column(unsigned( x/pixel_size)/ppbin),
      row   (unsigned(-y/pixel_size)/ppbin),
      ppb(ppbin) 
    {}
    virtual ~Asic() {}
  public:
    virtual void fill(Ami::DescImage& image) const = 0;
    virtual void fill(Ami::EntryImage& image,
		      double, double) const = 0;
    virtual void fill(Ami::EntryImage& image,
		      const uint16_t*  data) const = 0;
  public:
    virtual void boundary(unsigned& x0, unsigned& x1, 
			  unsigned& y0, unsigned& y1) const = 0;
  protected:
    unsigned column, row;
    unsigned ppb;
    uint16_t*  _sta[CsPad::MaxRowsPerASIC*CsPad::ColumnsPerASIC*2];
  };

#define AsicTemplate(classname,bi,ti,PPB)                               \
  class classname : public Asic {					\
  public:								\
    classname(double x, double y) : Asic(x,y,PPB) {}                    \
    void boundary(unsigned& dx0, unsigned& dx1,				\
		  unsigned& dy0, unsigned& dy1) const {			\
      FRAME_BOUNDS;							\
      dx0=x0; dx1=x1; dy0=y0; dy1=y1; }					\
    void fill(Ami::DescImage& image) const {				\
      FRAME_BOUNDS;							\
      image.add_frame(x0,y0,x1-x0+1,y1-y0+1);				\
    }									\
    void fill(Ami::EntryImage& image,                                   \
              double v0, double v1) const {                             \
      unsigned u0 = unsigned(v0);                                       \
      unsigned data = 0;                                                \
      ti                                                                \
    }									\
    void fill(Ami::EntryImage& image,					\
	      const uint16_t*  data) const { bi }                       \
  }

#define B1 { BIN_ITER1((*data)); }
#define T1 { BIN_ITER1(u0); }

#define CALC_X(a,b,c) (a+b)			    
#define CALC_Y(a,b,c) (a-c)			     
  AsicTemplate(  AsicD0B1, B1, T1, 1);
#undef CALC_X
#undef CALC_Y
#define CALC_X(a,b,c) (a+c)			    
#define CALC_Y(a,b,c) (a+b)			     
  AsicTemplate( AsicD90B1, B1, T1, 1);
#undef CALC_X
#undef CALC_Y
#define CALC_X(a,b,c) (a-b)			    
#define CALC_Y(a,b,c) (a+c)			     
  AsicTemplate(AsicD180B1, B1, T1, 1);
#undef CALC_X
#undef CALC_Y
#define CALC_X(a,b,c) (a-c)			    
#define CALC_Y(a,b,c) (a-b)			     
  AsicTemplate(AsicD270B1, B1, T1, 1);
#undef CALC_X
#undef CALC_Y

#undef B1
#undef T1
#undef AsicTemplate

  class AsicP : public Asic {
  public:
    AsicP(double x, double y, unsigned ppbin, 
          FILE* ped, FILE* status, FILE* gain, FILE* sigma) :
      Asic(x,y,ppbin)
    { // load offset-pedestal 
      size_t sz = 8 * 1024;
      char* linep = (char *)malloc(sz);
      char* pEnd;

      if (ped) {
        uint16_t* off = _off;
        for(unsigned col=0; col<CsPad::ColumnsPerASIC; col++) {
          getline(&linep, &sz, ped);
          *off++ = Offset - uint16_t(strtod(linep,&pEnd));
          for (unsigned row=1; row < 2*Pds::CsPad::MaxRowsPerASIC; row++)
            *off++ = Offset - uint16_t(strtod(pEnd, &pEnd));
        }
      }
      else
        memset(_off,0,sizeof(_off));

      if (status) {
        uint16_t*  off = _off;
        uint16_t** sta = _sta;
        for(unsigned col=0; col<CsPad::ColumnsPerASIC; col++) {
          getline(&linep, &sz, status);
          if (strtoul(linep,&pEnd,0)) *sta++ = off;
          off++;
          for (unsigned row=1; row < 2*Pds::CsPad::MaxRowsPerASIC; row++, off++)
            if (strtoul(pEnd,&pEnd,0)) *sta++ = off;
        }
      }
      else
        _sta[0] = 0;

      if (gain) {
        float* gn = _gn;
        for(unsigned col=0; col<CsPad::ColumnsPerASIC; col++) {
          getline(&linep, &sz, gain);
          *gn++ = strtod(linep,&pEnd);
          for (unsigned row=1; row < 2*Pds::CsPad::MaxRowsPerASIC; row++)
            *gn++ = strtod(pEnd,&pEnd);
        }
      }
      else {
        float* gn = _gn;
        for(unsigned col=0; col<CsPad::ColumnsPerASIC; col++) {
          for (unsigned row=0; row < 2*Pds::CsPad::MaxRowsPerASIC; row++)
            *gn++ = 1.;
        }
      }
      
      if (sigma) {
        float* sg = _sg;
        float* gn = _gn;
        for(unsigned col=0; col<CsPad::ColumnsPerASIC; col++) {
          getline(&linep, &sz, sigma);
          *sg++ = strtod(linep,&pEnd);
          for (unsigned row=1; row < 2*Pds::CsPad::MaxRowsPerASIC; row++)
            *sg++ = strtod(pEnd,&pEnd)*(*gn++);
        }
      }
      else {
        float* sg = _sg;
        for(unsigned col=0; col<CsPad::ColumnsPerASIC; col++) {
          for (unsigned row=0; row < 2*Pds::CsPad::MaxRowsPerASIC; row++)
            *sg++ = 0;
        }
      }
      
      if (linep) {
        free(linep);
      }
    }
  protected:
    uint16_t  _off[CsPad::MaxRowsPerASIC*CsPad::ColumnsPerASIC*2];
    float     _gn [CsPad::MaxRowsPerASIC*CsPad::ColumnsPerASIC*2];
    float     _sg [CsPad::MaxRowsPerASIC*CsPad::ColumnsPerASIC*2];
  };

  static uint16_t  off_no_ped[CsPad::MaxRowsPerASIC*CsPad::ColumnsPerASIC*2];

#define AsicTemplate(classname,bi,ti,PPB)                               \
  class classname : public AsicP {					\
  public:								\
    classname(double x, double y,                                       \
              FILE* p, FILE* s, FILE* g, FILE* r)                       \
      : AsicP(x,y,PPB,p,s,g,r) {}                                       \
    void boundary(unsigned& dx0, unsigned& dx1,				\
		  unsigned& dy0, unsigned& dy1) const {			\
      FRAME_BOUNDS;							\
      dx0=x0; dx1=x1; dy0=y0; dy1=y1; }					\
    void fill(Ami::DescImage& image) const {				\
      FRAME_BOUNDS;							\
      image.add_frame(x0,y0,x1-x0+1,y1-y0+1);				\
    }									\
    void fill(Ami::EntryImage& image,                                   \
              double v0, double v1) const {                             \
      unsigned data = 0;                                                \
      bool lsuppress  = image.desc().options()&CspadCalib::option_suppress_bad_pixels(); \
      uint16_t* zero = 0;                                               \
      const uint16_t* off = _off;                                       \
      const uint16_t* const * sta = lsuppress ? _sta : &zero;           \
      const float* rms = _sg;                                           \
      ti;                                                               \
    }									\
    void fill(Ami::EntryImage& image,					\
	      const uint16_t*  data) const {                            \
      bool lsuppress  = image.desc().options()&CspadCalib::option_suppress_bad_pixels(); \
      bool lcorrectfn = image.desc().options()&CspadCalib::option_correct_common_mode(); \
      bool lnopedestal= image.desc().options()&CspadCalib::option_no_pedestal(); \
      uint16_t* zero = 0;                                               \
      const uint16_t*  off = lnopedestal ? off_no_ped : _off;           \
      const uint16_t* const * sta = lsuppress ? _sta : &zero;           \
      const float* gn = _gn;                                            \
      double fn = lcorrectfn ? frameNoise(data,off,sta) : 0;            \
      bi;                                                               \
    }                                                                   \
  }

#define B1 { BIN_ITER1(sum1(data,off,sta,fn,gn)); }
#define T1 { BIN_ITER1(thr1(v0,v1,off,sta,rms)); }

#define CALC_X(a,b,c) (a+b)			    
#define CALC_Y(a,b,c) (a-c)			     
  AsicTemplate(  AsicD0B1P, B1, T1, 1);
#undef CALC_X
#undef CALC_Y
#define CALC_X(a,b,c) (a+c)			    
#define CALC_Y(a,b,c) (a+b)			     
  AsicTemplate( AsicD90B1P, B1, T1, 1);
#undef CALC_X
#undef CALC_Y
#define CALC_X(a,b,c) (a-b)			    
#define CALC_Y(a,b,c) (a+c)			     
  AsicTemplate(AsicD180B1P, B1, T1, 1);
#undef CALC_X
#undef CALC_Y
#define CALC_X(a,b,c) (a-c)			    
#define CALC_Y(a,b,c) (a-b)			     
  AsicTemplate(AsicD270B1P, B1, T1, 1);
#undef CALC_X
#undef CALC_Y

#undef B1
#undef T1
#undef AsicTemplate

  class TwoByTwo {
  public:
    TwoByTwo(double x, double y, unsigned ppb, Rotation r, 
	     const Ami::Cspad::TwoByTwoAlignment& a,
             FILE* f=0, FILE* s=0, FILE* g=0, FILE* rms=0) 
    {
      for(unsigned i=0; i<2; i++) {
	double tx(x), ty(y);
	_transform(tx,ty,a.xAsicOrigin[i<<1],a.yAsicOrigin[i<<1],r);
        if (f) {
          switch(r) {
          case D0  : asic[i] = new  AsicD0B1P  (tx,ty,f,s,g,rms); break;
          case D90 : asic[i] = new  AsicD90B1P (tx,ty,f,s,g,rms); break;
          case D180: asic[i] = new  AsicD180B1P(tx,ty,f,s,g,rms); break;
          case D270: asic[i] = new  AsicD270B1P(tx,ty,f,s,g,rms); break;
          default  : break;
          }
        }
        else {
          switch(r) {
          case D0  : asic[i] = new  AsicD0B1  (tx,ty); break;
          case D90 : asic[i] = new  AsicD90B1 (tx,ty); break;
          case D180: asic[i] = new  AsicD180B1(tx,ty); break;
          case D270: asic[i] = new  AsicD270B1(tx,ty); break;
          default  : break;
          }
        }
      }
    }
    ~TwoByTwo() {  for(unsigned i=0; i<2; i++) delete asic[i]; }
    void fill(Ami::DescImage& image,
              unsigned        mask) const
    {
      if (mask&1) asic[0]->fill(image);
      if (mask&2) asic[1]->fill(image);
    }
    void fill(Ami::EntryImage& image,
              unsigned mask,
              double v0, double v1) const
    {
      if (mask&1) asic[0]->fill(image,v0,v1);
      if (mask&2) asic[1]->fill(image,v0,v1);
    }
    void fill(Ami::EntryImage&           image,
	      const CspadElement&        element) const
    {
      asic[0]->fill(image,&element.pair[0][0].s0);
      asic[1]->fill(image,&element.pair[0][0].s1);
    }
  public:
    Asic* asic[2];
  };

  class Detector {
  public:
    Detector(const Src& src,
             FILE* f,    // offsets
             FILE* s,    // status
             FILE* g,    // gain
             FILE* rms,  // noise
             FILE* gm,   // geometry
             unsigned max_pixels) :
      _src   (src)
    {
      //  Determine layout : binning, origin
      double x,y;

      Ami::Cspad::TwoByTwoAlignment qalign = qalign_def[0].twobytwo(0);
      Rotation qrot = D0;
      if (gm) {
        qalign = Ami::Cspad::QuadAlignment::load(gm)->twobytwo(0);

        size_t sz=256;
        char* linep = (char *)malloc(sz);
        while(1) {
          if (getline(&linep, &sz, gm)==-1) break;
          if (linep[0]=='#') continue;
          qrot = Rotation(strtoul(linep,0,0));
          break;
        }
        free(linep);
      }
      //
      //  Create a default layout
      //
      _pixels = 2048-256;
      _ppb = 1;
      { const double frame = double(_pixels)*pixel_size;
	x =  0.5*frame;
	y = -0.5*frame;
      }
      mini = new TwoByTwo(x,y,_ppb,qrot,qalign);

      //
      //  Test extremes and narrow the focus
      //
      unsigned xmin(_pixels), xmax(0), ymin(_pixels), ymax(0);
      for(unsigned i=0; i<2; i++) {
        unsigned x0,x1,y0,y1;
        mini->asic[i&1]->boundary(x0,x1,y0,y1);
        if (x0<xmin) xmin=x0;
        if (x1>xmax) xmax=x1;
        if (y0<ymin) ymin=y0;
        if (y1>ymax) ymax=y1;
      }

      delete mini;

      int idx = xmax-xmin+1;
      int idy = ymax-ymin+1;
      int pixels = ((idx>idy) ? idx : idy);
      const int bin0 = 4;
      _ppb = 1;

      x += pixel_size*double(bin0*int(_ppb) - int(xmin));
      y -= pixel_size*double(bin0*int(_ppb) - int(ymin));

      _pixels = pixels + 2*bin0*_ppb;

      mini = new TwoByTwo(x,y,_ppb,qrot,qalign, f,s,g,rms);
    }
    ~Detector() { delete mini; }

    void fill(Ami::DescImage&    image,
	      Ami::FeatureCache& cache) const
    {
      char buff[64];
      _cache = &cache;
      const char* detname = DetInfo::name(static_cast<const DetInfo&>(_src));
      mini->fill(image, (1<<2)-1);
      for(unsigned a=0; a<4; a++) {
        sprintf(buff,"%s:Temp[%d]",detname,a);
        _feature[a] = cache.add(buff);
      }
#ifdef POST_INTEGRAL
      sprintf(buff,"%s:Cspad::Sum",detname);
      _feature[4] = cache.add(buff);
#endif
    }
    void fill(Ami::EntryImage& image,
	      const Xtc&       xtc) const
    {
      const CspadElement* elem = reinterpret_cast<const CspadElement*>(xtc.payload());
      mini->fill(image,*elem);
      for(int a=0; a<4; a++)
        _cache->cache(_feature[a],
                      CspadTemp::instance().getTemp(elem->sb_temp(a)));
#ifdef POST_INTEGRAL
      if (image.desc().options()&CspadCalib::option_post_integral()) {
        double s = 0;
        double p   = double(image.info(Ami::EntryImage::Pedestal));
        for(unsigned fn=0; fn<image.desc().nframes(); fn++) {
          int xlo(0), xhi(3000), ylo(0), yhi(3000);
          if (image.desc().xy_bounds(xlo, xhi, ylo, yhi, fn)) {
            for(int j=ylo; j<yhi; j++)
              for(int i=xlo; i<xhi; i++) {
                double v = double(image.content(i,j))-p;
                s += v;
              }
          }
        }
        _cache->cache(_feature[4],s);
      }
#endif
    }
    void fill(Ami::EntryImage& image, 
              double v0,double v1) const
    {
      mini->fill(image, (1<<2)-1, v0, v1);
    }
    unsigned ppb() const { return _ppb; }
    unsigned xpixels() { return _pixels; }
    unsigned ypixels() { return _pixels; }
  private:
    TwoByTwo* mini;
    const Src&  _src;
    mutable Ami::FeatureCache* _cache;
    mutable int _feature[4];
    unsigned _ppb;
    unsigned _pixels;
  };

  class CspadMiniPFF : public Ami::PeakFinderFn {
  public:
    CspadMiniPFF(FILE* gain,
                 FILE* rms,
                 const Detector*  detector,
                 Ami::DescImage& image) :
      _detector(*detector),
      _nbinsx  (image.nbinsx()),
      _nbinsy  (image.nbinsy()),
      _values  (new Ami::EntryImage(image))
    {
      image.pedcalib (true);
      image.gaincalib(gain!=0);
      image.rmscalib (rms !=0);
    }
    virtual ~CspadMiniPFF()
    {
      delete _values;
    }
  public:
    void     setup(double v0,
                   double v1)
    {
      const unsigned no_threshold = -1;
      for(unsigned k=0; k<_nbinsy; k++)
        for(unsigned j=0; j<_nbinsx; j++)
          _values->content(no_threshold,j,k);

      _detector.fill(*_values,v0,v1);
    }
    unsigned value(unsigned j, unsigned k) const
    {
      return _values->content(j,k);
    }
    Ami::PeakFinderFn* clone() const { return new CspadMiniPFF(*this); }
  private:
    CspadMiniPFF(const CspadMiniPFF& o) :
      _detector(o._detector),
      _nbinsx  (o._nbinsx),
      _nbinsy  (o._nbinsy),
      _values  (new Ami::EntryImage(o._values->desc())) {}
  private:
    const Detector&  _detector;
    unsigned         _nbinsx;
    unsigned         _nbinsy;
    Ami::EntryImage* _values;
  };
};

using namespace Ami;

static std::list<Pds::TypeId::Type> config_type_list()
{
  std::list<Pds::TypeId::Type> types;
  types.push_back(Pds::TypeId::Id_CspadConfig);
  types.push_back(Pds::TypeId::Id_Cspad2x2Config);
  return types;
}

CspadMiniHandler::CspadMiniHandler(const Pds::DetInfo& info, FeatureCache& features, unsigned max_pixels) :
  EventHandler(info, Pds::TypeId::Id_Cspad2x2Element, config_type_list()),
  _entry(0),
  _detector(0),
  _cache(features),
  _max_pixels(max_pixels),
  _options   (0)
{
  unsigned s = sizeof(CspadMiniGeometry::off_no_ped)/sizeof(uint16_t);
  for (unsigned i=0; i<s; i++) CspadMiniGeometry::off_no_ped[i] = (uint16_t)Offset;
}

CspadMiniHandler::~CspadMiniHandler()
{
  if (_detector)
    delete _detector;
  if (_entry)
    delete _entry;
}

unsigned CspadMiniHandler::nentries() const { return _entry ? 1 : 0; }

const Entry* CspadMiniHandler::entry(unsigned i) const { return i==0 ? _entry : 0; }

const Entry* CspadMiniHandler::hidden_entry(unsigned i) const { return 0; }

void CspadMiniHandler::reset() { _entry = 0; }

void CspadMiniHandler::_configure(Pds::TypeId type,const void* payload, const Pds::ClockTime& t)
{
  //
  //  Load pedestals
  //
  const int NameSize=128;
  char oname1[NameSize];
  char oname2[NameSize];

  sprintf(oname1,"ped.%08x.dat",info().phy());
  sprintf(oname2,"/reg/g/pcds/pds/cspadcalib/ped.%08x.dat",info().phy());
  FILE *f = fopen_dual(oname1, oname2, "pedestals");

  sprintf(oname1,"sta.%08x.dat",info().phy());
  sprintf(oname2,"/reg/g/pcds/pds/cspadcalib/sta.%08x.dat",info().phy());
  FILE *s = fopen_dual(oname1, oname2, "status map");

  sprintf(oname1,"gain.%08x.dat",info().phy());
  sprintf(oname2,"/reg/g/pcds/pds/cspadcalib/gain.%08x.dat",info().phy());
  FILE *g = fopen_dual(oname1, oname2, "gain map");

  sprintf(oname1,"res.%08x.dat",info().phy());
  sprintf(oname2,"/reg/g/pcds/pds/cspadcalib/res.%08x.dat",info().phy());
  FILE *rms = fopen_dual(oname1, oname2, "noise");

  sprintf(oname1,"geo.%08x.dat",info().phy());
  sprintf(oname2,"/reg/g/pcds/pds/cspadcalib/geo.%08x.dat",info().phy());
  FILE *gm = fopen_dual(oname1, oname2, "geometry");

  _create_entry( f,s,g,rms,gm, 
                 _detector, _entry, _max_pixels);

  Ami::PeakFinder::register_(info().phy(),   
                             new CspadMiniGeometry::CspadMiniPFF(g,rms,_detector,_entry->desc()));

  if (f ) fclose(f);
  if (s ) fclose(s);
  if (g ) fclose(g);
  if (rms) fclose(rms);
  if (gm) fclose(gm);
}

void CspadMiniHandler::_create_entry(FILE* f, FILE* s, FILE* g, FILE* rms, FILE* gm,
                                     CspadMiniGeometry::Detector*& detector,
                                     EntryImage*& entry, 
                                     unsigned max_pixels) 
{
  if (f ) rewind(f);
  if (s ) rewind(s);
  if (g ) rewind(g);
  if (rms) rewind(rms);
  if (gm) rewind(gm);

  if (detector)
    delete detector;

  detector = new CspadMiniGeometry::Detector(info(),f,s,g,rms,gm,max_pixels);

  if (entry) 
    delete entry;

  const unsigned ppb = detector->ppb();
  const DetInfo& det = static_cast<const DetInfo&>(info());
  DescImage desc(det, ChannelID::name(det,0), "photons",
                 detector->xpixels()/ppb, detector->ypixels()/ppb, 
                 ppb, ppb,
                 f!=0, g!=0, false);
  desc.set_scale(pixel_size*1e3,pixel_size*1e3);
    
  detector->fill(desc,_cache);

  entry = new EntryImage(desc);
  memset(entry->contents(),0,desc.nbinsx()*desc.nbinsy()*sizeof(unsigned));

  if (f)
    entry->info(Offset*ppb*ppb,EntryImage::Pedestal);
  else
    entry->info(0,EntryImage::Pedestal);
    
  entry->info(0,EntryImage::Normalization);
  entry->invalid();
}

void CspadMiniHandler::_calibrate(const void* payload, const Pds::ClockTime& t) {}
void CspadMiniHandler::_calibrate(Pds::TypeId::Type, const void* payload, const Pds::ClockTime& t) {}

void CspadMiniHandler::_event    (const void* payload, const Pds::ClockTime& t)
{
  const Xtc* xtc = reinterpret_cast<const Xtc*>(payload)-1;
  if (_entry) {
    unsigned o = _entry->desc().options();
    if (_options != o) {
      printf("CspadMiniHandler::event options %x -> %x\n", _options, o);
      _options = o;
    }

    _detector->fill(*_entry,*xtc);
    _entry->info(1,EntryImage::Normalization);
    _entry->valid(t);
  }
}

void CspadMiniHandler::_damaged() 
{
  if (_entry) 
    _entry->invalid(); 
}

