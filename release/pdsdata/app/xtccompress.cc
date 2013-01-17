//
//  Unofficial example of XTC compression
//
#include "pdsdata/camera/FrameV1.hh"
#include "pdsdata/compress/Camera_FrameV1.hh"

#include "pdsdata/cspad/ConfigV1.hh"
#include "pdsdata/cspad/ConfigV2.hh"
#include "pdsdata/cspad/ConfigV3.hh"
#include "pdsdata/cspad/ConfigV4.hh"
#include "pdsdata/cspad/ElementV1.hh"
#include "pdsdata/cspad/ElementV2.hh"
#include "pdsdata/compress/Cspad_ElementV2.hh"
#include "pdsdata/cspad/ElementIterator.hh"

#include "pdsdata/xtc/Dgram.hh"
#include "pdsdata/xtc/XtcFileIterator.hh"

#include <string.h>

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>

#include <map>
using std::map;

#include "xtccompress_test.icc"

using namespace Pds;

using Pds::Camera::FrameCompressor;

class myIter {
public:
  enum Status {Stop, Continue};
  myIter(bool extract) : _extract(extract), _aligned(true), _pwrite(0),
                         _outbuf(new char[0x1000000]),
                         _outsz (0x1000000),
                         _frameCmp(new FrameCompressor)
  {
  }
  ~myIter() 
  {
    delete _frameCmp;
    delete[] _outbuf;
  }  
private:
  //
  //  Iterate through the Xtc and compress, decompress, copy into new Xtc
  //
  void iterate(Xtc* root) {
    if (root->damage.value() & ( 1 << Damage::IncompleteContribution)) {
      return _write(root,root->extent);
    }
    
    int remaining = root->sizeofPayload();
    Xtc* xtc     = (Xtc*)root->payload();

    uint32_t* pwrite = _pwrite;
    _write(root, sizeof(Xtc));
    
    while(remaining > 0)
      {
	unsigned extent = xtc->extent;
	if(extent==0) {
	  printf("Breaking on zero extent\n");
	  break; // try to skip corrupt event
	}
	process(xtc);
	remaining -= extent;
	xtc        = (Xtc*)((char*)xtc+extent);
      }

    reinterpret_cast<Xtc*>(pwrite)->extent = (_pwrite-pwrite)*sizeof(uint32_t);
  }
  
  void process(Xtc* xtc) {

    //
    //  We're only interested in compressing/decompressing
    //
    switch (xtc->contains.id()) {
    case (TypeId::Id_Xtc):
      iterate(xtc);
      return;
    case (TypeId::Id_CspadConfig):
      _cache_config( xtc );
      break;
    default:
      break;
    }

    if (xtc->contains.compressed()) {

      if (_extract) {   // Anything to decompress?
        
        switch(xtc->contains.id()) {

        case (TypeId::Id_Frame) :
          switch (xtc->contains.compressed_version()) {
          case 1: {
            const Camera::CompressedFrameV1* pframe = 
              reinterpret_cast<const Camera::CompressedFrameV1*>(xtc->payload());
            if (_decompress(xtc, *pframe))
              return;
            else
              printf("decompress %x failed\n",xtc->contains.value());
            break;
          }
          default:
            printf("Compressed Id_Frame version %d unsupported\n",
                   xtc->contains.compressed_version());
            break;
          }
          break;

        case (TypeId::Id_CspadElement) : {
          switch (xtc->contains.compressed_version()) {
          case 2: {
            const CsPad::CompressedElementV2* pframe =
              reinterpret_cast<const CsPad::CompressedElementV2*>(xtc->payload());
            if (_decompress(xtc, *pframe))
              return;
            else
              printf("decompress %x failed\n",xtc->contains.value());
            break;
          }
          default:
            printf("Compressed Id_CspadElement version %d [%x] unsupported\n",
                   xtc->contains.compressed_version(),
                   xtc->contains.value());
            break;
          }
          break;
        }
          
        default:
          break;
        }
      }
    }
    else if (!_extract) {  // Anything to compress?

      switch (xtc->contains.id()) {

      case (TypeId::Id_Frame) : {
        if (!xtc->damage.value())
          switch (xtc->contains.version()) {
          case 1: {
            const Camera::FrameV1& pframe = *reinterpret_cast<const Camera::FrameV1*>(xtc->payload());
            if (_compress(xtc, pframe, CompressedPayload::None))
              return;
            break;
          }
          default:
            printf("Uncompressed Id_Frame version %d unsupported\n",
                   xtc->contains.compressed_version());
            break;
          }
        break;
      }

      case (TypeId::Id_CspadElement) : {
        if (!xtc->damage.value()) {
          CsPad::ElementIterator* iter = _lookup_iterator(xtc);
          if (_compress(xtc, iter, CompressedPayload::Hist16))
            return;
        }
        break;
      }

      default:
        break; 
      }
    }
    _write(xtc,xtc->extent);
  }

private:
  //
  //  Xtc headers are 32b aligned.  
  //  Compressed data is not.
  //  Enforce alignment during Xtc construction.
  //
  char* _new(ssize_t sz)
  {
    uint32_t* p = _pwrite;
    _pwrite += sz>>2;
    return (char*)p;
  }

  void _write(const void* p, ssize_t sz) 
  {
    if (!_aligned)
      perror("Writing 32b data alignment not guaranteed\n");

    const uint32_t* pread = (uint32_t*)p;
    if (_pwrite!=pread) {
      const uint32_t* const end = pread+(sz>>2);
      while(pread < end)
	*_pwrite++ = *pread++;
    }
    else
      _pwrite += sz>>2;
  }
  void _uwrite(const void* p, ssize_t sz) 
  {
    if (_aligned)
      perror("Writing 8b data when 32b alignment required\n");

    const uint8_t* pread = (uint8_t*)p;
    if (_upwrite!=pread) {
      const uint8_t* const end = pread+sz;
      while(pread < end)
	*_upwrite++ = *pread++;
    }
    else
      _upwrite += sz;
  }
  void _align_unlock()
  {
    _aligned = false;
    _upwrite = (uint8_t*)_pwrite;
  }
  void _align_lock()
  {
    _pwrite += (_upwrite - (uint8_t*)_pwrite +3)>>2;
    _aligned = true;
  }
  
public:
  void iterate(const Dgram* dg, uint32_t* pwrite) 
  {
    _pwrite = pwrite;
    _write(dg, sizeof(*dg)-sizeof(Xtc));
    iterate(const_cast<Xtc*>(&(dg->xtc)));
  }

private:
  //
  //  Compress Camera::FrameV1
  //
  bool _compress(const Xtc*                xtc,
                 const Camera::FrameV1&    frame,
                 CompressedPayload::Engine engine)
  {
    const CompressedPayload* pd = _compress(engine,
                                            frame.data(), 
                                            frame.data_size());

    if (!pd) 
      return false;

    // Copy the xtc header
    uint32_t* pwrite = _pwrite;
    { 
      Xtc nxtc( TypeId(xtc->contains.id(), xtc->contains.version(), true),
                xtc->src,
                xtc->damage );
      _write(&nxtc, sizeof(Xtc));
    }
    
    _align_unlock();
    
    Camera::CompressedFrameV1 cframe(frame);
    _uwrite(&cframe, sizeof(cframe)-sizeof(CompressedPayload));

    _uwrite(pd, sizeof(*pd)+pd->csize());

    _align_lock();

    //  Update the extent of the container
    reinterpret_cast<Xtc*>(pwrite)->extent = (_pwrite-pwrite)*sizeof(uint32_t);
    return true;
  }

  //
  //  Decompress Camera::FrameV1
  //
  bool _decompress(const Xtc* xtc,
                   const Camera::CompressedFrameV1& frame)
  {      
    // Copy the xtc header
    uint32_t* pwrite = _pwrite;
    { 
      Xtc nxtc(TypeId(xtc->contains.id(), xtc->contains.compressed_version()),
               xtc->src,
               xtc->damage);
      _write(&nxtc, sizeof(Xtc));
    }

    void* outbuf;
    _decompress(frame.pd(), outbuf);

    if (!outbuf)
      return false;

#if 1
    Camera::FrameV1 cframe(frame.width (),
			   frame.height(),
			   frame.depth (),
			   frame.offset());
    _write(&cframe, sizeof(cframe));
#else    
    new (_new(sizeof(Camera::FrameV1))) Camera::FrameV1(frame.width (),
							frame.height(),
							frame.depth (),
							frame.offset());
#endif
    _write(outbuf, frame.pd().dsize());

    //  Update the extent of the container
    reinterpret_cast<Xtc*>(pwrite)->extent = (_pwrite-pwrite)*sizeof(uint32_t);
    return true;
  }

  bool _compress(const Xtc* xtc, 
                 CsPad::ElementIterator* iter,
                 CompressedPayload::Engine engine) 
  {
    if (!iter) return false;

    // Copy the xtc header
    uint32_t* pwrite = _pwrite;
    { 
      Xtc nxtc(TypeId(xtc->contains.id(), xtc->contains.version(), true),
               xtc->src,
               xtc->damage);
      _write(&nxtc, sizeof(Xtc));
    }

    _align_unlock();

    bool result = true;

    const Pds::CsPad::ElementHeader* hdr;
    while( result && (hdr = iter->next()) ) {

      unsigned nsects = 0;
      const Pds::CsPad::Section* s;
      unsigned secnum;
      while( (s = iter->next(secnum)) )
        nsects++;

      _uwrite(hdr, sizeof(*hdr));

      const CompressedPayload* pd = _compress(engine,
                                              hdr+1, 
                                              nsects*Pds::CsPad::ColumnsPerASIC*Pds::CsPad::MaxRowsPerASIC*2*2);

      if (!pd) {
        result = false;
        break;
      }

      _uwrite(pd, sizeof(*pd)+pd->csize());
    
      uint32_t qw = iter->getQuadWord();
      _uwrite(&qw, sizeof(qw));
    }

    delete iter;

    if (!result) {
      //  undo
      return false;
    }

    _align_lock();

    //  Update the extent of the container
    reinterpret_cast<Xtc*>(pwrite)->extent = (_pwrite-pwrite)*sizeof(uint32_t);

    return true;
  }

  bool _decompress(const Xtc* xtc, 
                   const CsPad::CompressedElementV2& elem)
  {
    // Copy the xtc header
    uint32_t* pwrite = _pwrite;
    { 
      Xtc nxtc(TypeId(xtc->contains.id(),xtc->contains.compressed_version()),
               xtc->src,
               xtc->damage);
      _write(&nxtc, sizeof(Xtc));
    }

    const CsPad::CompressedElementV2* hdr = &elem;
    const char* end = xtc->payload()+xtc->sizeofPayload()-4;
    // Copy the quadrant headers and the section images
    while( ((const char*)hdr < end) ) {

      void*    outbuf;
      const CompressedPayload& pyl = hdr->pd();

      _decompress(pyl, outbuf);

      if (!outbuf) return false;

      _write(hdr, sizeof(*hdr)-sizeof(pyl));              // quadrant header
      _write(outbuf, pyl.dsize());                        // section images
      _write((const char*)pyl.cdata()+pyl.csize(), sizeof(uint32_t));  // 
      
      const uint8_t* next = (const uint8_t*)hdr->pd().cdata() + hdr->pd().csize() + sizeof(uint32_t);
      
      hdr = reinterpret_cast<const CsPad::CompressedElementV2*>( next );
    }

    //  Update the extent of the container
    reinterpret_cast<Xtc*>(pwrite)->extent = (_pwrite-pwrite)*sizeof(uint32_t);

    return true;
  }

  CsPad::ElementIterator* _lookup_iterator(const Xtc* xtc)
  {
#define CSPAD_VER(v) {                                                  \
      XtcMapKey key(*xtc, TypeId(TypeId::Id_CspadConfig, v));           \
      if (_xtcmap.find(key)!=_xtcmap.end())                             \
        return new CsPad::ElementIterator( *reinterpret_cast<const CsPad::ConfigV##v*>(_xtcmap[ key ]->payload()), *xtc ); \
    }

    CSPAD_VER(1);
    CSPAD_VER(2);
    CSPAD_VER(3);
    CSPAD_VER(4);
    return 0;

#undef CSPAD_VER
  }

private:

  //
  //  The real interface
  //
  const CompressedPayload* _compress(CompressedPayload::Engine engine,
                                     const void* inbuf,
                                     unsigned    insz)
  {
    bool     result;
    void*    outbuf;
    size_t   outsz;

    switch(engine) {
    case CompressedPayload::None:
      outbuf = const_cast<void*>(inbuf);
      outsz  = insz;
      result = true;
      break;
    case CompressedPayload::Hist16: {
      FrameCompressor::ImageParams params;
      params.width  = insz/2;
      params.height = 1;
      params.depth  = 2;
      
      result = _frameCmp->compress(inbuf,
                                   params, 
                                   outbuf, 
                                   outsz) == FrameCompressor::Success;
      break; }
    default:
      result = false;
      break;
    }

    if (!result) {
      return 0;
    }

    unsigned sz = sizeof(CompressedPayload)+outsz;
    if (sz > _outsz) {
      delete[] _outbuf;
      _outbuf = new char[_outsz = sz];
    }
    CompressedPayload* pd = new(_outbuf) CompressedPayload(engine, insz, outsz);
    memcpy(const_cast<void*>(pd->cdata()), outbuf, outsz);
    return pd;
  }

  void _decompress(const CompressedPayload& pyl,
                   void*&    outbuf)
  {
    bool result;

    switch(pyl.compressor()) {
    case CompressedPayload::None:
      outbuf = const_cast<void*>(pyl.cdata());
      result = true;
      break;
    case CompressedPayload::Hist16: {
      FrameCompressor::ImageParams params;
      params.width  = pyl.dsize()/2;
      params.height = 1;
      params.depth  = 2;

      int r = _frameCmp->decompress(pyl.cdata(),
                                     pyl.csize(),
                                     params, 
                                    outbuf);

      result = (r == FrameCompressor::Success);
      break; }
    default:
      result = false;
      break;
    }

    if (!result)
      outbuf = 0;
  }

private:
  //
  //  Cache a copy of this xtc
  //
  void _cache_config(const Xtc* xtc) {
    char* p = new char[xtc->extent];
    memcpy(p, xtc, xtc->extent);
    XtcMapKey key(*xtc);
    if (_xtcmap.find(key)!=_xtcmap.end())
      delete[] reinterpret_cast<const char*>(_xtcmap[key]);
    _xtcmap[ XtcMapKey(*xtc) ] = reinterpret_cast<const Xtc*>(p); 
  }
private:
  bool                                      _extract;
  bool                                      _aligned;
  uint32_t*                                 _pwrite;
  uint8_t*                                  _upwrite;

  char*                                     _outbuf;
  unsigned                                  _outsz;
  FrameCompressor*                          _frameCmp;

  map<XtcMapKey, const Xtc*>                _xtcmap;
};

void usage(char* progname) {
  fprintf(stderr,"Usage: %s -i <input filename> -o <output filename> [-x] [-h]\n  -x : extract (decompress)\n", progname);
}

int main(int argc, char* argv[]) {
  int c;
  char* inxtcname=0;
  char* outxtcname=0;
  bool extract=false;
  int parseErr = 0;
  unsigned nevents = -1;

  while ((c = getopt(argc, argv, "hxn:i:o:")) != -1) {
    switch (c) {
    case 'h':
      usage(argv[0]);
      exit(0);
    case 'x':
      extract = true;
      break;
    case 'n':
      nevents = atoi(optarg);
      break;
    case 'i':
      inxtcname = optarg;
      break;
    case 'o':
      outxtcname = optarg;
      break;
    default:
      parseErr++;
    }
  }
  
  if (!inxtcname || !outxtcname) {
    usage(argv[0]);
    exit(2);
  }

  int ifd = open(inxtcname, O_RDONLY | O_LARGEFILE);
  if (ifd < 0) {
    perror("Unable to open input file\n");
    exit(2);
  }

  FILE* ofd = fopen(outxtcname,"wx");
  if (ofd == 0) {
    perror("Unable to open output file\n");
    exit(2);
  }
  
  const unsigned MAX_DG_SIZE = 0x2000000;
  XtcFileIterator iter(ifd,MAX_DG_SIZE);
  Dgram* dg;

  uint32_t* obuff = new uint32_t[MAX_DG_SIZE>>2];

  myIter cmpiter(extract);

  while ((dg = iter.next())) {
    //    if (!dg->seq.isEvent())
    cmpiter.iterate(dg, obuff);
    
    const Dgram* odg = reinterpret_cast<const Dgram*>(obuff);
    fwrite(odg, sizeof(*odg) + odg->xtc.sizeofPayload(), 1, ofd);
    fflush(ofd);

    printf("%s transition: time 0x%x/0x%x, payloadSize %d (%d)\n",TransitionId::name(dg->seq.service()),
           dg->seq.stamp().fiducials(),dg->seq.stamp().ticks(), dg->xtc.sizeofPayload(), odg->xtc.sizeofPayload());

    if (dg->seq.isEvent())
      if (--nevents == 0)
        break;
  }
  
  close (ifd);
  fclose(ofd);
  return 0;
}

