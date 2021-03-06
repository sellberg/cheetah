#ifndef Ami_Aggregator_hh
#define Ami_Aggregator_hh

#include "ami/client/AbsClient.hh"
#include "ami/data/Cds.hh"
#include "ami/data/EntryList.hh"

class iovec;

namespace Ami {

  class BSocket;

  class Aggregator {
  public:
    Aggregator(AbsClient&);
    ~Aggregator();
  public:
    void connected       () ;
    void disconnected    () ;
    int  configure       (iovec*) ;
    int  configured      () ;
    void discovered      (const DiscoveryRx&,unsigned) ;
    int  read_description(Socket&,int,unsigned) ;
    int  read_payload    (Socket&,int,unsigned) ;
    bool svc             () const;
    void process         () ;
    void tmo             () ;
  public:
    void request_payload (const EntryList&);
  private:
    void _checkState     (const char*);
    void _checkState     (const char*, unsigned);
  private:
    AbsClient& _client;
    unsigned   _n;
    unsigned   _remaining;
    unsigned   _nsources;
    Cds        _cds;
    Cds        _ocds;
    unsigned   _niovload;
    iovec*     _iovload;
    iovec*     _iovdesc;
    BSocket*   _buffer;
    enum { Init, Connecting, Connected, 
           Discovering, Discovered, Configured, 
           Describing, Described, Processing } _state;
    double     _latest;
    unsigned   _current;
    EntryList  _request;
  };
};

#endif
