#include "Discovery.hh"

#include "ami/client/ClientManager.hh"
#include "ami/data/DescEntry.hh"
#include "ami/data/Discovery.hh"
#include "ami/service/ConnectionManager.hh"

static const Pds::DetInfo envInfo(0,Pds::DetInfo::NoDetector,0,Pds::DetInfo::Evr,0);
static const Pds::DetInfo noInfo (0,Pds::DetInfo::NoDetector,0,Pds::DetInfo::NoDevice,0);

static const int MaxConfigSize = 0x100000;

using namespace Ami::Python;

Discovery::Discovery(unsigned ppinterface,
		     unsigned interface,
		     unsigned serverGroup) :
  _interface  (interface),
  _serverGroup(serverGroup),
  _connect_mgr(new ConnectionManager(ppinterface)),
  _manager    (new ClientManager(interface,
                                 serverGroup, 
				 *_connect_mgr,
                                 *this))
{
  _manager->connect();
}

Discovery::~Discovery()
{
  delete _manager;
}

Ami::ClientManager* Discovery::allocate(Ami::AbsClient& c)
{
  return new ClientManager(_interface, 
			   _serverGroup, 
                           *_connect_mgr,
			   c);
}

void Discovery::connected       () { _manager->discover(); }

int Discovery:: configure       (iovec* iov) { return 0; }

int Discovery:: configured      () { return 0; }

void Discovery::discovered      (const DiscoveryRx& rx) 
{ 
  const Pds::DetInfo noInfo;
  const Ami::DescEntry* n;
  for(const Ami::DescEntry* e = rx.entries(); e < rx.end(); e = n) {
    n = reinterpret_cast<const Ami::DescEntry*>
      (reinterpret_cast<const char*>(e) + e->size());
    
    if (e->info() == noInfo)   // Skip unknown devices 
      continue;
  }
}

int Discovery::read_description(Socket& s,int len) {
  return len;
}

int Discovery::read_payload    (Socket& s,int) {
  return 0;
}

bool Discovery::svc             () const {
  return false;
}

void Discovery::process         () {
}

