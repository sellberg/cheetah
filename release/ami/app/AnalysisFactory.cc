#include "AnalysisFactory.hh"

#include "ami/app/EventFilter.hh"
#include "ami/data/Analysis.hh"
#include "ami/data/UserModule.hh"

#include "ami/data/FeatureCache.hh"
#include "ami/data/Message.hh"
#include "ami/data/ConfigureRequest.hh"
#include "ami/data/DescEntry.hh"
#include "ami/data/DescCache.hh"
#include "ami/data/Entry.hh"
#include "ami/data/EntryFactory.hh"
#include "ami/server/ServerManager.hh"

//#define DBUG

using namespace Ami;

AnalysisFactory::AnalysisFactory(std::vector<FeatureCache*>&  cache,
				 ServerManager& srv,
                                 UList&         user,
                                 EventFilter&   filter) :
  _srv       (srv),
  _cds       ("Analysis"),
  _ocds      ("Hidden"),
  _features  (cache),
  _user      (user),
  _user_cds  (user.size()),
  _filter    (filter),
  _waitingForConfigure(false)
{
  pthread_mutex_init(&_mutex, NULL);
  pthread_cond_init(&_condition, NULL);
  EntryFactory::source(*_features[PostAnalysis]);

  unsigned i=0;
  for(UList::iterator it=user.begin(); it!=user.end(); it++,i++) {
    _user_cds[i] = 0;
    (*it)->driver = this;
  }
}

AnalysisFactory::~AnalysisFactory()
{
  for(unsigned i=0; i<_user_cds.size(); i++)
    if (_user_cds[i])
      delete _user_cds[i];
}

std::vector<FeatureCache*>& AnalysisFactory::features() { return _features; }

Cds& AnalysisFactory::discovery() { return _cds; }
Cds& AnalysisFactory::hidden   () { return _ocds; }

//
//  The discovery cds has changed
//
void AnalysisFactory::discover(bool waitForConfigure) 
{
  pthread_mutex_lock(&_mutex);
  for(AnList::iterator it=_analyses.begin(); it!=_analyses.end(); it++) {
    Analysis* an = *it;
    delete an;
  }
  _analyses.clear();

  if (waitForConfigure) {
    if (_waitingForConfigure) {
      printf("Already waiting for configure?\n");
      *((char *)0) = 0;
    }
    _waitingForConfigure = true;
  }

  _features[PostAnalysis]->update();  // clear the update flag
  _srv.discover(); 

  if (waitForConfigure) {
    printf("AnalysisFactory: waiting for configuration...\n");
    for (;;) {
      if (! _waitingForConfigure) {
        break;
      }
      pthread_cond_wait(&_condition, &_mutex);
    }
    printf("AnalysisFactory: done waiting for configuration.\n");
  }

  pthread_mutex_unlock(&_mutex);
}

void AnalysisFactory::configure(unsigned       id,
				const Message& msg, 
				const char*    payload, 
				Cds&           cds)
{
  pthread_mutex_lock(&_mutex);

  //
  //  Segregate the entries belonging to the configuring client
  //
  AnList newlist;
  AnList oldList;
  for(AnList::iterator it=_analyses.begin(); it!=_analyses.end(); it++) {
    Analysis* a = *it;
    if (a->id() == id)
      oldList.push_back(a);
    else
      newlist.push_back(a);
  }
  _analyses = newlist;

  // create
  const char* const end = payload + msg.payload();
  while(payload < end) {
    const ConfigureRequest& req = *reinterpret_cast<const ConfigureRequest*>(payload);
    
    if (req.size()==0) {
      const uint32_t* p = reinterpret_cast<const uint32_t*>(&req);
      printf("AnalysisFactory::configure received corrupt request from id %d\n [%08x %08x %08x %08x %08x %08x]\n",
             id, p[0], p[1], p[2], p[3], p[4], p[5] );
      break;
    }

    if (req.source() == ConfigureRequest::User) {
      int i=0;
      for(UList::iterator it=_user.begin(); it!=_user.end(); it++,i++) {
        if (req.input() == i) {
	  _srv.register_key(i,id);
          Cds* ucds = _user_cds[i];
          if (ucds==0) {
            _user_cds[i] = ucds = new Cds((*it)->name());
            (*it)->clear();
            (*it)->create(*ucds);
          }
          ucds->mirror(cds);
          break;
        }
      }
    }
    else if (req.source() == ConfigureRequest::Filter) {
      unsigned o = *reinterpret_cast<const uint32_t*>(&req+1);
      _filter.enable(o);
    }
    else {
      const Cds* pcds = 0;
      char *reqType = NULL;
      switch(req.source()) {
        case ConfigureRequest::Discovery: pcds = &_cds; reqType = "Discovery"; break;
        case ConfigureRequest::Analysis : pcds = & cds; reqType = "Analysis"; break;
        case ConfigureRequest::Hidden   : pcds = &_ocds; reqType = "Hidden"; break;
        default:
          printf("AnalysisFactory::configure unknown source ConfigureRequest::%d\n",req.source());
          break;
      }
      if (!pcds) {
        printf("AnalysisFactory::configure failed (null pcds) for ConfigureRequest::%s\n", reqType);
        printf("\tinp %d  out %d  size %d\n",req.input(),req.output(),req.size());
        break;
      }
      const Cds& input_cds = *pcds;
      const Entry* input = input_cds.entry(req.input());
      if (!input) {
        printf("AnalysisFactory::configure failed (null input) for ConfigureRequest::%s\n", reqType);
        printf("\tinp %d  out %d  size %d\n",req.input(),req.output(),req.size());
        input_cds.showentries();
      }
      else if (req.state()==ConfigureRequest::Create) {
        //
        //  If the output entry already exists, don't need to create it.
        //
        bool lFound=false;
        for(AnList::iterator it=oldList.begin(); it!=oldList.end(); it++) {
          Analysis* a = *it;
          if (a->output().signature()==req.output() && a->valid()) {
#ifdef DBUG
            printf("Preserving analysis for %s [%d]\n",
                   a->output().name(),
                   a->output().signature());
#endif
            _analyses.push_back(a);
            oldList.remove(a);
            lFound=true;
            break;
          }
        }
        if (!lFound) {
          const char*  p     = reinterpret_cast<const char*>(&req+1);
          Analysis* a = new Analysis(id, *input, req.output(),
                                     cds, *_features[req.scalars()], p);
          _analyses.push_back(a);
#ifdef DBUG
          printf("Created analysis for %s [%d] features %d\n",
                 a->output().name(),
                 a->output().signature(),
		 req.scalars());
#endif
        }
      }
      else if (req.state()==ConfigureRequest::SetOpt) {
        unsigned o = *reinterpret_cast<const uint32_t*>(&req+1);
        printf("Set options %x on entry %p\n",o,input);
        const_cast<Ami::DescEntry&>(input->desc()).options(o);
      }
    }
    payload += req.size();
  }

  //
  //  Clean up any old client entries that were not requested
  //
  for(AnList::iterator it=oldList.begin(); it!=oldList.end(); it++) {
    Analysis* a = *it;
#ifdef DBUG
          printf("Removing analysis for %s [%d]\n",
                 a->output().name(),
                 a->output().signature());
#endif
    delete a;
  }

  cds.showentries();

  //
  //  Reconfigure Post Analyses whenever the PostAnalysis variable cache changes
  //    This guarantees that the post analyses are performed after the variables
  //      they depend upon are cached.
  //
  if (_features[PostAnalysis]->update()) {
    printf("PostAnalysis updated\n");
    _srv.discover_post();
  }

  if (_waitingForConfigure) {
    _waitingForConfigure = false;
    pthread_cond_signal(&_condition);
    printf("AnalysisFactory: configuration is done.\n");
  }

  pthread_mutex_unlock(&_mutex);
}

void AnalysisFactory::analyze  ()
{
  pthread_mutex_lock(&_mutex);
  for(AnList::iterator it=_analyses.begin(); it!=_analyses.end(); it++) {
    (*it)->analyze();
  }
  //
  //  Post-analyses go here
  //
  pthread_mutex_unlock(&_mutex);
}

void AnalysisFactory::remove(unsigned id)
{
  pthread_mutex_lock(&_mutex);
  AnList newlist;
  for(AnList::iterator it=_analyses.begin(); it!=_analyses.end(); it++) {
    Analysis* a = *it;
    if (a->id() == id) {
      delete a;
    }
    else
      newlist.push_back(a);
  }
  _analyses = newlist;
  pthread_mutex_unlock(&_mutex);
}

void AnalysisFactory::recreate(UserModule* user)
{
  int i=0;
  for(UList::iterator it=_user.begin(); it!=_user.end(); it++,i++)
    if (user == *it) {
      _user_cds[i] = 0;
      _srv.discover_key(i);
      break;
    }
}
