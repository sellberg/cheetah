#ifndef Pds_VPoll_hh
#define Pds_VPoll_hh

#include "ami/service/Routine.hh"

#include <poll.h>

namespace Ami {

  class Fd;
  class Task;
  class Socket;

  class VPoll : private Routine {
  public:
    VPoll(int timeout);
    ~VPoll();
  public:
    void manage  (Fd&);
    void unmanage(Fd&);
  private:
    int poll();
    void routine();
  protected:
    int  timeout() const;
    void timeout(int);
  private:
    virtual int processTmo() = 0;
  private:
    int        _timeout;
    Task*      _task;
    Socket*    _loopback;
    Fd*        _ofd;
    pollfd     _pfd[2];
  };

};

#endif
