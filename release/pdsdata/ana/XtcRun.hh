#ifndef Pds_XtcRun_hh
#define Pds_XtcRun_hh

#include "pdsdata/ana/XtcSlice.hh"

#include <list>
#include <vector>
#include <string>

namespace Pds
{
namespace Ana
{
extern bool _live;

class XtcRun {
public:
  XtcRun();
  ~XtcRun();

  void reset   (std::string fname);
  bool add_file(std::string fname);

  const char* base() const;
  unsigned    run_number() const;

  void   init();
  Result next(Pds::Dgram*& dg, int* piSlice = NULL, int64_t* pi64OffsetCur = NULL);

  int jump    (int calib, int jump, int& eventNum);

  int findTime(const char* sTime, int& iCalib, int& iEvent, bool& bExactMatch, bool& bOvertime);
  int findTime(uint32_t uSeconds, uint32_t uNanoseconds, int& iCalib, int& iEvent, bool& bExactMatch, bool& bOvertime);
  int getStartAndEndTime(ClockTime& start, ClockTime& end);
  int findNextFiducial
              (uint32_t uFiducialSearch, int iFidFromEvent, int& iCalib, int& iEvent);
  int numTotalEvent(int& iNumTotalEvent);
  int numCalib(int& iNumCalib);
  int numEventInCalib(int calib, int& iNumEvents);
  int curEventId(int& eventGlobal, int& calib, int& eventInCalib);

  int eventGlobalToCalib(int eventGlobal, int& calib, int& eventInCalib);
  int eventCalibToGlobal(int calib, int eventInCalib, int& eventGlobal);
  int nextCalibEventId(int iNumCalibAfter, int iNumEventAfter, bool bResetEventInCalib, int& eventGlobal, int& calib, int& eventInCalib);
  int nextGlobalEventId(int iNumEventAfter, int& eventGlobal, int& calib, int& eventInCalib);

  static void live_read(bool l);
  static void read_ahead(bool l);

private:
  int eventGlobalToSlice
              ( int iGlobalEvent, std::vector<int>& lSliceEvent );

  std::list<XtcSlice*>  _slices;
  std::list<XtcSlice*>  _doneSlices;
  std::string           _base;
  bool                  _startAndEndValid;
  int                   _iNumTotalEvent;
  int                   _iNumTotalCalib;
  std::vector<int>      _vNumEventsInCalib;
  ClockTime             _start;
  ClockTime             _end;
  int                   _iCurCalib;
  int                   _iCurEventGlobal;
  int                   _iCurCalibBaseEvt;
};

} // namespace Ana
} // namespace Pds

#endif // #ifndef Pds_XtcRun_hh
