#include "pdsdata/ana/XtcRun.hh"

namespace Pds
{
namespace Ana
{
bool _live=false;

void XtcRun::live_read(bool l) { _live=l; }

XtcRun::XtcRun() : _startAndEndValid(false),
  _iNumTotalEvent(0), _iNumTotalCalib(0) {}

static void deleteSlices(std::list<XtcSlice*>& slices) {
  for (std::list<XtcSlice*>::iterator it = slices.begin(); it != slices.end(); it++) {
    //printf("deleteSlices: deleting %p\n", *it);
    delete (*it);
  }
  slices.clear();
}


XtcRun::~XtcRun()
{
  deleteSlices(_slices);
  deleteSlices(_doneSlices);
}

void XtcRun::reset(std::string fname)
{
  deleteSlices(_slices);
  deleteSlices(_doneSlices);
  _slices.push_back(new XtcSlice(fname));
  _base = fname.substr(0,fname.find("-s"));
}

bool XtcRun::add_file(std::string fname)
{
  if (fname.compare(0,_base.size(),_base)!=0)
    return false;

  for(std::list<XtcSlice*>::iterator it=_slices.begin();
      it!=_slices.end(); it++)
    if ((*it)->add_file(fname))
      return true;
  _slices.push_back(new XtcSlice(fname));
  return true;
}

const char* XtcRun::base() const
{ return _base.c_str(); }

unsigned XtcRun::run_number() const
{ return atoi(_base.c_str()+_base.find("-r")+2); }


void XtcRun::init()
{
  for(std::list<XtcSlice*>::iterator it=_slices.begin();
      it!=_slices.end(); it++)
    (*it)->init();

  int iError;
  _iNumTotalEvent = 0;
  XtcSlice* firstSlice = NULL;
  XtcSlice* lastSlice = NULL;
  for(std::list<XtcSlice*>::iterator it=_slices.begin();
      it!=_slices.end(); it++) {
    int iNumSliceEvent = 0;
    iError = (*it)->numTotalEvent(iNumSliceEvent);
    if (iError != 0) printf( "XtcRun::init(): Query Slice Event# failed\n" );
    _iNumTotalEvent += iNumSliceEvent;

    if (firstSlice == NULL) {
      firstSlice = *it;
      iError = firstSlice->numCalib(_iNumTotalCalib);
      if (iError != 0) printf( "XtcRun::init(): Query Slice Calib# failed\n" );
    }
    lastSlice = *it;
  }

  _startAndEndValid = false;
  if (firstSlice == NULL || lastSlice == NULL) {
    return;
  }

  _vNumEventsInCalib.clear();
  _vNumEventsInCalib.push_back(0);
  for (int iCalib = 1; iCalib <= _iNumTotalCalib; ++iCalib)
  {
    int iCalibEventSum = 0;

    int iSlice = 0;
    for(std::list<XtcSlice*>::iterator
        it =  _slices.begin();
        it != _slices.end();
        it++, iSlice++)
    {
      int iNumEvent = -1;
      iError = (*it)->numEventInCalib(iCalib, iNumEvent);
      if ( iError != 0 )
      {
        printf( "XtcRun::init(): Query Event# for Calib# %d in Slice %d failed\n", iCalib, iSlice );
        continue;
      }

      iCalibEventSum += iNumEvent;
    }

    _vNumEventsInCalib.push_back(iCalibEventSum);
  }

  // Set _start
  uint32_t seconds, nanoseconds;
  int index = 1;
  iError = firstSlice->getTimeGlobal(index, seconds, nanoseconds);
  if (iError) {
    return;
  }
  _start = ClockTime(seconds, nanoseconds);

  // Set _end
  iError = lastSlice->numTotalEvent(index);
  if (iError) {
    return;
  }
  iError = lastSlice->getTimeGlobal(index, seconds, nanoseconds);
  if (iError) {
    return;
  }
  _end = ClockTime(seconds, nanoseconds);
  _startAndEndValid = true;

  _iCurCalib        = 0;
  _iCurEventGlobal  = 0;
  _iCurCalibBaseEvt = 0;
}

int XtcRun::numTotalEvent(int& iNumTotalEvent) {
  if (! _startAndEndValid) {
    iNumTotalEvent = -1;
    return 1;
  }
  iNumTotalEvent = _iNumTotalEvent;
  return 0;
}

int XtcRun::getStartAndEndTime(ClockTime& start, ClockTime& end) {
  if (! _startAndEndValid) {
    return 1;
  }
  start = _start;
  end = _end;
  return 0;
}

static bool hasNullSequence(std::list<XtcSlice*>::iterator it) {
  XtcSlice* slice = *it;
  const Sequence* sequence = &slice->hdr().seq;
  return (sequence == NULL);
}

Result XtcRun::next(Pds::Dgram*& dg, int* piSlice, int64_t* pi64OffsetCur)
{
  //
  //  Process L1A with lowest clock time first
  //
  Pds::ClockTime tmin(-1,-1);

  int iSlice      = 0;
  int iNextSlice  = 0;
  std::list<XtcSlice*>::iterator n  = _slices.begin();
  for(std::list<XtcSlice*>::iterator it = _slices.begin();
      it != _slices.end(); it++, iSlice++) {
    if (hasNullSequence(it)) {
      printf("XtcRun::next: %d: null sequence for slice=%p\n", iSlice, *it);
      continue;
    }
    if ((*it)->hdr().seq.service()==Pds::TransitionId::L1Accept &&
        tmin > (*it)->hdr().seq.clock())
    {
      tmin = (*(n = it))->hdr().seq.clock();
      iNextSlice  = iSlice;
    }
  }
  if (hasNullSequence(n)) {
    printf("XtcRun::next: no valid slices found, returning End\n");
    return End;
  }

  //
  //  On a transition, advance all slices
  //
  Pds::TransitionId::Value service = (*n)->hdr().seq.service();
  if (service == Pds::TransitionId::L1Accept)
    ++_iCurEventGlobal;
  else
  {
    if (service == Pds::TransitionId::BeginCalibCycle)
    {
      ++_iCurCalib;
      _iCurCalibBaseEvt = _iCurEventGlobal;
    }

    for(std::list<XtcSlice*>::iterator it = _slices.begin();
        it != _slices.end();) {
      if (it != n && (*it)->skip()==End) {
        _doneSlices.push_back(*it);
        std::list<XtcSlice*>::iterator ee = it++;
        _slices.erase(ee);
      }
      else
        it++;
    }
  }

  Result r = (*n)->next(dg, pi64OffsetCur);
  if (r == End) {
    _doneSlices.push_back(*n);
    _slices.erase(n);
    if (_slices.size())
      r = OK;
  }

  if (piSlice != NULL)
    *piSlice = iNextSlice;

  //printf("<%d>", iNextSlice);//!!debug
  return r;
}

int XtcRun::jump(int calib, int jump, int& eventNum)
{
  if (jump == 0)
  {
    int iErrorAll     = 0;
    int iEventSum = 0;
    int iSlice        = 0;
    for(std::list<XtcSlice*>::iterator
      it =  _slices.begin();
      it != _slices.end();
      it++, iSlice++)
    {
      int iEventSlice;
      int iError =
        (*it)->jump( calib, 0, iEventSlice );

      if ( iError == 0 )
        iEventSum += iEventSlice-1;
      else
        iErrorAll++;
    }

    if ( iErrorAll != 0 )
    {
      eventNum = -1;
      return 1;
    }

    ++iEventSum;
    eventNum          = iEventSum;
    _iCurCalib        = calib;
    _iCurEventGlobal  = eventNum;
    _iCurCalibBaseEvt = eventNum;
    return 0;
  }

  struct TSliceJumpInfo
  {
    int iNumEvent;
    int iLowerBound;
    int iUpperBound;
    int iCurrentCalib;
    int iCurrentIndex;
  } infoList[ _slices.size() ];

  int iNumEventTotal = 0;
  int iSlice         = 0;
  for(std::list<XtcSlice*>::iterator
    it =  _slices.begin();
    it != _slices.end();
    it++, iSlice++)
  {
    infoList[iSlice].iLowerBound = 1;

    int iError = (*it)->numEventInCalib(calib, infoList[iSlice].iNumEvent);
    if (iError != 0)
    {
      printf("XtcRun::jump(): Cannot get Event# for Slice %d in Calib# %d\n", iSlice, calib);
      return 2;
    }

    infoList[iSlice].iUpperBound    = infoList[iSlice].iNumEvent;
    infoList[iSlice].iCurrentCalib  = -1;
    infoList[iSlice].iCurrentIndex  = -1;
    iNumEventTotal                  += infoList[iSlice].iNumEvent;
  }
  if (jump < 0 || jump > iNumEventTotal)
  {
    printf("XtcRun::jump(): Invalid Event# %d (Max # = %d)\n", jump, iNumEventTotal);
    return 3;
  }

  int iJumpAvg  = (int) ( (jump+_slices.size()-1) / _slices.size());
  if (iJumpAvg < 1)
    iJumpAvg = 1;
  else if (iJumpAvg > infoList[0].iUpperBound )
    iJumpAvg = infoList[0].iUpperBound;

  iSlice = 0;
  for(std::list<XtcSlice*>::iterator
    it =  _slices.begin();
    it != _slices.end();
    it++, iSlice++)
  {
    int iTestIndex = ( ( iSlice == 0 ) ? iJumpAvg :
      (infoList[iSlice].iLowerBound + infoList[iSlice].iUpperBound)/2 );

    for (; infoList[iSlice].iLowerBound <= infoList[iSlice].iUpperBound;
      iTestIndex = (infoList[iSlice].iLowerBound + infoList[iSlice].iUpperBound)/2)
    {
      infoList[iSlice].iCurrentIndex = iTestIndex;

      uint32_t uSeconds, uNanoseconds;
      (*it)->getTime(calib, iTestIndex, uSeconds, uNanoseconds);

      int iLocalEventSum  = 0;
      int iSliceQuery     = 0;
      for(std::list<XtcSlice*>::iterator
        itQuery =  _slices.begin();
        itQuery != _slices.end();
        itQuery++, iSliceQuery++)
      {
        if (iSliceQuery == iSlice)
        {
          iLocalEventSum += iTestIndex-1;
          continue;
        }

        int   iCalibQuery       = -1;
        int   iEventQuery       = -1;
        bool  bExactMatchQuery  = false;
        bool  bOvertimeSlice    = false;
        int iError = (*itQuery)->findTime(uSeconds, uNanoseconds, iCalibQuery, iEventQuery, bExactMatchQuery, bOvertimeSlice);
        if (iError != 0)
        {
          if (bOvertimeSlice)
          {
            iCalibQuery = calib;
            iEventQuery = infoList[iSliceQuery].iNumEvent+1;
          }
          else // Abnormal error
            continue;
        }

        if (iCalibQuery != calib)
        {
          if (iCalibQuery < calib)
          {
            printf("XtcRun::jump(): Inconsistent Calib# in different Slices: %d [Slice %d] and %d [Slice %d]\n",
              calib, iSlice, iCalibQuery, iSliceQuery );
            return 4;
          }
          else
          {
            iCalibQuery = calib;
            iEventQuery = infoList[iSliceQuery].iNumEvent+1;
          }
        }

        infoList[iSliceQuery].iCurrentIndex = iEventQuery;

        iLocalEventSum += iEventQuery-1;
      } // itQuery loop

      ++iLocalEventSum;
      if (iLocalEventSum < jump)
      {
        infoList[iSlice].iLowerBound = infoList[iSlice].iCurrentIndex+1;
        for (int iSliceUpdate = iSlice+1; iSliceUpdate < (int) _slices.size(); iSliceUpdate++ )
        {
          if (infoList[iSliceUpdate].iCurrentIndex>infoList[iSliceUpdate].iLowerBound)
            infoList[iSliceUpdate].iLowerBound = infoList[iSliceUpdate].iCurrentIndex;
        }
      }
      else if (iLocalEventSum > jump)
      {
        infoList[iSlice].iUpperBound = infoList[iSlice].iCurrentIndex-1;
        for (int iSliceUpdate = iSlice+1; iSliceUpdate < (int) _slices.size(); iSliceUpdate++ )
        {
          if (infoList[iSliceUpdate].iCurrentIndex-1<infoList[iSliceUpdate].iUpperBound)
            infoList[iSliceUpdate].iUpperBound = infoList[iSliceUpdate].iCurrentIndex-1;
        }
      }
      else // (iLocalEventSum == jump) : we have found the correct event
      {
        int iErrorAll = 0;

        int iGlobalEventSum = 0;
        int iSliceUpdate = 0;
        for(std::list<XtcSlice*>::iterator
          itUpdate =  _slices.begin();
          itUpdate != _slices.end();
          itUpdate++, iSliceUpdate++)
        {
          int iEventSlice;
          int iError =
            (*itUpdate)->jump( calib, infoList[iSliceUpdate].iCurrentIndex, iEventSlice, true );

          if ( iError == 0 )
            iGlobalEventSum += iEventSlice-1;
          else
          {
            printf("XtcRun::jump(): Jump failed in Slice %d for Calib# %d Event# %d\n",
              iSliceUpdate, calib, infoList[iSliceUpdate].iCurrentIndex);
            iErrorAll++;
          }
        }

        if ( iErrorAll != 0 )
        {
          eventNum = -1;
          return 5;
        }

        ++iGlobalEventSum;
        eventNum          = iGlobalEventSum;
        _iCurCalib        = calib;
        _iCurEventGlobal  = eventNum - 1;
        _iCurCalibBaseEvt = eventNum - jump;

        printf("XtcRun::jump(): Jumped to Calib# %d Event# %d (Global# %d)\n",
          calib, iLocalEventSum, iGlobalEventSum); //!!debug

        return 0;
      }

    } // binary search loop
  } // for each slice

  return 6;
}

int XtcRun::findTime(const char* sTime, int& iCalib, int& iEvent, bool& bExactMatch, bool& bOvertime)
{
  uint32_t uSeconds = 0, uNanoseconds = 0;
  int iError = Index::convertTimeStringToSeconds( sTime, uSeconds, uNanoseconds );

  if (iError != 0)
    return 1;

  return findTime(uSeconds, uNanoseconds, iCalib, iEvent, bExactMatch, bOvertime);
}

int XtcRun::findTime(uint32_t uSeconds, uint32_t uNanoseconds, int& iCalib, int& iEvent, bool& bExactMatch, bool& bOvertime)
{
  iCalib      = -1;
  iEvent      = 0;
  bExactMatch = false;

  int iErrorAll     = 0;
  int iOvertimeSum  = 0;
  int iCalibMin     = 0;

  struct TSliceFindTimeInfo
  {
    int iCalib;
    int iEvent;
  } infoList[ _slices.size() ];

  int iSlice = 0;
  for(std::list<XtcSlice*>::iterator
    it =  _slices.begin();
    it != _slices.end();
    it++, iSlice++)
  {
    int   iCalibSlice = -1, iEventSlice = -1;
    bool  bExactMatchSlice  = false;
    bool  bOvertimeSlice    = false;
    int iError =
      (*it)->findTime( uSeconds, uNanoseconds, iCalibSlice, iEventSlice, bExactMatchSlice, bOvertimeSlice );

    if ( iError != 0 )
    {
      if ( bOvertimeSlice )
      {
        int iNumCalib = 0;
        int iError = (*it)->numCalib(iNumCalib);
        if (iError != 0)
        {
          printf("XtcRun::findTime(): Cannot get Calib# for Slice %d\n", iSlice);
          return 1;
        }
        iCalibSlice = iNumCalib; // set current calib to be the last calib

        int iNumEventsInCalib = 0;
        iError = (*it)->numEventInCalib(iCalibSlice, iNumEventsInCalib);
        if (iError != 0)
        {
          printf("XtcRun::findTime(): Cannot get Event# for Slice %d in Calib# %d\n", iSlice, iCalibSlice);
          return 1;
        }
        iEventSlice = iNumEventsInCalib+1; // set current event to be the last event+1

        printf("XtcRun::findTime(): Overtime for Slice %d, set Calib# %d Event# %d\n",
          iSlice, iCalibSlice, iEventSlice );

        ++iOvertimeSum;
      }
      else
      {
        printf("XtcRun::findTime(): Failed in Slice %d\n", iSlice );
        iErrorAll++;
        continue;
      }
    }

    infoList[iSlice].iCalib = iCalibSlice;
    infoList[iSlice].iEvent = iEventSlice;

    if ( iSlice == 0 || iCalibSlice < iCalibMin )
      iCalibMin = iCalibSlice;

    if (bExactMatchSlice)
      bExactMatch = true;
  }

  if ( iErrorAll != 0 )
    return 1;

  if (iOvertimeSum >= (int) _slices.size())
  {
    bOvertime = true;
    return 2;
  }

  /*
   * Special case:
   *   The event with the input timestamp may become the last event of some slice,
   *   so it will not be found in other slices. In those slices, the returned calib# = real calib# + 1
   *
   *   Also in those slices, all events in the real calib# are earlier than the timestamp,
   *   so we add up all these slice-based event# to compute the real event#
   */
  int iEventSum = 0;
  iSlice = 0;
  for(std::list<XtcSlice*>::iterator
    it =  _slices.begin();
    it != _slices.end();
    it++, iSlice++)
  {
    if (infoList[iSlice].iCalib != iCalibMin)
    {
      int iNumEvent = -1;
      int iError = (*it)->numEventInCalib(iCalibMin, iNumEvent);
      if (iError != 0)
      {
        printf("XtcRun::findTime(): Cannot get Event# for Slice %d in Calib# %d\n", iSlice, iCalibMin);
        return 3;
      }
      infoList[iSlice].iEvent = iNumEvent+1;
    }

    iEventSum += infoList[iSlice].iEvent-1;
  }

  iCalib = iCalibMin;
  iEvent = iEventSum + 1;
  return 0;
}

int XtcRun::eventGlobalToSlice( int iGlobalEvent, std::vector<int>& lSliceEvent )
{
  lSliceEvent.assign(_slices.size(), 0);
  struct TSliceJumpInfo
  {
    int iNumEvent;
    int iLowerBound;
    int iUpperBound;
    int iCurrentIndex;
  } infoList[ _slices.size() ];

  int iNumEventTotal = 0;
  int iSlice         = 0;
  for(std::list<XtcSlice*>::iterator
    it =  _slices.begin();
    it != _slices.end();
    it++, iSlice++)
  {
    infoList[iSlice].iLowerBound = 1;

    int iError = (*it)->numTotalEvent(infoList[iSlice].iNumEvent);
    if (iError != 0)
    {
      printf("XtcRun::jump(): Cannot get Event# for Slice %d\n", iSlice);
      return 1;
    }

    infoList[iSlice].iUpperBound   = infoList[iSlice].iNumEvent;
    infoList[iSlice].iCurrentIndex = -1;
    iNumEventTotal                 += infoList[iSlice].iNumEvent;
  }
  if (iGlobalEvent < 0 || iGlobalEvent > iNumEventTotal)
  {
    printf("XtcRun::eventGlobalToSlice(): Invalid Event# %d (Max # = %d)\n", iGlobalEvent, iNumEventTotal);
    return 2;
  }

  int iEventAvg  = (int) ( (iGlobalEvent+_slices.size()-1) / _slices.size());
  if (iEventAvg < 1)
    iEventAvg = 1;
  else if (iEventAvg > infoList[0].iUpperBound )
    iEventAvg = infoList[0].iUpperBound;

  //printf("XtcRun::eventGlobalToSlice(): Input global event# %d iEventAvg %d\n", iGlobalEvent, iEventAvg); //!!debug

  iSlice = 0;
  for(std::list<XtcSlice*>::iterator
    it =  _slices.begin();
    it != _slices.end();
    it++, iSlice++)
  {
    int iTestIndex = ( ( iSlice == 0 ) ? iEventAvg :
      (infoList[iSlice].iLowerBound + infoList[iSlice].iUpperBound)/2 );

    for (; infoList[iSlice].iLowerBound <= infoList[iSlice].iUpperBound;
      iTestIndex = (infoList[iSlice].iLowerBound + infoList[iSlice].iUpperBound)/2)
    {
      infoList[iSlice].iCurrentIndex = iTestIndex;

      //printf("XtcRun::eventGlobalToSlice(): Testing slice %d index %d\n", iSlice, iTestIndex); //!!debug

      uint32_t uSeconds, uNanoseconds;
      int iError = (*it)->getTimeGlobal(iTestIndex, uSeconds, uNanoseconds);
      if (iError != 0)
      {
        printf("XtcRun::eventGlobalToSlice(): Cannot get time for Slice# %d Event# %d\n", iSlice, iTestIndex);
        return 3;
      }

      int iLocalEventSum  = 0;
      int iSliceQuery     = 0;
      for(std::list<XtcSlice*>::iterator
        itQuery =  _slices.begin();
        itQuery != _slices.end();
        itQuery++, iSliceQuery++)
      {
        if (iSliceQuery == iSlice)
        {
          iLocalEventSum += iTestIndex-1;
          continue;
        }

        int   iSliceEventQuery  = -1;
        bool  bExactMatchQuery  = false;
        bool  bOvertimeSlice    = false;
        int iError = (*itQuery)->findTimeGlobal(uSeconds, uNanoseconds, iSliceEventQuery, bExactMatchQuery, bOvertimeSlice);
        if (iError != 0)
        {
          if (bOvertimeSlice)
          {
            iSliceEventQuery = infoList[iSliceQuery].iNumEvent+1;
          }
          else // Abnormal error
            continue;
        }

        //printf("XtcRun::eventGlobalToSlice(): slice %d has %d earlier L1s\n", iSliceQuery, iSliceEventQuery); //!!debug
        infoList[iSliceQuery].iCurrentIndex = iSliceEventQuery;
        iLocalEventSum += iSliceEventQuery-1;
      } // itQuery loop

      ++iLocalEventSum;

      //printf("XtcRun::eventGlobalToSlice(): resultant global event# %d\n", iLocalEventSum); //!!debug
      if (iLocalEventSum < iGlobalEvent)
      {
        infoList[iSlice].iLowerBound = infoList[iSlice].iCurrentIndex+1;
        for (int iSliceUpdate = iSlice+1; iSliceUpdate < (int) _slices.size(); iSliceUpdate++ )
        {
          if (infoList[iSliceUpdate].iCurrentIndex>infoList[iSliceUpdate].iLowerBound)
            infoList[iSliceUpdate].iLowerBound = infoList[iSliceUpdate].iCurrentIndex;
        }
      }
      else if (iLocalEventSum > iGlobalEvent)
      {
        infoList[iSlice].iUpperBound = infoList[iSlice].iCurrentIndex-1;
        for (int iSliceUpdate = iSlice+1; iSliceUpdate < (int) _slices.size(); iSliceUpdate++ )
        {
          if (infoList[iSliceUpdate].iCurrentIndex-1<infoList[iSliceUpdate].iUpperBound)
            infoList[iSliceUpdate].iUpperBound = infoList[iSliceUpdate].iCurrentIndex-1;
        }
      }
      else // (iLocalEventSum == iGlobalEvent) : we have found the correct event
      {
        for(int iSliceUpdate = 0;
            iSliceUpdate < (int) _slices.size();
            ++iSliceUpdate)
        {
          lSliceEvent[iSliceUpdate] = infoList[iSliceUpdate].iCurrentIndex;
          //printf("XtcRun::eventGlobalToSlice(): Convert Global Event# %d to Slice# %d Event# %d\n",
          //  iGlobalEvent, iSliceUpdate, lSliceEvent[iSliceUpdate]); //!!debug
        }
        return 0;
      }

    } // binary search loop
  } // for each slice

  return 4;
}

int XtcRun::findNextFiducial(uint32_t uFiducialSearch, int iFidFromEvent, int& iCalib, int& iEvent)
{
  std::vector<int> lSliceEvent;
  int iError = eventGlobalToSlice(iFidFromEvent, lSliceEvent);
  if (iError != 0)
  {
    printf("XtcRun::findNextFiducial(): Failed to convert Global Event# %d to slice events\n", iFidFromEvent);
    return 1;
  }

  int      iSlice       = 0;
  uint32_t uSeconds     = 0;
  uint32_t uNanoseconds = 0;
  for(std::list<XtcSlice*>::iterator
    it =  _slices.begin();
    it != _slices.end();
    it++, iSlice++)
  {
    XtcSlice* pSlice      = *it;
    int       iSliceEvent = -1;
    int iError            = pSlice->findNextFiducial(uFiducialSearch, lSliceEvent[iSlice], iSliceEvent);

    if (iError == 0) // The event with the input timestamp has been found
    {
      uint32_t  uSecondsCur     = 0;
      uint32_t  uNanosecondsCur = 0;
      int iError = pSlice->getTimeGlobal(iSliceEvent, uSecondsCur, uNanosecondsCur);
      if (iError != 0)
      {
        printf("XtcRun::findNextFiducial(): Cannot get time for Slice# %d\n", iSlice);
        return 2;
      }

      if
      (
        (uSeconds == 0 && uNanoseconds == 0) ||
         uSecondsCur < uSeconds ||
        (uSecondsCur == uSeconds && uNanosecondsCur < uNanoseconds)
      )
      {
        uSeconds     = uSecondsCur;
        uNanoseconds = uNanosecondsCur;
      }
    }
  }

  if (uSeconds == 0 && uNanoseconds == 0)
    return 3;

  bool bExactMatch = false;
  bool bOvertime   = false;
  iError = findTime(uSeconds, uNanoseconds, iCalib, iEvent, bExactMatch, bOvertime);
  if (iError != 0)
  {
    printf("XtcRun::findNextFiducial(): Failed to find event with the searched timestamp of input fiducial\n");
    return 4;
  }
  return 0;
}

int XtcRun::numCalib(int& iNumCalib)
{
  if (!_startAndEndValid)
  {
    iNumCalib = -1;
    return 1;
  }
  iNumCalib = _iNumTotalCalib;
  return 0;
}

int XtcRun::numEventInCalib(int calib, int& iNumEvents)
{
  if (!_startAndEndValid)
  {
    iNumEvents = -1;
    return 1;
  }

  if (calib <= 0 || calib > _iNumTotalCalib)
  {
    iNumEvents = 0;
    return 2;
  }

  iNumEvents =  _vNumEventsInCalib[calib];
  return 0;
}

int XtcRun::curEventId(int& eventGlobal, int& calib, int& eventInCalib)
{
  eventGlobal   = _iCurEventGlobal;
  calib         = _iCurCalib;
  eventInCalib  = _iCurEventGlobal - _iCurCalibBaseEvt;
  return 0;
}

int XtcRun::eventGlobalToCalib(int eventGlobal, int& calib, int& eventInCalib)
{
  calib = -1;
  eventInCalib = -1;

  if (eventGlobal <= 0 || eventGlobal > _iNumTotalEvent)
    return 1;

  int eventRemaining = eventGlobal;
  int iCalibTest = 1;
  for (; iCalibTest <= _iNumTotalCalib; ++iCalibTest)
  {
    if (eventRemaining <= _vNumEventsInCalib[iCalibTest])
    {
      calib = iCalibTest;
      eventInCalib = eventRemaining;
      return 0;
    }

    eventRemaining -= _vNumEventsInCalib[iCalibTest];
  }

  return 2;
}

int XtcRun::eventCalibToGlobal(int calib, int eventInCalib, int& eventGlobal)
{
  eventGlobal = -1;
  if (calib <= 0 || calib > _iNumTotalCalib)
    return 1;

  if (eventInCalib <= 0 || eventInCalib > _vNumEventsInCalib[calib])
    return 2;

  int eventSum = eventInCalib;
  for (int iCalibTest = 1; iCalibTest < calib; ++iCalibTest)
    eventSum += _vNumEventsInCalib[iCalibTest];

  eventGlobal = eventSum;
  return 0;
}

int XtcRun::nextCalibEventId(int iNumCalibAfter, int iNumEventAfter, bool bResetEventInCalib, int& eventGlobal, int& calib, int& eventInCalib)
{
  int iCalibFinal         = _iCurCalib + iNumCalibAfter;
  int iEventInCalibFinal;
  if (bResetEventInCalib && iNumCalibAfter != 0)
    iEventInCalibFinal = iNumEventAfter;
  else
    iEventInCalibFinal = _iCurEventGlobal - _iCurCalibBaseEvt + iNumEventAfter;
  if (iCalibFinal <= 0 || iCalibFinal > _iNumTotalCalib)
    return 1;
  if (iEventInCalibFinal <= 0 || iEventInCalibFinal > _vNumEventsInCalib[iCalibFinal])
    return 2;

  int iEventGlobalFinal = -1;
  int iError = eventCalibToGlobal(iCalibFinal, iEventInCalibFinal, iEventGlobalFinal);
  if (iError != 0)
    return 3;

  eventGlobal   = iEventGlobalFinal;
  calib         = iCalibFinal;
  eventInCalib  = iEventInCalibFinal;
  return 0;
}

int XtcRun::nextGlobalEventId(int iNumEventAfter, int& eventGlobal, int& calib, int& eventInCalib)
{
  int iEventGlobalFinal = _iCurEventGlobal + iNumEventAfter;
  if (eventGlobal <= 0 || eventGlobal > _iNumTotalEvent)
    return 1;

  int iCalibFinal         = -1;
  int iEventInCalibFinal  = -1;
  int iError = eventGlobalToCalib(iEventGlobalFinal, iCalibFinal, iEventInCalibFinal);
  if (iError != 0)
    return 2;

  eventGlobal   = iEventGlobalFinal;
  calib         = iCalibFinal;
  eventInCalib  = iEventInCalibFinal;
  return 0;
}

} // namespace Ana
} // namespace Pds

