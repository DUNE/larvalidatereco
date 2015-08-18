#include "LArPID.h"

ClassImp(LArPID)

  void LArPID::Clear()
  {
    EventEventRun =        kUnassigned;
    EventID       =        kUnassigned;
    EventSubRun   =        kUnassigned;

    NTracks=kUnassigned;
  }

  LArPID::LArPID():kUnassigned(0xDEADBEEF){this->Clear();}
  LArPID::~LArPID(){}

