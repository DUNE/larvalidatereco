#include "LArPID.h"

ClassImp(LArPID)

  void LArPID::Clear()
  {
    EventEventRun =        kUnassigned;
    EventID       =        kUnassigned;
    EventSubRun   =        kUnassigned;

    EigenValues  .clear();
    EigenVectors .clear();
    Covariance   .clear();
    PCAHitsSpacePoints .clear();

    NTracks=kUnassigned;
  }

  LArPID::LArPID():kUnassigned(0xDEADBEEF){this->Clear();}
  LArPID::~LArPID(){}

