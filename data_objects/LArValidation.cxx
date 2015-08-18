#include "LArValidation.h"

ClassImp(LArValidation)

  void LArValidation::Clear()
  {
    EventEventRun =        kUnassigned;
    EventID       =        kUnassigned;
    EventSubRun   =        kUnassigned;

    NTracks=kUnassigned;
    TrackTrajID.clear();
    TrackCompleteness.clear();
    TrackPurity.clear();
    TrackStitched.clear();
  }

  LArValidation::LArValidation():kUnassigned(0xDEADBEEF){this->Clear();}
  LArValidation::~LArValidation(){}

