#include "LArTracking.h"

ClassImp(LArTracking)

  void LArTracking::Clear()
  {

  NTrajMC=kUnassigned;
  TrajIDMC.clear();
  TrajPointsMC.clear();

  NTracks=kUnassigned;
  TrackID.clear();
  TrackPoints.clear();

  }

  LArTracking::LArTracking():kUnassigned(0xDEADBEEF){this->Clear();}
  LArTracking::~LArTracking(){}

