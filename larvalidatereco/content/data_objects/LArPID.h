
#ifndef __LAR_PID__
#define __LAR_PID__

#include "TObject.h"

///
/// The LArAnalysis class definintion.
/// This object will constitute a single TBranch.
///
class LArPID : public TObject{

///
/// Abbreviations
///
public:

  const int               kUnassigned;           ///< unassigned parameters initialised to this

  ///
  /// Event Bookeeping Information
  ///
  int                     EventEventRun;         ///< run number
  int                     EventID;               ///< event ID
  int                     EventSubRun;           ///< sub run number

  int NTracks;
  std::vector<int> TrackID;
  ///
  /// Reset parameter values and vector containers
  ///
  virtual void Clear();

  ///
  /// Default constructor/destructor
  ///
  LArPID();
  virtual ~LArPID();
  ClassDef(LArPID,1);
};

#endif