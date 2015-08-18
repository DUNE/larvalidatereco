
#ifndef __LAR_VALIDATION__
#define __LAR_VALIDATION__

#include "TObject.h"

///
/// The LArAnalysis class definintion.
/// This object will constitute a single TBranch.
///
class LArValidation : public TObject{

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
  std::vector<int> TrackTrajID;
  std::vector<double> TrackCompleteness;
  std::vector<double> TrackPurity;
  std::vector<int> TrackStitched;
  ///
  /// Reset parameter values and vector containers
  ///
  virtual void Clear();

  ///
  /// Default constructor/destructor
  ///
  LArValidation();
  virtual ~LArValidation();
  ClassDef(LArValidation,2);
};

#endif
