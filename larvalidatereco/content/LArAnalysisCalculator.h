#ifndef LAR_ANALYSIS_CALCULATOR_H
#define LAR_ANALYSIS_CALCULATOR_H

#include "data_objects/LArAnalysis.h"
#include "larvalidatereco/framework/CalculatorBase.h"

namespace lar_valrec{

  ///Class to fill the basic DST information in the LArAnalysis IO object
  class LArAnalysisCalculator : public CalculatorBase{

  public:

    LArAnalysisCalculator(){}

    ///Implement Calculate method declared in CalculatorBase. This calls all the private
    ///methods which actually fill the various bits of output.
    virtual void Calculate(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper);

    virtual ~LArAnalysisCalculator(){}

private:

    ///Fill global metadata (run/event number, electric field, trigger times...)
    static void FillEventMetadata(LArAnalysis* outputPtr,const EventHelper& evHelper);

    ///Fill Monte Carlo trajectory data
    static void FillEventMCTraj(LArAnalysis* outputPtr,const EventHelper& evHelper);

    ///Fill Monte Carlo vertex data
    static void FillEventMCVertices(LArAnalysis* outputPtr,const EventHelper& evHelper);

    ///Fill reconstructed tracks
    static void FillEventRecoTracks(LArAnalysis* outputPtr,const EventHelper& evHelper);
    
    //Fill reconstructed showers
    static void FillEventRecoShowers(LArAnalysis* outputPtr,const EventHelper& evHelper);

    ///Fill reconstructed clusters
    static void FillEventRecoClusters(LArAnalysis* outputPtr,const EventHelper& evHelper);

    ///Fill reconstructed hits
    static void FillEventHits(LArAnalysis* outputPtr,const EventHelper& evHelper);

    ///Helper method to calculate track length inside the detector from the track points
    static double GetTrackLength(std::vector<TVector3>& points, int TrackIndex=-1);

    ///Helper method which uses the geometry service to figure out if a coordinate
    ///is in an active region of the detector.
    static bool IsInActiveRegion(const TVector3& position);
  };
}

#endif //#ifndef LAR_ANALYSIS_CALCULATOR_H
