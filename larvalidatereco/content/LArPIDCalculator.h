#ifndef LAR_PID_CALCULATOR_H
#define LAR_PID_CALCULATOR_H

#include "data_objects/LArPID.h"
#include "larvalidatereco/framework/CalculatorBase.h"
#include "AnalysisAlg/CalorimetryAlg.h"
#include "TGraph2D.h"
#include "TVector3.h"

namespace lar_valrec{

  ///Class to fill the basic DST information in the LArAnalysis IO object
  class LArPIDCalculator : public CalculatorBase{

  public:

    LArPIDCalculator(){}

    ///Implement Calculate method declared in CalculatorBase. This calls all the private
    ///methods which actually fill the various bits of output.
    virtual void Calculate(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper);

    //static void FillClusterPCA(LArPID* outputPtr,const EventHelper& evHelper);

    virtual ~LArPIDCalculator(){}

    bool IsInActiveRegion(const TVector3& position);

    void FillEventPID(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper);

    double CalcDistSqPointLine(const TVector3& point, const TVector3& linePoint, const TVector3& lineDirection);

    TVector3 CalcNearestPointOnLine(const TVector3& point);

    double CalcResRangeFraction(const TVector3& point, const TVector3& trackStart, const TVector3& trackEnd);
 
    //int FitTrack(TGraph2D *gHitsSpacepoints);
    int FitTrack(TH1D *hHitsSpacepoints);

    double GetTrackLength(std::vector<TVector3>& points,int trackIndex);

    void FillEventMetadata(LArPID* outputPtr,const EventHelper& evHelper);

    void FillEventMCTraj(LArPID* outputPtr,const EventHelper& evHelper);

    void FillEventTracksAndHits(LArPID* outputPtr,const EventHelper& evHelper);

  private:

    bool trackFitMade;

    TVector3 fittedTrackPoint;
    TVector3 fittedTrackVector;

  };

  //This function calculates the sum of squared residuals between the track hits and a straight line fitted to the track.
  //This is the quantity minimised during the fit of the track.
  //This function has to be outside the class so that it can be passed to SetFCN().
  void CalcSumSqResidual(Int_t &npar, Double_t *, Double_t &f, Double_t *par, Int_t iflag);

}

#endif //#ifndef LAR_PID_CALCULATOR_H
