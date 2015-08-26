#ifndef LAR_PID_CALCULATOR_H
#define LAR_PID_CALCULATOR_H

#include "data_objects/LArPID.h"
#include "larvalidatereco/framework/CalculatorBase.h"
#include "AnalysisAlg/CalorimetryAlg.h"

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
 
    double GetTrackLength(std::vector<TVector3>& points,int trackIndex);

    void FillEventMetadata(LArPID* outputPtr,const EventHelper& evHelper);

    void FillEventMCTraj(LArPID* outputPtr,const EventHelper& evHelper);

    void FillEventdEdx(LArPID* outputPtr,const EventHelper& evHelper);

  private:
  };
}

#endif //#ifndef LAR_PID_CALCULATOR_H
