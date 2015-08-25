#ifndef LAR_TRACKING_CALCULATOR_H
#define LAR_TRACKING_CALCULATOR_H

#include "data_objects/LArTracking.h"
#include "larvalidatereco/framework/CalculatorBase.h"

namespace lar_valrec{

  ///Class to fill the basic DST information in the LArAnalysis IO object
  class LArTrackingCalculator : public CalculatorBase{

  public:

    LArTrackingCalculator(){}

    ///Implement Calculate method declared in CalculatorBase. This calls all the private
    ///methods which actually fill the various bits of output.
    virtual void Calculate(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper);

    virtual ~LArTrackingCalculator(){}
  };
}

#endif //#ifndef LAR_TRACKING_CALCULATOR_H
