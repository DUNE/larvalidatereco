#ifndef LAR_PID_CALCULATOR_H
#define LAR_PID_CALCULATOR_H

#include "data_objects/LArPID.h"
#include "../framework/CalculatorBase.h"

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

private:
  };
}

#endif //#ifndef LAR_PID_CALCULATOR_H
