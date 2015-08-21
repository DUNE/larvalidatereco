#ifndef LAR_VALIDATION_CALCULATOR_H
#define LAR_VALIDATION_CALCULATOR_H

#include "/data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArValidation.h"
#include "/data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/framework/CalculatorBase.h"

namespace lar_valrec{

  ///Calculates validation variables for event, for now
  ///mainly completeness and purity of reconstructed tracks
  class LArValidationCalculator : public CalculatorBase{

  public:

    LArValidationCalculator(){}

    ///Implement Calculate method declared in CalculatorBase. This calls all the private
    ///methods which actually fill the various bits of output.
    virtual void Calculate(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper);

    virtual ~LArValidationCalculator(){}

private:
  };
}

#endif //#ifndef LAR_VALIDATION_CALCULATOR_H
