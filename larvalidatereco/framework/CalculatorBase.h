#ifndef CALCULATOR_BASE_H
#define CALCULATOR_BASE_H


namespace lar_valrec{

  class VarHelper;
  class EventHelper;

  ///Base class to allow VarHelper to store concrete calculators for the various outputs
  class CalculatorBase{

  protected:

    CalculatorBase(){}

    virtual ~CalculatorBase(){}

  public:

    ///Will be overriden in derived class to actually fill the elements of the tObjectPointer.
    ///Note this will always require a downcast of the TObject* to the concrete output class.
    virtual void Calculate(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper)=0;

  };
}
  
#endif //#ifndef CALCULATOR_BASE_H
