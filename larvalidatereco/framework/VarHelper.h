/**
 *  @file   larvalidatereco/framework/VarHelper.h
 *
 *  @brief  Helper class to register outputs and invoke their calculation
 *
 */

#ifndef VARHELPER_H
#define VARHELPER_H 1

// Framework Includes
#include "larvalidatereco/interface/PandoraTypedefs.h"

#include "TTree.h"
#include "CalculatorBase.h"

namespace lar_valrec
{

  class EventHelper;

/**
 *  @brief  VarHelper class
 */
class VarHelper
{
public:
  
  ///Initialize output tree
  VarHelper(TTree* outputTree);
  
  ///Register an output type to the framework.
  ///Templated since we need to know the concrete class to create a TBranch
  template <class T, class U> static void RegisterOutput(std::string name){
    GetOutputs()[name]=new T;
    GetCalculators()[name]=new U;
  }
  
  ///Tell the code that an output will actually be put into the tree and calculated.
  ///Outputs will be calculated in the same order as they are added so calculations
  ///can use previous results.
  void AddOutput(std::string name);
  
  ///Called each event; invokes calculation of each output and filling of the tree
  void CalculateVars(const EventHelper& eh);

  ///Get list of outputs registered to the framework (read or write)
  static std::map<std::string,TObject*>& GetOutputs();
  
  ~VarHelper();
  
 private:
  
  ///Output tree
  TTree* fTree;

  ///Get list of calculators registered to the framework (read or write)
  static std::map<std::string,CalculatorBase*>& GetCalculators();

  ///Ordered list of outputs to actually be calculated
  std::vector<std::string> fOutputPipeline;
};

} // namespace lar_valrec


#endif // #ifndef VARHELPER_H
