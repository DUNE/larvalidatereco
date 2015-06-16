#include "VarHelper.h"
#include "TROOT.h"
#include "TCollection.h"

namespace lar_valrec{

  VarHelper::VarHelper(TTree* outputTree) : fTree(outputTree){}

  std::map<std::string,TObject* >& VarHelper::GetOutputs(){
    //Initialize static container on first use
    static std::map<std::string,TObject* > outputs;
    return outputs;
  }

  std::map<std::string,CalculatorBase*>& VarHelper::GetCalculators(){
    //Initialize static container on first use
    static std::map<std::string,CalculatorBase*> calculators;
    return calculators;
  }

  void VarHelper::CalculateVars(const EventHelper& eh){
    for(auto varIter=fOutputPipeline.begin();varIter!=fOutputPipeline.end();++varIter){
      GetCalculators()[*varIter]->Calculate(GetOutputs()[*varIter],*this,eh);
   }

    fTree->Fill();
}

  void VarHelper::AddOutput(std::string name){
    fOutputPipeline.push_back(name);
    TObject*& objPtr=GetOutputs()[name];
    fTree->Branch(name.c_str(),GetOutputs()[name]->ClassName(),(void**)(&objPtr));
  }
  
VarHelper::~VarHelper(){}

}
