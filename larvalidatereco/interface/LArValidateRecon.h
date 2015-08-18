/**
 *  @file   larvalidatereco/interface/LArValidateRecon.h
 *
 *  @brief  Guts of an analyser module to produce validation/DST information from a reconstructed event
 *
 */

#ifndef LAR_VALIDATERECON_H
#define LAR_VALIDATERECON_H 1

// Framework Includes
#include "art/Framework/Core/EDAnalyzer.h"

#include "EventHelper.h"
#include "larvalidatereco/framework/VarHelper.h"
#include "TH1F.h"
#include "TTree.h"

namespace lar_valrec
{

  class VarHelper;
  class EventHelper;

/**
 *  @brief  LArValidateRecon class
 */
class LArValidateRecon : public art::EDAnalyzer
{
public:
    LArValidateRecon(fhicl::ParameterSet const &pset);

    virtual ~LArValidateRecon();

    virtual void beginJob();
    //    virtual void endJob();

    void analyze(const art::Event &evt);

 private:

    ///Output tree
    TTree* fTree;

    ///Object to calculate outputs and write tree
    VarHelper* fVarHelper;
    
    ///Object to read in data from event
    EventHelper* fEventHelper;

    ///Input fcl parameters
    fhicl::ParameterSet fPset;

    ///Name of output tree (from fcl)
    std::string fTreeName;

    ///List of variables to calculate (from fcl)
    std::vector<std::string> fVariables;
};

} // namespace lar_valrec

#endif // #ifndef LAR_ValidateRecon_H
