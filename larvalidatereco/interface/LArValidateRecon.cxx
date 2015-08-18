// Framework includes
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/search_path.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/TimeService.h"

#include "SimpleTypesAndConstants/RawTypes.h" // raw::TDCtick_t
#include "SimulationBase/MCTruth.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/PFParticle.h"

// System includes
#include <iostream>
#include <limits>
#include <algorithm>

#include "LArValidateRecon.h"
#include "LArPandoraEventHelper.h"

namespace lar_valrec {

  LArValidateRecon::LArValidateRecon(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset), fPset(pset)
{
  fVariables = pset.get<std::vector<std::string>>("Variables");
  fTreeName = pset.get<std::string>("TreeName");
}

void LArValidateRecon::analyze(const art::Event &evt)
{
  fEventHelper->LoadEvent(evt);
  fVarHelper->CalculateVars(*fEventHelper);
}

LArValidateRecon::~LArValidateRecon()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArValidateRecon::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>(fTreeName.c_str(),fTreeName.c_str());
  fEventHelper=new LArPandoraEventHelper(fPset);
  fVarHelper=new VarHelper(fTree);
  for(auto var=fVariables.begin();var!=fVariables.end();++var){
    fVarHelper->AddOutput(*var);
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

}
