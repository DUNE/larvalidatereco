/**
 *  @file   larvalidaterecon/interface/LArValidateReco_module.cc
 *
 *  @brief  Provide proper module filename for LArValidateRecon class
 *
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "PandoraTypedefs.h"
#include "LArPandoraInterface/LArPandoraHelper.h"
#include "Utilities/AssociationUtil.h"
// Local includes

// std includes
#include <string>

#include "TTree.h"
#include "TFile.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_valrec
{

/**
 *  @brief  LArValRec class
 */
class DumpGENIEParticles : public art::EDAnalyzer
{
public: 
  DumpGENIEParticles(fhicl::ParameterSet const &pset);
  
  virtual ~DumpGENIEParticles(){}
  
  virtual void beginJob();

  virtual void endJob();
  
  virtual void analyze(const art::Event &evt);

private:

  std::string fGENIEModuleName;
  TFile* fOutFile;
  TTree* fTree;
  const simb::MCTruth* fBranchPtr;
};

DEFINE_ART_MODULE(DumpGENIEParticles)

} // namespace lar_valrec

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows


namespace lar_valrec {

  DumpGENIEParticles::DumpGENIEParticles(fhicl::ParameterSet const &pset) :
    art::EDAnalyzer(pset)
  {
  }

//------------------------------------------------------------------------------------------------------------------------------------------

void DumpGENIEParticles::beginJob()
{
  fOutFile=new TFile("particles.root","RECREATE");
  fTree=new TTree("Particles","Particles");
  fBranchPtr=0;
  fTree->Branch("MCTruth",&fBranchPtr);
}

void DumpGENIEParticles::endJob()
{
  fOutFile->Write();
  fOutFile->Close();
}


void DumpGENIEParticles::analyze(const art::Event& evt)
{
  art::Handle<std::vector<simb::MCTruth> > events;
  evt.getByLabel("generator",events);
  
  for(auto eIter=events->begin();eIter!=events->end();++eIter){
    for(auto iPart=0;iPart<eIter->NParticles();++iPart){
      const simb::MCParticle& part(eIter->GetParticle(iPart));

      //Need to change particle to a primary with no mother. Only
      //way is to build a new one from scratch.
      simb::MCParticle outPart(0,               //trackId
			       part.PdgCode(),  //PDG code
			       "primary",       //process producing particle
			       -1);             //TrackId of mother (i.e. no mother)

      for(unsigned int i=0;i<part.NumberTrajectoryPoints();++i){
	outPart.AddTrajectoryPoint(part.Position(i),part.Momentum(i));
      }
      simb::MCTruth mcTruth;
      mcTruth.Add(outPart);
      mcTruth.SetOrigin(simb::kSingleParticle);
      fBranchPtr=&mcTruth;
      fTree->Fill();
    }
  }
}

} // namespace lar_valrec
