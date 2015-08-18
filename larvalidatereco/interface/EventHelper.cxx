/**
 *  @file   larvalidatereco/interface/EventHelper.cxx
 *
 *  @brief  Base producer module for recon validation information
 *
 */

#include "EventHelper.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"

namespace lar_valrec
{

  void EventHelper::LoadEvent(const art::Event& evt){

    art::ServiceHandle<util::LArProperties> larProp;
    art::ServiceHandle<util::DetectorProperties> detProp;
    fRun=evt.run();
    fSubRun=evt.subRun();           ///< sub run number
    fEventID=evt.id().event();

    fEField=larProp->Efield();
    fEventT0=detProp->TriggerOffset();

    //Work for derived class to do.
    this->LoadEvent_(evt);
  }

} // namespace lar_valrec

