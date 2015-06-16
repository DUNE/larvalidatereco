#include "TDatabasePDG.h"

#include "Geometry/TPCGeo.h"
#include "Geometry/Geometry.h"
#include "larvalidatereco/framework/VarHelper.h"
#include "larvalidatereco/interface/EventHelper.h"
#include "LArPIDCalculator.h"

#include "TGraph.h"
#include "TFile.h"

namespace lar_valrec{

  class LArPID_register{
  public:
    LArPID_register(){
        VarHelper::RegisterOutput<LArPID,LArPIDCalculator>("PID");
    }
  };
  static LArPID_register __LArPIDRegister;

  void LArPIDCalculator::Calculate(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper);

  static void LArPIDCalculator::FillClusterPCA(LArPID* outputPtr,const EventHelper& evHelper);
