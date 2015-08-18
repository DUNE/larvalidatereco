/**
 *  @file   larvalidaterecon/interface/LArValidateReco_module.cc
 *
 *  @brief  Provide proper module filename for LArValidateRecon class
 *
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"

// Local includes
#include "LArValidateRecon.h"

// std includes
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_valrec
{

/**
 *  @brief  LArValRec class
 */
class LArValRec : public LArValidateRecon
{
public: 
    LArValRec(fhicl::ParameterSet const &pset);

    virtual ~LArValRec();

private:

};

DEFINE_ART_MODULE(LArValRec)

} // namespace lar_valrec

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

// LArSoft includes
#include "Geometry/Geometry.h"

namespace lar_valrec {

LArValRec::LArValRec(fhicl::ParameterSet const &pset) : LArValidateRecon(pset)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArValRec::~LArValRec()
{
}

} // namespace lar_valrec
