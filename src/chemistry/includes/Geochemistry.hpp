#ifndef __Geochemistry_hpp__
#define __Geochemistry_hpp__

#include "Teuchos_RCP.hpp"
#include "Chemistry_PK.hpp"

#include "Species.hpp"
#include "PrimarySpecies.hpp"
#include "SecondarySpecies.hpp"
#include "Speciation.hpp"
#include "MineralReaction.hpp"
#include "GasExchange.hpp"
#include "IonExchange.hpp"
#include "SurfaceComplexation.hpp"


// Chemistry Process Kernel Interface

class Geochemistry : public Chemistry_PK {

public:
  Geochemistry(Teuchos::RCP<Chemistry_State> CS_);
  ~Geochemistry();



private:

  std::vector<PrimarySpecies> primarySpecies_;
  std::vector<PrimarySpecies> secondarySpecies_;
  std::vector<ActivityCoefficient> activityCoefficients_;

  Speciation speciation_;
  std::vector<MineralReaction> mineralReactions_;
  std::vector<GasExchange> gasReactions_;
  std::vector<IonExchange> ionExchangeReactions_;
  std::vector<SurfaceComplexation> surfaceComplexationReactions_;

};

#endif
