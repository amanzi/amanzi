#ifndef __Geochemistry_hpp__
#define __Geochemistry_hpp__

#include "Species.hpp"
#include "SecondarySpecies.hpp"
#include "Speciation.hpp"
#include "MineralReaction.hpp"
#include "GasExchange.hpp"
#include "IonExchange.hpp"
#include "SurfaceComplexation.hpp"


// Driver class for evalating geochemical related processes at a
// single computational node

class Geochemistry : {

public:
  Geochemistry();
  ~Geochemistry();

  evaluateChemistry(std::vector<double> current_primary_species);

private:

  std::vector<double> primarySpecies_;
  std::vector<SecondarySpecies> secondarySpecies_;
  std::vector<ActivityCoefficient*> activityCoefficients_;

  Speciation speciation_;
  std::vector<MineralReaction*> mineralReactions_;
  std::vector<GasExchange*> gasReactions_;
  std::vector<IonExchange*> ionExchangeReactions_;
  std::vector<SurfaceComplexation*> surfaceComplexationReactions_;

};

#endif
