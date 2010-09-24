#ifndef __Geochemistry_hpp__
#define __Geochemistry_hpp__

#include "Species.hpp"
#include "AqueousEquilibriumComplex.hpp"
//#include "MineralReaction.hpp"
//#include "GasExchange.hpp"
//#include "IonExchange.hpp"
//#include "SurfaceComplexation.hpp"
#include "Block.hpp"
#include "LU.hpp"

#include <vector>

// Driver class for evalating geochemical related processes at a
// single computational node

class Geochemistry {

public:
  Geochemistry();
  ~Geochemistry();

  // inheriting classes setup the species, etc
  virtual void setup(std::vector<double> *total);

};

#endif // __Geochemistry_hpp__
