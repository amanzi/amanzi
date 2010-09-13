#ifndef __Speciation_Process_hpp__
#define __Speciation_Process_hpp__

#include "Teuchos_RCP.hpp"
#include "State.hpp"
#include "Chemistry_State.hpp"


// class for speciation process

class SpeciationProcess : public ReactionProcess {
public:
  SpeciationProcess();

  ~SpeciationProcess();

  void basisSwitch(void);
  void directSubstitution(void);
  std::vector<SecondarySpecies> speciate(std::vector<PrimarySpecies>& primary);

private:


};

#endif
