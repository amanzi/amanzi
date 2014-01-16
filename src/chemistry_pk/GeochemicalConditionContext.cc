/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi Chemistry

License: see COPYRIGHT
Author: Jeffrey Johnson

This is a point of contact for the chemistry engine exposed by Alquimia to 
the rest of Amanzi--it provides the ability to enforce geochemical conditions 
and to integrate reactions given a chemical configuration.

 ------------------------------------------------------------------------- */

#include "GeochemicalConditionContext.hh"

namespace Amanzi {
namespace AmanziChemistry {

// This subclass of Amanzi::Function provides an interface by which a geochemical condition can be 
// enforced on a given species.
class GeochemicalConcentrationFunction: public Function {

 public:

  // Constructs a GeochemicalConcentrationFunction that reports species concentrations 
  // computed by the given GeochemicalConditionContext object.
  GeochemicalConcentrationFunction(GeochemicalConditionContext* context, double* data, int speciesIndex):
    context_(context), data_(data), index_(speciesIndex)
  {
  }

  // Destructor.
  ~GeochemicalConcentrationFunction();

  // Overridden methods.
  Function* clone() const
  {
    return new GeochemicalConcentrationFunction(context_, data_, index_);
  }

  double operator() (const double* xt) const
  {
    context_->ComputeConcentrations(t);
    return data_[index_];
  }

 private:

  GeochemicalConditionContext* context_;
  double data_;
  int index_;
};

GeochemicalConditionContext::GeochemicalConditionContext(Teuchos::RCP<Chemistry_Engine> chem_engine, 
                                                         const std::string& geochem_condition):
  chem_engine_(chem_engine_), condition_(geochem_condition), functions_(chem_engine->NumPrimarySpecies())
{
}

GeochemicalConditionContext::~GeochemicalConditionContext()
{
}

Teuchos::RCP<Function> 
GeochemicalConditionContext::speciesFunction(const std::string& species)
{
  // Find the index of the given species within the list kept in the chemistry engine.
  int speciesIndex = -1;
  std::vector<std::string> speciesNames;
  chem_engine_->GetPrimarySpeciesNames(speciesNames);
  for (size_t i = 0; i < speciesNames.size(); ++i)
  {
    if (species == speciesNames[i])
    {
      speciesIndex = (int)i;
      break;
    }
  }
  // FIXME: Check for speciesIndex == -1.
  functions_[speciesIndex] = Teuchos::RCP<Function>(*this, speciesIndex);
  return functions_[speciesIndex];
}

void GeochemicalConditionContext::computeConcentrations(double t)
{
}

}  // namespace AmanziTransport
}  // namespace Amanzi
