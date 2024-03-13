/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry

  Base class for colloids.
*/

#ifndef AMANZI_CHEMISTRY_COLLOID_HH_
#define AMANZI_CHEMISTRY_COLLOID_HH_

#include <vector>

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class Colloid {
 public:
  Colloid();
  Colloid(int id, const std::string& name, const Teuchos::ParameterList& plist)
    : id_(id), name_(name){};
  virtual ~Colloid(){};

  // returns sorbed concentration
  virtual void Evaluate(const Species& primary_species){};
  virtual void EvaluateDerivative(const Species& primary_species){};

  // access
  std::string name() const { return name_; }

 protected:
  int id_;
  std::string name_;
};

} // namespace AmanziChemistry
} // namespace Amanzi

#endif
