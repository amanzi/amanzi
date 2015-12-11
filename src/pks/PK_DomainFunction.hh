/*
  Process Kernel

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "MultiFunction.hh"
#include "UniqueMeshFunction.hh"

namespace Amanzi {

class PK_DomainFunction : public Functions::UniqueMeshFunction {
 public:
  PK_DomainFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : 
      UniqueMeshFunction(mesh),
      finalized_(false) {};

  virtual void Define(const std::vector<std::string>& regions,
                      const Teuchos::RCP<const MultiFunction>& f,
                      int action, int submodel);

  virtual void Define(const std::string& region,
                      const Teuchos::RCP<const MultiFunction>& f,
                      int action, int submodel);

  // source term on time interval (t0, t1]
  virtual void Compute(double t0, double t1, Teuchos::RCP<const Epetra_Vector> weight);

  // a place keeper
  void Finalize() {};
 
  // iterator methods
  typedef std::map<int, double>::const_iterator Iterator;
  Iterator begin() const { return value_.begin(); }
  Iterator end() const { return value_.end(); }
  std::map<int, double>::size_type size() { return value_.size(); }

  // extract internal information
  int CollectActionsList();

 private:
  void ComputeDensity_(double t0, double t1);
  void ComputeIntegral_(double t0, double t1);
  void ComputeIntegral_(double t0, double t1, double* weight);
  
 protected:
  std::map<int, double> value_;
  bool finalized_;

  std::vector<int> actions_;
  std::vector<int> submodel_;
};

}  // namespace Amanzi

#endif

