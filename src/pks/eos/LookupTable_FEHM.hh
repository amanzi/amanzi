/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Tabulated equations of state.

  The underlying assumption is that all intersections of saturation
  line with non-uniform mesh are of type 6, i.e. through mesh vertices.
*/

#ifndef AMANZI_LOOKUP_TABLE_FEHM_HH_
#define AMANZI_LOOKUP_TABLE_FEHM_HH_

#include "Teuchos_ParameterList.hpp"

#include "LookupTable.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class LookupTable_FEHM : public LookupTable {
 public:
  LookupTable_FEHM(Teuchos::ParameterList& plist);
  ~LookupTable_FEHM() {};

  virtual double Function(double T, double p, int* ierr) override;
  virtual double DFunctionDT(double T, double p, int* ierr) override;
  virtual double DFunctionDp(double T, double p, int* ierr) override;
  virtual int Location(double T, double p, int* ierr) override;

 private:
  double Function_(double T, double p, int* ierr, const std::vector<std::vector<double>>& F);

  std::string ReadBlock_(std::ifstream& ifs, int nP, int nT, bool flag);
  void ReadBlockSat_(std::ifstream& ifs, int nS, std::vector<std::vector<double>>& satF);

  void ComputeDownwardMap_(int nP, int nT, int nS);

 private:
  std::string field_;

  std::vector<std::vector<int>> map_;
  std::vector<double> satT_, satP_;
  std::vector<std::vector<double>> satFl_, satFg_;

  std::vector<std::vector<double>> dFdP_, dFdT_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
