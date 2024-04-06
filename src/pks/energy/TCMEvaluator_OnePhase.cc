/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Energy

  Interface for a thermal conductivity model with one phase.
*/

#include "dbc.hh"
#include "EOSFactory.hh"
#include "EOS_ThermalConductivity.hh"
#include "EOS_Utils.hh"

#include "TCMEvaluator_OnePhase.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor.
****************************************************************** */
TCMEvaluator_OnePhase::TCMEvaluator_OnePhase(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                                             Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(
      make_pair(plist_.get<std::string>("thermal conductivity key"), Tags::DEFAULT));
  }
  std::string domain = Keys::getDomain(my_keys_[0].first);

  temperature_key_ = plist_.get<std::string>("temperature key", Keys::getKey(domain, "temperature"));
  dependencies_.insert(std::make_pair(temperature_key_, Tags::DEFAULT));

  porosity_key_ = plist_.get<std::string>("porosity key", Keys::getKey(domain, "porosity"));
  dependencies_.insert(std::make_pair(porosity_key_, Tags::DEFAULT));

  // EOS for thermal conductivity
  std::vector<std::vector<std::string>> region_list;
  AmanziEOS::EOSFactory<AmanziEOS::EOS_ThermalConductivity> eos_fac;

  for (auto lcv = plist.begin(); lcv != plist.end(); ++lcv) {
    std::string name = lcv->first;
    if (plist.isSublist(name)) {
      Teuchos::ParameterList sublist = plist.sublist(name);
      region_list.push_back(sublist.get<Teuchos::Array<std::string>>("regions").toVector());

      tc_.push_back(eos_fac.Create(sublist));
      k_rock_.push_back(sublist.get<double>("thermal conductivity of rock"));
    }
  }

  partition_ = Teuchos::rcp(new Functions::MeshPartition());
  partition_->Initialize(mesh, AmanziMesh::Entity_kind::CELL, region_list, -1);
  partition_->Verify();
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
TCMEvaluator_OnePhase::TCMEvaluator_OnePhase(const TCMEvaluator_OnePhase& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    tc_(other.tc_),
    temperature_key_(other.temperature_key_){};


/* ******************************************************************
* TBW.
****************************************************************** */
Teuchos::RCP<Evaluator>
TCMEvaluator_OnePhase::Clone() const
{
  return Teuchos::rcp(new TCMEvaluator_OnePhase(*this));
}


/* ******************************************************************
* Evaluator body.
****************************************************************** */
void
TCMEvaluator_OnePhase::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  // pull out the dependencies
  const auto& temp_c = *S.Get<CompositeVector>(temperature_key_).ViewComponent("cell");
  const auto& poro_c = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  Epetra_MultiVector& result_c = *results[0]->ViewComponent("cell");

  int ierr(0);
  int ncomp = results[0]->size("cell", false);
  for (int i = 0; i != ncomp; ++i) {
    int id = (*partition_)[i];
    double k_liq = tc_[id]->ThermalConductivity(temp_c[0][i], 0.0);
    ierr = std::max(ierr, tc_[id]->error_code());

    double phi = poro_c[0][i];
    result_c[0][i] = phi * k_liq + (1.0 - phi) * k_rock_[id];
  }
  AmanziEOS::ErrorAnalysis(S.Get<CompositeVector>(temperature_key_).Comm(), ierr, tc_[0]->error_msg());
}


/* ******************************************************************
* Evaluator of derivarives.
****************************************************************** */
void
TCMEvaluator_OnePhase::EvaluatePartialDerivative_(const State& S,
                                                  const Key& wrt_key,
                                                  const Tag& wrt_tag,
                                                  const std::vector<CompositeVector*>& results)
{
  for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
    const auto& temp_v = *S.Get<CompositeVector>(temperature_key_).ViewComponent(*comp);
    const auto& poro_v = *S.Get<CompositeVector>(porosity_key_).ViewComponent(*comp);
    Epetra_MultiVector& result_v = *results[0]->ViewComponent(*comp);
    int ncells = results[0]->size(*comp);

    if (wrt_key == porosity_key_) {
      for (int i = 0; i != ncells; ++i) {
        int id = (*partition_)[i];
        double k_liq = tc_[id]->ThermalConductivity(temp_v[0][i], 0.0);
        result_v[0][i] = k_liq - k_rock_[id];
      }
    } else if (wrt_key == temperature_key_) {
      for (int i = 0; i != ncells; ++i) {
        int id = (*partition_)[i];
        double k_liq = tc_[id]->DThermalConductivityDT(temp_v[0][i], 0.0);
        result_v[0][i] = poro_v[0][i] * k_liq;
      }
    }
  }
}


/* ******************************************************************
* Compatibility check is not needed at this level.
****************************************************************** */
void
TCMEvaluator_OnePhase::EnsureCompatibility_Units_(State& S)
{
  S.GetRecordSetW(my_keys_[0].first).set_units("W/m/K");
}


} // namespace Energy
} // namespace Amanzi
