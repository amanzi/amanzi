/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Energy

  The internal energy model evaluator simply calls the IEM with
  the correct arguments.
*/

#include "EOS_Utils.hh"

#include "IEMEvaluator.hh"
#include "IEMFactory.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructors.
****************************************************************** */
IEMEvaluator::IEMEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  /*
  AMANZI_ASSERT(plist_.isSublist("IEM parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("IEM parameters");
  IEMFactory fac;
  iem_ = fac.CreateIEM(sublist);
  */

  InitializeFromPlist_();
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
IEMEvaluator::IEMEvaluator(const IEMEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    temperature_key_(other.temperature_key_),
    pressure_key_(other.pressure_key_),
    iem_(other.iem_){};


/* ******************************************************************
* Copy assinement
****************************************************************** */
Teuchos::RCP<Evaluator>
IEMEvaluator::Clone() const
{
  return Teuchos::rcp(new IEMEvaluator(*this));
}


/* ******************************************************************
* Itialization.
****************************************************************** */
void
IEMEvaluator::InitializeFromPlist_()
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(
      std::make_pair(plist_.get<std::string>("internal energy key"), Tags::DEFAULT));
  }

  // Set up my dependencies.
  tag_ = Tags::DEFAULT;
  domain_ = Keys::getDomain(my_keys_[0].first);
  std::string prefix = Keys::getDomainPrefix(my_keys_[0].first);

  temperature_key_ = plist_.get<std::string>("temperature key", prefix + "temperature");
  pressure_key_ = plist_.get<std::string>("pressure key", prefix + "pressure");
  dependencies_.insert(std::make_pair(temperature_key_, tag_));
  dependencies_.insert(std::make_pair(pressure_key_, tag_));
}


/* ******************************************************************
* Evaluation body.
****************************************************************** */
void
IEMEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  if (iem_ == Teuchos::null) { CreateIEMPartition_(S.GetMesh(domain_), plist_); }

  auto temp = S.GetPtr<CompositeVector>(temperature_key_, tag_);
  auto pres = S.GetPtr<CompositeVector>(pressure_key_, tag_);

  for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
    auto kind = AmanziMesh::createEntityKind(*comp);
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
    const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
    Epetra_MultiVector& result_v = *results[0]->ViewComponent(*comp);

    int ierr(0);
    std::string msg;

    int ncomp = results[0]->size(*comp);
    for (int i = 0; i != ncomp; ++i) {
      auto id = MyModel_(kind, i);
      auto model = iem_->second[(*iem_->first)[id]];
      result_v[0][i] = model->InternalEnergy(temp_v[0][i], pres_v[0][i]);

      if (model->error_code() > 0) {
        ierr = 1;
        msg = model->error_msg();
      }
    }
    AmanziEOS::ErrorAnalysis(temp->Comm(), ierr, msg);
  }
}


/* ******************************************************************
* Evaluation of field derivatives.
****************************************************************** */
void
IEMEvaluator::EvaluatePartialDerivative_(const State& S,
                                         const Key& wrt_key,
                                         const Tag& wrt_tag,
                                         const std::vector<CompositeVector*>& results)
{
  if (iem_ == Teuchos::null) { CreateIEMPartition_(S.GetMesh(domain_), plist_); }

  auto temp = S.GetPtr<CompositeVector>(temperature_key_, tag_);
  auto pres = S.GetPtr<CompositeVector>(pressure_key_, tag_);

  for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
    auto kind = AmanziMesh::createEntityKind(*comp);
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
    const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
    Epetra_MultiVector& result_v = *results[0]->ViewComponent(*comp);

    int ncomp = results[0]->size(*comp);

    if (wrt_key == temperature_key_) {
      for (int i = 0; i != ncomp; ++i) {
        auto id = MyModel_(kind, i);
        result_v[0][i] =
          iem_->second[(*iem_->first)[id]]->DInternalEnergyDT(temp_v[0][i], pres_v[0][i]);
      }
    } else if (wrt_key == pressure_key_) {
      for (int i = 0; i != ncomp; ++i) {
        auto id = MyModel_(kind, i);
        result_v[0][i] =
          iem_->second[(*iem_->first)[id]]->DInternalEnergyDp(temp_v[0][i], pres_v[0][i]);
      }
    }
  }
}


/* ******************************************************************
* Compatibility check is not needed.
****************************************************************** */
void
IEMEvaluator::EnsureCompatibility_Units_(State& S)
{
  S.GetRecordSetW(my_keys_[0].first).set_units("J/mol");
}


/* ******************************************************************
* Create partition
****************************************************************** */
void
IEMEvaluator::CreateIEMPartition_(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                  const Teuchos::ParameterList& plist)
{
  mesh_ = mesh;

  std::vector<Teuchos::RCP<IEM>> iem_list;
  std::vector<std::vector<std::string>> region_list;

  IEMFactory fac;
  const Teuchos::ParameterList& tmp = plist.sublist("IEM parameters");

  for (auto lcv = tmp.begin(); lcv != tmp.end(); ++lcv) {
    std::string name = lcv->first;
    if (tmp.isSublist(name)) {
      const Teuchos::ParameterList& aux = tmp.sublist(name);
      region_list.push_back(aux.get<Teuchos::Array<std::string>>("regions").toVector());

      Teuchos::ParameterList model_list = aux.sublist("IEM parameters");
      iem_list.push_back(fac.CreateIEM(model_list));
    } else {
      AMANZI_ASSERT(0);
    }
  }

  auto partition = Teuchos::rcp(new Functions::MeshPartition());
  partition->Initialize(mesh, AmanziMesh::Entity_kind::CELL, region_list, -1);
  partition->Verify();

  iem_ = Teuchos::rcp(new IEMPartition(partition, iem_list));
}


/* ******************************************************************
* Evaluation at a point
****************************************************************** */
double
IEMEvaluator::EvaluateFieldSingle(int c, double T, double p)
{
  return iem_->second[(*iem_->first)[c]]->InternalEnergy(T, p);
}


/* ******************************************************************
* Return model id
****************************************************************** */
AmanziMesh::Entity_ID
IEMEvaluator::MyModel_(AmanziMesh::Entity_kind kind, AmanziMesh::Entity_ID id)
{
  if (kind == AmanziMesh::Entity_kind::CELL) {
    return id;
  } else if (kind == AmanziMesh::Entity_kind::BOUNDARY_FACE) {
    auto cells = mesh_->getFaceCells(id, AmanziMesh::Parallel_kind::ALL);
    return cells[0];
  }
  return -1;
}

} // namespace Energy
} // namespace Amanzi
