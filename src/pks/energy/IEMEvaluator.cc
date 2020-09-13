/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  The internal energy model evaluator simply calls the IEM with 
  the correct arguments.
*/

#include "IEMEvaluator.hh"
#include "IEMFactory.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructors.
****************************************************************** */
IEMEvaluator::IEMEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
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
IEMEvaluator::IEMEvaluator(const IEMEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    iem_(other.iem_),
    temp_key_(other.temp_key_) {};


/* ******************************************************************
* Copy assinement
****************************************************************** */
Teuchos::RCP<FieldEvaluator> IEMEvaluator::Clone() const {
  return Teuchos::rcp(new IEMEvaluator(*this));
}


/* ******************************************************************
* Itialization.
****************************************************************** */
void IEMEvaluator::InitializeFromPlist_()
{
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("internal energy key");
  }

  // Set up my dependencies.
  domain_ = Keys::getDomain(my_key_);
  std::string prefix = Keys::getDomainPrefix(my_key_);

  temp_key_ = plist_.get<std::string>("temperature key", prefix + "temperature");
  dependencies_.insert(temp_key_);
}


/* ******************************************************************
* Evaluation body.
****************************************************************** */
void IEMEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result)
{
  if (iem_ == Teuchos::null) {
    CreateIEMPartition_(S->GetMesh(domain_), plist_);
  }

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);

  for (auto comp = result->begin(); comp != result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp);

    int ncomp = result->size(*comp);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = iem_->second[(*iem_->first)[i]]->InternalEnergy(temp_v[0][i]);
    }
  }
}


/* ******************************************************************
* Evaluation of field derivatives.
****************************************************************** */
void IEMEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  if (iem_ == Teuchos::null) {
    CreateIEMPartition_(S->GetMesh(domain_), plist_);
  }

  AMANZI_ASSERT(wrt_key == temp_key_);
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);

  for (auto comp = result->begin(); comp != result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp);

    int ncomp = result->size(*comp);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = iem_->second[(*iem_->first)[i]]->DInternalEnergyDT(temp_v[0][i]);
    }
  }
}


/* ******************************************************************
* Create partition
****************************************************************** */
void IEMEvaluator::CreateIEMPartition_(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    const Teuchos::ParameterList& plist)
{
  std::vector<Teuchos::RCP<IEM> > iem_list;
  std::vector<std::vector<std::string> > region_list;

  IEMFactory fac;
  const Teuchos::ParameterList& tmp = plist.sublist("IEM parameters");

  for (auto lcv = tmp.begin(); lcv != tmp.end(); ++lcv) {
    std::string name = lcv->first;
    if (tmp.isSublist(name)) {
      const Teuchos::ParameterList& aux = tmp.sublist(name);
      region_list.push_back(aux.get<Teuchos::Array<std::string> >("regions").toVector());

      Teuchos::ParameterList model_list = aux.sublist("IEM parameters");
      iem_list.push_back(fac.CreateIEM(model_list));
    } else {
      AMANZI_ASSERT(0);
    }
  }

  auto partition = Teuchos::rcp(new Functions::MeshPartition());
  partition->Initialize(mesh, AmanziMesh::CELL, region_list, -1);
  partition->Verify();

  iem_ = Teuchos::rcp(new IEMPartition(partition, iem_list));
}


/* ******************************************************************
* Evaluation at a point
****************************************************************** */
double IEMEvaluator::EvaluateFieldSingle(int c, double T)
{
  return iem_->second[(*iem_->first)[c]]->InternalEnergy(T);
}

}  // namespace Energy
}  // namespace Amanzi
