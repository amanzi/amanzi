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
    temperature_key_(other.temperature_key_),
    pressure_key_(other.pressure_key_) {};


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

  temperature_key_ = plist_.get<std::string>("temperature key", prefix + "temperature");
  pressure_key_ = plist_.get<std::string>("pressure key", prefix + "pressure");
  dependencies_.insert(temperature_key_);
  dependencies_.insert(pressure_key_);
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

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pressure_key_);

  for (auto comp = result->begin(); comp != result->end(); ++comp) {
    auto kind = AmanziMesh::entity_kind(*comp);
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
    const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp);

    int ncomp = result->size(*comp);
    for (int i = 0; i != ncomp; ++i) {
      auto id = MyModel_(kind, i);
      result_v[0][i] = iem_->second[(*iem_->first)[id]]->InternalEnergy(temp_v[0][i], pres_v[0][i]);
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

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pressure_key_);

  for (auto comp = result->begin(); comp != result->end(); ++comp) {
    auto kind = AmanziMesh::entity_kind(*comp);
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
    const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp);

    int ncomp = result->size(*comp);

    if (wrt_key == temperature_key_) {
      for (int i = 0; i != ncomp; ++i) {
        auto id = MyModel_(kind, i);
        result_v[0][i] = iem_->second[(*iem_->first)[id]]->DInternalEnergyDT(temp_v[0][i], pres_v[0][i]);
      }
    } else if (wrt_key == pressure_key_) {
      for (int i = 0; i != ncomp; ++i) {
        auto id = MyModel_(kind, i);
        result_v[0][i] = iem_->second[(*iem_->first)[id]]->DInternalEnergyDp(temp_v[0][i], pres_v[0][i]);
      }
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
  mesh_ = mesh;

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
double IEMEvaluator::EvaluateFieldSingle(int c, double T, double p)
{
  return iem_->second[(*iem_->first)[c]]->InternalEnergy(T, p);
}


/* ******************************************************************
* Return model id
****************************************************************** */
AmanziMesh::Entity_ID IEMEvaluator::MyModel_(
    AmanziMesh::Entity_kind kind, AmanziMesh::Entity_ID id)
{
  if (kind == AmanziMesh::CELL) {
    return id;
  } else if (kind == AmanziMesh::BOUNDARY_FACE) {
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(id, AmanziMesh::Parallel_type::ALL, &cells);
    return cells[0];
  }
  return -1;
}

}  // namespace Energy
}  // namespace Amanzi
