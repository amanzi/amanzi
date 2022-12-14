/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef AMANZI_OBSERVABLE_LINE_SEGMENT_AQUEOUS_HH
#define AMANZI_OBSERVABLE_LINE_SEGMENT_AQUEOUS_HH

#include "LimiterCell.hh"
#include "ReconstructionCellLinear.hh"

#include "ObservableAmanzi.hh"
#include "ObservableAqueous.hh"
#include "ObservableLineSegment.hh"


namespace Amanzi {

class ObservableLineSegmentAqueous : public ObservableAqueous, public ObservableLineSegment {
 public:
  ObservableLineSegmentAqueous(std::string variable,
                               std::string region,
                               std::string functional,
                               Teuchos::ParameterList& plist,
                               Teuchos::ParameterList& units_plist,
                               Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  virtual void
  ComputeObservation(State& S, double* value, double* volume, std::string& unit, double dt);
  virtual int ComputeRegionSize();
  void InterpolatedValues(State& S,
                          std::string& var,
                          std::string& interpolation,
                          AmanziMesh::Entity_ID_List& ids,
                          std::vector<AmanziGeometry::Point>& line_pnts,
                          std::vector<double>& values);
};


ObservableLineSegmentAqueous::ObservableLineSegmentAqueous(
  std::string variable,
  std::string region,
  std::string functional,
  Teuchos::ParameterList& plist,
  Teuchos::ParameterList& units_plist,
  Teuchos::RCP<const AmanziMesh::Mesh> mesh)
  : Observable(variable, region, functional, plist, units_plist, mesh),
    ObservableAqueous(variable, region, functional, plist, units_plist, mesh),
    ObservableLineSegment(variable, region, functional, plist, units_plist, mesh){};


int
ObservableLineSegmentAqueous::ComputeRegionSize()
{
  return ObservableLineSegment::ComputeRegionSize();
}


void
ObservableLineSegmentAqueous::ComputeObservation(State& S,
                                                 double* value,
                                                 double* volume,
                                                 std::string& unit,
                                                 double dt)
{
  Errors::Message msg;
  int dim = mesh_->space_dimension();

  std::vector<double> values(region_size_);
  double weight_corr = 1e-15;

  InterpolatedValues(S, variable_, interpolation_, entity_ids_, line_points_, values);

  *value = 0.0;
  *volume = 0.0;
  unit = "m";

  if (weighting_ == "none") {
    for (int i = 0; i < region_size_; i++) {
      *value += values[i] * lofs_[i];
      *volume += lofs_[i];
    }
  } else if (weighting_ == "flux norm") {
    if (S.HasRecord("darcy_velocity")) {
      const auto& darcy_vel = *S.Get<CompositeVector>("darcy_velocity").ViewComponent("cell");
      for (int i = 0; i < region_size_; i++) {
        int c = entity_ids_[i];

        double norm = 0.0;
        for (int j = 0; j < dim; j++) norm += darcy_vel[j][c] * darcy_vel[j][c];
        norm = sqrt(norm) + weight_corr;

        *value += values[i] * lofs_[i] * norm;
        *volume += lofs_[i] * norm;
      }
    }
  }
}


void
ObservableLineSegmentAqueous::InterpolatedValues(State& S,
                                                 std::string& var,
                                                 std::string& interpolation,
                                                 AmanziMesh::Entity_ID_List& ids,
                                                 std::vector<AmanziGeometry::Point>& line_pnts,
                                                 std::vector<double>& values)
{
  Teuchos::RCP<const Epetra_MultiVector> vector;
  Teuchos::RCP<const CompositeVector> cv;

  if (var == "hydraulic head") {
    if (!S.HasRecord("hydraulic_head")) {
      Errors::Message msg;
      msg << "InterpolatedValue: field hydraulic_head doesn't exist in state";
      Exceptions::amanzi_throw(msg);
    }
    cv = S.GetPtr<CompositeVector>("hydraulic_head", Tags::DEFAULT);
    vector = cv->ViewComponent("cell", true);
  } else {
    if (!S.HasRecord(var)) {
      Errors::Message msg;
      msg << "InterpolatedValue: field " << var << " doesn't exist in state";
      Exceptions::amanzi_throw(msg);
    }
    cv = S.GetPtr<CompositeVector>(var, Tags::DEFAULT);
    vector = cv->ViewComponent("cell", true);
  }

  if (interpolation == "linear") {
    Teuchos::ParameterList plist;
    auto lifting = Teuchos::rcp(new Operators::ReconstructionCellLinear(mesh_));

    cv->ScatterMasterToGhosted();

    lifting->Init(plist);
    lifting->Compute(ids, vector, 0);

    if (limiter_) {
      plist.set<std::string>("limiter", "Kuzmin");
      std::vector<int> bc_model;
      std::vector<double> bc_value;

      Operators::LimiterCell limiter(mesh_);
      limiter.Init(plist);
      limiter.ApplyLimiter(ids, vector, 0, lifting, bc_model, bc_value);
    }

    for (int i = 0; i < ids.size(); i++) {
      int c = ids[i];
      values[i] = lifting->getValue(c, line_pnts[i]);
    }
  } else if (interpolation == "constant") {
    for (int i = 0; i < ids.size(); i++) {
      int c = ids[i];
      values[i] = (*vector)[0][c];
    }
  } else {
    Errors::Message msg;
    msg << "InterpolatedValue: unknown interpolation method " << interpolation;
    Exceptions::amanzi_throw(msg);
  }
}

} // namespace Amanzi

#endif
