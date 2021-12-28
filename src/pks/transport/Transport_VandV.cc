/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <vector>

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "errors.hh"

#include "Transport_PK.hh"

namespace Amanzi {
namespace Transport {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ****************************************************************
* Construct default state for unit tests.
**************************************************************** */
void Transport_PK::CreateDefaultState(
    Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int ncomponents) 
{
  std::string name("state"); 
  S_->Require<double>("const_fluid_density", Tags::DEFAULT, name);

  if (!S_->HasData(saturation_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(saturation_liquid_key_, Tags::DEFAULT, name)
      .SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  
  if (!S_->HasData(prev_saturation_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(prev_saturation_liquid_key_, Tags::DEFAULT, name)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S_->HasData(darcy_flux_key_)) {
    S_->Require<CV_t, CVS_t>(darcy_flux_key_, Tags::DEFAULT, name)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("face", AmanziMesh::FACE, 1);
  }
  
  if (!S_->HasData(tcc_key_)) {
    S_->Require<CV_t, CVS_t>(tcc_key_, Tags::DEFAULT, name, component_names_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, ncomponents);
  }

  // initialize fields
  S_->Setup();
 
  // set popular default values
  S_->GetW<double>("const_fluid_density", name) = 1000.0;
  S_->GetRecordW("const_fluid_density", name).set_initialized();

  S_->GetW<CV_t>(saturation_liquid_key_, name).PutScalar(1.0);
  S_->GetRecordW(saturation_liquid_key_, name).set_initialized();

  S_->GetW<CV_t>(prev_saturation_liquid_key_, name).PutScalar(1.0);
  S_->GetRecordW(prev_saturation_liquid_key_, name).set_initialized();

  S_->GetW<CV_t>(tcc_key_, name).PutScalar(0.0);
  S_->GetRecordW(tcc_key_, name).set_initialized();

  S_->GetW<CV_t>(darcy_flux_key_, name).PutScalar(0.0);
  S_->GetRecordW(darcy_flux_key_, name).set_initialized();

  S_->InitializeFields();
}


/* *******************************************************************
* Routine verifies that the velocity field is divergence free                 
******************************************************************* */
void Transport_PK::Policy(Teuchos::Ptr<State> S)
{
  if (mesh_->get_comm()->NumProc() > 1) {
    if (!S->Get<CV_t>(tcc_key_).Ghosted()) {
      Errors::Message msg;
      msg << "Field \"total component concentration\" has no ghost values."
          << " Transport PK is giving up.\n";
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* *******************************************************************
* Calculates extrema of specified solutes and print them.
******************************************************************* */
void Transport_PK::VV_PrintSoluteExtrema(
    const Epetra_MultiVector& tcc_next, double dT_MPC, const std::string& mesh_id)
{
  const auto& darcy_flux = *S_->Get<CV_t>(darcy_flux_key_).ViewComponent("face", true);

  int num_components = tcc_next.NumVectors();
  double tccmin_vec[num_components];
  double tccmax_vec[num_components];

  tcc_next.MinValue(tccmin_vec);
  tcc_next.MaxValue(tccmax_vec);

  for (int n = 0; n < runtime_solutes_.size(); n++) {
    int i = FindComponentNumber(runtime_solutes_[n]);
    double tccmin, tccmax;
    tcc_next.Comm().MinAll(&(tccmin_vec[i]), &tccmin, 1);  // find the global extrema
    tcc_next.Comm().MaxAll(&(tccmax_vec[i]), &tccmax, 1); 

    int nregions = runtime_regions_.size();
    double solute_flux(0.0);
    bool flag(false);

    for (int k = 0; k < nregions; k++) {
      if (mesh_->valid_set_name(runtime_regions_[k], AmanziMesh::FACE)) {
        flag = true;
        AmanziMesh::Entity_ID_List block;
        mesh_->get_set_entities(runtime_regions_[k], AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED, &block);
        int nblock = block.size();

        for (int m = 0; m < nblock; m++) {
          int f = block[m];

          Amanzi::AmanziMesh::Entity_ID_List cells;
          mesh_->face_get_cells(f, Amanzi::AmanziMesh::Parallel_type::ALL, &cells);
          int dir, c = cells[0];

          mesh_->face_normal(f, false, c, &dir);
          double u = darcy_flux[0][f] * dir;
          if (u > 0) solute_flux += u * tcc_next[i][c];
        }
      }
    }
    solute_flux *= units_.concentration_factor();

    double tmp = solute_flux;
    mesh_->get_comm()->SumAll(&tmp, &solute_flux, 1);

    *vo_->os() << runtime_solutes_[n] << mesh_id << ": min=" << units_.OutputConcentration(tccmin) 
               << " max=" << units_.OutputConcentration(tccmax);
    if (flag) *vo_->os() << ", flux=" << solute_flux << " mol/s";

    // old capability
    mass_solutes_exact_[i] += VV_SoluteVolumeChangePerSecond(i) * dT_MPC;
    double mass_solute(0.0);
    for (int c = 0; c < ncells_owned; c++) {
      double vol = mesh_->cell_volume(c);
      mass_solute += (*ws)[0][c] * (*phi)[0][c] * tcc_next[i][c] * vol;
    }
    mass_solute *= units_.concentration_factor();

    double tmp1 = mass_solute, tmp2 = mass_solutes_exact_[i], mass_exact;
    mesh_->get_comm()->SumAll(&tmp1, &mass_solute, 1);
    mesh_->get_comm()->SumAll(&tmp2, &mass_exact, 1);

    *vo_->os() << ", total=" << mass_solute << " mol" << std::endl;
  }
}


/* *******************************************************************
* Fancy output of limiter statistics
******************************************************************* */
void Transport_PK::VV_PrintLimiterStatistics()
{
  if (vo_->getVerbLevel() > Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    for (int n = 0; n < runtime_solutes_.size(); n++) {
      std::string& name = runtime_solutes_[n];

      if (FindComponentNumber(name) == current_component_) {
        const Epetra_Vector& limiter = *limiter_->limiter();
        double vmin, vavg, vmax;

        limiter.MinValue(&vmin);
        limiter.MaxValue(&vmax);
        limiter.MeanValue(&vavg);

        *vo_->os() << name << ": limiter min/avg/max: " 
                   << vmin << " " << vavg<< " " << vmax << std::endl;
      }
    }
  }
}


/********************************************************************
* Check completeness of influx boundary conditions.                        
****************************************************************** */
void Transport_PK::VV_CheckInfluxBC() const
{
  const auto& darcy_flux = *S_->Get<CV_t>(darcy_flux_key_).ViewComponent("face", true);

  int number_components = tcc->ViewComponent("cell")->NumVectors();
  std::vector<int> influx_face(nfaces_wghost);

  for (int i = 0; i < number_components; i++) {
    influx_face.assign(nfaces_wghost, 0);

    for (int m = 0; m < bcs_.size(); m++) {
      std::vector<int>& tcc_index = bcs_[m]->tcc_index();
      int ncomp = tcc_index.size();

      for (int k = 0; k < ncomp; k++) {
        if (i == tcc_index[k]) {
          for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
            int f = it->first;
            influx_face[f] = 1;
          }
        }
      }
    }

    for (int m = 0; m < bcs_.size(); m++) {
      std::vector<int>& tcc_index = bcs_[m]->tcc_index();
      int ncomp = tcc_index.size();

      for (int k = 0; k < ncomp; k++) {
        if (i == tcc_index[k]) {
          for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
            int f = it->first;
            if (darcy_flux[0][f] < 0 && influx_face[f] == 0) {
              Errors::Message msg;
              msg << "No influx boundary condition has been found for component " << i << ".\n";
              Exceptions::amanzi_throw(msg);
            }
          }
        }
      }
    }
  }
}


/* *******************************************************************
 * Check that global extrema diminished                          
 ****************************************************************** */
void Transport_PK::VV_CheckGEDproperty(Epetra_MultiVector& tracer) const
{
  int i, num_components = tracer.NumVectors();
  double tr_min[num_components];
  double tr_max[num_components];

  tracer.MinValue(tr_min);
  tracer.MaxValue(tr_max);

  for (i = 0; i < num_components; i++) {
    if (tr_min[i] < 0) {
      std::cout << "Transport_PK: concentration violates GED property" << std::endl;
      std::cout << "    Make an Amanzi ticket or turn off internal transport tests" << std::endl;
      std::cout << "    MyPID = " << MyPID << std::endl;
      std::cout << "    component = " << i << std::endl;
      std::cout << "    time = " << t_physics_ << std::endl;
      std::cout << "    min/max values = " << tr_min[i] << " " << tr_max[i] << std::endl;

      Errors::Message msg;
      msg << "Concentration violates GED property." << "\n";
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* ******************************************************************
* Check that the tracer is between 0 and 1.                        
****************************************************************** */
void Transport_PK::VV_CheckTracerBounds(Epetra_MultiVector& tracer,
                                        int component,
                                        double lower_bound,
                                        double upper_bound,
                                        double tol) const
{
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    double value = tracer[component][c];
    if (value < lower_bound - tol || value > upper_bound + tol) {
      std::cout << "Transport_PK: tracer violates bounds" << std::endl;
      std::cout << "    Make an Amanzi ticket or turn off internal transport tests" << std::endl;
      std::cout << "    MyPID = " << MyPID << std::endl;
      std::cout << "    component = " << component << std::endl;
      std::cout << "    simulation time = " << t_physics_ << std::endl;
      std::cout << "      cell = " << c << std::endl;
      std::cout << "      center = " << mesh_->cell_centroid(c) << std::endl;
      std::cout << "      value (old) = " << tcc_prev[component][c] << std::endl;
      std::cout << "      value (new) = " << value << std::endl;

      Errors::Message msg;
      msg << "Tracer violates bounds." << "\n";
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* ******************************************************************
* Calculate change of tracer volume per second due to boundary flux.
* This is the simplified version (lipnikov@lanl.gov).
****************************************************************** */
double Transport_PK::VV_SoluteVolumeChangePerSecond(int idx_tracer)
{
  const auto& darcy_flux = *S_->Get<CV_t>(darcy_flux_key_).ViewComponent("face", true);

  double volume = 0.0;

  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (tcc_index[i] == idx_tracer) {
        for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
          int f = it->first;

          std::vector<double>& values = it->second; 

          if (downwind_cells_[f].size() > 0) {
            int c2 = downwind_cells_[f][0];

            if (f < nfaces_owned && c2 >= 0) {
              double u = fabs(darcy_flux[0][f]);
              volume += u * values[i];
            }
          }
        }
      }
    }
  }
  return volume;
}


/* *******************************************************************
* Error estimate uses analytic function and solution.
* ***************************************************************** */
void Transport_PK::CalculateLpErrors(
    AnalyticFunction f, double t, Epetra_Vector* sol, double* L1, double* L2)
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  *L1 = *L2 = 0.0;
  for (int c = 0; c < ncells; c++) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    double d = (*sol)[c] - f(xc, t);

    double volume = mesh_->cell_volume(c);
    *L1 += fabs(d) * volume;
    *L2 += d * d * volume;
  }

  double tmp_out[2], tmp_in[2] = {*L1, *L2};
  mesh_->get_comm()->SumAll(tmp_in, tmp_out, 2);
  *L1 = tmp_out[0];
  *L2 = sqrt(tmp_out[1]);
}

}  // namespace Transport
}  // namespace Amanzi

