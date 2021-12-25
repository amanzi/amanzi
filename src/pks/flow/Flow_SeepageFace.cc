/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1) 
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "GMVMesh.hh"
#include "Mesh.hh"
#include "Mesh_Algorithms.hh"
#include "OperatorDefs.hh"
#include "State.hh"

#include "Flow_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Due to round-off errors, we have to use tolerances.
****************************************************************** */
void Flow_PK::SeepageFacePFloTran(const CompositeVector& u, int* nseepage, double* area_seepage)
{
  const auto& flux = *S_->Get<CompositeVector>(darcy_flux_key_, Tags::DEFAULT).ViewComponent("face", true);
  const auto& u_cell = *u.ViewComponent("cell");

  std::vector<int>& bc_model = op_bc_->bc_model();
  std::vector<double>& bc_value = op_bc_->bc_value();
  std::vector<double>& bc_mixed = op_bc_->bc_mixed();

  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "seepage" && 
        bcs_[i]->seepage_model() == "PFloTran") {
      double ref_pressure = bcs_[i]->ref_pressure();
      double tol = ref_pressure * 1e-14;
      double flux_threshold = bcs_[i]->seepage_flux_threshold();

      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        double face_value = BoundaryFaceValue(f, u);

        if (face_value < ref_pressure - tol) {
          bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
          bc_value[f] = it->second[0] * flux_units_;
        } else {
          int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh_, f);
          if (u_cell[0][c] < face_value) {
            bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
            bc_value[f] = it->second[0] * flux_units_;
          } else {
            bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
            bc_value[f] = ref_pressure;
          }
        }
        if (flux[0][f] < it->second[0] * flux_units_ * flux_threshold) {
          bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
          bc_value[f] = it->second[0] * flux_units_;
        }
      }
    }
  }

  // calculate area of the seepage face
  *nseepage = 0;
  *area_seepage = 0.0;

  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "seepage") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
          (*nseepage)++;
          (*area_seepage) += mesh_->face_area(f);
        }
      }
    }
  }
}


/* ******************************************************************
* Model II.
****************************************************************** */
void Flow_PK::SeepageFaceFACT(const CompositeVector& u, int* nseepage, double* area_seepage)
{
  std::vector<int>& bc_model = op_bc_->bc_model();
  std::vector<double>& bc_value = op_bc_->bc_value();
  std::vector<double>& bc_mixed = op_bc_->bc_mixed();

  *nseepage = 0;
  *area_seepage = 0.0;

  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "seepage" && 
        bcs_[i]->seepage_model() == "FACT") {
      double ref_pressure = bcs_[i]->ref_pressure();

      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        double face_value = BoundaryFaceValue(f, u);

        double pcreg = -FLOW_BC_SEEPAGE_FACE_REGULARIZATION;
        double influx = it->second[0] * flux_units_;
        double pcmin = 3 * pcreg / 2;
        double pcmax = pcreg / 2;

        double pc = face_value - ref_pressure;
        if (pc < pcmin) {
          bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
          bc_value[f] = influx;
        } else if (pc < pcmax) {
          bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
          double a = 2 * (pcreg - pc) / pcreg;
          double q = (7 - 2 * a - a * a) / 8;
          bc_value[f] = q * influx; 
        } else {
          double I = influx / pcreg;
          bc_model[f] = Operators::OPERATOR_BC_MIXED;
          bc_value[f] = -I * ref_pressure;
          bc_mixed[f] = I;  // Impedance I should be positive.
          (*nseepage)++;
          (*area_seepage) += mesh_->face_area(f);
        }
      }
    }
  }
}

}  // namespace Flow
}  // namespace Amanzi

