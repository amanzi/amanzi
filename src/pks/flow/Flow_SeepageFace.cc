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
#include "mfd3d.hh"
#include "OperatorDefs.hh"
#include "State.hh"

#include "Flow_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Due to round-off errors, we have to use tolerances.
****************************************************************** */
bool Flow_PK::SeepageFacePFloTran(const CompositeVector& u, int* nseepage, double* area_seepage)
{
  const Epetra_MultiVector& flux = *S_->GetFieldData("darcy_flux")->ViewComponent("face", true);
  const Epetra_MultiVector& u_cell = *u.ViewComponent("cell");
  double ref_pressure = bc_seepage->reference_pressure();
  double tol = ref_pressure * 1e-14;

  bool flag(true);
  FlowBoundaryFunction::Iterator bc;
  for (bc = bc_seepage->begin(); bc != bc_seepage->end(); ++bc) {
    int f = bc->first;

    if (bc_submodel[f] & FLOW_BC_SUBMODEL_SEEPAGE_PFLOTRAN) {
      double face_value = BoundaryFaceValue(f, u);

      if (face_value < ref_pressure - tol) {
        bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value[f] = bc->second * rainfall_factor[f] * flux_units_;
      } else {
        int c = BoundaryFaceGetCell(f);
        if (u_cell[0][c] < face_value) {
          bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
          bc_value[f] = bc->second * rainfall_factor[f] * flux_units_;
        } else {
          bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value[f] = ref_pressure;
        }
      }
      if (flux[0][f] < 0.0) {
        bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value[f] = bc->second * rainfall_factor[f] * flux_units_;
      }
    } else {
      flag = false;
    }
  }

  // calculate area of the seepage face
  *nseepage = 0;
  *area_seepage = 0.0;

  for (bc = bc_seepage->begin(); bc != bc_seepage->end(); ++bc) {
    int f = bc->first;
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
      (*nseepage)++;
      (*area_seepage) += mesh_->face_area(f);
    }
  }

  return flag;
}


/* ******************************************************************
* Model II.
****************************************************************** */
bool Flow_PK::SeepageFaceFACT(const CompositeVector& u, int* nseepage, double* area_seepage)
{
  const Epetra_MultiVector& u_cell = *u.ViewComponent("cell");
  double ref_pressure = bc_seepage->reference_pressure();

  *nseepage = 0;
  *area_seepage = 0.0;

  bool flag(true);
  FlowBoundaryFunction::Iterator bc;
  for (bc = bc_seepage->begin(); bc != bc_seepage->end(); ++bc) {
    int f = bc->first;

    if (bc_submodel[f] & FLOW_BC_SUBMODEL_SEEPAGE_FACT) {
      double face_value = BoundaryFaceValue(f, u);

      double pcreg = -FLOW_BC_SEEPAGE_FACE_REGULARIZATION;
      double influx = bc->second * rainfall_factor[f] * flux_units_;
      double pcmin = 3 * pcreg / 2;
      double pcmax = pcreg / 2;

      double pc = face_value - ref_pressure;
      if (pc < pcmin) {
        bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value[f] = influx;
      } else if (pc < pcmax) {
        bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
        double f = 2 * (pcreg - pc) / pcreg;
        double q = (7 - 2 * f - f * f) / 8;
        bc_value[f] = q * influx; 
      } else {
        double I = influx / pcreg;
        bc_model[f] = Operators::OPERATOR_BC_MIXED;
        bc_value[f] = -I * ref_pressure;
        bc_mixed[f] = I;  // Impedance I should be positive.
        (*nseepage)++;
        (*area_seepage) += mesh_->face_area(f);
      }
    } else {
      flag = false;
    }
  }
  return flag;
}

}  // namespace Flow
}  // namespace Amanzi

