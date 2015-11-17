/*
  This is the transport component of Amanzi. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Teuchos_RCP.hpp"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#include "mfd3d_diffusion.hh"
#include "nlfv.hh"
#include "Tensor.hh"
#include "PreconditionerFactory.hh"

#include "TransportDefs.hh"
#include "Transport_PK.hh"

namespace Amanzi {
namespace Transport {

/* *******************************************************************
* Calculate a dispersive tensor the from Darcy fluxes. The flux is
* assumed to be scaled by face area.
******************************************************************* */
void Transport_PK::CalculateDispersionTensor_(
    const Epetra_MultiVector& darcy_flux, 
    const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation)
{
  D_.resize(ncells_owned);
  for (int c = 0; c < ncells_owned; c++) D_[c].Init(dim, 2);

  AmanziGeometry::Point velocity(dim), omega(dim);
  AmanziMesh::Entity_ID_List faces;
  WhetStone::MFD3D_Diffusion mfd3d(mesh_);

  for (int mb = 0; mb < mat_properties_.size(); mb++) {
    Teuchos::RCP<MaterialProperties> spec = mat_properties_[mb]; 

    std::vector<AmanziMesh::Entity_ID> block;
    for (int r = 0; r < (spec->regions).size(); r++) {
      std::string region = (spec->regions)[r];
      mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

      AmanziMesh::Entity_ID_List::iterator c;
      if (spec->model == TRANSPORT_DISPERSIVITY_MODEL_SCALAR) {
        for (c = block.begin(); c != block.end(); c++) {
          D_[*c].PutScalar(0.0); 
          for (int i = 0; i < dim; i++) {
            D_[*c](i, i) = spec->alphaLH;
          }
          D_[*c] *= porosity[0][*c] * saturation[0][*c];
        }
      // Isotropic dispersitivity model.
      // Darcy velocity is reconstructed for its normal components.
      } else if (spec->model == TRANSPORT_DISPERSIVITY_MODEL_BEAR) {
        for (c = block.begin(); c != block.end(); c++) {
          D_[*c].PutScalar(0.0); 

          mesh_->cell_get_faces(*c, &faces);
          int nfaces = faces.size();

          std::vector<double> flux(nfaces);
          for (int n = 0; n < nfaces; n++) flux[n] = darcy_flux[0][faces[n]];
          mfd3d.RecoverGradient_MassMatrix(*c, flux, velocity);
          velocity /= porosity[0][*c];  // pore velocity

          double vel_norm = norm(velocity);
          if (vel_norm == 0.0) continue;

          double anisotropy = (spec->alphaLH - spec->alphaTH) / vel_norm;
          for (int i = 0; i < dim; i++) {
            D_[*c](i, i) = spec->alphaTH * vel_norm;
            for (int j = i; j < dim; j++) {
              double s = anisotropy * velocity[i] * velocity[j];
              D_[*c](j, i) = D_[*c](i, j) += s;
            }
          }

          D_[*c] *= porosity[0][*c] * saturation[0][*c];
        }
      // Anisotropic dispersitivity model.
      // This model assumes that space is 3D which is checked in Transport_IO.cc.
      // and that permeability tensor is diagonal. 
      } else if (spec->model == TRANSPORT_DISPERSIVITY_MODEL_BURNETT_FRIND ||
                 spec->model == TRANSPORT_DISPERSIVITY_MODEL_LICHTNER_KELKAR_ROBINSON) {
        const Epetra_MultiVector& perm = *S_->GetFieldData("permeability")->ViewComponent("cell");

        for (c = block.begin(); c != block.end(); c++) {
          D_[*c].PutScalar(0.0); 

          // Reconstruct Darcy velocity for its normal components.
          mesh_->cell_get_faces(*c, &faces);
          int nfaces = faces.size();

          std::vector<double> flux(nfaces);
          for (int n = 0; n < nfaces; n++) flux[n] = darcy_flux[0][faces[n]];
          mfd3d.RecoverGradient_MassMatrix(*c, flux, velocity);
          velocity /= porosity[0][*c];  // pore velocity

          double vel_norm = norm(velocity);
          if (vel_norm == 0.0) continue;

          double a1, a2, a3;
          if (spec->model == TRANSPORT_DISPERSIVITY_MODEL_BURNETT_FRIND) {
            a1 = spec->alphaTV * vel_norm;
            a2 = (spec->alphaLH - spec->alphaTV) / vel_norm;
            a3 = (spec->alphaTH - spec->alphaTV) / vel_norm;

            for (int i = 0; i < dim; i++) {
              D_[*c](i, i) = a1;
              for (int j = i; j < dim; j++) {
                D_[*c](i, j) += a2 * velocity[i] * velocity[j];
                D_[*c](j, i) = D_[*c](i, j);
              }
            }
            D_[*c](0, 0) += a3 * velocity[1] * velocity[1];
            D_[*c](1, 1) += a3 * velocity[0] * velocity[0];

            D_[*c](0, 1) -= a3 * velocity[0] * velocity[1];
            D_[*c](1, 0) -= a3 * velocity[0] * velocity[1];

          } else if (spec->model == TRANSPORT_DISPERSIVITY_MODEL_LICHTNER_KELKAR_ROBINSON) {
            int k = axi_symmetry_[*c];
            double theta = velocity[k] / vel_norm;  // cosine of angle theta
            double theta2 = theta * theta; 

            // define direction orthogonal to symmetry axis
            omega = velocity * (-theta / vel_norm);
            omega[k] += 1.0;

            // we use formula (46) of Lichtner, Water Res. Research, 38 (2002)
            double alphaL = spec->alphaLH + theta2 * (spec->alphaLV - spec->alphaLH);  
            a1 = spec->alphaTH * vel_norm;
            a2 = (alphaL - spec->alphaTH) / vel_norm;
            a3 = (spec->alphaTV - spec->alphaTH) * vel_norm;

            for (int i = 0; i < dim; i++) {
              D_[*c](i, i) = a1;
              for (int j = i; j < dim; j++) {
                D_[*c](i, j) += a2 * velocity[i] * velocity[j] + a3 * omega[i] * omega[j];
                D_[*c](j, i) = D_[*c](i, j);
              }
            }
          }

          D_[*c] *= porosity[0][*c] * saturation[0][*c];
        }
      }
    }
  }
}


/* *******************************************************************
* Calculate diffusion tensor and add it to the dispersion tensor.
******************************************************************* */
void Transport_PK::CalculateDiffusionTensor_(
    double md, int phase, 
    const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation)
{
  if (D_.size() == 0) {
    D_.resize(ncells_owned);
    for (int c = 0; c < ncells_owned; c++) D_[c].Init(dim, 1);
  }

  for (int mb = 0; mb < mat_properties_.size(); mb++) {
    Teuchos::RCP<MaterialProperties> spec = mat_properties_[mb]; 

    std::vector<AmanziMesh::Entity_ID> block;
    for (int r = 0; r < (spec->regions).size(); r++) {
      std::string region = (spec->regions)[r];
      mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

      AmanziMesh::Entity_ID_List::iterator c;
      if (phase == TRANSPORT_PHASE_LIQUID) {
        for (c = block.begin(); c != block.end(); c++) {
          D_[*c] += md * spec->tau[phase] * porosity[0][*c] * saturation[0][*c];
        }
      } else if (phase == TRANSPORT_PHASE_GAS) {
        for (c = block.begin(); c != block.end(); c++) {
          D_[*c] += md * spec->tau[phase] * porosity[0][*c] * (1.0 - saturation[0][*c]);
        }
      }
    }
  }
}


/* ******************************************************************
* Check all phases for the given name.
****************************************************************** */
int Transport_PK::FindDiffusionValue(const std::string& tcc_name, double* md, int* phase)
{
  for (int i = 0; i < TRANSPORT_NUMBER_PHASES; i++) {
    if (diffusion_phase_[i] == Teuchos::null) continue;
    int ok = diffusion_phase_[i]->FindDiffusionValue(tcc_name, md);
    if (ok == 0) {
      *phase = i;
      return 0;
    }
  }

  *md = 0.0;
  *phase = -1;
  return -1;
}


/* ******************************************************************
*  Find direction of axi-symmetry.                                               
****************************************************************** */
void Transport_PK::CalculateAxiSymmetryDirection()
{
  axi_symmetry_.resize(ncells_owned, -1);
  const Epetra_MultiVector& perm = *S_->GetFieldData("permeability")->ViewComponent("cell");

  for (int mb = 0; mb < mat_properties_.size(); mb++) {
    Teuchos::RCP<MaterialProperties> spec = mat_properties_[mb];

    std::vector<AmanziMesh::Entity_ID> block;
    for (int r = 0; r < (spec->regions).size(); r++) {
      std::string region = (spec->regions)[r];
      mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

      if (spec->model == TRANSPORT_DISPERSIVITY_MODEL_BURNETT_FRIND ||
          spec->model == TRANSPORT_DISPERSIVITY_MODEL_LICHTNER_KELKAR_ROBINSON) {
        AmanziMesh::Entity_ID_List::iterator c;
        for (c = block.begin(); c != block.end(); c++) {
          int k;
          if (perm[0][*c] != perm[1][*c] && perm[1][*c] == perm[2][*c]) {
            k = 0;
          } else if (perm[1][*c] != perm[2][*c] && perm[2][*c] == perm[0][*c]) {
            k = 1;
          } else if (perm[2][*c] != perm[0][*c] && perm[0][*c] == perm[1][*c]) {
            k = 2;
          } else {
            Errors::Message msg;
            msg << "Transport PK: Dispersivity model \"Burnett-Frind\" or " 
                << "\"Lichtner-Kelkar_Robinson\" can be applied only to an axi-symmetric materials.\n";
            Exceptions::amanzi_throw(msg);  
          }
          axi_symmetry_[*c] = k;
        }
      }
    }
  }
}

}  // namespace Transport
}  // namespace Amanzi



