// PDE_DiffusionMFD: elliptic operators using the MFD family of discretizations.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_MFD_TRACER_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_MFD_TRACER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "exceptions.hh"
#include "Tensor.hh"
#include "Point.hh"
#include "CompositeVector.hh"
#include "DenseMatrix.hh"

#include "BCs.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionMFD.hh"


namespace Amanzi {
namespace Operators {

class Local_Surface {
  public:
  Local_Surface(int surface_id) : surface_id_(surface_id) {};
  Local_Surface(int surface_id,
                AmanziGeometry::Point sur_cntr) :
    surface_id_(surface_id), sur_cntr_(sur_cntr)  {};
  Local_Surface(int surface_id,
                AmanziGeometry::Point sur_cntr,
                AmanziGeometry::Point sur_norm) :
    surface_id_(surface_id), sur_cntr_(sur_cntr), sur_norm_(sur_norm) {};
 
  protected:
    int surface_id_;
    AmanziGeometry::Point sur_cntr_;
    AmanziGeometry::Point sur_norm_;
    std::vector< std::vector<AmanziGeometry::Point> > surface_pnt_;
    std::vector< std::vector<int> > edges_intersected_;
    std::vector< std::vector<double> > edges_coef;
};

  
class PDE_DiffusionMFD_Tracer : public virtual PDE_DiffusionMFD {
 public:
 PDE_DiffusionMFD_Tracer(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<Operator>& global_op) :
   PDE_Diffusion(global_op),
   PDE_DiffusionMFD(plist,global_op)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_TRACER;
    InitDiffusion_(plist);
  }

  PDE_DiffusionMFD_Tracer(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
    PDE_Diffusion(mesh),
    PDE_DiffusionMFD(plist, mesh)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_TRACER;
    InitDiffusion_(plist);
  }

  PDE_DiffusionMFD_Tracer(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
    PDE_Diffusion(mesh),
    PDE_DiffusionMFD(plist, mesh)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_TRACER;
    InitDiffusion_(plist);
  }

  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u) override;

 protected:
  void InitDiffusion_(Teuchos::ParameterList& plist);
  void CellSurfaceInterception_(int c, const AmanziGeometry::Point& cntr, const AmanziGeometry::Point& norm);
  bool LineLineIntersect_(const AmanziGeometry::Point& p1,
                          const AmanziGeometry::Point& p2,
                          const AmanziGeometry::Point& p3,
                          const AmanziGeometry::Point& p4,
                          double EPS,
                          AmanziGeometry::Point& int_p1,
                          AmanziGeometry::Point& int_p2,
                          double *mu_1, double *mu_2);
  
  std::vector< std::vector<Local_Surface> > cell_surfaces_;

};



}
}

#endif
