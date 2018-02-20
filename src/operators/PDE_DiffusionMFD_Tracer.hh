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

struct Intercetion_pnt{
  AmanziGeometry::Point pnt_;
  std::vector<int> v_ids_;
  int edge_id_;
  double inter_coef_;
};

class Local_Surface {
  public:
  Local_Surface(int surface_id) : surface_id_(surface_id), nfaces_(0) {};
  Local_Surface(int surface_id,
                AmanziGeometry::Point sur_norm) :
    surface_id_(surface_id),  sur_norm_(sur_norm), nfaces_(0)  {};
  Local_Surface(int surface_id,
                AmanziGeometry::Point sur_norm,
                double d) :
    surface_id_(surface_id), sur_norm_(sur_norm), sur_d_(d), nfaces_(0) {};

  std::vector< std::vector<AmanziGeometry::Point> >& surface_pnt(){return surface_pnt_;}
  std::vector< std::vector<int> >& v_ids(){return v_ids_;}
  std::vector< std::vector<int> >& edge_ids(){return edge_ids_;}
  std::vector< std::vector<double> >& inter_coef() {return inter_coef_;}
  int num_faces() {return nfaces_;}
  int surface_id(){ return surface_id_;}
  
  void Add_Face( std::vector<AmanziGeometry::Point>& vert_pnt, std::vector<int>& v_ids, std::vector<double>& wgts){
    surface_pnt_.push_back(vert_pnt);
    v_ids_.push_back(v_ids);
    inter_coef_.push_back(wgts);
    nfaces_++;
  }

  void Add_Face( std::vector<AmanziGeometry::Point>& vert_pnt, std::vector<int>& v_ids, std::vector<double>& wgts, std::vector<int>& edge ){
    surface_pnt_.push_back(vert_pnt);
    v_ids_.push_back(v_ids);
    inter_coef_.push_back(wgts);
    edge_ids_.push_back(edge);
    nfaces_++;
  }

  void print(){
    int nfaces = surface_pnt_.size();
    for (int i=0;i<nfaces;i++){
      std::cout<<"points: "<<surface_pnt_[i][0]<<" -- "<<surface_pnt_[i][1]<<"\n";
      std::cout<<"ids: "<<v_ids_[i][0]<<" "<<v_ids_[i][1]<<" -- "<<v_ids_[i][2]<<" "<<v_ids_[i][3]<<"\n";
      std::cout<<"weights: "<<inter_coef_[i][0]<<" -- "<<inter_coef_[i][1]<<"\n";
      std::cout<<"edges: "<<edge_ids_[i][0]<<" -- "<<edge_ids_[i][1]<<"\n";
      std::cout<<"\n";
    }
  }
 
  protected:
    int surface_id_;
    int nfaces_;
    AmanziGeometry::Point sur_norm_;
    double sur_d_;
    std::vector< std::vector<AmanziGeometry::Point> > surface_pnt_;
    std::vector< std::vector<int> > v_ids_;
    std::vector< std::vector<int> > edge_ids_;
    std::vector< std::vector<double> > inter_coef_;
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
                              const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& surf_presence,
                              const Teuchos::Ptr<const CompositeVector>& surface_parm);
  
  virtual void ApplyBCs(const Epetra_MultiVector& marker, bool primary, bool eliminate);
  void OutputGMV_Surface(int surf_id, const Epetra_MultiVector& data, std::string filename);

 protected:
  void InitDiffusion_(Teuchos::ParameterList& plist);
  void ApplyBCs_Nodal_(const Epetra_MultiVector& marker,
                       const Teuchos::Ptr<BCs>& bc_f,
                       const Teuchos::Ptr<BCs>& bc_v,
                       bool primary, bool eliminate);
  void CellSurfaceInterception_(int c, int surf_id, const AmanziGeometry::Point& norm, double surf_d);
  bool LineLineIntersect_(const AmanziGeometry::Point& p1,
                          const AmanziGeometry::Point& p2,
                          const AmanziGeometry::Point& p3,
                          const AmanziGeometry::Point& p4,
                          double EPS,
                          AmanziGeometry::Point& int_p1,
                          AmanziGeometry::Point& int_p2,
                          double *mu_1, double *mu_2);
  
  std::vector< std::vector< Teuchos::RCP<Local_Surface> > > cell_surfaces_;

};



}
}

#endif
