/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  An elevation evaluator getting values from the volumetric mesh.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
#include <boost/algorithm/string/predicate.hpp>

#include "Mesh.hh"
#include "Mesh_MSTK.hh"
#include "Point.hh"
#include "elevation_evaluator_column.hh"

namespace Amanzi {
namespace Flow {

ElevationEvaluatorColumn::ElevationEvaluatorColumn(Teuchos::ParameterList& plist) :
  ElevationEvaluator(plist) {
  
  Key domain = Keys::getDomain(my_keys_[0]);

  if (domain == "surface_star") {
    for (int c =0; c < 4; c++){
      std::stringstream name;
      name << "column_"<< c;
      base_por_key_ = Keys::readKey(plist_, name.str(), "base porosity", "base_porosity");
      dependencies_.insert(base_por_key_);
    }
  } else {
    Errors::Message msg("ElevationEvaluatorColumn: this evaluator should be used for columnar meshes only.");
    Exceptions::amanzi_throw(msg);
  }  
};

ElevationEvaluatorColumn::ElevationEvaluatorColumn(const ElevationEvaluatorColumn& other) :
  ElevationEvaluator(other),
  base_por_key_(other.base_por_key_)
{};

Teuchos::RCP<FieldEvaluator>
ElevationEvaluatorColumn::Clone() const {
  return Teuchos::rcp(new ElevationEvaluatorColumn(*this));
}

void ElevationEvaluatorColumn::EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {

  Teuchos::Ptr<CompositeVector> elev = results[0];
  Teuchos::Ptr<CompositeVector> slope = results[1];
  Epetra_MultiVector& elev_c = *elev->ViewComponent("cell", false);
  Epetra_MultiVector& slope_c = *slope->ViewComponent("cell", false);
 
  
  // Get the elevation and slope values from the domain mesh.
  Key domain = Keys::getDomain(my_keys_[0]);

  auto surface_mesh =
    Teuchos::rcp_static_cast<const AmanziMesh::Mesh_MSTK>(S->GetMesh(domain));
  
  
  // Set the elevation on cells by getting the corresponding face and its
  // centroid.
  int ncells = elev_c.MyLength();
  std::vector<AmanziGeometry::Point> my_centroid; 
  for (int c=0; c!=ncells; ++c) {
    std::stringstream my_name;
    int id = S->GetMesh("surface_star")->cell_map(false).GID(c);
    my_name << "column_" << id;
    
    int nfaces = S->GetMesh(my_name.str())->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
    
    std::vector<AmanziGeometry::Point> coord; 
    S->GetMesh(my_name.str())->face_get_coordinates(nfaces-1, &coord);
    
    
    elev_c[0][c] = coord[0][2];
    if (domain == "surface_star"){
    
      my_centroid.push_back(coord[0]);
      AmanziGeometry::Entity_ID_List nadj_cellids;
      S->GetMesh("surface_star")->cell_get_face_adj_cells(id, AmanziMesh::USED, &nadj_cellids);
      int nface_pcell = S->GetMesh("surface_star")->cell_get_num_faces(id);
      
      int ngb_cells = nadj_cellids.size();
      std::vector<AmanziGeometry::Point> ngb_centroids(ngb_cells);
      std::cout<<"---------ELEVATION EVAL COLUMN----------"<<c<<" "<<id<<"\n";
      
      //get the neighboring cell's centroids
      for(unsigned i=0; i<ngb_cells; i++){
	std::stringstream col_name;
	col_name << "column_" <<nadj_cellids[i];
	std::vector<AmanziGeometry::Point> crd; 
	S->GetMesh(col_name.str())->face_get_coordinates(nfaces-1, &crd);
	ngb_centroids[i] = crd[0];
	std::cout<<"Centroid: "<< my_name.str()<<" "<<col_name.str()<<" "<<my_centroid[c]<<" "<<nadj_cellids[i]<<" "<<crd[0]<< " \n";
      }

      
      
      std::vector<AmanziGeometry::Point> Normal;
      AmanziGeometry::Point N, PQ, PR, Nor_avg(3);
      
      if (ngb_cells >1){
	for (int i=0; i <ngb_cells-1; i++){
	  PQ = my_centroid[c] - ngb_centroids[i];
	  PR = my_centroid[c] - ngb_centroids[i+1];
	  N = PQ^PR;
          if (N[2] < 0)
            N *= -1.; // all normals upward
	  Normal.push_back(N);
	}

	AmanziGeometry::Point fnor = S->GetMesh(my_name.str())->face_normal(nfaces-1);
	std::cout<<"Normal: "<<my_name.str()<<" "<<N<<" "<<fnor<<std::endl;
	Nor_avg = (nface_pcell - Normal.size()) * fnor; 
	for (int i=0; i <Normal.size(); i++)
	  Nor_avg += Normal[i];
	
    	std::cout<<"Neighbor normal: "<<fnor<<" "<<Nor_avg<<" "<<std::endl;
	
	Nor_avg /= nface_pcell;
	std::cout<<"Avg normal: "<<Nor_avg<<" "<<std::endl;
	slope_c[0][c] = (std::sqrt(std::pow(Nor_avg[0],2) + std::pow(Nor_avg[1],2)))/ std::abs(Nor_avg[2]);
        
      }
      else{
	PQ = my_centroid[c] - ngb_centroids[0];
	slope_c[0][c] = std::abs(PQ[2]) / (std::sqrt(std::pow(PQ[0],2) + std::pow(PQ[1],2)));
        
      }
      std::cout<<"SLOPE: "<<c<<" "<<slope_c[0][c]<<"\n";
    }
    else{
      slope_c[0][c]=0; // if domain is surface_column_*, slope is zero.
    }

  } 

  // Set the elevation on faces by getting the corresponding face and its
  // centroid.
  /* if (elev->HasComponent("face")) {
    Epetra_MultiVector& elev_f = *elev->ViewComponent("face", false);
    int nfaces = elev_f.MyLength();
    for (int f=0; f!=nfaces; ++f) {
        // Note that a surface face is a surface mesh's face.
      AmanziMesh::Entity_ID domain_face =
        surface_mesh->entity_get_parent(AmanziMesh::FACE, f);
      AmanziGeometry::Point x = domain_mesh->face_centroid(domain_face);
      
        elev_f[0][f] = x[2];
    }
    }*/

}

} //namespace
} //namespace
