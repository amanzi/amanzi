/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 * Author: Daniil Svyatsky, Xu Chonggang
 *
 * DOCUMENT ME:

 *
 *
 * ------------------------------------------------------------------------- */


#include "fates_pk.hh"

namespace Amanzi {
namespace Vegetation {


FATES_PK::FATES_PK(Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution):
  PK(pk_tree, global_list,  S, solution),
  PK_Physical_Default(pk_tree, global_list,  S, solution)
{
  
  // set up additional primary variables -- this is very hacky...
  // -- surface energy source
  domain_surf_ = plist_->get<std::string>("surface domain name", "surface");

  // set up additional primary variables -- this is very hacky...
  Teuchos::ParameterList& FElist = S->FEList();
  // -- transpiration
  // trans_key_ = Keys::readKey(*plist_, domain_, "transpiration", "transpiration");
  // Teuchos::ParameterList& trans_sublist =
  //     FElist.sublist(trans_key_);
  // trans_sublist.set("field evaluator type", "primary variable");

  // -- shortwave incoming shading
  // shaded_sw_key_ = Keys::readKey(*plist_, domain_surf_, "shaded shortwave radiation", "shaded_shortwave_radiation");
  // Teuchos::ParameterList& sw_sublist =
  //     FElist.sublist(shaded_sw_key_);
  // sw_sublist.set("field evaluator type", "primary variable");

  // // -- lai
  // total_lai_key_ = Keys::readKey(*plist_, domain_surf_, "total leaf area index", "total_leaf_area_index");
  // Teuchos::ParameterList& lai_sublist =
  //     FElist.sublist(total_lai_key_);
  // lai_sublist.set("field evaluator type", "primary variable");  

  //timestep size
  dt_ = plist_->get<double>("max time step", 1.e99);
  patchno_ = plist_->get<int>("number of patches", 1);
  nlevdecomp_ = plist_->get<int>("number of decomposition levels", 1);
  

}

void FATES_PK::Setup(const Teuchos::Ptr<State>& S){

  PK_Physical_Default::Setup(S);

  dt_ = plist_->get<double>("initial time step", 1.);

  // my mesh is the subsurface mesh, but we need the surface mesh, index by column, as well
  mesh_surf_ = S->GetMesh("surface");
  //  soil_part_name_ = plist_->get<std::string>("soil partition name");

  

  
  // requirements: primary variable
  S->RequireField(key_, name_)->SetMesh(mesh_surf_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  
}

  
void FATES_PK::Initialize(const Teuchos::Ptr<State>& S){

  PK_Physical_Default::Initialize(S);

  int iulog=10;
  int masterproc = 0;

  CFI_cdesc_t fatesdesc, clmdesc;
  int retval;
 
  
  char *fates_file = "parameter_files/fates_params_c052019.nc";
  char *clm_file = "parameter_files/clm_params_c180301.nc";

  retval = CFI_establish(&fatesdesc, fates_file, CFI_attribute_other, CFI_type_char, strlen(fates_file), 0, NULL);
  retval = CFI_establish(&clmdesc, clm_file, CFI_attribute_other, CFI_type_char, strlen(clm_file), 0, NULL); 
  fatesreadparameters(&masterproc, &clmdesc, &fatesdesc);


  ncells_owned_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

}


bool FATES_PK::AdvanceStep(double t_old, double t_new, bool reinit){
  

}

  

void FATES_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S){


}


// // helper function for pushing field to column
// void FATES_PK::FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
//         Teuchos::Ptr<Epetra_SerialDenseVector> col_vec, bool copy) {
//   if (col_vec == Teuchos::null) {
//     col_vec = Teuchos::ptr(new Epetra_SerialDenseVector(ncells_per_col_));
//   }

//   ColIterator col_iter(*mesh_, mesh_surf_->entity_get_parent(AmanziMesh::CELL, col), ncells_per_col_);
//   for (std::size_t i=0; i!=col_iter.size(); ++i) {
//     (*col_vec)[i] = vec[col_iter[i]];
//   }
// }

// // helper function for collecting column dz and depth
// void FATES_PK::ColDepthDz_(AmanziMesh::Entity_ID col,
//                             Teuchos::Ptr<Epetra_SerialDenseVector> depth,
//                             Teuchos::Ptr<Epetra_SerialDenseVector> dz) {
//   AmanziMesh::Entity_ID f_above = mesh_surf_->entity_get_parent(AmanziMesh::CELL, col);
//   ColIterator col_iter(*mesh_, f_above, ncells_per_col_);

//   AmanziGeometry::Point surf_centroid = mesh_->face_centroid(f_above);
//   AmanziGeometry::Point neg_z(3);
//   neg_z.set(0.,0.,-1);

//   for (std::size_t i=0; i!=col_iter.size(); ++i) {
//     // depth centroid
//     (*depth)[i] = surf_centroid[2] - mesh_->cell_centroid(col_iter[i])[2];

//     // dz
//     // -- find face_below
//     AmanziMesh::Entity_ID_List faces;
//     std::vector<int> dirs;
//     mesh_->cell_get_faces_and_dirs(col_iter[i], &faces, &dirs);

//     // -- mimics implementation of build_columns() in Mesh
//     double mindp = 999.0;
//     AmanziMesh::Entity_ID f_below = -1;
//     for (std::size_t j=0; j!=faces.size(); ++j) {
//       AmanziGeometry::Point normal = mesh_->face_normal(faces[j]);
//       if (dirs[j] == -1) normal *= -1;
//       normal /= AmanziGeometry::norm(normal);

//       double dp = -normal * neg_z;
//       if (dp < mindp) {
//         mindp = dp;
//         f_below = faces[j];
//       }
//     }

//     // -- fill the val
//     (*dz)[i] = mesh_->face_centroid(f_above)[2] - mesh_->face_centroid(f_below)[2];
//     AMANZI_ASSERT( (*dz)[i] > 0. );
//     f_above = f_below;
//   }

// }



  
}
}
