/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#include <vector>

// Amanzi
#include "MeshAlgorithms.hh"
#include "MFD3D_Diffusion.hh"
#include "nlfv.hh"
#include "ParallelCommunication.hh"

#include "Op_Face_CellBndFace.hh"
#include "OperatorDefs.hh"
#include "Operator_CellBndFace.hh"
#include "PDE_DiffusionNLFVwithBndFaces.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization
* 
* Purpose: Initializes the NLFV diffusion operator with boundary faces.
*          Sets up operator schemas, creates local and global operators,
*          and configures nonlinear coefficient treatment and Newton 
*          correction options.
*
* Key tasks:
*   1. Define operator schemas (data layout for DOFs and stencils)
*   2. Create or validate global operator with cell and boundary_face DOFs
*   3. Create local face-based operators for matrix assembly
*   4. Configure upwind treatment for nonlinear coefficients
*   5. Set up Newton correction for nonlinear problems
*   6. Store mesh dimension information
****************************************************************** */
void
PDE_DiffusionNLFVwithBndFaces::Init_(Teuchos::ParameterList& plist)
{
  // Define stencil for the FV diffusion method.
  // Local operator schema specifies:
  //   - BASE_FACE: matrices are assembled on faces  CHECK!!!
  //   - DOFS_CELL: primary unknowns are cell-centered
  //   - DOFS_BNDFACE: boundary values are explicit DOFs
  local_op_schema_ =
    OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_BNDFACE;

  // create or check the existing Operator
  if (global_op_ == Teuchos::null) {
    // constructor was given a mesh - need to create the global operator
    // Global schema defines where degrees of freedom live (cells + boundary faces)
    global_op_schema_ = OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_BNDFACE;

    // build the Composite Vector Space (CVS) from the global schema
    // CVS defines the data structure for solution vectors
    // Ghosted = true means we include ghost cells for parallel communication
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);
    cvs->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);              // 1 DOF per cell
    cvs->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1); // 1 DOF per boundary face

    // Create specialized operator for cell-boundary_face coupling
    global_op_ = Teuchos::rcp(new Operator_CellBndFace(cvs, plist, global_op_schema_));

  } else {
    // constructor was given an Operator - extract schema and mesh from it
    global_op_schema_ = global_op_->schema();
    mesh_ = global_op_->DomainMap().Mesh();
  }

  // create the local Op and register it with the global Operator
  // Local operator stores face-based 2x2 matrices connecting adjacent cells/boundary faces
  std::string name = "Diffusion: FACE_CELLBNDFACE";
  local_op_ = Teuchos::rcp(new Op_Face_CellBndFace(name, mesh_));
  global_op_->OpPushBack(local_op_);

  // Configure upwind treatment for nonlinear scalar coefficients (e.g., relative permeability)
  // Upwinding determines how to average coefficients across faces based on flux direction
  std::string uwname = plist.get<std::string>("nonlinear coefficient", "upwind: face");
  little_k_ = OPERATOR_LITTLE_K_UPWIND;
  if (uwname == "none") {
    little_k_ = OPERATOR_LITTLE_K_NONE;  // No upwinding - use simple averaging
  }

  // Configure Newton correction terms for nonlinear problems
  // When permeability depends on solution (e.g., k(p)), we need Jacobian corrections
  std::string jacobian = plist.get<std::string>("Newton correction", "none");
  if (jacobian == "none") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_NONE;  // Linear problem
  } else if (jacobian == "approximate Jacobian") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE;

    // Create additional operator to store Jacobian correction terms
    name = "Diffusion: FACE_CELLBNDFACE Jacobian terms";
    jac_op_ = Teuchos::rcp(new Op_Face_CellBndFace(name, mesh_));

    global_op_->OpPushBack(jac_op_);
  } else if (jacobian == "true Jacobian") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_TRUE;
    Errors::Message msg;
    msg << "PDE_DiffusionNLFVwithBndFaces: \"true Jacobian\" not supported -- maybe you mean "
           "\"approximate Jacobian\"?";
    Exceptions::amanzi_throw(msg);
  } else {
    Errors::Message msg;
    msg << "PDE_DiffusionNLFVwithBndFaces: invalid parameter \"" << jacobian
        << "\" for option \"Newton correction\" -- valid are: \"none\", \"approximate Jacobian\"";
    Exceptions::amanzi_throw(msg);
  }

  // Store mesh dimensions
  // dim_ = space dimension (2 or 3)
  // manifold_dim_ = manifold dimension (for surface/fracture meshes embedded in higher dim space)
  dim_ = mesh_->getSpaceDimension();
  manifold_dim_ = mesh_->getManifoldDimension();
}


/* ******************************************************************
* Setup methods: krel and dkdp must be called after calling a
* setup with K absolute
****************************************************************** */
void
PDE_DiffusionNLFVwithBndFaces::SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                                    const Teuchos::RCP<const CompositeVector>& dkdp)
{
  stencil_initialized_ = false;

  k_ = k;
  dkdp_ = dkdp;

  if (k_ != Teuchos::null) {
    if (little_k_ == OPERATOR_LITTLE_K_UPWIND) {
      AMANZI_ASSERT(k_->HasComponent("face"));
    }
  }
  // if (dkdp_ != Teuchos::null) AMANZI_ASSERT(dkdp_->HasComponent("cell"));
}


/* ******************************************************************
* Compute harmonic averaging points (function of geometry and tensor)
* and the positive decomposition of face conormals. The face-based
* data from left and right cells are ordered by the global cells ids.
*
* Purpose: Builds the geometric stencils needed for NLFV discretization.
*          This is a preprocessing step that depends only on mesh geometry
*          and diffusion tensor K, not on the solution.
*
* Algorithm Overview:
*   PART 1: Allocate storage for stencil data
*   PART 2: Compute Harmonic Averaging Points (HAPs) for all faces
*   PART 3: Compute positive decomposition of conormals for interior faces
*   PART 4: Compute positive decomposition for boundary faces
*   PART 5: Parallel communication of stencil data
*
* Key Concepts:
*   - HAP: A special point on each face where flux interpolation is centered
*          Computed to ensure monotonicity and positivity of the scheme
*   - Conormal: K·n where K is diffusion tensor and n is face normal
*   - Positive decomposition: Expresses conormal as positive linear combination
*          of vectors pointing to neighboring cells/faces (ensures M-matrix)
*   - Stencil: For each face, identifies which neighboring cells contribute
*          to the flux, along with their weights
*
* Data Structures Created:
*   stencil_data_: CompositeVector containing:
*     - "hap": Coordinates of harmonic averaging points [dim_ per face]
*     - "gamma": Weighting factor (1-hap_weight) for flux interpolation [1 per face]
*     - "weight": Decomposition weights [2*dim_ per face, dim_ for each side]
*     - "flux_data": Precomputed flux coefficients [2*dim_ per face]
*   stencil_faces_[2*dim_]: For each face, lists which faces are in its stencil
*   stencil_cells_[2*dim_]: For each face, lists which cells are in its stencil
*
* Storage pattern for each face f:
*   Indices [0, dim_-1]: Stencil for "left" cell (lower global ID)
*   Indices [dim_, 2*dim_-1]: Stencil for "right" cell or boundary
****************************************************************** */
void
PDE_DiffusionNLFVwithBndFaces::InitStencils_()
{
  // ==================== PART 1: ALLOCATE PERSISTENT MEMORY ====================
  
  // Create storage for all stencil-related data
  // This data persists across multiple UpdateMatrices() calls
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_)
    ->SetGhosted(true)  // Include ghost cells for parallel computing
    ->AddComponent("hap", AmanziMesh::Entity_kind::FACE, dim_)           // HAP coordinates
    ->AddComponent("gamma", AmanziMesh::Entity_kind::FACE, 1)            // Interpolation weight
    ->AddComponent("weight", AmanziMesh::Entity_kind::FACE, 2 * dim_)   // Decomposition weights
    ->AddComponent("flux_data", AmanziMesh::Entity_kind::FACE, 2 * dim_); // Flux coefficients
  stencil_data_ = Teuchos::rcp(new CompositeVector(cvs));

  // Get references to the components for easy access
  Epetra_MultiVector& hap = *stencil_data_->ViewComponent("hap", true);
  Epetra_MultiVector& gamma = *stencil_data_->ViewComponent("gamma", true);
  Epetra_MultiVector& weight = *stencil_data_->ViewComponent("weight", true);

  stencil_data_->PutScalarMasterAndGhosted(0.0);

  // Allocate integer arrays to store stencil topology
  // For each face, we store up to 2*dim_ stencil entries (dim_ per side)
  // In 3D: up to 6 entries (3 for left cell, 3 for right cell/boundary)
  stencil_faces_.resize(2 * dim_);
  stencil_cells_.resize(2 * dim_);
  for (int i = 0; i < 2 * dim_; ++i) {
    stencil_faces_[i] =
      Teuchos::rcp(new Epetra_IntVector(mesh_->getMap(AmanziMesh::Entity_kind::FACE, true)));
    stencil_cells_[i] =
      Teuchos::rcp(new Epetra_IntVector(mesh_->getMap(AmanziMesh::Entity_kind::FACE, true)));

    stencil_faces_[i]->PutValue(0);
    stencil_cells_[i]->PutValue(0);
  }

  // ==================== PREPARE DIFFUSION TENSOR DATA ====================
  
  // allocate temporary memory for distributed tensor
  // We need to communicate tensor K to ghost cells for stencil computation
  CompositeVectorSpace cvs_tmp;
  cvs_tmp.SetMesh(mesh_)->SetGhosted(true)->AddComponent(
    "tensor", AmanziMesh::Entity_kind::CELL, dim_ * dim_);
  Teuchos::RCP<CompositeVector> cv_tmp = Teuchos::rcp(new CompositeVector(cvs_tmp));
  Epetra_MultiVector& Ktmp = *cv_tmp->ViewComponent("tensor", true);

  // instantiate variables to access supporting tools
  WhetStone::NLFV nlfv(mesh_);           // NLFV geometric algorithms
  WhetStone::MFD3D_Diffusion mfd3d(mesh_); // MFD discretization tools (if needed)

  // Convert diffusion tensor from tensor format to vector format for storage
  // For each cell, flatten K matrix into a vector [K11, K12, K13, K21, K22, ...]
  WhetStone::DenseVector data(dim_ * dim_);
  for (int c = 0; c < ncells_owned; ++c) {
    if (K_ != Teuchos::null) {
      WhetStone::TensorToVector((*K_)[c], data);  // Use provided tensor
    } else {
      // Default to identity tensor if none provided
      WhetStone::Tensor Kc(dim_, 1);
      Kc(0, 0) = 1.0;
      WhetStone::TensorToVector(Kc, data);
    }

    for (int i = 0; i < dim_ * dim_; ++i) {
      Ktmp[i][c] = data(i);
    }
  }
  // Communicate tensor data to ghost cells (needed for faces on processor boundaries)
  cv_tmp->ScatterMasterToGhosted();

  // ==================== PART 2: COMPUTE HARMONIC AVERAGING POINTS ====================
  
  // For each interior face, compute the HAP based on geometry and diffusion tensor
  // HAP is located between cell centroids such that the discrete flux is consistent
  int c1, c2;
  double hap_weight;
  WhetStone::Tensor T(dim_, 2);
  AmanziGeometry::Point Kn1(dim_), Kn2(dim_), p(dim_);

  for (int f = 0; f < nfaces_owned; f++) {
    auto cells = mesh_->getFaceCells(f);
    int ncells = cells.size();

    if (ncells == 2) {
      // Interior face: compute HAP between two cells
      const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
      OrderCellsByGlobalId_(cells, c1, c2);  // Ensure consistent ordering

      // Compute conormals Kn1 = K1·n and Kn2 = K2·n for both cells
      for (int i = 0; i < dim_ * dim_; ++i) data(i) = Ktmp[i][c1];
      VectorToTensor(data, T);
      Kn1 = T * normal;

      for (int i = 0; i < dim_ * dim_; ++i) data(i) = Ktmp[i][c2];
      VectorToTensor(data, T);
      Kn2 = T * normal;

      // Compute HAP location p and weight using NLFV algorithm
      // hap_weight determines how close HAP is to cell c1 (0) vs c2 (1)
      nlfv.HarmonicAveragingPoint(f, c1, c2, Kn1, Kn2, p, hap_weight);
    } else {
      // Boundary face: HAP is at the face centroid
      p = mesh_->getFaceCentroid(f);
      hap_weight = 0.0;
    }

    // Store HAP coordinates
    for (int i = 0; i < dim_; ++i) hap[i][f] = p[i];
    // gamma = (1 - hap_weight) is used later for flux interpolation
    gamma[0][f] = 1.0 - hap_weight;
  }

  // Communicate HAP data to ghost faces
  stencil_data_->ScatterMasterToGhosted("hap");
  stencil_data_->ScatterMasterToGhosted("gamma");

  // ==================== PART 3: POSITIVE DECOMPOSITION FOR INTERIOR FACES ====================
  
  // For each cell, decompose the conormal of each adjacent face as a positive
  // linear combination of vectors pointing to neighboring cells/faces
  // This ensures the scheme produces an M-matrix (diagonal dominant with positive diagonal)
  
  AmanziGeometry::Point conormal(dim_), v(dim_);
  std::vector<AmanziGeometry::Point> tau;  // Candidate vectors for decomposition

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);

    // Get all faces of this cell and their orientations
    const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
    int nfaces = faces.size();

    // Get diffusion tensor for this cell
    WhetStone::Tensor Kc(mesh_->getSpaceDimension(), 1);
    Kc(0, 0) = 1.0;
    if (K_.get() ) Kc = (*K_)[c];

    // Build candidate vectors tau: from cell centroid to HAP of each face
    // These form the basis for positive decomposition
    tau.clear();
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      for (int i = 0; i < dim_; ++i) v[i] = hap[i][f] - xc[i];
      tau.push_back(v);
    }

    // For each face of this cell, decompose its conormal
    int ierr, ids[3];   // ids stores which tau vectors are selected
    double ws[3];       // ws stores the weights for the decomposition
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
      // Conormal = K·n with correct orientation for this cell
      conormal = ((*K_)[c] * normal) * dirs[n];

      // Find positive decomposition: conormal = sum(ws[i] * tau[ids[i]])
      // where ws[i] >= 0 for all i
      ierr = nlfv.PositiveDecomposition(n, tau, conormal, manifold_dim_, ws, ids);
      AMANZI_ASSERT(ierr == 0);

      // Store decomposition weights and stencil connectivity
      auto cells = mesh_->getFaceCells(f);
      OrderCellsByGlobalId_(cells, c1, c2);
      // k = 0 for left cell, dim_ for right cell (storage offset)
      int k = (c == c1) ? 0 : dim_;

      for (int i = 0; i < dim_; i++) {
        weight[k + i][f] = ws[i];                    // Store weight
        (*stencil_faces_[k + i])[f] = faces[ids[i]]; // Which face in stencil
        // Which cell is across that face from current cell c
        (*stencil_cells_[k + i])[f] = AmanziMesh::getFaceAdjacentCell(*mesh_, c, faces[ids[i]]);
      }
    }
  }

  // ==================== PART 4: POSITIVE DECOMPOSITION FOR BOUNDARY FACES ====================
  
  // Boundary faces require special treatment because they have only one adjacent cell
  // We build vectors from the boundary HAP to interior cell HAP and the cell centroid
  
  for (int f = 0; f < nfaces_owned; f++) {
    auto cells = mesh_->getFaceCells(f);
    int ncells = cells.size();
    if (ncells == 1) {
      // This is a boundary face
      int c = cells[0];
      const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
      const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
      int nfaces = faces.size();

      // Get diffusion tensor for the adjacent cell
      WhetStone::Tensor Kc(mesh_->getSpaceDimension(), 1);
      Kc(0, 0) = 1.0;
      if (K_.get() ) Kc = (*K_)[c];
      
      int face_itself;  // Index of boundary face f within cell c's face list

      // Build candidate vectors tau
      // For boundary face itself: vector from HAP to cell centroid (inward pointing)
      // For other faces: vector from boundary HAP to their HAPs
      tau.clear();
      for (int n = 0; n < nfaces; n++) {
        int g = faces[n];
        if (g == f) {
          face_itself = n;
          // Inward vector: from boundary HAP to cell centroid
          for (int i = 0; i < dim_; ++i) v[i] = xc[i] - hap[i][f];
        } else {
          // Vector to neighboring face's HAP
          for (int i = 0; i < dim_; ++i) v[i] = hap[i][g] - hap[i][f];
        }
        tau.push_back(v);
      }

      // Compute conormal for boundary face (negative orientation, pointing outward from cell)
      const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
      conormal = -(Kc * normal) * dirs[face_itself];

      int ierr, ids[3];
      double ws[3];

      // Positive decomposition for boundary face
      ierr = nlfv.PositiveDecomposition(face_itself, tau, conormal, manifold_dim_, ws, ids);
      AMANZI_ASSERT(ierr == 0);

      // Store decomposition data
      // First component (i=0) corresponds to cell itself (tau vector points to cell center)
      weight[dim_][f] = ws[0];
      (*stencil_faces_[dim_])[f] = -100;  // Special marker: stencil includes cell, not a face
      (*stencil_cells_[dim_])[f] = c;     // The adjacent cell
      
      // Remaining components (i=1 to dim_-1) correspond to neighboring faces
      for (int i = 1; i < dim_; i++) {
        weight[dim_ + i][f] = ws[i];
        (*stencil_faces_[dim_ + i])[f] = faces[ids[i]];
        (*stencil_cells_[dim_ + i])[f] = AmanziMesh::getFaceAdjacentCell(*mesh_, c, faces[ids[i]]);
      }
    }
  }

  // ==================== PART 5: PARALLEL COMMUNICATION OF STENCIL DATA ====================
  
  // In parallel computations, faces on processor boundaries need stencil information
  // from neighboring processors. We gather/scatter to synchronize this data.
  
  // Synchronize weights across processors
  stencil_data_->GatherGhostedToMaster("weight");
  stencil_data_->ScatterMasterToGhosted("weight");

  // Synchronize stencil topology (which faces and cells are in each stencil)
  ParallelCommunication pp(mesh_);
  for (int i = 0; i < 2 * dim_; ++i) {
    // Combine contributions from ghost to master entities
    pp.CombineGhostEntity2MasterEntity(
      AmanziMesh::Entity_kind::FACE, *stencil_faces_[i], (Epetra_CombineMode)Add);
    pp.CombineGhostEntity2MasterEntity(
      AmanziMesh::Entity_kind::FACE, *stencil_cells_[i], (Epetra_CombineMode)Add);

    // Copy master entity data back to ghost entities
    pp.CopyMasterEntity2GhostEntity(AmanziMesh::Entity_kind::FACE, *stencil_faces_[i]);
    pp.CopyMasterEntity2GhostEntity(AmanziMesh::Entity_kind::FACE, *stencil_cells_[i]);
  }


  // std::cout<<"stencil 9: ";
  // for (int n=0; n<6; n++) std::cout<<(*stencil_faces_[n])[9]<<" ";
  // std::cout<<"\n";

  // Mark stencils as initialized - this expensive computation only done once
  stencil_initialized_ = true;
}


/* ******************************************************************
* Populate face-based 2x2 matrices on interior faces 
* We avoid round-off operations since the stencils
* already incorporate them.
****************************************************************** */
void
PDE_DiffusionNLFVwithBndFaces::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                              const Teuchos::Ptr<const CompositeVector>& u)
{
  if (!stencil_initialized_) InitStencils_();
  if (k_ != Teuchos::null) k_->ScatterMasterToGhosted("face");

  u->ScatterMasterToGhosted("cell");

  Epetra_MultiVector& hap_gamma = *stencil_data_->ViewComponent("gamma", true);
  Epetra_MultiVector& weight = *stencil_data_->ViewComponent("weight", true);
  Epetra_MultiVector& flux_data = *stencil_data_->ViewComponent("flux_data", true);

  // allocate auxiliary matrix structure
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_)->SetGhosted(true)->AddComponent("face", AmanziMesh::Entity_kind::FACE, 2);
  CompositeVector matrix_cv(cvs), sideflux_cv(cvs);

  Epetra_MultiVector& matrix = *matrix_cv.ViewComponent("face", true);
  Epetra_MultiVector& sideflux = *sideflux_cv.ViewComponent("face", true);

  // calculate one-sides flux corrections. Since a flux stencil can
  // occupy (dim_ + 1) cells, we need parallel communications.
  OneSidedFluxCorrections_(1, *u, sideflux_cv);
  //std::cout<<"sideflux\n"<<sideflux<<"\n";

  //OneSidedNeumannCorrections_(*u, sideflux_cv);
  sideflux_cv.GatherGhostedToMaster();
  sideflux_cv.ScatterMasterToGhosted();

  // un-rolling little-k data
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  if (k_ != Teuchos::null) k_face = k_->ViewComponent("face");

  // split each stencil between different local matrices
  int c1, c2, c3, c4, k1, k2;


  matrix_cv.PutScalarMasterAndGhosted(0.0);
  flux_data.PutScalar(0.0);

  for (int c = 0; c < ncells_owned; ++c) {
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      auto cells = mesh_->getFaceCells(f);
      int ncells = cells.size();

      OrderCellsByGlobalId_(cells, c1, c2);
      k1 = (c1 == c) ? 0 : 1;
      k2 = k1 * dim_;

      // calculate little_k on the current face
      double kf(1.0);
      if (k_face.get() ) kf = (*k_face)[0][f];

      // Calculate solution-dependent weigths using corrections to the
      // two-point flux. Note mu does not depend on one-sided flux.
      double gamma, g1, g2, gg(-1.0), w1, w2(0.0), tpfa, mu(1.0), ntpfa1, ntpfa2;

      gamma = hap_gamma[0][f];
      if (ncells == 2) {
        w1 = weight[0][f] * gamma;
        w2 = weight[dim_][f] * (1.0 - gamma);

        g1 = sideflux[0][f];
        g2 = sideflux[1][f];
        gg = g1 * g2;


        g1 = fabs(g1);
        g2 = fabs(g2);
        mu = (g1 + g2 == 0.0) ? 0.5 : g2 / (g1 + g2);
        
        tpfa = mu * w1 + (1.0 - mu) * w2;

        //if (f==9) std::cout<<"matrix f9 "<<"c"<<c<<" mu "<<mu<<" w1 "<<w1<<" w2 "<<w2<<" kf "<<kf<<"\n";
        //if (f==9) std::cout<<"matrix f9 "<<"c"<<c<<" "<<matrix[0][9]<<" "<<matrix[1][9]<<"\n";
        
        matrix[k1][f] += kf * tpfa;
        flux_data[k2][f] = kf * tpfa;

        //if (f==9) std::cout<<"matrix f9 "<<"c"<<c<<" "<<matrix[0][9]<<" "<<matrix[1][9]<<"\n";

      } else {
        g1 = sideflux[0][f];
        g2 = sideflux[1][f];
        mu = (g1 + g2 == 0.0) ? 0.5 : g2 / (g1 + g2);
        w1 = weight[0][f];
        w2 = weight[dim_][f];

        double tc1, tc2;

        NLTPFAContributions_(f, tc1, tc2);

        ntpfa1 = mu * w1 + (1.0 - mu) * w2 + mu * tc1;
        ntpfa2 = mu * w1 + (1.0 - mu) * w2 + (1.0 - mu) * tc2;

        gg = 0.;
        matrix[0][f] += kf * ntpfa1;
        matrix[1][f] += kf * ntpfa2;
        flux_data[0][f] = kf * ntpfa1;
        flux_data[dim_][f] = kf * ntpfa2;

        //}
      }

      // remaining terms of one-sided flux in cell c. Now we need
      // to select mu depending on the one-sided flux.
      if (gg < 0.0) {
        if (c1 != c) mu = 1.0 - mu;
        for (int i = 1; i < dim_; i++) {
          int f1 = (*stencil_faces_[i + k2])[f];
          auto cells_tmp = mesh_->getFaceCells(f1);

          gamma = hap_gamma[0][f1];
          OrderCellsByGlobalId_(cells_tmp, c3, c4);

          //if (f1==9) std::cout<<"Aux matrix c="<<c<<" f="<<f<<" stencil "<<(*stencil_faces_[1 + k2])[f]<<" "<<(*stencil_faces_[2 + k2])[f] <<"\n";
          
          k1 = 0;
          if (c3 != c) {
            gamma = 1.0 - gamma;
            k1 = 1;
          }

          double tmp = ncells * weight[i + k2][f] * gamma * mu;
          matrix[k1][f1] += kf * tmp;
          flux_data[k2 + i][f] = kf * tmp;
          
          if (cells_tmp.size() == 1) {
            int kb = 1 - k1;
            matrix[kb][f1] += kf * tmp;
          }

        }
      }
    }
  }

  //std::cout<<"matrix 9 " <<matrix[0][9]<<" "<<matrix[1][9]<<"\n";
  //exit(0);
  
  stencil_data_->GatherGhostedToMaster("flux_data");
  matrix_cv.GatherGhostedToMaster();

  // populate local matrices
  for (int f = 0; f < nfaces_owned; ++f) {
    auto cells = mesh_->getFaceCells(f);
    int ncells = cells.size();

    WhetStone::DenseMatrix Aface(2, 2);

    if (ncells == 2) {
      k1 = OrderCellsByGlobalId_(cells, c3, c4);
      k2 = 1 - k1;
      Aface(0, 0) = matrix[k1][f];
      Aface(0, 1) = -matrix[k1][f];

      Aface(1, 0) = -matrix[k2][f];
      Aface(1, 1) = matrix[k2][f];
    } else {
      Aface(0, 0) = matrix[0][f];
      Aface(0, 1) = -matrix[1][f];
      Aface(1, 1) = flux_data[dim_][f];
      Aface(1, 0) = -flux_data[0][f];
    }

    //    if (ncells==2) std::cout<<"Aface "<<f<<" "<<ncells<<"\n"<<Aface<<"\n";
    
    local_op_->matrices[f] = Aface;
  }

  //exit(0);
}


/* ******************************************************************
* Modify operator by adding an upwind approximation of the Newton
* correction term.
****************************************************************** */
void
PDE_DiffusionNLFVwithBndFaces::UpdateMatricesNewtonCorrection(
  const Teuchos::Ptr<const CompositeVector>& flux,
  const Teuchos::Ptr<const CompositeVector>& u,
  double scalar_factor)
{
  // ignore correction if no flux provided.
  if (flux == Teuchos::null) return;

  // Correction is zero for linear problems
  if (k_ == Teuchos::null || dkdp_ == Teuchos::null) return;

  if (k_->HasComponent("face") ) k_->ScatterMasterToGhosted("face");
  if (dkdp_->HasComponent("face") ) dkdp_->ScatterMasterToGhosted("face");

  // Correction is not required
  if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_NONE) return;

  // only works on upwinded methods
  if (little_k_ == OPERATOR_UPWIND_NONE) return;

  const Epetra_MultiVector& kf = *k_->ViewComponent("face");
  const Epetra_MultiVector& dkdp_f = *dkdp_->ViewComponent("face");
  const Epetra_MultiVector& flux_f = *flux->ViewComponent("face");

  // populate the local matrices
  for (int f = 0; f < nfaces_owned; f++) {
    auto cells = mesh_->getFaceCells(f);
    int ncells = cells.size();
    WhetStone::DenseMatrix Aface(2, 2);
    Aface.PutScalar(0.0);

    // We use the upwind discretization of the generalized flux.
    double v = std::abs(kf[0][f]) > 0.0 ? flux_f[0][f] * dkdp_f[0][f] / kf[0][f] : 0.0;
    double vmod = std::abs(v);

    // prototype for future limiters
    vmod *= scalar_factor;

    // We use the upwind discretization of the generalized flux.
    int i, dir, c1;
    c1 = cells[0];
    mesh_->getFaceNormal(f, c1, &dir);
    i = (v * dir >= 0.0) ? 0 : 1;

    if (ncells == 2) {
      Aface(i, i) = vmod;
      Aface(1 - i, i) = -vmod;
    } else if (i == 0) {
      Aface(0, 0) = vmod;
    }

    jac_op_->matrices[f] = Aface;
  }
}


/* ******************************************************************
* Modify operator by adding an upwind approximation of the Newton
* correction term.
****************************************************************** */
void
PDE_DiffusionNLFVwithBndFaces::UpdateMatricesNewtonCorrection(
  const Teuchos::Ptr<const CompositeVector>& flux,
  const Teuchos::Ptr<const CompositeVector>& u,
  const Teuchos::Ptr<const CompositeVector>& factor)
{
  // ignore correction if no flux provided.
  if (flux == Teuchos::null) return;

  // Correction is zero for linear problems
  if (k_ == Teuchos::null || dkdp_ == Teuchos::null) return;

  // Correction is not required
  if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_NONE) return;

  // only works on upwinded methods
  if (little_k_ == OPERATOR_UPWIND_NONE) return;

  const Epetra_MultiVector& kf = *k_->ViewComponent("face");
  const Epetra_MultiVector& dkdp_f = *dkdp_->ViewComponent("face");
  const Epetra_MultiVector& flux_f = *flux->ViewComponent("face");
  const Epetra_MultiVector& factor_cell = *factor->ViewComponent("cell");

  // populate the local matrices
  double v, vmod;

  for (int f = 0; f < nfaces_owned; f++) {
    auto cells = mesh_->getFaceCells(f);
    int ncells = cells.size();
    WhetStone::DenseMatrix Aface(2, 2);
    Aface.PutScalar(0.0);

    // We use the upwind discretization of the generalized flux.
    v = std::abs(kf[0][f]) > 0.0 ? flux_f[0][f] * dkdp_f[0][f] / kf[0][f] : 0.0;
    vmod = std::abs(v);

    double scalar_factor = 0.;
    for (int j = 0; j < ncells; j++) scalar_factor += factor_cell[0][cells[j]];
    scalar_factor *= 1.0 / ncells;

    // prototype for future limiters
    vmod *= scalar_factor;

    // We use the upwind discretization of the generalized flux.
    int i, dir, c1;
    c1 = cells[0];
    mesh_->getFaceNormal(f, c1, &dir);
    i = (v * dir >= 0.0) ? 0 : 1;

    if (ncells == 2) {
      Aface(i, i) = vmod;
      Aface(1 - i, i) = -vmod;
    } else if (i == 0) {
      Aface(0, 0) = vmod;
    }

    jac_op_->matrices[f] = Aface;
  }
}


/* ******************************************************************
* Calculate one-sided fluxes (i0=0) or flux corrections (i0=1).
****************************************************************** */
void
PDE_DiffusionNLFVwithBndFaces::OneSidedFluxCorrections_(int i0,
                                                        const CompositeVector& u,
                                                        CompositeVector& flux_cv)
{
  // un-rolling composite vectors
  const Epetra_MultiVector& uc = *u.ViewComponent("cell", true);
  const Epetra_MultiVector& ubnd = *u.ViewComponent("boundary_face", true);
  Epetra_MultiVector& flux = *flux_cv.ViewComponent("face", true);

  // std::cout<<uc<<"\n";
  // std::cout<<ubnd<<"\n";
  
  Epetra_MultiVector& hap_gamma = *stencil_data_->ViewComponent("gamma", true);
  Epetra_MultiVector& weight = *stencil_data_->ViewComponent("weight", true);

  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  // un-rolling little-k data
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  if (k_ != Teuchos::null) k_face = k_->ViewComponent("face");

  int c1, c2, c3, k1, k2;
  double gamma, tmp;

  flux_cv.PutScalarMasterAndGhosted(0.0);

  for (int c = 0; c < ncells_owned; ++c) {
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      auto cells = mesh_->getFaceCells(f);

      OrderCellsByGlobalId_(cells, c1, c2);
      k1 = (c1 == c) ? 0 : 1;
      k2 = k1 * dim_;

      // scalar (nonlinear) coefficient
      double kf(1.0);
      if (k_face.get() ) kf = (*k_face)[0][f];

      double sideflux(0.0);
      double neumanflux(0.0);
      
      for (int i = i0; i < dim_; ++i) {
        int f1 = (*stencil_faces_[i + k2])[f];
        c3 = (*stencil_cells_[i + k2])[f];
        if (c3 >= 0) {
          auto cells_tmp = mesh_->getFaceCells(f1);
          OrderCellsByGlobalId_(cells_tmp, c1, c2);

          gamma = hap_gamma[0][f1];
          if (c1 != c) gamma = 1.0 - gamma;

          tmp = weight[i + k2][f] * gamma;
          sideflux += tmp * (uc[0][c] - uc[0][c3]);
          //std::cout<<"f "<<f<<" c "<<c<<" c3 "<<c3<<" sideflux "<<sideflux<<"\n";
        } else if (bc_model[f1] == OPERATOR_BC_DIRICHLET) {
          tmp = weight[i + k2][f];
          int bf = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false)
                     .LID(mesh_->getMap(AmanziMesh::Entity_kind::FACE, false).GID(f1));

          ubnd[0][bf] = bc_value[f1];
          
          sideflux += tmp * (uc[0][c] - ubnd[0][bf]); // * dir;
        } else if (bc_model[f1] == OPERATOR_BC_NEUMANN) {
          tmp = weight[i + k2][f];
          int bf = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false)
                     .LID(mesh_->getMap(AmanziMesh::Entity_kind::FACE, false).GID(f1));
          //sideflux += tmp * (uc[0][c] - ubnd[0][bf]); //* dir;
          neumanflux += tmp * bc_value[f1] * mesh_->getFaceArea(f1);

        }
      }

      //if (f==9) std::cout<<"k1 "<<k1<<" sideflux "<<sideflux<<"\n";
      
      flux[k1][f] = kf * sideflux;
    }
  }
}

void
PDE_DiffusionNLFVwithBndFaces::OneSidedNeumannCorrections_(const CompositeVector& u,
                                                           CompositeVector& flux_cv)
{
  // un-rolling composite vectors
  const Epetra_MultiVector& uc = *u.ViewComponent("cell", true);
  const Epetra_MultiVector& ubnd = *u.ViewComponent("boundary_face", true);
  Epetra_MultiVector& flux = *flux_cv.ViewComponent("face", true);

  Epetra_MultiVector& hap_gamma = *stencil_data_->ViewComponent("gamma", true);
  Epetra_MultiVector& weight = *stencil_data_->ViewComponent("weight", true);

  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();
  
  // un-rolling little-k data
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  if (k_ != Teuchos::null) k_face = k_->ViewComponent("face");

  int c1, c2, c3, k1, k2;
  double gamma, tmp;
  std::vector<int> dirs;

  for (int f = 0; f < nfaces_owned; f++) {
    if ((bc_model[f] == OPERATOR_BC_NEUMANN) || (bc_model[f] == OPERATOR_BC_DIRICHLET)) {
      auto cells = mesh_->getFaceCells(f);
      AMANZI_ASSERT(cells.size() == 1);

      k1 = 1;
      k2 = dim_;

      // scalar (nonlinear) coefficient
      double kf(1.0);
      if (k_face.get() ) kf = (*k_face)[0][f];

      for (k1 = 0; k1 < 2; k1++) {
        k2 = k1 * dim_;
        double sideflux(0.0);
        for (int i = 1; i < dim_; ++i) {
          int f1 = (*stencil_faces_[i + k2])[f];
          c3 = (*stencil_cells_[i + k2])[f];
          if (f1 >= 0) {
            if (c3 >= 0) {
              auto cells_tmp = mesh_->getFaceCells(f1);
              OrderCellsByGlobalId_(cells_tmp, c1, c2);

              gamma = hap_gamma[0][f1];

              double uf = (1 - gamma) * uc[0][c1] + gamma * uc[0][c2];
              sideflux += weight[i + k2][f] * uf;

            } else if (bc_model[f1] == OPERATOR_BC_DIRICHLET) {
              tmp = weight[i + k2][f];
              int bf = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false)
                         .LID(mesh_->getMap(AmanziMesh::Entity_kind::FACE, false).GID(f1));
              sideflux += tmp * ubnd[0][bf];

            } else if (bc_model[f1] == OPERATOR_BC_NEUMANN) {
              tmp = weight[i + k2][f];
              int bf = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false)
                         .LID(mesh_->getMap(AmanziMesh::Entity_kind::FACE, false).GID(f1));
              sideflux += tmp * ubnd[0][bf];
            }
          }
        }
        flux[k1][f] = kf * sideflux;
      }
    }
  }
  // flux_cv.GatherGhostedToMaster();
  // flux_cv.ScatterMasterToGhosted();
}


/* ******************************************************************
* Calculate one-sided fluxes (i0=0) or flux corrections (i0=1).
****************************************************************** */
void
PDE_DiffusionNLFVwithBndFaces::OneSidedWeightFluxes_(int i0,
                                                     const CompositeVector& u,
                                                     CompositeVector& flux_cv)
{
  // un-rolling composite vectors
  const Epetra_MultiVector& uc = *u.ViewComponent("cell", true);
  const Epetra_MultiVector& ubnd = *u.ViewComponent("boundary_face", true);
  Epetra_MultiVector& flux = *flux_cv.ViewComponent("face", true);

  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();
  
  Epetra_MultiVector& flux_data = *stencil_data_->ViewComponent("flux_data", true);

  int c1, c2, c3, k1, k2;

  flux_cv.PutScalarMasterAndGhosted(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      auto cells = mesh_->getFaceCells(f);

      if (cells.size() > 1) {
        OrderCellsByGlobalId_(cells, c1, c2);
        k1 = (c1 == c) ? 0 : 1;
        k2 = k1 * dim_;

        double sideflux(0.0);
        double neumann_flux(0.0);
        for (int i = i0; i < dim_; ++i) {
          c3 = (*stencil_cells_[i + k2])[f];
          int f1 = (*stencil_faces_[i + k2])[f];

          if (c3 >= 0) {
            sideflux += flux_data[i + k2][f] * (uc[0][c] - uc[0][c3]);
          } else if (bc_model[f1] == OPERATOR_BC_DIRICHLET) {
            int bf = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false)
                       .LID(mesh_->getMap(AmanziMesh::Entity_kind::FACE, false).GID(f1));
            
            ubnd[0][bf] = bc_value[f1];
            
            sideflux += flux_data[i + k2][f] * (uc[0][c] - ubnd[0][bf]);
          } else if (bc_model[f1] == OPERATOR_BC_NEUMANN) {
            int bf = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false)
                       .LID(mesh_->getMap(AmanziMesh::Entity_kind::FACE, false).GID(f1));

            //sideflux += flux_data[i + k2][f] * (uc[0][c] - ubnd[0][bf]);
            neumann_flux += flux_data[i + k2][f] * bc_value[f1] * mesh_->getFaceArea(f1);
          }
        }
        flux[k1][f] = sideflux + neumann_flux;
      } else if (cells.size() == 1) {
        int bf = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false)
                   .LID(mesh_->getMap(AmanziMesh::Entity_kind::FACE, false).GID(f));


        flux[0][f] = flux_data[0][f] * uc[0][c] - flux_data[dim_][f] * ubnd[0][bf];
      }
    }
  }

  flux_cv.GatherGhostedToMaster();
  flux_cv.ScatterMasterToGhosted();
}


/* ******************************************************************
* Matrix-based implementation of boundary conditions.
****************************************************************** */
void
PDE_DiffusionNLFVwithBndFaces::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  Epetra_MultiVector& rhs_bnd = *global_op_->rhs()->ViewComponent("boundary_face");
  Epetra_MultiVector& hap = *stencil_data_->ViewComponent("hap", true);

  // un-rolling little-k data
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  if (k_ != Teuchos::null) k_face = k_->ViewComponent("face");

  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] != OPERATOR_BC_NONE) {
      auto cells = mesh_->getFaceCells(f);
      int c = cells[0];
      int bf = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false)
                 .LID(mesh_->getMap(AmanziMesh::Entity_kind::FACE, false).GID(f));

      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        WhetStone::DenseMatrix& Aface = local_op_->matrices[f];
        local_op_->matrices_shadow[f] = Aface;

        Aface(1, 0) = 0.;
        // Aface(1,1) = 1e-10;
        rhs_bnd[0][bf] = Aface(1, 1) * bc_value[f];
      } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
        WhetStone::DenseMatrix& Aface = local_op_->matrices[f];
        local_op_->matrices_shadow[f] = Aface;

        WhetStone::Tensor Kc(dim_, 1);
        Kc(0, 0) = 1.0;
        if (K_.get() ) Kc = (*K_)[c];

        const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
        AmanziGeometry::Point v(dim_);

        for (int i = 0; i < dim_; ++i) v[i] = hap[i][f] - xc[i];

        double ub = bc_value[f] * mesh_->getFaceArea(f);
        rhs_bnd[0][bf] -= ub;
      }
    }
  }

  return;
}


/* ******************************************************************
* Calculate flux using cell-centered data.
* **************************************************************** */
void
PDE_DiffusionNLFVwithBndFaces::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                                          const Teuchos::Ptr<CompositeVector>& flux)
{
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_)->SetGhosted(true)->AddComponent("face", AmanziMesh::Entity_kind::FACE, 2);
  CompositeVector wgt_sideflux_cv(cvs);

  Epetra_MultiVector& wgt_sideflux = *wgt_sideflux_cv.ViewComponent("face", true);
  Epetra_MultiVector& flux_data = *flux->ViewComponent("face", true);

  u->ScatterMasterToGhosted("cell");
  OneSidedWeightFluxes_(0, *u, wgt_sideflux_cv);

  int c1, c2, dir;

  for (int f = 0; f < nfaces_owned; ++f) {
    if ((bc_model[f] == OPERATOR_BC_DIRICHLET)) {
      auto cells = mesh_->getFaceCells(f);
      mesh_->getFaceNormal(f, cells[0], &dir);
      flux_data[0][f] = wgt_sideflux[0][f] * dir;
    } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
      //flux_data[0][f] = bc_value[f] * mesh_->getFaceArea(f);
      auto cells = mesh_->getFaceCells(f);
      mesh_->getFaceNormal(f, cells[0], &dir);
      flux_data[0][f] = wgt_sideflux[0][f] * dir;
    } else if (bc_model[f] == OPERATOR_BC_NONE) {
      auto cells = mesh_->getFaceCells(f);
      OrderCellsByGlobalId_(cells, c1, c2);
      mesh_->getFaceNormal(f, c1, &dir);

      double wg1 = wgt_sideflux[0][f];
      double wg2 = wgt_sideflux[1][f];

      // if (cells.size() == 2) flux_data[0][f] = 0.5 * (wg1 - wg2) * dir;
      // else flux_data[0][f] = dir * wg1; 
      flux_data[0][f] = 0.5 * (wg1 - wg2) * dir;  
    }
  }
}


/* ******************************************************************
* Order cells by their global ids. Returns 1 if cells were swapped.
****************************************************************** */
int
PDE_DiffusionNLFVwithBndFaces::OrderCellsByGlobalId_(const AmanziMesh::cEntity_ID_View& cells,
                                                     int& c1,
                                                     int& c2)
{
  c1 = cells[0];
  c2 = -1;

  int ncells = cells.size();
  if (ncells == 1) return 0;

  c2 = cells[1];
  if (mesh_->getMap(AmanziMesh::Entity_kind::CELL, true).GID(c1) >
      mesh_->getMap(AmanziMesh::Entity_kind::CELL, true).GID(c2)) {
    int c(c1);
    c1 = c2;
    c2 = c;
    return 1;
  }

  return 0;
}


/* ******************************************************************
* TBW
****************************************************************** */
int
PDE_DiffusionNLFVwithBndFaces::NLTPFAContributions_(int f, double& tc1, double& tc2)
{
  int c, c1, c2, c3, f1;
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();

  //AMANZI_ASSERT(bc_model[f] == OPERATOR_BC_NEUMANN);
  auto cells = mesh_->getFaceCells(f);
  c = cells[0];
  tc1 = 0.;
  tc2 = 0.;

  Epetra_MultiVector& hap_gamma = *stencil_data_->ViewComponent("gamma", true);
  Epetra_MultiVector& weight = *stencil_data_->ViewComponent("weight", true);

  for (int i = 1; i < dim_; i++) {
    c3 = (*stencil_cells_[i])[f];
    f1 = (*stencil_faces_[i])[f];
    if (c3 >= 0) {
      auto cells_tmp = mesh_->getFaceCells(f1);
      OrderCellsByGlobalId_(cells_tmp, c1, c2);

      double gamma = hap_gamma[0][f1];
      if (c1 != c) gamma = 1.0 - gamma;
      tc1 += weight[i][f];

    } else if (bc_model[f1] == OPERATOR_BC_DIRICHLET) {
      tc1 += weight[i][f];
    } else if (bc_model[f1] == OPERATOR_BC_NEUMANN) {
      tc1 += weight[i][f];
    }
  }

  for (int i = 1; i < dim_; i++) {
    c3 = (*stencil_cells_[dim_ + i])[f];
    f1 = (*stencil_faces_[dim_ + i])[f];
    if (c3 >= 0) {
      auto cells_tmp = mesh_->getFaceCells(f1);
      OrderCellsByGlobalId_(cells_tmp, c1, c2);

      double gamma = hap_gamma[0][f1];
      if (c1 != c) gamma = 1.0 - gamma;
      tc2 += weight[i + dim_][f];

    } else if (bc_model[f1] == OPERATOR_BC_DIRICHLET) {
      tc2 += weight[i + dim_][f];
    } else if (bc_model[f1] == OPERATOR_BC_NEUMANN) {
      tc2 += weight[i + dim_][f];
    }
  }

  return 0;
}

} // namespace Operators
} // namespace Amanzi
