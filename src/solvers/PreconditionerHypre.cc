/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (coonet@ornl.gov)
*/

//! Hypre based preconditioners include Algebraic MultiGrid and global ILU
#include "HYPRE_parcsr_ls.h"
#include "Ifpack_Hypre.h"
#include "Teuchos_RCP.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "VerboseObject.hh"
#include "PreconditionerHypre.hh"
#include "Epetra_IntVector.h"

namespace Amanzi {
namespace AmanziSolvers {

void
PreconditionerHypre::copy_matrix_()
{
  Teuchos::RCP<const Epetra_CrsMatrix> Matrix =
    Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(A_);
  if (Matrix.is_null())
    throw std::runtime_error("Hypre<MatrixType>: Unsupported matrix "
                             "configuration: Tpetra::CrsMatrix required");

  std::vector<int> new_indices(Matrix->MaxNumEntries());
  for (int i = 0; i < Matrix->NumMyRows(); i++) {
    int numEntries;
    int* indices;
    double* values;
    Matrix->ExtractMyRowView(i, numEntries, values, indices);
    for (int j = 0; j < numEntries; j++) {
      new_indices[j] = GloballyContiguousColMap_->GID(indices[j]);
    }
    int GlobalRow[1];
    GlobalRow[0] = GloballyContiguousRowMap_->GID(i);
    HYPRE_IJMatrixSetValues(HypreA_, 1, &numEntries, GlobalRow, new_indices.data(), values);
  }
  HYPRE_IJMatrixAssemble(HypreA_);
  HYPRE_IJMatrixGetObject(HypreA_, (void**)&ParMatrix_);
}


/* ******************************************************************
 * Apply the preconditioner.
 ****************************************************************** */
int
PreconditionerHypre::ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const

{
  hypre_Vector* XLocal_ = hypre_ParVectorLocalVector(XVec_);
  hypre_Vector* YLocal_ = hypre_ParVectorLocalVector(YVec_);

  double* XValues = const_cast<double*>(&v[0]);
  double* YValues = const_cast<double*>(&hv[0]);

  double* XTemp = XLocal_->data;
  XLocal_->data = XValues;
  double* YTemp = YLocal_->data;
  YLocal_->data = YValues;

  HYPRE_ParVectorSetConstantValues(ParY_, 0.0);
  if (method_type_ == Boomer) {
    // Use the solver methods
    HYPRE_BoomerAMGSolve(method_, ParMatrix_, ParX_, ParY_);
  } else if (method_type_ == ILU) {
    // Apply the preconditioner
    HYPRE_ILUSolve(method_, ParMatrix_, ParX_, ParY_);
  } else if (method_type_ == MGR) {
    HYPRE_MGRSolve(method_, ParMatrix_, ParX_, ParY_);
  } else if (method_type_ == AMS) {
    HYPRE_AMSSolve(method_, ParMatrix_, ParX_, ParY_);
  }

  XLocal_->data = XTemp;
  YLocal_->data = YTemp;

  return 0;
}

/* ******************************************************************
 * Initialize the preconditioner.
 ****************************************************************** */
void
PreconditionerHypre::set_inverse_parameters(Teuchos::ParameterList& list)
{
  plist_ = list;

  std::string vo_name = this->name() + " (" + plist_.get<std::string>("method") + ")";
  vo_ = Teuchos::rcp(new VerboseObject(vo_name, plist_));
}

/* ******************************************************************
* Initialize AMG preconditioner.
****************************************************************** */
void
PreconditionerHypre::InitBoomer_()
{
#ifdef HAVE_HYPRE

  // check for old input spec and error
  if (plist_.isParameter("number of cycles")) {
    Errors::Message msg("\"boomer amg\" ParameterList uses old style, \"number "
                        "of cycles\".  Please update to the new style using "
                        "\"cycle applications\" and \"smoother sweeps\"");
    Exceptions::amanzi_throw(msg);
  }

  // verbosity
  int vlevel_int = 0;
  if (vo_->getVerbLevel() == Teuchos::VERB_EXTREME) {
    vlevel_int = 3;
  }

  if (plist_.isParameter("verbosity")) vlevel_int = plist_.get<int>("verbosity");
  HYPRE_BoomerAMGSetPrintLevel(method_, vlevel_int);

  HYPRE_BoomerAMGSetTol(method_, plist_.get<double>("tolerance", 0.0));
  HYPRE_BoomerAMGSetMaxIter(method_, plist_.get<int>("cycle applications", 5));
  HYPRE_BoomerAMGSetCoarsenType(method_, plist_.get<int>("coarsen type", 0));
  HYPRE_BoomerAMGSetStrongThreshold(method_, plist_.get<double>("strong threshold", 0.5));
  HYPRE_BoomerAMGSetCycleType(method_, plist_.get<int>("cycle type", 1));
  HYPRE_BoomerAMGSetNumSweeps(method_, plist_.get<int>("smoother sweeps", 3));

  if (plist_.isParameter("relaxation type down") && plist_.isParameter("relaxation type up")) {
    HYPRE_BoomerAMGSetCycleRelaxType(method_, plist_.get<int>("relaxation type down"), 1);
    HYPRE_BoomerAMGSetCycleRelaxType(method_, plist_.get<int>("relaxation type up"), 2);
  } else if (plist_.isParameter("relaxation type")) {
    HYPRE_BoomerAMGSetRelaxType(method_, plist_.get<int>("relaxation type"));
  } else {
    // use Hypre's defaults
  }

  if (plist_.isParameter("coarsening type"))
    HYPRE_BoomerAMGSetCoarsenType(method_, plist_.get<int>("coarsening type"));
  if (plist_.isParameter("interpolation type"))
    HYPRE_BoomerAMGSetInterpType(method_, plist_.get<int>("interpolation type"));
  if (plist_.isParameter("relaxation order"))
    HYPRE_BoomerAMGSetRelaxOrder(method_, plist_.get<int>("relaxation order"));
  if (plist_.isParameter("max multigrid levels"))
    HYPRE_BoomerAMGSetMaxLevels(method_, plist_.get<int>("max multigrid levels"));
  if (plist_.isParameter("max coarse size"))
    HYPRE_BoomerAMGSetMaxCoarseSize(method_, plist_.get<int>("max coarse size"));

  if (plist_.get<bool>("use block indices", false)) {
    num_blocks_ = plist_.get<int>("number of unique block indices");
    HYPRE_BoomerAMGSetNumFunctions(method_, num_blocks_);
    // Block indices is an array of integers, indicating what unknowns are
    // coarsened as a system.
    block_indices_ = plist_.get<Teuchos::RCP<std::vector<int>>>("block indices");
  }

  if (block_indices_.get()) {
    // must NEW the index array EVERY time, as it gets freed every time by
    // Hypre on cleanup.  This is pretty stupid, but we can't reuse it.  --etc
    //
    // Note this is not a memory leak -- this gets freed the second time
    // IfpHypre_::Compute() gets called (for every call but the last) and when
    // IfpHypre_ gets destroyed (for the last call).
    int* indices = new int[block_indices_->size()];
    for (int i = 0; i != block_indices_->size(); ++i) { indices[i] = (*block_indices_)[i]; }
    HYPRE_BoomerAMGSetDofFunc(method_, indices);
  }

  if (plist_.isParameter("number of functions")) {
    if (num_blocks_ > 0) {
      Errors::Message msg("Hypre (BoomerAMG) cannot be given both \"use block "
                          "indices\" and \"number of functions\" options as "
                          "these are two ways of specifying the same thing.");
      Exceptions::amanzi_throw(msg);
    }

    // num_funcs > 1 tells BoomerAMG to use the automatic "systems of
    // PDEs" code.  Note that, to use this approach, unknowns must be
    // ordered with DoF fastest varying (i.e. not the native
    // Epetra_MultiVector order).  By default, it uses the "unknown"
    // approach in which each equation is coarsened and interpolated
    // independently.  Comments below are taken from Allison Baker's
    // email to the PETSc mailing list, 25 Apr 2007, as these features
    // of BoomerAMG are not documented very well.  Here we ignore her
    // option 2, as she warns it is inefficient and likely not useful.
    // http://lists.mcs.anl.gov/pipermail/petsc-users/2007-April/001487.html

    int num_funcs = plist_.get<int>("number of functions");
    HYPRE_BoomerAMGSetNumFunctions(method_, num_funcs);

    // additional options
    if (num_funcs > 1) {
      // HYPRE_BOOMERAMGSetNodal(solver, int nodal ) tells AMG to coarsen such
      // that each variable has the same coarse grid - sometimes this is more
      // "physical" for a particular problem. The value chosen here for nodal
      // determines how strength of connection is determined between the
      // coupled system.  I suggest setting nodal = 1, which uses a Frobenius
      // norm.  This does NOT tell AMG to use nodal relaxation.
      if (plist_.isParameter("nodal strength of connection norm")) {
        int nodal = plist_.get<int>("nodal strength of connection norm", 0);
        HYPRE_BoomerAMGSetNodal(method_, nodal);
      }

      // You can additionally do nodal relaxation via the schwarz
      // smoother option in hypre.
      //   HYPRE_BoomerAMGSetSmoothType(solver, 6);
      //   HYPRE_BoomerAMGSetDomainType(solver, 1);
      //   HYPRE_BoomerAMGSetOverlap(solver, 0);
      //   HYPRE_BoomerAMGSetSmoothNumLevels(solver, num_levels);
      // Set num_levels to number of levels on which you want nodal
      // smoothing, i.e. 1=just the fine grid, 2= fine grid and the grid
      // below, etc.  I find that doing nodal relaxation on just the finest
      // level is generally sufficient.)  Note that the interpolation scheme
      // used will be the same as in the unknown approach - so this is what
      // we call a hybrid systems method.
      //
      if (plist_.isParameter("nodal relaxation levels")) {
        int num_levels = plist_.get<int>("nodal relaxation levels");

        // I believe this works, but needs testing -- we do not pop previous
        // settings, and instead just call the function twice. --ETC
        HYPRE_BoomerAMGSetSmoothType(method_, 6);
        HYPRE_BoomerAMGSetDomainType(method_, 1);
        HYPRE_BoomerAMGSetOverlap(method_, 0);
        HYPRE_BoomerAMGSetSmoothNumLevels(method_, num_levels);
        HYPRE_BoomerAMGSetSchwarzUseNonSymm(method_,
                                            1); // should provide an option for non-sym

        // Note that if num_levels > 1, you MUST also do nodal coarsening (to
        // maintain the nodes on coarser grids).
        if (num_levels > 1) {
          int nodal = plist_.get<int>("nodal strength of connection norm", 1);
          HYPRE_BoomerAMGSetNodal(method_, nodal);
        }
      }
    }
  }

#else
  Errors::Message msg("Hypre (BoomerAMG) is not available in this installation "
                      "of Amanzi.  To use Hypre, please reconfigure.");
  Exceptions::amanzi_throw(msg);
#endif
}


/* ******************************************************************
* Initialize modified ILU preconditioner.
****************************************************************** */
void
PreconditionerHypre::InitILU_()
{
#ifdef HAVE_HYPRE

  if (plist_.isParameter("verbosity"))
    HYPRE_ILUSetPrintLevel(method_, plist_.get<int>("verbosity"));

  if (plist_.isParameter("ilu(k) fill level"))
    HYPRE_ILUSetLevelOfFill(method_, plist_.get<int>("ilu(k) fill level"));
  HYPRE_ILUSetTol(method_, 0.0);

  // Default is 1 (suggested preconditioner value) for now FIXME
  // if (plist_.isParameter("maximum number of iterations"))
  //  HYPRE_ILUSetMaxIter(method_, plist_.get<int>("maximum number of iterations"));
  HYPRE_ILUSetMaxIter(method_, 1);

#else
  Errors::Message msg("Hypre (ILU) is not available in this installation of Amanzi.  To use "
                      "Hypre, please reconfigure.");
  Exceptions::amanzi_throw(msg);
#endif
}

/* ******************************************************************
* Initialize MRG preconditioner.
****************************************************************** */
void
PreconditionerHypre::InitMGR_()
{
#ifdef HAVE_HYPRE
  method_type_ = MGR;

  // Set up block data with information about coarse indexes for reduction.
  // Here, the user specifies the number of reduction levels, as well as the
  // coarse nodes for each level of the reduction. These coarse nodes are
  // indexed by their index in the block of unknowns. This is used internally
  // to tag the appropriate indexes of the linear system matrix as coarse nodes.
  // HYPRE_MGRSetCpointsByBlock(
  //   method_, block_size, max_num_levels, *num_block_coarse_points, **block_coarse_indexes);

  // Prescribe a subset of nodes to be kept as coarse nodes until the coarsest level.
  // These nodes are transferred onto the coarsest grid of the BoomerAMG coarse grid solver.
  // HYPRE_MGRSetReservedCoarseNodes(method_, reserved_coarse_size, *reserved_coarse_nodes);

  // Set points not prescribed as C points to be fixed as F points for intermediate levels.
  // Setting this to 1 uses the user input to define the C/F splitting. Otherwise, a
  // BoomerAMG coarsening routine is used to determine the C/F splitting for intermediate levels.
  // HYPRE_MGRSetNonCpointsToFpoints(method_, nonCptToFptFlag);

  // This function sets the BoomerAMG solver to be used for the solve on the coarse grid.
  // The user can define their own BoomerAMG solver with their preferred options and pass
  // this to the MGR solver. Otherwise, an internal BoomerAMG solver is used as the coarse
  // grid solver instead.
  // HYPRE_MGRSetCoarseSolver(
  //   method_, coarse_grid_solver_solve, coarse_grid_solver_setup, coarse_grid_solver);

#else
  Errors::Message msg("Hypre (MRG) is not available in this installation of Amanzi.  To use "
                      "Hypre, please reconfigure.");
  Exceptions::amanzi_throw(msg);
#endif
}


/* ******************************************************************
* Initialize specialized AMS preconditioner.
****************************************************************** */
void
PreconditionerHypre::InitAMS_()
{
#ifdef HAVE_HYPRE
  method_type_ = AMS;

  if (plist_.isParameter("verbosity"))
    HYPRE_AMSSetPrintLevel(method_, plist_.get<int>("verbosity"));

  // additional data: parallel vector with coordinates and gradient matrix
  // G_ * xyz_ should be the valid operation
  Teuchos::RCP<Epetra_MultiVector> xyz;
  Teuchos::RCP<Epetra_CrsMatrix> G;

  if (plist_.isParameter("discrete gradient operator") &&
      plist_.isType<Teuchos::RCP<Epetra_CrsMatrix>>("discrete gradient operator")) {
    G = plist_.get<Teuchos::RCP<Epetra_CrsMatrix>>("discrete gradient operator");
    SetDiscreteGradient_(G);
  }

  if (plist_.isParameter("graph coordinates") &&
      plist_.isType<Teuchos::RCP<Epetra_MultiVector>>("graph coordinates")) {
    xyz = plist_.get<Teuchos::RCP<Epetra_MultiVector>>("graph coordinates");
    SetCoordinates_(xyz);
  }

  if (!xyz.get() || !G.get()) {
    Errors::Message msg("Hypre (AMS) needs additional data: gradient matrix and graph coordinates");
    Exceptions::amanzi_throw(msg);
  }

  // PList must go first
  //Teuchos::ParameterList tmp;
  //tmp.set<std::string>("hypre: Preconditioner", "AMS")
  //  .set<std::string>("hypre: SolveOrPrecondition", "Preconditioner");
  //tmp.sublist("Coordinates").set<Teuchos::RCP<Epetra_MultiVector>>("Coordinates", xyz);
  //tmp.sublist("Operators").set<Teuchos::RCP<const Epetra_CrsMatrix>>("G", G);
  //IfpHypre_->SetParameters(tmp);

  HYPRE_AMSSetAlphaAMGOptions(method_, 6, 0, 6, 0.25, 0, 0);
  HYPRE_AMSSetBetaAMGOptions(method_, 6, 0, 6, 0.25, 0, 0);

  HYPRE_AMSSetTol(method_, plist_.get<double>("tolerance", 0.0));
  HYPRE_AMSSetPrintLevel(method_, plist_.get<int>("verbosity", 0));
  HYPRE_AMSSetMaxIter(method_, plist_.get<int>("cycle applications", 5));
  HYPRE_AMSSetCycleType(method_, plist_.get<int>("cycle type", 1));
  HYPRE_AMSSetAlphaAMGCoarseRelaxType(method_, plist_.get<int>("coarse level relaxation", 9));
  HYPRE_AMSSetBetaAMGCoarseRelaxType(method_, plist_.get<int>("coarse level relaxation", 9));

#else
  Errors::Message msg("Hypre (AMS) is not available in this installation of Amanzi. To use Hypre, "
                      "please reconfigure.");
  Exceptions::amanzi_throw(msg);
#endif
}

/* ******************************************************************
* Amanzi preconditioner initialization
****************************************************************** */
void
PreconditionerHypre::InitializeInverse()
{
  // must be row matrix
  A_ = h_;

  assert(A_->RowMatrixRowMap().SameAs(A_->OperatorRangeMap()));

  if (A_->RowMatrixRowMap().LinearMap()) {
    GloballyContiguousRowMap_ = Teuchos::rcpFromRef(A_->RowMatrixRowMap());
    GloballyContiguousColMap_ = Teuchos::rcpFromRef(A_->RowMatrixColMap());
  } else {
    assert(false);
  }

  MPI_Comm comm = dynamic_cast<const Epetra_MpiComm*>(&A_->Comm())->GetMpiComm();

  // Next create vectors that will be used when ApplyInverse() is called
  int ilower = GloballyContiguousRowMap_->MinMyGID();
  int iupper = GloballyContiguousRowMap_->MaxMyGID();
  // X in AX = Y
  HYPRE_IJVectorCreate(comm, ilower, iupper, &XHypre_);
  HYPRE_IJVectorSetObjectType(XHypre_, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(XHypre_);
  HYPRE_IJVectorAssemble(XHypre_);
  HYPRE_IJVectorGetObject(XHypre_, (void**)&ParX_);
  XVec_ = Teuchos::rcp((hypre_ParVector*)hypre_IJVectorObject(((hypre_IJVector*)XHypre_)), false);

  // Y in AX = Y
  HYPRE_IJVectorCreate(comm, ilower, iupper, &YHypre_);
  HYPRE_IJVectorSetObjectType(YHypre_, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(YHypre_);
  HYPRE_IJVectorAssemble(YHypre_);
  HYPRE_IJVectorGetObject(YHypre_, (void**)&ParY_);
  YVec_ = Teuchos::rcp((hypre_ParVector*)hypre_IJVectorObject(((hypre_IJVector*)YHypre_)), false);

  // A in AX = Y
  HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &HypreA_);

  std::string method_name = plist_.get<std::string>("method");
  if (method_name == "boomer amg") {
    method_type_ = Boomer;
    HYPRE_BoomerAMGCreate(&method_);
    InitBoomer_();
  } else if (method_name == "ILU") {
    method_type_ = ILU;
    HYPRE_ILUCreate(&method_);
    InitILU_();
  } else if (method_name == "ams") {
    method_type_ = AMS;
    HYPRE_AMSCreate(&method_);
    InitAMS_();
  } else if (method_name == "MGR") {
    method_type_ = MGR;
    HYPRE_MGRCreate(&method_);
    InitMGR_();
  } else {
    Errors::Message msg;
    msg << "PreconditionerHypre: unknown method name \"" << method_name << "\"";
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
 * Rebuild the preconditioner using the given matrix A.
 ****************************************************************** */
void
PreconditionerHypre::ComputeInverse()
{
#ifdef HAVE_HYPRE
  HYPRE_IJMatrixSetObjectType(HypreA_, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(HypreA_);
  copy_matrix_();
  if (method_type_ == Boomer) {
    // Hypre Setup must be called after matrix has values
    HYPRE_BoomerAMGSetup(method_, ParMatrix_, ParX_, ParY_);
  } else if (method_type_ == ILU) {
    HYPRE_ILUSetup(method_, ParMatrix_, ParX_, ParY_);
  } else if (method_type_ == MGR) {
    HYPRE_MGRSetup(method_, ParMatrix_, ParX_, ParY_);
  } else if (method_type_ == AMS) {
    HYPRE_AMSSetup(method_, ParMatrix_, ParX_, ParY_);
  }
#endif
}

int
PreconditionerHypre::SetCoordinates_(Teuchos::RCP<Epetra_MultiVector> coords)
{
  if (!G_.is_null() && !G_->DomainMap().SameAs(coords->Map()))
    throw std::runtime_error("Ifpack_Hypre: Node map mismatch: G->DomainMap() and coords");

  if (method_type_ != AMS) return 0;

  double* xPtr;
  double* yPtr;
  double* zPtr;

  IFPACK_CHK_ERR(((*coords)(0))->ExtractView(&xPtr));
  IFPACK_CHK_ERR(((*coords)(1))->ExtractView(&yPtr));
  IFPACK_CHK_ERR(((*coords)(2))->ExtractView(&zPtr));

  MPI_Comm comm = dynamic_cast<const Epetra_MpiComm*>(&A_->Comm())->GetMpiComm();

  int NumEntries = coords->MyLength();
  int* indices = GloballyContiguousNodeRowMap_->MyGlobalElements();

  int ilower = GloballyContiguousNodeRowMap_->MinMyGID();
  int iupper = GloballyContiguousNodeRowMap_->MaxMyGID();

  assert(NumEntries == iupper - ilower + 1);

  HYPRE_IJVectorCreate(comm, ilower, iupper, &xHypre_);
  HYPRE_IJVectorSetObjectType(xHypre_, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(xHypre_);
  HYPRE_IJVectorSetValues(xHypre_, NumEntries, indices, xPtr);
  HYPRE_IJVectorAssemble(xHypre_);
  HYPRE_IJVectorGetObject(xHypre_, (void**)&xPar_);

  HYPRE_IJVectorCreate(comm, ilower, iupper, &yHypre_);
  HYPRE_IJVectorSetObjectType(yHypre_, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(yHypre_);
  HYPRE_IJVectorSetValues(yHypre_, NumEntries, indices, yPtr);
  HYPRE_IJVectorAssemble(yHypre_);
  HYPRE_IJVectorGetObject(yHypre_, (void**)&yPar_);

  HYPRE_IJVectorCreate(comm, ilower, iupper, &zHypre_);
  HYPRE_IJVectorSetObjectType(zHypre_, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(zHypre_);
  HYPRE_IJVectorSetValues(zHypre_, NumEntries, indices, zPtr);
  HYPRE_IJVectorAssemble(zHypre_);
  HYPRE_IJVectorGetObject(zHypre_, (void**)&zPar_);

  HYPRE_AMSSetCoordinateVectors(method_, xPar_, yPar_, zPar_);

  return 0;

} //SetCoordinate

int
PreconditionerHypre::SetDiscreteGradient_(Teuchos::RCP<const Epetra_CrsMatrix> G)
{
  // Sanity check
  if (!A_->RowMatrixRowMap().SameAs(G->RowMap()))
    throw std::runtime_error("Ifpack_Hypre: Edge map mismatch: A and discrete gradient");

  // Get the maps for the nodes (assuming the edge map from A is OK);
  GloballyContiguousNodeRowMap_ = Teuchos::rcp(new Epetra_Map(
    G->DomainMap().NumGlobalElements(), G->DomainMap().NumMyElements(), 0, A_->Comm()));
  Teuchos::RCP<const Epetra_RowMatrix> Grow = Teuchos::rcp_dynamic_cast<const Epetra_RowMatrix>(G);
  GloballyContiguousNodeColMap_ = MakeContiguousColumnMap_(Grow);

  // Start building G
  MPI_Comm comm = dynamic_cast<const Epetra_MpiComm*>(&A_->Comm())->GetMpiComm();
  int ilower = GloballyContiguousRowMap_->MinMyGID();
  int iupper = GloballyContiguousRowMap_->MaxMyGID();
  int jlower = GloballyContiguousNodeRowMap_->MinMyGID();
  int jupper = GloballyContiguousNodeRowMap_->MaxMyGID();
  HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, &HypreG_);
  HYPRE_IJMatrixSetObjectType(HypreG_, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(HypreG_);

  std::vector<int> new_indices(G->MaxNumEntries());
  for (int i = 0; i < G->NumMyRows(); i++) {
    int numEntries;
    double* values;
    int* indices;
    G->ExtractMyRowView(i, numEntries, values, indices);
    for (int j = 0; j < numEntries; j++) {
      new_indices[j] = GloballyContiguousNodeColMap_->GID(indices[j]);
    }
    int GlobalRow[1];
    GlobalRow[0] = GloballyContiguousRowMap_->GID(i);
    HYPRE_IJMatrixSetValues(HypreG_, 1, &numEntries, GlobalRow, new_indices.data(), values);
  }
  HYPRE_IJMatrixAssemble(HypreG_);
  HYPRE_IJMatrixGetObject(HypreG_, (void**)&ParMatrixG_);

  HYPRE_AMSSetDiscreteGradient(method_, ParMatrixG_);

  return 0;
} //SetDiscreteGradient()

Teuchos::RCP<const Epetra_Map>
PreconditionerHypre::MakeContiguousColumnMap_(Teuchos::RCP<const Epetra_RowMatrix>& MatrixRow) const
{
  // Must create GloballyContiguous DomainMap (which is a permutation of Matrix_'s
  // DomainMap) and the corresponding permuted ColumnMap.
  //   Epetra_GID  --------->   LID   ----------> HYPRE_GID
  //           via DomainMap.LID()       via GloballyContiguousDomainMap.GID()
  Teuchos::RCP<const Epetra_CrsMatrix> Matrix =
    Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(MatrixRow);
  if (Matrix.is_null())
    throw std::runtime_error(
      "Ifpack_Hypre: Unsupported matrix configuration: Epetra_CrsMatrix required");
  const Epetra_Map& DomainMap = Matrix->DomainMap();
  const Epetra_Map& ColumnMap = Matrix->ColMap();
  const Epetra_Import* importer = Matrix->Importer();

  if (DomainMap.LinearMap()) {
    // If the domain map is linear, then we can just use the column map as is.
    return Teuchos::rcpFromRef(ColumnMap);
  } else {
    // The domain map isn't linear, so we need a new domain map
    Teuchos::RCP<Epetra_Map> ContiguousDomainMap = Teuchos::rcp(
      new Epetra_Map(DomainMap.NumGlobalElements(), DomainMap.NumMyElements(), 0, A_->Comm()));
    if (importer) {
      // If there's an importer then we can use it to get a new column map
      Epetra_IntVector MyGIDsHYPRE(View, DomainMap, ContiguousDomainMap->MyGlobalElements());

      // import the HYPRE GIDs
      Epetra_IntVector ColGIDsHYPRE(ColumnMap);
      ColGIDsHYPRE.Import(MyGIDsHYPRE, *importer, Insert);

      // Make a HYPRE numbering-based column map.
      return Teuchos::rcp(new Epetra_Map(
        ColumnMap.NumGlobalElements(), ColGIDsHYPRE.MyLength(), &ColGIDsHYPRE[0], 0, A_->Comm()));
    } else {
      // The problem has matching domain/column maps, and somehow the domain map isn't linear, so just use the new domain map
      return Teuchos::rcp(new Epetra_Map(ColumnMap.NumGlobalElements(),
                                         ColumnMap.NumMyElements(),
                                         ContiguousDomainMap->MyGlobalElements(),
                                         0,
                                         A_->Comm()));
    }
  }
}

} // namespace AmanziSolvers
} // namespace Amanzi
