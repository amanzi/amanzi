/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (coonet@ornl.gov)
*/

//! Hypre based preconditioners include Algebraic MultiGrid and global ILU

#include "Teuchos_RCP.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "VerboseObject.hh"
#include "PreconditionerHypre.hh"

#include "cuda_decl.h"

namespace Amanzi {

namespace AmanziSolvers {

bool PreconditionerHypre::inited = false; 

void PreconditionerHypre::copy_matrix_(){
  nvtxRangePush("HP: copy");

  Teuchos::RCP<const Matrix_type> Matrix = Teuchos::rcp_dynamic_cast<const Matrix_type>(h_row);
  if(Matrix.is_null()) 
    throw std::runtime_error("Hypre<MatrixType>: Unsupported matrix configuration: Tpetra::CrsMatrix required");

  LO nrows = Matrix->getLocalNumRows();

  #if 1
  Kokkos::View<LO*,Kokkos::DefaultExecutionSpace> colsperrow("ColsPerRow",nrows);
  auto rowPtrs = Matrix->getCrsGraph()->getLocalGraphDevice().row_map;
  Kokkos::parallel_for(nrows,
  KOKKOS_LAMBDA(const int i){
    colsperrow(i) = rowPtrs[i+1]-rowPtrs[i];
  });

  auto rowindices = Matrix->getRowMap()->getMyGlobalIndices();

  Kokkos::View<GO *, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace> rid("",rowindices.size());
  Kokkos::deep_copy(rid,rowindices);

  auto values = Matrix->getLocalValuesDevice(Tpetra::Access::ReadOnly); 
  auto colindices = Matrix->getCrsGraph()->getLocalGraphDevice().entries;
  decltype(colindices) nci("",colindices.size()); 
  auto ge = GloballyContiguousColMap_->getMyGlobalIndices(); 

  Kokkos::View<GO*,Kokkos::DefaultExecutionSpace> dge("",ge.size());
  Kokkos::deep_copy(dge,ge);  
  Kokkos::parallel_for(
    "", 
    colindices.size(), 
    KOKKOS_LAMBDA(const int i){
      nci[i] = dge(colindices[i]);
     }
  );
  #else 
  Kokkos::View<LO*,Kokkos::HostSpace> colsperrow("ColsPerRow",nrows);
  auto rowPtrs = Matrix->getCrsGraph()->getLocalGraphHost().row_map;
  for(int i = 0 ; i < nrows; ++i){
    colsperrow(i) = rowPtrs[i+1]-rowPtrs[i];
  }
  auto rowindices = Matrix->getRowMap()->getMyGlobalIndices();

  Kokkos::View<GO *, Kokkos::LayoutLeft, Kokkos::HostSpace> rid("",rowindices.size());
  Kokkos::deep_copy(rid,rowindices);

  auto values = Matrix->getLocalValuesHost(Tpetra::Access::ReadOnly); 
  auto colindices = Matrix->getCrsGraph()->getLocalGraphHost().entries;
  decltype(colindices) nci("",colindices.size()); 
  auto ge = GloballyContiguousColMap_->getMyGlobalIndices(); 

  Kokkos::View<GO*,Kokkos::HostSpace> dge("",ge.size());
  Kokkos::deep_copy(dge,ge);  
  for(int i = 0 ; i < colindices.size() ; ++i){
      nci[i] = dge(colindices[i]);
  }
  #endif 

  HYPRE_IJMatrixSetValues(HypreA_,nrows,colsperrow.data(),
                                         rid.data(),nci.data(),values.data());
  HYPRE_IJMatrixAssemble(HypreA_);
  HYPRE_IJMatrixGetObject(HypreA_, (void**)&ParMatrix_);
  nvtxRangePop();
}

/* ******************************************************************
* Apply the preconditioner.
****************************************************************** */
int PreconditionerHypre::applyInverse(const Vector_type& v, Vector_type& hv) const
{
  nvtxRangePush("HP: apply");
  assert(v.getNumVectors() == 1);
  assert(hv.getNumVectors() == 1);
  assert(&v != &hv);

  hypre_Vector *XLocal_ = hypre_ParVectorLocalVector(XVec_);
  hypre_Vector *YLocal_ = hypre_ParVectorLocalVector(YVec_);

  double_type * XValues = const_cast<double_type*>(v.getLocalViewDevice(Tpetra::Access::ReadOnly).data());
  double_type * YValues = const_cast<double_type*>(hv.getLocalViewDevice(Tpetra::Access::ReadWrite).data());

  double_type *XTemp = XLocal_->data;
  XLocal_->data = XValues;
  double_type *YTemp = YLocal_->data;
  YLocal_->data = YValues;
  
  HYPRE_ParVectorSetConstantValues(ParY_, 0.0); 
  if(PrecondType == Boomer){
    HYPRE_BoomerAMGSolve(HyprePrecond_,ParMatrix_ ,ParX_, ParY_);
  } else if (PrecondType == Euclid){
    HYPRE_EuclidSolve(HyprePrecond_,ParMatrix_ ,ParX_, ParY_);
  }

  XLocal_->data = XTemp;
  YLocal_->data = YTemp;
  nvtxRangePop();
  return 0;
}


/* ******************************************************************
* Initialize the preconditioner.
****************************************************************** */
void PreconditionerHypre::set_inverse_parameters(Teuchos::ParameterList& list)
{
  plist_ = list;

  std::string vo_name = this->name()+" ("+plist_.get<std::string>("method")+")";
  vo_ = Teuchos::rcp(new VerboseObject(vo_name, plist_));
}

void PreconditionerHypre::InitBoomer_()
{
  nvtxRangePush("HP: initBoomer");
#ifdef HAVE_IFPACK2_HYPRE

  // check for old input spec and error
  if (plist_.isParameter("number of cycles")) {
    Errors::Message msg("\"boomer amg\" ParameterList uses old style, \"number of cycles\".  Please update to the new style using \"cycle applications\" and \"smoother sweeps\"");
    Exceptions::amanzi_throw(msg);
  }

  // verbosity
  int vlevel_int = 1;
  std::string vlevel = plist_.sublist("verbose object").get<std::string>("verbosity level", "medium");
  if (vlevel == "high") {
    vlevel_int = 2;
  } else if (vlevel == "extreme") {
    vlevel_int = 3;
  }

  if (plist_.isParameter("verbosity")) vlevel_int = plist_.get<int>("verbosity");
  HYPRE_BoomerAMGSetPrintLevel(HyprePrecond_,vlevel_int);
  HYPRE_BoomerAMGSetTol(HyprePrecond_, plist_.get<double>("tolerance", 0.0));
  HYPRE_BoomerAMGSetMaxIter(HyprePrecond_, plist_.get<int>("cycle applications", 5));
  HYPRE_BoomerAMGSetCoarsenType(HyprePrecond_, plist_.get<int>("coarsen type", 0));
  HYPRE_BoomerAMGSetStrongThreshold(HyprePrecond_, plist_.get<double>("strong threshold", 0.5));
  HYPRE_BoomerAMGSetCycleType(HyprePrecond_, plist_.get<int>("cycle type", 1));
  HYPRE_BoomerAMGSetNumSweeps(HyprePrecond_, plist_.get<int>("smoother sweeps", 3));

  if (plist_.isParameter("relaxation type down") && plist_.isParameter("relaxation type up")) {
    HYPRE_BoomerAMGSetCycleRelaxType(HyprePrecond_, plist_.get<int>("relaxation type down"), 1);
    HYPRE_BoomerAMGSetCycleRelaxType(HyprePrecond_, plist_.get<int>("relaxation type up"), 2);
  } else if (plist_.isParameter("relaxation type")) {
    HYPRE_BoomerAMGSetRelaxType(HyprePrecond_, plist_.get<int>("relaxation type")); 
  } else { 
    // use Hypre's defaults
  }

  if (plist_.isParameter("coarsening type")) 
    HYPRE_BoomerAMGSetCoarsenType(HyprePrecond_, plist_.get<int>("coarsening type"));
  if (plist_.isParameter("interpolation type")) 
    HYPRE_BoomerAMGSetInterpType(HyprePrecond_, plist_.get<int>("interpolation type"));
  if (plist_.isParameter("relaxation order")) 
   HYPRE_BoomerAMGSetRelaxOrder(HyprePrecond_, plist_.get<int>("relaxation order"));
  if (plist_.isParameter("max multigrid levels"))
    HYPRE_BoomerAMGSetMaxLevels(HyprePrecond_, plist_.get<int>("max multigrid levels"));
  if (plist_.isParameter("max coarse size"))
    HYPRE_BoomerAMGSetMaxCoarseSize(HyprePrecond_, plist_.get<int>("max coarse size"));


  if (block_indices_.get()) {
    // must NEW the index array EVERY time, as it gets freed every time by
    // Hypre on cleanup.  This is pretty stupid, but we can't reuse it.  --etc
    //
    // Note this is not a memory leak -- this gets freed the second time
    // IfpHypre_::Compute() gets called (for every call but the last) and when
    // IfpHypre_ gets destroyed (for the last call).
    int* indices = new int[block_indices_->size()];
    for (int i=0; i!=block_indices_->size(); ++i) {
      indices[i] = (*block_indices_)[i];
    }
    HYPRE_BoomerAMGSetDofFunc(HyprePrecond_, indices);
  }

  if (plist_.get<bool>("use block indices", false)) {
    num_blocks_ = plist_.get<int>("number of unique block indices");
    HYPRE_BoomerAMGSetNumFunctions(HyprePrecond_, num_blocks_);
    // Block indices is an array of integers, indicating what unknowns are
    // coarsened as a system.
    block_indices_ = plist_.get<Teuchos::RCP<std::vector<int> > >("block indices");
  }

  if (plist_.isParameter("number of functions")) {
    if (num_blocks_ > 0) {
      Errors::Message msg("Hypre (BoomerAMG) cannot be given both \"use block indices\" and \"number of functions\" options as these are two ways of specifying the same thing.");
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
    HYPRE_BoomerAMGSetNumFunctions(HyprePrecond_, num_funcs);

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
        HYPRE_BoomerAMGSetNodal(HyprePrecond_, nodal);
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
        HYPRE_BoomerAMGSetSmoothType(HyprePrecond_, 6);
        HYPRE_BoomerAMGSetDomainType(HyprePrecond_, 1);
        HYPRE_BoomerAMGSetOverlap(HyprePrecond_, 0);
        HYPRE_BoomerAMGSetSmoothNumLevels(HyprePrecond_, num_levels);
        HYPRE_BoomerAMGSetSchwarzUseNonSymm(HyprePrecond_, 1); // should provide an option for non-sym

        // Note that if num_levels > 1, you MUST also do nodal coarsening (to maintain the nodes on coarser grids).
        if (num_levels > 1) {
          int nodal = plist_.get<int>("nodal strength of connection norm", 1);
          HYPRE_BoomerAMGSetNodal(HyprePrecond_, nodal);
        }
      }
    }
  }

#else
  Errors::Message msg("Hypre (BoomerAMG) is not available in this installation of Amanzi.  To use Hypre, please reconfigure.");
  Exceptions::amanzi_throw(msg);
#endif
  nvtxRangePop();
}


/* ******************************************************************
 * Initialize the preconditioner.
 ****************************************************************** */
void PreconditionerHypre::InitEuclid_()
{
  Init_();

#ifdef HAVE_IFPACK2_HYPRE

  HYPRE_EuclidSetStats(HyprePrecond_,plist_.get<int>("verbosity"));

  HYPRE_EuclidSetLevel(HyprePrecond_, plist_.get<int>("ilu(k) fill level"));

  if (plist_.isParameter("rescale rows")) {
    bool rescale_rows = plist_.get<bool>("rescale rows");
    HYPRE_EuclidSetRowScale(HyprePrecond_,rescale_rows ? 1 : 0);
  }

  if (plist_.isParameter("ilut drop tolerance"))
    HYPRE_EuclidSetILUT(HyprePrecond_, plist_.get<double>("ilut drop tolerance"));
#else
  Errors::Message msg("Hypre (Euclid) is not available in this installation of Amanzi.  To use Hypre, please reconfigure.");
  Exceptions::amanzi_throw(msg);
#endif
}


void PreconditionerHypre::initializeInverse()
{
  nvtxRangePush("HP: initInverse");
  int count =0; 

  // must be row matrix
  h_row = h_;
  //IfpHypre_ = Teuchos::rcp(new Ifpack2::Hypre<RowMatrix_type>(h_row));

  if (h_row->getRowMap()->isContiguous()) {
    GloballyContiguousRowMap_ = h_row->getRowMap();
    GloballyContiguousColMap_ = h_row->getColMap();
  } else {  
    assert(false); 
  }

  MPI_Comm comm = * (Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(h_row->getRowMap()->getComm())->getRawMpiComm());

  // Next create vectors that will be used when ApplyInverse() is called
  GO ilower = GloballyContiguousRowMap_->getMinGlobalIndex();
  GO iupper = GloballyContiguousRowMap_->getMaxGlobalIndex();
  // X in AX = Y
  HYPRE_IJVectorCreate(comm, ilower, iupper, &XHypre_);
  HYPRE_IJVectorSetObjectType(XHypre_, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(XHypre_);
  HYPRE_IJVectorAssemble(XHypre_);
  HYPRE_IJVectorGetObject(XHypre_, (void**) &ParX_);
  XVec_ = Teuchos::rcp((hypre_ParVector *) hypre_IJVectorObject(((hypre_IJVector *) XHypre_)),false);

  // Y in AX = Y
  HYPRE_IJVectorCreate(comm, ilower, iupper, &YHypre_);
  HYPRE_IJVectorSetObjectType(YHypre_, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(YHypre_);
  HYPRE_IJVectorAssemble(YHypre_);
  HYPRE_IJVectorGetObject(YHypre_, (void**) &ParY_);
  YVec_ = Teuchos::rcp((hypre_ParVector *) hypre_IJVectorObject(((hypre_IJVector *) YHypre_)),false);

  std::string method_name = plist_.get<std::string>("method");
  if (method_name == "boomer amg") {
    PrecondType = Boomer; 
  } else if (method_name == "euclid") {
    PrecondType = Euclid; 
  } else {
    Errors::Message msg;
    msg << "PreconditionerHypre: unknown method name \"" << method_name << "\"";
    Exceptions::amanzi_throw(msg);
  }

}

Teuchos::RCP<const Map_type> 
PreconditionerHypre::make_contiguous_(Teuchos::RCP<const RowMatrix_type> &Matrix){
  using import_type     = Tpetra::Import<LO,GO,NT>;
  using go_vector_type  = Tpetra::Vector<GO,LO,GO,NT>;

  // Must create GloballyContiguous DomainMap (which is a permutation of Matrix_'s
  // DomainMap) and the corresponding permuted ColumnMap.
  //   Epetra_GID  --------->   LID   ----------> HYPRE_GID
  //           via DomainMap.LID()       via GloballyContiguousDomainMap.GID()
  if(Matrix.is_null()) 
    throw std::runtime_error("Hypre<MatrixType>: Unsupported matrix configuration: Tpetra::CrsMatrix required");
  Teuchos::RCP<const Map_type> DomainMap = Matrix->getDomainMap();
  Teuchos::RCP<const Map_type> ColumnMap = Matrix->getColMap();
  //Teuchos::RCP<const import_type> importer = Matrix->getGraph()->getImporter();

  if(DomainMap->isContiguous() ) {
    // If the domain map is linear, then we can just use the column map as is.
    return ColumnMap;
  }
  else {
    assert(false); 
    return ColumnMap; 
  }  
  nvtxRangePop();
}

/* ******************************************************************
* Rebuild the preconditioner using the given matrix A.
****************************************************************** */
void PreconditionerHypre::computeInverse()
{
  nvtxRangePush("HP: compute");
#ifdef HAVE_IFPACK2_HYPRE
  MPI_Comm comm = * (Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(h_row->getRowMap()->getComm())->getRawMpiComm());
  GO ilower = GloballyContiguousRowMap_->getMinGlobalIndex();
  GO iupper = GloballyContiguousRowMap_->getMaxGlobalIndex();
  HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &HypreA_);
  HYPRE_IJMatrixSetObjectType(HypreA_, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(HypreA_);
  copy_matrix_();
  if(PrecondType == Boomer){
    HYPRE_BoomerAMGCreate(&HyprePrecond_);
    // can only set parameter when the preconditioner is created  
    InitBoomer_();
    // Hypre Setup must be called after matrix has values
    HYPRE_BoomerAMGSetup(HyprePrecond_, ParMatrix_, ParX_, ParY_);
  } else if (PrecondType == Euclid){
    HYPRE_EuclidCreate(comm,&HyprePrecond_); 
    InitEuclid_(); 
    HYPRE_EuclidSetup(HyprePrecond_, ParMatrix_, ParX_, ParY_); 
  }
#endif
  nvtxRangePop();
}

}  // namespace AmanziSolvers
}  // namespace Amanzi
