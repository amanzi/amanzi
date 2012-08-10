#include <RichardSolver.H>
#include <RICHARDSOLVER_F.H>

#include <Utility.H>


#undef BUILD_DENSE_J
#define BUILD_DENSE_J 1
#undef PETSC_3_2
#define PETSC_3_2 1

static bool write_matrix = false;

// FIXME: Should be user-settable at runtime
static int pressure_maxorder = 4;
//static Real errfd = 1.e-10;
static Real errfd = 1.e-8;

#if defined(BUILD_DENSE_J)
static int max_number_nonzero_J_cols = 2048;
#else
static int max_number_nonzero_J_cols = 1 + (pressure_maxorder-1)*(2*BL_SPACEDIM); // 3.3 seems to need this
#endif

static int richard_max_ls_iterations = 10;
static Real richard_min_ls_factor = 1.e-8;
static Real richard_ls_acceptance_factor = 2;
static Real richard_ls_reduction_factor = 0.1;
static int richard_monitor_line_search = 1;

#if defined(PETSC_3_2)
PetscErrorCode FormLineSearch(SNES snes,void* user,Vec X,Vec F,Vec G,Vec Y,Vec W,PetscReal fnorm,PetscReal xnorm,
                              PetscReal *ynorm,PetscReal *gnorm,PetscBool  *flag)
{
  PetscErrorCode ierr;
  /*
      Input parameters
	snes 	- nonlinear context
	lsctx 	- optional user-defined context for line search
	x 	- current iterate
	f 	- residual evaluated at x
	y 	- search direction
	fnorm 	- 2-norm of f

      Output parameters
	g 	- residual evaluated at new iterate y
	w 	- new iterate
	gnorm 	- 2-norm of g
	ynorm 	- 2-norm of search length
	flag 	- set to PETSC_TRUE if the line search succeeds; PETSC_FALSE on failure. 
   */
  PetscScalar mone=-1.0;
  *flag=PETSC_TRUE;
  ierr=VecNorm(Y,NORM_2,ynorm);
  ierr=VecWAXPY(W,mone,Y,X); /* W = -Y + X */
  ierr=SNESComputeFunction(snes,W,G);CHKPETSC(ierr);
  ierr=VecNorm(G,NORM_2,gnorm);CHKPETSC(ierr);

#else
    PetscErrorCode FormLineSearch(SNESLineSearch linesearch, void * user)
    {
        
        Vec  X, Y, F, W, G;
        
        SNES snes;
        PetscErrorCode ierr;
        PetscReal fnorm, xnorm, *ynorm, *gnorm;
        PetscBool *flag;
        PetscScalar mone=-1.0;

        PetscFunctionBegin;
        
        ierr = SNESLineSearchGetSNES(linesearch, &snes);CHKPETSC(ierr);
        
        ierr = SNESLineSearchSetSuccess(linesearch, PETSC_TRUE);CHKPETSC(ierr);

        // X (old), F = f(X), Y = dX, W = X + Y new solution, G = f(W)
        ierr = SNESLineSearchGetVecs(linesearch, &X, &Y, &F, &W, &G);CHKPETSC(ierr);

        PetscErrorCode (**func)(SNES,Vec,Vec,void*);
        void *fctx;
        ierr = SNESGetFunction(snes,PETSC_NULL,func,&fctx);
        std::cout << "********************************** doing f1" << std::endl;
        (*func)(snes,X,F,fctx);
        std::cout << "********************************** doing f2" << std::endl;
        (*func)(snes,W,G,fctx);
        std::cout << "********************************** doing f2a" << std::endl;

        ierr = SNESLineSearchComputeNorms(linesearch);CHKPETSC(ierr);
        ierr = SNESLineSearchGetNorms(linesearch,&xnorm,&fnorm,ynorm);CHKPETSC(ierr);
#endif
        ierr=VecNorm(G,NORM_2,gnorm);CHKPETSC(ierr);

  RichardSolver* rs = static_cast<RichardSolver*>(user);

  std::cout << "********************************** in my ls" << std::endl;

  std::string tag = "       Newton step: ";
  std::string tag_ls = "  line-search:  ";

  bool norm_acceptable = *gnorm < fnorm * richard_ls_acceptance_factor;
  int ls_iterations = 0;
  Real ls_factor = 1;
  bool finished = norm_acceptable 
    || ls_iterations > richard_max_ls_iterations
    || ls_factor <= richard_min_ls_factor;

  while (!finished) 
    {
      ls_factor *= richard_ls_reduction_factor;
      if (ls_factor < richard_min_ls_factor) {
          ls_factor = richard_min_ls_factor;
      }
      ierr=VecWAXPY(W,mone*ls_factor,Y,X); /* W = -Y + X */

      PetscErrorCode (**func)(SNES,Vec,Vec,void*);
      void *fctx;
      ierr = SNESGetFunction(snes,PETSC_NULL,func,&fctx);
      std::cout << "********************************** doing f3" << std::endl;
      (*func)(snes,W,G,fctx);
      std::cout << "********************************** doing f4" << std::endl;
      ierr=VecNorm(G,NORM_2,gnorm);CHKPETSC(ierr);
      norm_acceptable = *gnorm < fnorm * richard_ls_acceptance_factor;

      if (ls_factor < 1 
	  && richard_monitor_line_search 
	  && ParallelDescriptor::IOProcessor())
	{
          std::cout << tag << tag_ls
                    << "iter=" << ls_iterations
                    << ", step length=" << ls_factor
                    << ", init residual norm=" << fnorm
                    << ", new residual norm=" << *gnorm << '\n';
	}

      finished = norm_acceptable 
	|| ls_iterations > richard_max_ls_iterations
	|| ls_factor <= richard_min_ls_factor;      
      ls_iterations++;
    }

  if (ls_iterations > richard_max_ls_iterations) 
    {
      std::string reason = "Solution rejected.  Linear system solved, but ls_iterations too large";
      *flag = PETSC_FALSE;
      if (ParallelDescriptor::IOProcessor() && richard_monitor_line_search) {
	std::cout << tag << tag_ls << reason << std::endl;
      }
    }
  else if (ls_factor <= richard_min_ls_factor) {
    std::string reason = "Solution rejected.  Linear system solved, but ls_factor too small";
    *flag = PETSC_FALSE;
    if (ParallelDescriptor::IOProcessor() && richard_monitor_line_search) {
      std::cout << tag << tag_ls << reason << std::endl;
    }
  }
  else {
    *flag = PETSC_TRUE;
    if (ls_factor == 1) {
        std::string reason = "Full linear step accepted";
        if (ParallelDescriptor::IOProcessor() && richard_monitor_line_search) {
            std::cout << tag << tag_ls << reason << std::endl;
        }
    }
    else {
        // Set update to the one actually used
        ierr=VecScale(Y,ls_factor);
    }
  }
  return ierr;
}

PetscErrorCode RichardRes_DpDt(SNES snes,Vec x,Vec f,void *dummy)
{
    PetscErrorCode ierr; 
    RichardSolver* rs = static_cast<RichardSolver*>(dummy);
    if (!rs) {
        BoxLib::Abort("Bad cast in RichardRes_DpDt");
    }
    Layout& layout = PMAmr::GetLayout();
    MFTower& xMFT = rs->GetPressure();
    MFTower& fMFT = rs->GetResidual();

    ierr = layout.VecToMFTower(xMFT,x,0); CHKPETSC(ierr);
    ierr = layout.VecToMFTower(fMFT,f,0); CHKPETSC(ierr);

    Real time = rs->Time();
    Real dt = rs->Dt();
    rs->DpDtResidual(fMFT,xMFT,time,dt);

    ierr = layout.MFTowerToVec(x,xMFT,0); CHKPETSC(ierr);
    ierr = layout.MFTowerToVec(f,fMFT,0); CHKPETSC(ierr);

    return 0;
}

PetscErrorCode RichardRes_DpDt_Alt(SNES snes,Vec x,Vec f,void *dummy)
{
    PetscErrorCode ierr; 
    RichardSolver* rs = static_cast<RichardSolver*>(dummy);
    if (!rs) {
        BoxLib::Abort("Bad cast in RichardRes_DpDt");
    }
    Layout& layout = PMAmr::GetLayout();
    MFTower& xMFT = rs->GetPressure();
    MFTower& fMFT = rs->GetResidual();

    ierr = layout.VecToMFTower(xMFT,x,0); CHKPETSC(ierr);
    ierr = layout.VecToMFTower(fMFT,f,0); CHKPETSC(ierr);

    Real time = rs->Time();
    Real dt = rs->Dt();
    rs->DpDtResidual_Alt(fMFT,xMFT,time,dt);

    ierr = layout.MFTowerToVec(x,xMFT,0); CHKPETSC(ierr);
    ierr = layout.MFTowerToVec(f,fMFT,0); CHKPETSC(ierr);

    return 0;
}

static RichardSolver* static_rs_ptr = 0;

void
RichardSolver::SetTheRichardSolver(RichardSolver* ptr) 
{
    static_rs_ptr = ptr;
}
#if defined(PETSC_3_2)
#include <private/matimpl.h>
#else
#include <petsc-private/matimpl.h> 
#endif
PetscErrorCode  RichardMatFDColoringApply(Mat J,MatFDColoring coloring,Vec x1,MatStructure *flag,void *sctx)
{
  PetscErrorCode (*f)(void*,Vec,Vec,void*) = (PetscErrorCode (*)(void*,Vec,Vec,void *))coloring->f;
  PetscErrorCode ierr;
  PetscInt       k,start,end,l,row,col,srow,**vscaleforrow,m1,m2;
  PetscScalar    dx,*y,*xx,*w3_array;
  PetscScalar    *vscale_array, *sigma_array;
  PetscReal      epsilon = coloring->error_rel,umin = coloring->umin,unorm; 
  Vec            w1=coloring->w1,w2=coloring->w2,w3;
  void           *fctx = coloring->fctx;
  PetscBool      flg = PETSC_FALSE;
  PetscInt       ctype=coloring->ctype,N,col_start=0,col_end=0;
  Vec            x1_tmp;

  PetscFunctionBegin;    
  PetscValidHeaderSpecific(J,MAT_CLASSID,1);
  PetscValidHeaderSpecific(coloring,MAT_FDCOLORING_CLASSID,2);
  PetscValidHeaderSpecific(x1,VEC_CLASSID,3);
  if (!f) SETERRQ(((PetscObject)J)->comm,PETSC_ERR_ARG_WRONGSTATE,"Must call MatFDColoringSetFunction()");

  ierr = PetscLogEventBegin(MAT_FDColoringApply,coloring,J,x1,0);CHKPETSC(ierr);
  ierr = MatSetUnfactored(J);CHKPETSC(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-mat_fd_coloring_dont_rezero",&flg,PETSC_NULL);CHKPETSC(ierr);
  if (flg) {
    ierr = PetscInfo(coloring,"Not calling MatZeroEntries()\n");CHKPETSC(ierr);
  } else {
    PetscBool  assembled;
    ierr = MatAssembled(J,&assembled);CHKPETSC(ierr);
    if (assembled) {
      ierr = MatZeroEntries(J);CHKPETSC(ierr);
    }
  }

  x1_tmp = x1; 
  if (!coloring->vscale){ 
    ierr = VecDuplicate(x1_tmp,&coloring->vscale);CHKPETSC(ierr);
  }
    
  /*
    This is a horrible, horrible, hack. See DMMGComputeJacobian_Multigrid() it inproperly sets
    coloring->F for the coarser grids from the finest
  */
  if (coloring->F) {
    ierr = VecGetLocalSize(coloring->F,&m1);CHKPETSC(ierr);
    ierr = VecGetLocalSize(w1,&m2);CHKPETSC(ierr);
    if (m1 != m2) {  
      coloring->F = 0; 
      }    
    }   


  RichardSolver* rs = static_rs_ptr;
  BL_ASSERT(rs);
  
  Layout& layout = PMAmr::GetLayout();
  int nLevs = layout.NumLevels();
  PArray<MultiFab> pcap_params(nLevs,PArrayNoManage);
  for (int lev=0; lev<nLevs; ++lev) {
      pcap_params.set(lev, rs->GetPMlevel(lev).PCapParams());
  }
  MFTower PCapParamsMFT(layout,pcap_params);
  Vec& SigmaV = rs->GetPCapParamsV();
  ierr = layout.MFTowerToVec(SigmaV,PCapParamsMFT,2); CHKPETSC(ierr);
  ierr = VecGetOwnershipRange(w1,&start,&end);CHKPETSC(ierr); /* OwnershipRange is used by ghosted x! */
      
  /* Set w1 = F(x1) */
  if (coloring->F) {
    w1          = coloring->F; /* use already computed value of function */
    coloring->F = 0; 
  } else {
    ierr = PetscLogEventBegin(MAT_FDColoringFunction,0,0,0,0);CHKPETSC(ierr);
    ierr = (*f)(sctx,x1_tmp,w1,fctx);CHKPETSC(ierr);
    ierr = PetscLogEventEnd(MAT_FDColoringFunction,0,0,0,0);CHKPETSC(ierr);
  }
      
  if (!coloring->w3) {
    ierr = VecDuplicate(x1_tmp,&coloring->w3);CHKPETSC(ierr);
    ierr = PetscLogObjectParent(coloring,coloring->w3);CHKPETSC(ierr);
  }
  w3 = coloring->w3;

    /* Compute all the local scale factors, including ghost points */
  ierr = VecGetLocalSize(x1_tmp,&N);CHKPETSC(ierr);
  ierr = VecGetArray(x1_tmp,&xx);CHKPETSC(ierr);
  ierr = VecGetArray(SigmaV,&sigma_array);CHKPETSC(ierr);
  ierr = VecGetArray(coloring->vscale,&vscale_array);CHKPETSC(ierr);
  if (ctype == IS_COLORING_GHOSTED){
    col_start = 0; col_end = N;
  } else if (ctype == IS_COLORING_GLOBAL){
    xx = xx - start;
    sigma_array = sigma_array - start;
    vscale_array = vscale_array - start;
    col_start = start; col_end = N + start;
  }
  for (col=col_start; col<col_end; col++) { 
      vscale_array[col] = (PetscScalar)sigma_array[col] / epsilon;
  } 
  if (ctype == IS_COLORING_GLOBAL)  vscale_array = vscale_array + start;      
  ierr = VecRestoreArray(coloring->vscale,&vscale_array);CHKPETSC(ierr);
  ierr = VecRestoreArray(SigmaV,&sigma_array);CHKPETSC(ierr);

  if (ctype == IS_COLORING_GLOBAL){
      ierr = VecGhostUpdateBegin(coloring->vscale,INSERT_VALUES,SCATTER_FORWARD);CHKPETSC(ierr);
      ierr = VecGhostUpdateEnd(coloring->vscale,INSERT_VALUES,SCATTER_FORWARD);CHKPETSC(ierr);
  }
  
  if (coloring->vscaleforrow) {
    vscaleforrow = coloring->vscaleforrow;
  } else SETERRQ(((PetscObject)J)->comm,PETSC_ERR_ARG_NULL,"Null Object: coloring->vscaleforrow");

  /*
    Loop over each color
  */
  ierr = VecGetArray(coloring->vscale,&vscale_array);CHKPETSC(ierr);
  for (k=0; k<coloring->ncolors; k++) { 
    coloring->currentcolor = k;
    ierr = VecCopy(x1_tmp,w3);CHKPETSC(ierr);
    ierr = VecGetArray(w3,&w3_array);CHKPETSC(ierr);
    if (ctype == IS_COLORING_GLOBAL) w3_array = w3_array - start;
    /*
      Loop over each column associated with color 
      adding the perturbation to the vector w3.
    */
    for (l=0; l<coloring->ncolumns[k]; l++) {
      col = coloring->columns[k][l];    /* local column of the matrix we are probing for */
      w3_array[col] += 1/vscale_array[col];
    } 
    if (ctype == IS_COLORING_GLOBAL) w3_array = w3_array + start;
    ierr = VecRestoreArray(w3,&w3_array);CHKPETSC(ierr);

    /*
      Evaluate function at w3 = x1 + dx (here dx is a vector of perturbations)
                           w2 = F(x1 + dx) - F(x1)
    */
    ierr = PetscLogEventBegin(MAT_FDColoringFunction,0,0,0,0);CHKPETSC(ierr);
    ierr = (*f)(sctx,w3,w2,fctx);CHKPETSC(ierr);        
    ierr = PetscLogEventEnd(MAT_FDColoringFunction,0,0,0,0);CHKPETSC(ierr);
    ierr = VecAXPY(w2,-1.0,w1);CHKPETSC(ierr); 
        
    /*
      Loop over rows of vector, putting results into Jacobian matrix
    */
    ierr = VecGetArray(w2,&y);CHKPETSC(ierr);
    for (l=0; l<coloring->nrows[k]; l++) {
      row    = coloring->rows[k][l];             /* local row index */
      col    = coloring->columnsforrow[k][l];    /* global column index */
      y[row] *= vscale_array[vscaleforrow[k][l]];
      srow   = row + start;
      ierr   = MatSetValues(J,1,&srow,1,&col,y+row,INSERT_VALUES);CHKPETSC(ierr);
    }
    ierr = VecRestoreArray(w2,&y);CHKPETSC(ierr);
  } /* endof for each color */
  if (ctype == IS_COLORING_GLOBAL) xx = xx + start; 
  ierr = VecRestoreArray(coloring->vscale,&vscale_array);CHKPETSC(ierr);
  ierr = VecRestoreArray(x1_tmp,&xx);CHKPETSC(ierr);
   
  coloring->currentcolor = -1;
  ierr  = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKPETSC(ierr);
  ierr  = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKPETSC(ierr);
  ierr = PetscLogEventEnd(MAT_FDColoringApply,coloring,J,x1,0);CHKPETSC(ierr);

  flg  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-mat_null_space_test",&flg,PETSC_NULL);CHKPETSC(ierr);
  if (flg) {
    ierr = MatNullSpaceTest(J->nullsp,J,PETSC_NULL);CHKPETSC(ierr);
  }

  if (write_matrix) {
      PetscViewer viewer;
#if defined(BUILD_DENSE_J)
      std::cout << "Writing/checking dense matrix" << std::endl;

      std::string fname = "junk_dense.bin";
      
      std::string oname = "junk_sparse.bin";
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,oname.c_str(),FILE_MODE_READ,&viewer);
      Mat A;
      MatCreate(PETSC_COMM_WORLD,&A);
      MatSetType(A,MATSEQAIJ);
      MatLoad(A,viewer);
      PetscViewerDestroy(&viewer);

      int rstart, rend;
      ierr = MatGetOwnershipRange(J,&rstart,&rend);CHKPETSC(ierr);
      int nrows = 0;
      int Jncols, Ancols;
      const PetscInt *Jcols, *Acols;
      const PetscScalar *Jvals, *Avals;
      PetscReal dtol = 1.e-20;
      for (int row=rstart; row<rend; row++){
          ierr = MatGetRow(J,row,&Jncols,&Jcols,&Jvals);CHKPETSC(ierr);
          ierr = MatGetRow(A,row,&Ancols,&Acols,&Avals);CHKPETSC(ierr);
          for (int j=0; j<Jncols; j++){
              PetscScalar Jval = Jvals[j];
              int Aptr = -1;
              for (int jj=0; jj<Ancols; ++jj) {
                  if (Acols[jj] == Jcols[j]) {
                      Aptr = jj;
                  }
              }
              if (Aptr<0 && std::abs(Jval)>dtol) {
                  std::cout << "Row: " << row << " missing col: " << Jcols[j] << std::endl;
              }
              else {
                  PetscScalar Aval = Avals[Aptr];
                  if (std::abs(Jval-Aval) >= dtol) {
                      std::cout << "Row: " << row << " missing col: " << Jcols[j] 
                                << " vals: " << Jval  << " " << Aval << std::endl;
                  }
                  else {
                      std::cout << "Row: " << row << " col: " << Jcols[j] << " diff: " 
                                << Jval - Aval << " mag: " << std::abs(Jval) << std::endl;
                  }
              }
          }
          ierr = MatRestoreRow(J,row,&Jncols,&Jcols,&Jvals);CHKPETSC(ierr);
          ierr = MatRestoreRow(A,row,&Ancols,&Acols,&Avals);CHKPETSC(ierr);
      }
#else
      std::cout << "Writing sparse matrix" << std::endl;
      std::string fname = "junk_sparse.bin";
#endif
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,fname.c_str(),FILE_MODE_WRITE,&viewer);
      MatView(J,viewer);
      PetscViewerDestroy(&viewer);

      BoxLib::Abort();
  }
  PetscFunctionReturn(0);
}

PetscErrorCode RichardComputeJacobianColor(SNES snes,Vec x1,Mat *J,Mat *B,MatStructure *flag,void *ctx)
{
  MatFDColoring  color = (MatFDColoring) ctx;
  PetscErrorCode ierr;
  Vec            f;
  PetscErrorCode (*ff)(void),(*fd)(void);

  PetscFunctionBegin;
  PetscValidHeaderSpecific(color,MAT_FDCOLORING_CLASSID,6);
  *flag = SAME_NONZERO_PATTERN;
  ierr  = SNESGetFunction(snes,&f,(PetscErrorCode (**)(SNES,Vec,Vec,void*))&ff,0);CHKPETSC(ierr);
  ierr  = MatFDColoringGetFunction(color,&fd,PETSC_NULL);CHKPETSC(ierr);
  if (fd == ff) { /* reuse function value computed in SNES */
    ierr  = MatFDColoringSetF(color,f);CHKPETSC(ierr);
  }
  ierr  = RichardMatFDColoringApply(*B,color,x1,flag,snes);CHKPETSC(ierr);
  if (*J != *B) {
    ierr = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);CHKPETSC(ierr);
    ierr = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);CHKPETSC(ierr);
  }
  PetscFunctionReturn(0);
}

RichardSolver::RichardSolver(PMAmr& _pm_amr)
  : pm_amr(_pm_amr)
{
  mftfp = new MFTFillPatch(PMAmr::GetLayout());
  nLevs = pm_amr.finestLevel() + 1;
  pm.resize(nLevs,PArrayNoManage);
  for (int lev = 0; lev < nLevs; lev++)  {
    pm.set(lev,dynamic_cast<PorousMedia*>(&pm_amr.getLevel(lev)));
  }
      
  Layout& layout = PMAmr::GetLayout();
  BL_ASSERT(nLevs == layout.NumLevels());

  // Manipulate data already allocated in state into structures more convenient for ml solves
  RhoSatOld = 0;
  RhoSatNew = 0;

  // ResetRhoSat();
  PArray<MultiFab> lambda(nLevs,PArrayNoManage);
  PArray<MultiFab> porosity(nLevs,PArrayNoManage);
  PArray<MultiFab> P_new(nLevs,PArrayNoManage);
  PArray<MultiFab> utmp(nLevs,PArrayNoManage);
  
  for (int lev=0; lev<nLevs; ++lev) {
    lambda.set(lev,pm[lev].LambdaCC_Curr());
    porosity.set(lev,pm[lev].Porosity());
    P_new.set(lev,&(pm[lev].get_new_data(Press_Type)));
  }

  Lambda = new MFTower(layout,lambda);
  Porosity = new MFTower(layout,porosity);
  Pnew = new MFTower(layout,P_new);
      
  // Because kpedge has a weird boxarray, we maintain our own copy....should fix this
  ktmp.resize(BL_SPACEDIM);
  ctmp.resize(BL_SPACEDIM);
  Rhs = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1);
  
  Kappa.resize(BL_SPACEDIM,PArrayManage);
  RichardCoefs.resize(BL_SPACEDIM,PArrayManage);
  DarcyVelocity.resize(BL_SPACEDIM,PArrayNoManage);
  for (int d=0; d<BL_SPACEDIM; ++d) {
    ktmp[d].resize(nLevs,PArrayManage);
    ctmp[d].resize(nLevs,PArrayManage);
    for (int lev=0; lev<nLevs; ++lev) {
      BoxArray ba = BoxArray(lambda[lev].boxArray()).surroundingNodes(d);
      ktmp[d].set(lev, new MultiFab(ba,1,0)); ktmp[d][lev].copy(pm[lev].KappaEC()[d]);
      ctmp[d].set(lev, new MultiFab(ba,1,0));
      utmp.set(lev,&(pm[lev].UMac_Curr()[d]));
    }
    Kappa.set(d,new MFTower(layout,ktmp[d]));
    DarcyVelocity.set(d, new MFTower(layout,utmp));
    RichardCoefs.set(d, new MFTower(layout,ctmp[d]));
    utmp.clear();
  }

  // Make some temporary space
  RhoSatStar = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1);
  AlphaStar = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1);
  PStar = new MFTower(layout,IndexType(IntVect::TheZeroVector()),1,1);

  PetscErrorCode ierr;       
  int n = layout.NumberOfLocalNodeIds();
  int N = layout.NumberOfGlobalNodeIds();
  MPI_Comm comm = ParallelDescriptor::Communicator();
  ierr = VecCreateMPI(comm,n,N,&RhsV); CHKPETSC(ierr);
  ierr = VecDuplicate(RhsV,&SolnV); CHKPETSC(ierr);
  ierr = VecDuplicate(RhsV,&PCapParamsV); CHKPETSC(ierr);

  const BCRec& pressure_bc = pm[0].get_desc_lst()[Press_Type].getBC(0);
  mftfp->BuildStencil(pressure_bc, pressure_maxorder);

  gravity.resize(BL_SPACEDIM,0);
  gravity[BL_SPACEDIM-1] = PorousMedia::getGravity();
  density = PorousMedia::Density();
  
  int d_nz = 1 + 2*BL_SPACEDIM; // Estmated number of nonzero local columns of J
  int o_nz = 0; // Estimated number of nonzero nonlocal (off-diagonal) columns of J
  Mat Jac;
#if defined(PETSC_3_2)
  ierr = MatCreateMPIAIJ(comm, n, n, N, N, d_nz, PETSC_NULL, o_nz, PETSC_NULL, &Jac); CHKPETSC(ierr);
#else
  ierr = MatCreate(comm, &Jac); CHKPETSC(ierr);
  ierr = MatSetSizes(Jac,n,n,N,N);  CHKPETSC(ierr);
  ierr = MatSetFromOptions(Jac); CHKPETSC(ierr);
  d_nz = max_number_nonzero_J_cols;
  ierr = MatSeqAIJSetPreallocation(Jac, d_nz*d_nz, PETSC_NULL); CHKPETSC(ierr);
  ierr = MatMPIAIJSetPreallocation(Jac, d_nz, PETSC_NULL, o_nz, PETSC_NULL); CHKPETSC(ierr);
#endif  
  BuildOpSkel(Jac);
  
  matfdcoloring = 0;
  ierr = SNESCreate(comm,&snes); CHKPETSC(ierr);
  ierr = SNESSetFunction(snes,RhsV,RichardRes_DpDt,(void*)(this)); CHKPETSC(ierr);
  ierr = MatGetColoring(Jac,MATCOLORINGSL,&iscoloring); CHKPETSC(ierr);
  ierr = MatFDColoringCreate(Jac,iscoloring,&matfdcoloring); CHKPETSC(ierr);
  ierr = MatFDColoringSetFunction(matfdcoloring,
				  (PetscErrorCode (*)(void))RichardRes_DpDt,
				  (void*)(this)); CHKPETSC(ierr);
  ierr = SNESSetJacobian(snes,Jac,Jac,RichardComputeJacobianColor,matfdcoloring);CHKPETSC(ierr); 
  ierr = MatFDColoringSetParameters(matfdcoloring,errfd,PETSC_DEFAULT);CHKPETSC(ierr);
#if defined(PETSC_3_2)
  ierr = SNESLineSearchSet(snes,FormLineSearch,(void *)(this));CHKPETSC(ierr);
#else
  ierr = SNESGetSNESLineSearch(snes, &ls);CHKPETSC(ierr);


  Vec X, Y, F, W, G;
  ierr = VecDuplicate(RhsV,&X); CHKPETSC(ierr);
  ierr = VecDuplicate(RhsV,&Y); CHKPETSC(ierr);
  ierr = VecDuplicate(RhsV,&F); CHKPETSC(ierr);
  ierr = VecDuplicate(RhsV,&W); CHKPETSC(ierr);
  ierr = VecDuplicate(RhsV,&G); CHKPETSC(ierr);

  ierr = VecSet(G,2);
  PetscReal *gnorm;
  ierr=VecNorm(G,NORM_2,gnorm);CHKPETSC(ierr);

  ierr = SNESLineSearchSetVecs(ls,X,F,Y,W,G);
  ierr = SNESLineSearchSetType(ls, SNESLINESEARCHSHELL);CHKPETSC(ierr);
  ierr = SNESLineSearchShellSetUserFunc(ls,FormLineSearch,(void*)(this));CHKPETSC(ierr);
#endif

  // TODO: Add deletes/destroys for all structures allocated/created here

}

void
RichardSolver::ResetRhoSat()
{
  PArray<MultiFab> S_new(nLevs,PArrayNoManage);
  PArray<MultiFab> S_old(nLevs,PArrayNoManage);
  
  for (int lev=0; lev<nLevs; ++lev) {
    S_new.set(lev,&(pm[lev].get_new_data(State_Type)));
    S_old.set(lev,&(pm[lev].get_old_data(State_Type)));
  }

  Layout& layout = PMAmr::GetLayout();
  delete RhoSatOld; RhoSatOld = new MFTower(layout,S_old);
  delete RhoSatNew; RhoSatNew = new MFTower(layout,S_new);

  for (int lev=0; lev<nLevs; ++lev) {
      MultiFab::Copy((*PStar)[lev],pm[lev].get_new_data(Press_Type),0,0,1,0);
      MultiFab::Copy((*RhoSatStar)[lev],(*RhoSatNew)[lev],0,0,1,0);
      pm[lev].calc_richard_alpha(&((*AlphaStar)[lev]),(*RhoSatStar)[lev]);
  }
}

void
RichardSolver::SetTime(Real time) 
{
  mytime = time;
}

Real
RichardSolver::Time() const 
{
  return mytime;
}

void 
RichardSolver::SetDt(Real dt) 
{
  mydt = dt;
}

Real
RichardSolver::Dt() const 
{
  return mydt;
}

void
RichardSolver::BuildOpSkel(Mat& J)
{
  Layout& layout = PMAmr::GetLayout();
  int num_rows = 1;
  int rows[1]; // At the moment, only set one row at a time
  Array<Real> vals;
  Array<int> cols;
  
  const Array<Geometry>& geomArray = layout.GeomArray();
  const Array<BoxArray>& gridArray = layout.GridArray();
  const Array<IntVect>& refRatio = layout.RefRatio();
  const PArray<Layout::MultiNodeFab>& nodes = layout.Nodes();
  const PArray<Layout::MultiIntFab>& nodeIds = layout.NodeIds();
  const Array<BoxArray>& bndryCells = layout.BndryCells();
  const Array<Array<IVSMap> >& growCellStencil = mftfp->GrowCellStencil();
  int nLevs = layout.NumLevels();
  
  PetscErrorCode ierr;
  int num_nbrs_reg = 2*BL_SPACEDIM+1;
  Layout::IntFab reg_neighbors;
  std::set<int> neighbors;
  typedef BaseFab<std::set<int> > ISetFab;
  typedef FabArray<ISetFab> MultiSetFab;
  PArray<MultiSetFab> crseContribs(nLevs,PArrayManage);
  
  int myproc = ParallelDescriptor::MyProc();
  int numprocs = ParallelDescriptor::NProcs();
  
  for (int lev=nLevs-1; lev>=0; --lev) 
    {
      const Array<IVSMap>& growCellStencilLev = growCellStencil[lev];
      const Layout::MultiNodeFab& nodeLev = nodes[lev];
      const Layout::MultiIntFab& nodeIdsLev = nodeIds[lev];

      Layout::MultiIntFab crseIds; // coarse cell ids surrounding fine grid, distributed to align with fine patches
      crseContribs.set(lev,new MultiSetFab);
      if (lev>0) {
	BoxArray bacg = BoxArray(gridArray[lev]).coarsen(refRatio[lev-1]).grow(1);
	crseIds.define(bacg,1,0,Fab_allocate);
            
	const Layout::MultiIntFab& crseIds_orig = nodeIds[lev-1]; // Need the following to crse cells through periodic boundary
	BoxArray gcba = BoxArray(crseIds_orig.boxArray()).grow(crseIds_orig.nGrow());
	Layout::MultiIntFab tmp(gcba,1,0);
	for (MFIter mfi(crseIds_orig); mfi.isValid(); ++mfi) {
	  tmp[mfi].copy(crseIds_orig[mfi]); // NOTE: Assumes grow cells already filled
	}
	crseIds.copy(tmp); // Parallel copy

	crseContribs[lev].define(bacg,1,0,Fab_allocate);
      }

      std::map<IntVect,std::set<int>,IntVect::Compare> stencil;
      if (lev<nLevs-1) {
	// Pack up the crseContribs for a parallel copy
	const BoxArray& ba = gridArray[lev];
	MultiSetFab& crseContribsFine = crseContribs[lev+1];
	const DistributionMapping& dm = crseContribsFine.DistributionMap();
	std::map<int,Array<int> > ccArrays;
	for (MFIter mfi(crseContribsFine); mfi.isValid(); ++mfi) {
	  const ISetFab& ccFab = crseContribsFine[mfi];
	  const Box& vbox = mfi.validbox();
	  std::vector< std::pair<int,Box> > isects = ba.intersections(vbox);
	  for (int i=0; i<isects.size(); ++i) {
	    int dst_proc = dm[isects[i].first];

            // HACK  This was originally written for parallel, but when I tried it in serial, the entire 
            // crseContribs structure was ignored!!  For now, set this up as a communication, even if 
            // serial...probably an easy logic issue to clear up....famous last words...
	    if (1 || dst_proc != myproc) {
	      for (IntVect iv(vbox.smallEnd()), iEnd=vbox.bigEnd(); iv<=iEnd; vbox.next(iv))
		{
		  const std::set<int>& ids = ccFab(iv,0);
		  int thisSize = ids.size();
		  if (thisSize) {
		    Array<int>& ints = ccArrays[dst_proc];
		    int old_cc_size = ints.size();
		    int delta_cc = BL_SPACEDIM + 1 + ids.size();
		    int new_cc_size = old_cc_size + delta_cc;

		    ints.resize(new_cc_size);
		    for (int d=0; d<BL_SPACEDIM; ++d) {
		      ints[old_cc_size+d] = iv[d];
		    }
		    ints[old_cc_size+BL_SPACEDIM] = ids.size();
		    int cnt=0;
		    for (std::set<int>::const_iterator it=ids.begin(), End=ids.end(); it!=End; ++it, ++cnt) {
		      ints[old_cc_size+BL_SPACEDIM+1+cnt] = *it;
		    }
		  }
		}
	    }
	  }
	}

	int total_num_to_send = 0;
	Array<int> sends(numprocs,0);
	Array<int> soffsets(numprocs,0);
	for (int i=0; i<numprocs; ++i) {
	  sends[i] = ccArrays[i].size();
	  total_num_to_send += sends[i];
	  if (i>0) {
	    soffsets[i] = soffsets[i-1] + ccArrays[i-1].size();
	  }
	}
	Array<int> sbuf(total_num_to_send);
	for (int i=0; i<numprocs; ++i) {
	  for (int j=0; j<ccArrays[i].size(); ++j) {
	    sbuf[soffsets[i] + j] = ccArrays[i][j];
	  }
	}

	Array<int> recvs(numprocs);
	BL_MPI_REQUIRE( MPI_Alltoall(sends.dataPtr(),
				     1,
				     ParallelDescriptor::Mpi_typemap<int>::type(),
				     recvs.dataPtr(),
				     1,
				     ParallelDescriptor::Mpi_typemap<int>::type(),
				     ParallelDescriptor::Communicator()) );
            
	int total_num_to_recv = 0;
	Array<int> roffsets(numprocs,0);
	for (int i=0; i<numprocs; ++i) {
	  total_num_to_recv += recvs[i];
	  if (i>0) {
	    roffsets[i] = roffsets[i-1] + recvs[i-1];
	  }
	}
	Array<int> rbuf(total_num_to_recv);
	BL_MPI_REQUIRE( MPI_Alltoallv(total_num_to_send == 0 ? 0 : sbuf.dataPtr(),
				      sends.dataPtr(),
				      soffsets.dataPtr(),
				      ParallelDescriptor::Mpi_typemap<int>::type(),
				      total_num_to_recv == 0 ? 0 : rbuf.dataPtr(),
				      recvs.dataPtr(),
				      roffsets.dataPtr(),
				      ParallelDescriptor::Mpi_typemap<int>::type(),
				      ParallelDescriptor::Communicator()) );
            
	for (int i=0; i<numprocs; ++i) {
	  int jcnt = roffsets[i];
	  while (jcnt < roffsets[i] + recvs[i]) {
	    IntVect iv(&(rbuf[jcnt]));
	    int size = rbuf[jcnt+BL_SPACEDIM];
	    std::set<int>& iset = stencil[iv];
	    for (int k=0; k<size; ++k) {
	      iset.insert(rbuf[jcnt+BL_SPACEDIM+1+k]);
	    }
	    jcnt += BL_SPACEDIM+1+size;
	  }
	}
      }

      for (MFIter mfi(nodeLev); mfi.isValid(); ++mfi) {
	const Layout::NodeFab& nodeFab = nodeLev[mfi];
	const Layout::IntFab& nodeIdFab = nodeIdsLev[mfi];
	const Layout::IntFab* crseIdFab = (lev>0  ?  &(crseIds[mfi])  : 0);
	const Box& vbox = mfi.validbox();
	Box gbox = Box(vbox).grow(1);

	for (IntVect iv(vbox.smallEnd()), iEnd=vbox.bigEnd(); iv<=iEnd; vbox.next(iv))
	  {
	    const Node& nC = nodeFab(iv,0);
	    if (nC.type==Node::VALID) {
	      rows[0] = nodeIdFab(iv,0);
	      neighbors.clear();

	      std::map<IntVect,std::set<int>,IntVect::Compare>::const_iterator sit=stencil.find(iv);
	      if (sit!=stencil.end()) {
		const std::set<int>& iset = sit->second;
		neighbors.insert(iset.begin(),iset.end());
	      }
	      neighbors.insert(rows[0]);

	      for (int d=0; d<BL_SPACEDIM; ++d) {
		for (int pm = -1; pm<2; pm+=2) {
		  std::set<int> nd;
		  IntVect ivA = iv  +  pm * BoxLib::BASISV(d);
		  IVScit it=growCellStencilLev[d].find(ivA);
		  if (it!=growCellStencilLev[d].end()) {
		    const Stencil& s = it->second;
		    for (Stencil::const_iterator it=s.begin(), End=s.end(); it!=End; ++it) {
		      const Node& node = it->first;
		      const IntVect& ivs = node.iv;
		      int slev = node.level;
		      if (slev==lev) {
			BL_ASSERT(nodeIdFab.box().contains(ivs));
			int idx = nodeIdFab(ivs,0);
			if (ivs != iv && idx>=0) { // idx<0 is where we hold the Dirichlet data, and iv was already added above
			  nd.insert(idx);
			}
		      }
		      else if (slev==lev-1) {
			BL_ASSERT(crseIdFab);
			BL_ASSERT(crseIdFab->box().contains(ivs));
			nd.insert((*crseIdFab)(ivs,0));
		      }
		      else {
			std::cout << "stencil: " << s << std::endl;
			BoxLib::Abort("Bad stencil");
		      }
		    }

		    // contribute to coarse cell stencil, if appropriate
		    const Node& offcenter_node = nodeFab(ivA,0);
		    if (offcenter_node.type==Node::VALID  &&  offcenter_node.level==lev-1) {
		      crseContribs[lev][mfi](offcenter_node.iv,0).insert(rows[0]);
		      crseContribs[lev][mfi](offcenter_node.iv,0).insert(nd.begin(),nd.end());
		    }
		  }
		  else {
		    int idx = nodeIdFab(ivA,0);
		    if (idx>=0) { // idx<0 is a covered cell
		      neighbors.insert(idx);
		    }
		  }

		  // Merge this arm into full set
		  neighbors.insert(nd.begin(),nd.end());

		}
	      }

#if defined(BUILD_DENSE_J)
	      int num_cols = layout.NumberOfGlobalNodeIds();
	      cols.resize(num_cols);
	      vals.resize(num_cols,0);
	      int cnt = 0;
              for (int i=0; i<num_cols; ++i) cols[i] = i;
#else
	      int num_cols = neighbors.size();
	      cols.resize(num_cols);
	      vals.resize(num_cols,0);
	      int cnt = 0;
	      for (std::set<int>::const_iterator it=neighbors.begin(), End=neighbors.end(); it!=End; ++it) {
		cols[cnt++] = *it;
	      }
#endif
	      ierr = MatSetValues(J,num_rows,rows,num_cols,cols.dataPtr(),vals.dataPtr(),INSERT_VALUES); CHKPETSC(ierr);
	    }
	  }
      }
    }

  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKPETSC(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKPETSC(ierr);
}

void
RichardSolver::CenterToEdgeUpwind(PArray<MFTower>&       mfte,
				  MFTower&               mftc,
				  const PArray<MFTower>& sgn,
				  int                    nComp,
                                  const BCRec&           bc) const
{
    for (int lev=0; lev<nLevs; ++lev) {
        MultiFab& clev = mftc[lev];
        BL_ASSERT(nComp<=clev.nComp());
        const Box& domain = GeomArray()[lev].Domain();
        for (MFIter mfi(clev); mfi.isValid(); ++mfi) {
            FArrayBox& cfab = clev[mfi];
            const Box& vccbox = mfi.validbox();
            for (int d=0; d<BL_SPACEDIM; ++d) {            
                FArrayBox& efab = mfte[d][lev][mfi];
                const FArrayBox& sgnfab = sgn[d][lev][mfi];
                BL_ASSERT(nComp<=efab.nComp());
                BL_ASSERT(nComp<=sgnfab.nComp());
                BL_ASSERT(Box(vccbox).surroundingNodes(d).contains(efab.box()));
                BL_ASSERT(Box(vccbox).surroundingNodes(d).contains(sgnfab.box()));
                
                int dir_bc_lo = bc.lo()[d]==EXT_DIR && vccbox.smallEnd()[d]==domain.smallEnd()[d];
                int dir_bc_hi = bc.hi()[d]==EXT_DIR && vccbox.bigEnd()[d]==domain.bigEnd()[d];
                
                FORT_RS_CTE_UPW(efab.dataPtr(), ARLIM(efab.loVect()), ARLIM(efab.hiVect()),
				cfab.dataPtr(), ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
				sgnfab.dataPtr(),ARLIM(sgnfab.loVect()), ARLIM(sgnfab.hiVect()),
				vccbox.loVect(), vccbox.hiVect(), &d, &nComp, &dir_bc_lo, &dir_bc_hi);
            }
        }
    }
}

void 
RichardSolver::XmultYZ(MFTower&       X,
		       const MFTower& Y,
		       const MFTower& Z,
		       int            sCompY,
		       int            sCompZ,
		       int            dComp,
		       int            nComp,
		       int            nGrow)
{
  Layout& layout = PMAmr::GetLayout(); 
  BL_ASSERT(layout.IsCompatible(X));
  BL_ASSERT(layout.IsCompatible(Y));
  BL_ASSERT(layout.IsCompatible(Z));
  BL_ASSERT(X.NComp()>=dComp+nComp);
  BL_ASSERT(Y.NComp()>=sCompY+nComp);
  BL_ASSERT(Z.NComp()>=sCompZ+nComp);
  BL_ASSERT(X.NGrow()>=nGrow);
  BL_ASSERT(Y.NGrow()>=nGrow);
  BL_ASSERT(Z.NGrow()>=nGrow);
  const Array<BoxArray>& gridArray = GridArray();
  const Array<Geometry>& geomArray = GeomArray();
  const Array<IntVect>& refRatio = RefRatio();

  FArrayBox tfabY, tfabZ;
  for (int lev=0; lev<nLevs; ++lev)
    {
      MultiFab& Xlev = X[lev];
      const MultiFab& Ylev = Y[lev];
      const MultiFab& Zlev = Z[lev];
      BoxArray fba;
      if (lev<nLevs-1) {
	fba = BoxArray(X[lev+1].boxArray()).coarsen(refRatio[lev]);
      }

      for (MFIter mfi(Y[lev]); mfi.isValid(); ++mfi)
        {
	  const Box& vbox = mfi.validbox();
	  Box gbox = Box(vbox).grow(nGrow);
	  FArrayBox& Xfab = Xlev[mfi];
	  const FArrayBox& Yfab = Ylev[mfi];
	  const FArrayBox& Zfab = Zlev[mfi];

	  tfabY.resize(gbox,nComp);
	  tfabZ.resize(gbox,nComp);
	  tfabY.copy(Yfab,sCompY,0,nComp);
	  tfabZ.copy(Zfab,sCompZ,0,nComp);

	  // Zero out parts of Y,Z covered by fine grid (to ensure valid data)
	  if (lev<nLevs-1) {
	    std::vector< std::pair<int,Box> > isects = fba.intersections(gbox);
	    for (int i = 0; i < isects.size(); i++)
	      {
		tfabY.setVal(0,isects[i].second,0,nComp);
		tfabZ.setVal(0,isects[i].second,0,nComp);
	      }
	  }

	  FORT_RS_XMULTYZ(Xfab.dataPtr(dComp),ARLIM(Xfab.loVect()), ARLIM(Xfab.hiVect()),
			  tfabY.dataPtr(),ARLIM(Yfab.loVect()), ARLIM(Yfab.hiVect()),
			  tfabZ.dataPtr(),ARLIM(Zfab.loVect()), ARLIM(Zfab.hiVect()),
			  gbox.loVect(), gbox.hiVect(), &nComp);
        }
    }
}



//
// Compute mfte[dir][comp] = Grad(mftc[comp]) + a[dir][comp]
//
void
RichardSolver::CCtoECgradAdd(PArray<MFTower>& mfte,
			     const MFTower&   mftc,
			     const FArrayBox& a,
			     int              sComp,
			     int              dComp,
			     int              nComp) const
{
  const Layout& layout = PMAmr::GetLayout();
  for (int d=0; d<BL_SPACEDIM; ++d) {            
    BL_ASSERT(layout.IsCompatible(mfte[d]));
  }
  BL_ASSERT(layout.IsCompatible(mftc));

  const Array<Geometry>& geomArray = layout.GeomArray();
  for (int lev=0; lev<nLevs; ++lev) {
    const MultiFab& mfc = mftc[lev];
    BL_ASSERT(mfc.nGrow()>=1);
    BL_ASSERT(sComp+nComp<=mfc.nComp());
    const Real* dx = geomArray[lev].CellSize();

    for (MFIter mfi(mfc); mfi.isValid(); ++mfi) {
      const FArrayBox& cfab = mfc[mfi];
      const Box& vcbox = mfi.validbox();
            
      for (int d=0; d<BL_SPACEDIM; ++d) {            
	FArrayBox& efab = mfte[d][lev][mfi];
	BL_ASSERT(dComp+nComp<=efab.nComp());
	BL_ASSERT(Box(vcbox).surroundingNodes(d).contains(efab.box()));
	efab.setVal(0);
	FORT_RS_GXPA(efab.dataPtr(dComp),ARLIM(efab.loVect()), ARLIM(efab.hiVect()),
		     cfab.dataPtr(sComp),ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
		     vcbox.loVect(),vcbox.hiVect(),dx,a.dataPtr(d),&d,&nComp);
      }
    }
  }
}

void
RichardSolver::FillPatch(MFTower& mft,
			 int sComp,
			 int nComp,
			 bool do_piecewise_constant)
{
  mftfp->FillGrowCells(mft,sComp,nComp,do_piecewise_constant);
}

void 
RichardSolver::SetInflowVelocity(PArray<MFTower>& velocity,
				 Real             time)
{
  const Layout& layout = PMAmr::GetLayout();
  for (int d=0; d<BL_SPACEDIM; ++d) {            
    BL_ASSERT(layout.IsCompatible(velocity[d]));
  }

  const Array<Geometry>& geomArray = layout.GeomArray();

  FArrayBox inflow;
  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation face = oitr();
    int dir = face.coordDir();
    for (int lev=0; lev<nLevs; ++lev) {
      MultiFab& uld = velocity[dir][lev];
      if (pm[lev].get_inflow_velocity(face,inflow,time)) {
	int shift = ( face.isHigh() ? -1 : +1 );
	inflow.shiftHalf(dir,shift);
	for (MFIter mfi(uld); mfi.isValid(); ++mfi) {
	  FArrayBox& u = uld[mfi];
	  Box ovlp = inflow.box() & u.box();
	  if (ovlp.ok()) {
	    u.copy(inflow);
	  }
	}
      }
    }
  }
}

void 
RichardSolver::DivU(MFTower& divu,
		    MFTower& rhoSat,
		    MFTower& pressure,
		    Real     time)
{
    // Coming in, the grow cells of pressure area assumed to contained valid data.  On
    // Dirichlet boundaries, the grow cell will hold the dirichlet value to apply at
    // the cell wall.  Note that since we use "calcInvPressure" to fill rho.sat, these
    // are then values on the wall as well.  As a result, the lambda values computed
    // with rho.sat are evaluated at the wall as well.  

    // We use the FillPatch operation to set pressure values in the grow cells using 
    // polynomial extrapolation, and will then use these p values only for the puposes
    // of evaluating the pressure gradient on cell faces via a simple centered difference.

    // Assumes lev=0 here corresponds to Amr.level=0, sets dirichlet values of rho.sat and
    // lambda on dirichlet pressure faces
    const Layout& layout = PMAmr::GetLayout();
    for (int lev=0; lev<nLevs; ++lev) {
        pm[lev].calcInvPressure(rhoSat[lev],pressure[lev]); // FIXME: Writes/reads only to comp=0, does 1 grow
        pm[lev].calcLambda(&(GetLambda()[lev]),rhoSat[lev]); // FIXME: Writes/reads only to comp=0, does 1 grow
    }

    int nComp = 1;
    Box abox(IntVect::TheZeroVector(),(nComp-1)*BoxLib::BASISV(0));
    FArrayBox a(abox,BL_SPACEDIM); // Make a funny box for a to simplify passing to Fortran
    for (int d=0; d<BL_SPACEDIM; ++d) {
        Real* ap = a.dataPtr(d);
        for (int n=0; n<nComp; ++n) {
            ap[n] = GetDensity()[n] * GetGravity()[d];
        }
    }

    // Convert grow cells of pressure into extrapolated values so that from here on out,
    // the values are only used to compute gradients at faces.
    bool do_piecewise_constant = false;
    //bool do_piecewise_constant = true;
    FillPatch(pressure,0,nComp,do_piecewise_constant);

    // Get  -(Grad(p) + rho.g)
    CCtoECgradAdd(GetDarcyVelocity(),pressure,a);

    // Get edge-centered lambda (= krel/mu) based on the sign of -(Grad(p) + rho.g)
    const BCRec& pressure_bc = pm[0].get_desc_lst()[Press_Type].getBC(0);
    CenterToEdgeUpwind(GetRichardCoefs(),
		       GetLambda(),
		       GetDarcyVelocity(),nComp,pressure_bc);

    // Get Darcy flux = - lambda * kappa * (Grad(p) + rho.g)
    for (int d=0; d<BL_SPACEDIM; ++d) {
        XmultYZ(GetDarcyVelocity()[d],
		GetRichardCoefs()[d],
		GetKappa()[d]);
    }    

    // Overwrite fluxes at boundary with boundary conditions
    SetInflowVelocity(GetDarcyVelocity(),time);

    // Average down fluxes
    for (int d=0; d<BL_SPACEDIM; ++d) {
        MFTower::AverageDown(GetDarcyVelocity()[d]);
    }    

    // Get the divergence of the Darcy Flux
    MFTower::ECtoCCdiv(divu,
		       GetDarcyVelocity());
}

void
RichardSolver::DpDtResidual(MFTower& residual,
			    MFTower& pressure,
			    Real     time,
			    Real     dt)
{
  DivU(residual,
       GetRhoSatNp1(),
       pressure,
       time);

  const Layout& layout = PMAmr::GetLayout();
  const Array<BoxArray>& gridArray = layout.GridArray();
  const Array<IntVect>& refRatio = layout.RefRatio();

  int nComp = 1;
  for (int lev=0; lev<nLevs; ++lev)
    {
      MultiFab& Rlev = residual[lev];
      for (MFIter mfi(Rlev); mfi.isValid(); ++mfi) {
	const Box& vbox = mfi.validbox();
	FArrayBox& Res = Rlev[mfi];
	const FArrayBox& RSn = GetRhoSatN()[lev][mfi];
	const FArrayBox& RSnp1 = GetRhoSatNp1()[lev][mfi];
	const FArrayBox& Poros = GetPorosity()[lev][mfi];
	
	FORT_RS_PDOTRES(Res.dataPtr(),ARLIM(Res.loVect()), ARLIM(Res.hiVect()),
			RSn.dataPtr(),ARLIM(RSn.loVect()), ARLIM(RSn.hiVect()),
			RSnp1.dataPtr(),ARLIM(RSnp1.loVect()), ARLIM(RSnp1.hiVect()),
			Poros.dataPtr(),ARLIM(Poros.loVect()), ARLIM(Poros.hiVect()),
		        GetDensity().dataPtr(), &dt,
			vbox.loVect(), vbox.hiVect(), &nComp);
        }
    }
}

void
RichardSolver::DpDtResidual_Alt(MFTower& residual,
                                MFTower& pressure,
                                Real     time,
                                Real     dt)
{
  DivU(residual,
       GetRhoSatNp1(),
       pressure,
       time);

  const Layout& layout = PMAmr::GetLayout();
  const Array<BoxArray>& gridArray = layout.GridArray();
  const Array<IntVect>& refRatio = layout.RefRatio();

  int nComp = 1;
  for (int lev=0; lev<nLevs; ++lev)
    {
      MultiFab& Rlev = residual[lev];
      for (MFIter mfi(Rlev); mfi.isValid(); ++mfi) {
	const Box& vbox = mfi.validbox();
	FArrayBox& Res = Rlev[mfi];
	const FArrayBox& RSn = GetRhoSatN()[lev][mfi];
        const FArrayBox& Pnp1 = pressure[lev][mfi];
	const FArrayBox& RSstar = GetRhoSatStar()[lev][mfi];
	const FArrayBox& Pstar = GetPStar()[lev][mfi];
	const FArrayBox& Astar = GetAlphaStar()[lev][mfi];
	const FArrayBox& Poros = GetPorosity()[lev][mfi];
	
	FORT_RS_PDOTALT(Res.dataPtr(),ARLIM(Res.loVect()), ARLIM(Res.hiVect()),
                        RSn.dataPtr(),ARLIM(RSn.loVect()), ARLIM(RSn.hiVect()),
                        Pnp1.dataPtr(),ARLIM(Pnp1.loVect()), ARLIM(Pnp1.hiVect()),
                        RSstar.dataPtr(),ARLIM(RSstar.loVect()), ARLIM(RSstar.hiVect()),
                        Pstar.dataPtr(),ARLIM(Pstar.loVect()), ARLIM(Pstar.hiVect()),
                        Astar.dataPtr(),ARLIM(Astar.loVect()), ARLIM(Astar.hiVect()),
                        Poros.dataPtr(),ARLIM(Poros.loVect()), ARLIM(Poros.hiVect()),
                        GetDensity().dataPtr(), &dt,
                        vbox.loVect(), vbox.hiVect(), &nComp);
        }
    }
}
