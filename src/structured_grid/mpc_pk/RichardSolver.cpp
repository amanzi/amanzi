#include <RichardSolver.H>
#include <RICHARDSOLVER_F.H>

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
  PArray<MultiFab> utmp(BL_SPACEDIM,PArrayNoManage);
  
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

  PetscErrorCode ierr;       
  int n = layout.NumberOfLocalNodeIds();
  int N = layout.NumberOfGlobalNodeIds();
  MPI_Comm comm = ParallelDescriptor::Communicator();
  ierr = VecCreateMPI(comm,n,N,&RhsV); CHKPETSC(ierr);
  ierr = VecCreateMPI(comm,n,N,&SolnV); CHKPETSC(ierr);
  
  int pressure_maxorder = 4;
  const BCRec& pressure_bc = pm[0].get_desc_lst()[Press_Type].getBC(0);
  mftfp->BuildStencil(pressure_bc, pressure_maxorder);

  gravity.resize(BL_SPACEDIM,0);
  gravity[BL_SPACEDIM-1] = PorousMedia::getGravity();
  density = PorousMedia::Density();
  
  int d_nz = 1 + 2*BL_SPACEDIM; // Estmated number of nonzero local columns of J
  int o_nz = 0; // Estimated number of nonzero nonlocal (off-diagonal) columns of J
  Mat Jac;
  ierr = MatCreateMPIAIJ(comm, n, n, N, N, d_nz, PETSC_NULL, o_nz, PETSC_NULL, &Jac); CHKPETSC(ierr);
  
  BuildOpSkel(Jac);
  
  matfdcoloring = 0;
  ierr = SNESCreate(comm,&snes); CHKPETSC(ierr);
  ierr = SNESSetFunction(snes,RhsV,RichardRes_DpDt,(void*)(this)); CHKPETSC(ierr);
  ierr = MatGetColoring(Jac,MATCOLORINGSL,&iscoloring); CHKPETSC(ierr);
  ierr = MatFDColoringCreate(Jac,iscoloring,&matfdcoloring); CHKPETSC(ierr);
  ierr = MatFDColoringSetFunction(matfdcoloring,
				  (PetscErrorCode (*)(void))RichardRes_DpDt,
				  (void*)(this)); CHKPETSC(ierr);
  ierr = MatFDColoringSetFromOptions(matfdcoloring); CHKPETSC(ierr);
  ierr = SNESSetJacobian(snes,Jac,Jac,SNESDefaultComputeJacobianColor,matfdcoloring);CHKPETSC(ierr); 
  ierr = SNESSetFromOptions(snes);CHKPETSC(ierr);

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
  std::string name = BoxLib::Concatenate("stuff",myproc,1);
  
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
	    if (dst_proc != myproc) {
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

	      int num_cols = neighbors.size();
	      cols.resize(num_cols);
	      vals.resize(num_cols,0);
	      int cnt = 0;
	      for (std::set<int>::const_iterator it=neighbors.begin(), End=neighbors.end(); it!=End; ++it) {
		cols[cnt++] = *it;
	      }
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
				  int                    nComp) const
{
    for (int lev=0; lev<nLevs; ++lev) {
        MultiFab& clev = mftc[lev];
        BL_ASSERT(nComp<=clev.nComp());
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
                
                FORT_RS_CTE_UPW(efab.dataPtr(), ARLIM(efab.loVect()), ARLIM(efab.hiVect()),
				cfab.dataPtr(), ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
				sgnfab.dataPtr(),ARLIM(sgnfab.loVect()), ARLIM(sgnfab.hiVect()),
				vccbox.loVect(), vccbox.hiVect(), &d, &nComp);
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
			 bool do_piecewise_constant)
{
  (*mftfp)(mft,do_piecewise_constant);
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
    // FIXME: Need to update dirichlet data in pressure grow cells

    // fill grow cells of pressure
    FillPatch(pressure);

    // Assumes lev=0 here corresponds to Amr.level=0
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
            ap[n] = - GetDensity()[n] * GetGravity()[d];
        }
    }

    // Get  -Grad(p) - rho.g
    CCtoECgradAdd(GetDarcyVelocity(),pressure,a);

    // Get edge-centered lambda (= krel/mu) based on the sign of (Grad(p) - rho.g)
    CenterToEdgeUpwind(GetRichardCoefs(),
		       GetLambda(),
		       GetDarcyVelocity(),nComp);

    // Get Darcy flux = lambda * kappa * (Grad(p) - rho.g)
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
