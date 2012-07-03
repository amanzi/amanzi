#include <Richard.H>
#include <RICHARD_F.H>

#include <Utility.H>

#ifdef BL_USE_PETSC
#include <petscmat.h>
#define CHKPETSC(n) CHKERRABORT(ParallelDescriptor::Communicator(),n);
#endif

#include <ostream>
std::ostream& operator<< (std::ostream&  os, const MultiFab& mf)
{
  int myproc = ParallelDescriptor::MyProc();
    os << "MultiFab" << '\n';
    for (int i=0; i<mf.size(); ++i) {
      const FArrayBox& fab = mf[i];
      if (mf.DistributionMap()[i]==myproc) {
	const Box& box = fab.box();
	os << "Box: " << box << '\n';
	for (IntVect iv(box.smallEnd()); iv<=box.bigEnd(); box.next(iv))
	  {
	    os << iv << " ";
	    for (int n=0; n<fab.nComp(); ++n) {
	      os << fab(iv,n) << " ";
	    }
	    os << '\n';
	  } 
      }
      ParallelDescriptor::Barrier();
    }
    return os;
}


CalcCoefficients::CalcCoefficients(MFTFillPatch& _mftfp,
                                   PMAmr&        _pmamr)
    : mftfp(_mftfp),
      pmamr(_pmamr),
      layout(_mftfp.GetLayout()),
      nLevs(_mftfp.GetLayout().NumLevels())
{
}


void
CalcCoefficients::operator()(PArray<MFTower>&       coefficients,
                             MFTower&               pressure,
                             MFTower&               saturation,
                             MFTower&               lambda)
{
    BL_ASSERT(layout.IsCompatible(pressure));
    BL_ASSERT(pressure.NGrow()>=1);
    BL_ASSERT(lambda.NGrow()>=pressure.NGrow());
    BL_ASSERT(saturation.NGrow()>=pressure.NGrow());

    mftfp.FillGrowCells(pressure);

    // Assumes lev=0 here corresponds to Amr.level=0
    for (int lev=0; lev<nLevs; ++lev) {
        PorousMedia* pmp = dynamic_cast<PorousMedia*>(&(pmamr.getLevel(lev)));
        if (!pmp) {
            BoxLib::Abort("Bad cast in CalcCoefficients::operator()");
        }
        MultiFab& pressureLev = pressure[lev];
        MultiFab& saturationLev = saturation[lev];
        MultiFab& lambdaLev = lambda[lev];
        pmp->calcInvPressure(saturationLev,pressureLev); // FIXME: Writes/reads only to comp=0
        pmp->calcLambda(&lambdaLev,saturationLev); // FIXME: Writes/reads only to comp=0
    }
}

RichardContext::RichardContext(PMAmr&           _pmamr,
                               MFTFillPatch&     _mftfp,
                               MFTower&         _pressure_old,
                               MFTower&         _saturation,
                               MFTower&         _lambda,
                               PArray<MFTower>& _darcyVelocity,
                               const BCRec&     pressure_bc,
                               int              pressure_maxorder)
    : pmamr(_pmamr), 
      mftfp(_mftfp),
      pressure_old(_pressure_old),
      saturation(_saturation),
      lambda(_lambda),
      darcyVelocity(_darcyVelocity),
      density(PorousMedia::Density()),
      gravity(BL_SPACEDIM,0)
{
    gravity[BL_SPACEDIM-1] = PorousMedia::getGravity();
    mftfp.BuildStencil(pressure_bc, pressure_maxorder);
    calcCoefs = new CalcCoefficients(mftfp,pmamr);
}

RichardOp::RichardOp(RichardContext& _richardContext)
    : richardContext(_richardContext),
      pmamr(_richardContext.PMAMR()),
      coefs(BL_SPACEDIM,PArrayManage),
      bc_initialized(false)
{
    Array<IndexType> itype(BL_SPACEDIM);
    const Layout& layout = GetLayout(); 
    for (int d=0; d<BL_SPACEDIM; ++d) {
        coefs.set(d, new MFTower(layout,MFTower::EC[d]));
    }
}

RichardOp::~RichardOp()
{
}

void
RichardOp::CenterToEdgeUpwind(PArray<MFTower>&       mfte,
                              MFTower&               mftc,
                              const PArray<MFTower>& sgn,
                              int                    nComp) const
{
    int nLevs = richardContext.GetLayout().NumLevels();
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
                
                FORT_RICHARD_CTE_UPW(efab.dataPtr(), ARLIM(efab.loVect()), ARLIM(efab.hiVect()),
                                     cfab.dataPtr(), ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
                                     sgnfab.dataPtr(),ARLIM(sgnfab.loVect()), ARLIM(sgnfab.hiVect()),
                                     vccbox.loVect(), vccbox.hiVect(), &d, &nComp);
            }
        }
    }
}

void 
RichardOp::YmultX(MFTower&           Y,
                  const MFTower&     X,
                  int                sComp,
                  int                dComp,
                  int                nComp,
                  int                nGrow) const
{
    const Layout& layout = GetLayout(); 
    BL_ASSERT(layout.IsCompatible(X));
    BL_ASSERT(layout.IsCompatible(Y));
    BL_ASSERT(X.NComp()>=sComp+nComp);
    BL_ASSERT(Y.NComp()>=dComp+nComp);
    BL_ASSERT(X.NGrow()>=nGrow);
    BL_ASSERT(Y.NGrow()>=nGrow);
    const Array<BoxArray>& gridArray = GridArray();
    const Array<Geometry>& geomArray = GeomArray();
    const Array<IntVect>& refRatio = RefRatio();

    FArrayBox fab;
    int nLevs = layout.NumLevels();
    for (int lev=0; lev<nLevs; ++lev)
    {
        MultiFab& Ylev = Y[lev];
        const MultiFab& Xlev = X[lev];
        const BoxArray& ba = Xlev.boxArray();
        BoxArray fba;
        if (lev<nLevs-1) {
            fba = BoxArray(ba).coarsen(refRatio[lev]);
        }

        for (MFIter mfi(Y[lev]); mfi.isValid(); ++mfi)
        {
            const Box& vbox = mfi.validbox();
            Box gbox = Box(vbox).grow(nGrow);
            FArrayBox& Yfab = Ylev[mfi];
            const FArrayBox& Xfab = Xlev[mfi];

            fab.resize(gbox,nComp);
            fab.copy(Xfab,sComp,0,nComp);

            // Zero out parts of X covered by fine grid (to ensure valid data)
            if (lev<nLevs-1) {
                std::vector< std::pair<int,Box> > isects = fba.intersections(gbox);
                for (int i = 0; i < isects.size(); i++)
                {
                    fab.setVal(0,isects[i].second,0,nComp);
                }
            }

            FORT_RICHARD_YMULTX(Yfab.dataPtr(sComp),ARLIM(Yfab.loVect()), ARLIM(Yfab.hiVect()),
                                Xfab.dataPtr(dComp),ARLIM(Xfab.loVect()), ARLIM(Xfab.hiVect()),
                                gbox.loVect(), gbox.hiVect(), &nComp);
        }
    }
}



//
// Compute mfte[dir][comp] = Grad(mftc[comp]) + a[dir][comp]
//
void
RichardOp::CCtoECgradAdd(PArray<MFTower>& mfte,
                         const MFTower&   mftc,
                         const FArrayBox& a,
                         int              sComp,
                         int              dComp,
                         int              nComp) const
{
    const Layout& layout = richardContext.GetLayout();
    for (int d=0; d<BL_SPACEDIM; ++d) {            
        BL_ASSERT(layout.IsCompatible(mfte[d]));
    }
    BL_ASSERT(layout.IsCompatible(mftc));

    const Array<Geometry>& geomArray = layout.GeomArray();
    int nLevs = layout.NumLevels();
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
                FORT_RICHARD_GXPA(efab.dataPtr(dComp),ARLIM(efab.loVect()), ARLIM(efab.hiVect()),
                                  cfab.dataPtr(sComp),ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
                                  vcbox.loVect(),vcbox.hiVect(),dx,a.dataPtr(d),&d,&nComp);
            }
        }
    }
}

void
RichardContext::FillPatch(MFTower& mft,
                          bool     do_piecewise_constant)
{
    mftfp(mft,do_piecewise_constant);
}

void
RichardOp::FillPatch(MFTower& mft,
                     bool do_piecewise_constant)
{
    richardContext.FillPatch(mft,do_piecewise_constant);
}

void
RichardOp::Residual(MFTower&   residual,
                    MFTower&   pressure,
                    const Real dt)
{
    // Grab refs to temporary data
    MFTower& pressure_old = richardContext.PressureOld();
    MFTower& saturation = richardContext.Saturation();
    MFTower& lambda = richardContext.Lambda();
    PArray<MFTower>& velocity = richardContext.DarcyVelocity();
    CalcCoefficients& calcCoefs = richardContext.CalcCoefs();
    int nComp = 1;
    int sComp = 0; 
    int dComp = 0;
    Box abox(IntVect::TheZeroVector(),(nComp-1)*BoxLib::BASISV(0));
    FArrayBox a(abox,BL_SPACEDIM); // Make a funny box for a to simplify passing to Fortran
    for (int d=0; d<BL_SPACEDIM; ++d) {
        Real* ap = a.dataPtr(d);
        for (int n=0; n<nComp; ++n) {
            ap[n] = - richardContext.Density()[n] * richardContext.Gravity()[d];
        }
    }
    
    // Update cell-centered coefficients
    calcCoefs(coefs,pressure,saturation,lambda);

    // Get  Grad(p) - rho.g
    FillPatch(pressure);
    CCtoECgradAdd(velocity,pressure,a,sComp,dComp,nComp);

    // Get edge-centered coefficients based on the sign of (Grad(p) - rho.g)
    CenterToEdgeUpwind(coefs,lambda,velocity,nComp);

    // Get Darcy flux = H * (Grad(p) - rho.g). average down
    for (int d=0; d<BL_SPACEDIM; ++d) {
        BL_ASSERT(velocity[d].NComp()>dComp);
        YmultX(velocity[d],coefs[d]);
        MFTower::AverageDown(velocity[d]);
    }

    // Get the divergence of the Darcy Flux
    MFTower::ECtoCCdiv(residual,velocity);
}

void
RichardOp::BuildOpSkel()
{
#ifdef BL_USE_PETSC
    Layout& layout = richardContext.GetLayout();
    const MFTFillPatch& fp = richardContext.GetFillPatch();

    Mat& J = layout.Jacobian();
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
    const Array<Array<IVSMap> >& growCellStencil = fp.GrowCellStencil();
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
#endif
}

void
RichardOp::Write(const std::string& fileName) const
{
    std::string FullPath = fileName;
    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
        FullPath += '/';
    
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
            BoxLib::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();
    
    const Layout& layout = richardContext.GetLayout();
    int nLevs = layout.NumLevels();
    for (int lev=0; lev<nLevs; ++lev) {
        std::string FullPath = fileName;
        if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
        {
            FullPath += '/';
        }
        std::string LevelPath = FullPath + "RichardOp_Level_";
        LevelPath = BoxLib::Concatenate(LevelPath,lev,1);
        
        if (ParallelDescriptor::IOProcessor()) {
            if (!BoxLib::UtilCreateDirectory(LevelPath, 0755)) {
                BoxLib::CreateDirectoryFailed(LevelPath);
            }
        }
        ParallelDescriptor::Barrier();

        LevelPath += '/';
        for (int d=0; d<BL_SPACEDIM; ++d) {
            std::string DirLevelPathCoefs = LevelPath + "coefs_";
            DirLevelPathCoefs = BoxLib::Concatenate(DirLevelPathCoefs,d,1);
            VisMF::Write(coefs[d][lev],DirLevelPathCoefs);
        }
        std::string LevelPathPressureOld = LevelPath + "pressure_old";
        VisMF::Write(richardContext.PressureOld()[lev],LevelPathPressureOld);
        
        std::string LevelPathLambda = LevelPath + "lambda";
        VisMF::Write(richardContext.Lambda()[lev],LevelPathLambda);
        
        std::string LevelPathSaturation = LevelPath + "saturation";
        VisMF::Write(richardContext.Saturation()[lev],LevelPathSaturation);
        
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::string MiscPath = FullPath + "miscData.dat";
        std::ofstream os; os.open(MiscPath.c_str());
        os << "Number of Levels: " << nLevs << '\n';
        os << "Refinement Ratios: ";
        for (int lev=1; lev<nLevs; ++lev) {
            os << layout.RefRatio()[lev-1] << " ";
        }
        os << '\n';
        os << "Grid Arrays: " << '\n';
        for (int lev=0; lev<nLevs; ++lev) {
            os << layout.GridArray()[lev] << '\n';
        }
        os.close();
    }
}
