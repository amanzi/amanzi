#include <Layout.H>

std::ostream& operator<< (std::ostream&  os, const Node& node)
{
    os << "Node: IntVect=" << node.iv << ", level=" << node.level << ", grid=" << node.grid << ", type="; 
    if (node.type==Node::INIT)
        os << "INIT";
    else if (node.type==Node::COVERED)
        os << "COVERED";
    else
        os << "VALID";
    if (os.fail())
        BoxLib::Error("operator<<(std::ostream&,Node&) failed");
    return os;
}

static
BoxArray
GetBndryCells (const BoxArray& ba,
               const IntVect&  ngrow,
               const Geometry& geom,
               const DistributionMapping* dm=0) // last arg, get local boundary
{
    //
    // First get list of all ghost cells.
    //
    BoxList gcells, bcells;

    for (int i = 0; i < ba.size(); ++i) {
	gcells.join(BoxLib::boxDiff(BoxLib::grow(ba[i],ngrow),ba[i]));
    }

    BoxArray lba;
    if (dm) {
        BoxList sbl;
        for (int i=0; i<ba.size(); ++i) {
            if ((*dm)[i]==ParallelDescriptor::MyProc()) {
                sbl.push_back(ba[i]);
            }
        }
        lba = BoxArray(sbl);
    }
    else {
        lba = ba;
    }
    //
    // Now strip out intersections with original BoxArray (or, local boxes)
    //
    for (BoxList::const_iterator it = gcells.begin(); it != gcells.end(); ++it)
    {
        std::vector< std::pair<int,Box> > isects = lba.intersections(*it);

        if (isects.empty())
            bcells.push_back(*it);
        else
        {
            //
            // Collect all the intersection pieces.
            //
            BoxList pieces;
            for (int i = 0; i < isects.size(); i++)
                pieces.push_back(isects[i].second);
            BoxList leftover = BoxLib::complementIn(*it,pieces);
            bcells.catenate(leftover);
        }
    }
    //
    // Now strip out overlaps.
    //
    gcells.clear();
    gcells = BoxLib::removeOverlap(bcells);
    bcells.clear();

    if (geom.isAnyPeriodic())
    {
        Array<IntVect> pshifts(27);

        const Box domain = geom.Domain();

        for (BoxList::const_iterator it = gcells.begin(); it != gcells.end(); ++it)
        {
            if (!domain.contains(*it))
            {
                //
                // Add in periodic ghost cells shifted to valid region.
                //
                geom.periodicShift(domain, *it, pshifts);

                for (int i = 0; i < pshifts.size(); i++)
                {
                    const Box shftbox = *it + pshifts[i];

                    const Box ovlp = domain & shftbox;
                    BoxList bl = BoxLib::complementIn(ovlp,BoxList(lba));
                    bcells.catenate(bl);
                }
            }
        }

        gcells.catenate(bcells);
    }

    return BoxArray(gcells);
}

MFTower::MFTower()
{
}

MFTower::MFTower(const Array<BoxArray>& bat,
                 const Array<IntVect>&  ratio,
                 int                    nComp,
                 int                    nGrow)
{
    ref_ratio.clear();
    mft.clear();
    define(bat,ratio,nComp,nGrow);
}

MFTower::MFTower(const PArray<MultiFab>& pamf,
                 const Array<IntVect>&   ratio)
{
    int nLevs = pamf.size();
    mft.resize(nLevs, PArrayNoManage);
    ref_ratio.resize(nLevs-1);
    for (int lev=0; lev<nLevs; ++lev)
    {
        if (lev<nLevs-1) {
            ref_ratio[lev] = ratio[lev];
        }
        mft.set(lev,&(const_cast<MultiFab&>(pamf[lev]))); // careful...
    }
}

void
MFTower::define(const Array<BoxArray>& bat,
                const Array<IntVect>&  ratio,
                int                    nComp,
                int                    nGrow)
{
    int nLevs = bat.size();
    BL_ASSERT(nLevs == ratio.size()+1);
    mft.resize(bat.size(), PArrayManage);
    ref_ratio.resize(nLevs);
    for (int lev=0; lev<bat.size(); ++lev)
    {
        if (lev<nLevs-1) {
            ref_ratio[lev] = ratio[lev];
        }
        mft.set(lev,new MultiFab(bat[lev],nComp,nGrow));
    }
}

void
MFTower::AXPY(const MFTower& rhs,
              Real           p)
{
    BL_ASSERT(IsCompatible(rhs));
    int nLevs = mft.size();
    FArrayBox fab;
    for (int lev=0; lev<nLevs; ++lev)
    {
        for (MFIter mfi(mft[lev]); mfi.isValid(); ++mfi)
        {
            Box vbox = mfi.validbox();
            fab.resize(vbox,1);
            fab.copy(rhs[lev][mfi]);
            if (lev<nLevs-1) {
                BoxArray cba = BoxArray(mft[lev-1].boxArray()).refine(refRatio(lev));
                std::vector< std::pair<int,Box> > isects = cba.intersections(vbox);
                for (int i = 0; i < isects.size(); i++)
                {
                    Box ovlp = isects[i].second & vbox;
                    if (ovlp.ok()) {
                        fab.setVal(0);
                    }
                }
            }

            if (p!=1) {
                fab.mult(p);
            }
            mft[lev][mfi].plus(fab);
        }
    }
}

Real
MFTower::norm() const // currently only max norm supported
{
    int nLevs = mft.size();
    FArrayBox fab;
    Real norm = 0;
    for (int lev=0; lev<nLevs; ++lev)
    {
        for (MFIter mfi(mft[lev]); mfi.isValid(); ++mfi)
        {
            Box vbox = mfi.validbox();
            fab.resize(vbox,1);
            fab.copy(mft[lev][mfi]);
            if (lev<nLevs-1) {
                BoxArray cfba = BoxArray(mft[lev+1].boxArray()).coarsen(refRatio(lev));
                std::vector< std::pair<int,Box> > isects = cfba.intersections(vbox);
                for (int i = 0; i < isects.size(); i++)
                {
                    Box ovlp = isects[i].second & vbox;
                    if (ovlp.ok()) {
                        fab.setVal(0);
                    }
                }
            }
            norm = std::max(norm, fab.norm(2,0,1));
        }
    }
    ParallelDescriptor::ReduceRealMax(norm);
    return norm;
}

bool
MFTower::IsCompatible(const MFTower& rhs) const
{
    int nLevs = mft.size();
    bool isok = nLevs==rhs.size();
    for (int lev=0; lev<nLevs; ++lev)
    {
        isok &= rhs[lev].boxArray()==rhs[lev].boxArray();
        if (lev < nLevs-1) {
            isok &= rhs.refRatio(lev)==refRatio(lev);
        }
    }
    return isok;
}

#ifdef BL_USE_PETSC
Layout::Layout(Amr* _parent)
    : initialized(false),
      nGrow(1),
      parent(_parent)
{
}

Layout::~Layout()
{
    Clear();
}

void
Layout::SetParent(Amr* new_parent)
{
    if (parent) {
        Clear();
    }
    parent = new_parent;
}

void
Layout::DestroyPetscStructures()
{
    for (int i=0; i<mats_I_created.size(); ++i) {
        PetscErrorCode ierr = MatDestroy(mats_I_created[i]); CHKPETSC(ierr);
    }
    mats_I_created.clear();
    for (int i=0; i<vecs_I_created.size(); ++i) {
        PetscErrorCode ierr = VecDestroy(vecs_I_created[i]); CHKPETSC(ierr);
    }
    vecs_I_created.clear();
}

void
Layout::Clear()
{
    nodes.clear();
    nodeIds.clear();
    DestroyPetscStructures();
    initialized = false;
}

#include <VisMF.H>
void 
Write(const Layout::MultiIntFab& imf,
      const std::string& filename,
      int comp,
      int nGrow) 
{
    BL_ASSERT(comp<=imf.nComp() && comp>=0);
    BL_ASSERT(nGrow<=imf.nGrow());
    MultiFab mf(imf.boxArray(), 1, imf.nGrow());
    for (MFIter mfi(imf); mfi.isValid(); ++mfi)
    {
        const int* iptr = imf[mfi].dataPtr(comp);
        Real* rptr = mf[mfi].dataPtr();
        int numpts = imf[mfi].box().numPts();
        for (int i=0; i<numpts; ++i) {
            rptr[i] = Real(iptr[i]);
        }
    }
    VisMF::Write(mf,filename);
}

void 
Layout::Rebuild()
{
    nGrow = 1;
    BL_ASSERT(parent);

    nLevs = parent->finestLevel() + 1;

    gridArray.resize(nLevs);
    geomArray.resize(nLevs);
    refRatio.resize(nLevs-1);
    for (int i=0; i<nLevs; ++i) {
        gridArray[i] = parent->boxArray(i);
        geomArray[i] = parent->Geom(i);
        if (i<nLevs-1) {
            refRatio[i] = parent->refRatio(i);
        }
    }

    // Declare here so we can use again below
    Array<BoxArray> bndryCells(nLevs);

    Clear();
    nodes.resize(nLevs,PArrayManage);
    nodeIds.resize(nLevs,PArrayManage);
    int cnt = 0; // local node counter
    for (int lev=0; lev<nLevs; ++lev)
    {
        nodes.set(lev,new MultiNodeFab(gridArray[lev],1,nGrow));
        nodeIds.set(lev,new MultiIntFab(gridArray[lev],1,nGrow));

        for (MFIter fai(nodes[lev]); fai.isValid(); ++fai)
        {
            NodeFab& ifab = nodes[lev][fai];
            const Box box = ifab.box() & gridArray[lev][fai.index()];
            for (IntVect iv=box.smallEnd(); iv<=box.bigEnd(); box.next(iv))
                ifab(iv,0) = Node(iv,lev,fai.index(),Node::VALID);
        }
            
        if (lev > 0)
        {
            const Box rangeBox = Box(IntVect::TheZeroVector(),
                                     refRatio[lev-1] - IntVect::TheUnitVector());

            bndryCells[lev] = GetBndryCells(nodes[lev].boxArray(),refRatio[lev-1],geomArray[lev]);

            for (MFIter fai(nodes[lev]); fai.isValid(); ++fai)
            {
                const Box box = Box(fai.validbox()).grow(refRatio[lev-1]) & gridArray[lev][fai.index()];
                NodeFab& ifab = nodes[lev][fai];
                std::vector< std::pair<int,Box> > isects = bndryCells[lev].intersections(box);
                for (int i = 0; i < isects.size(); i++)
                {
                    Box co = isects[i].second & fai.validbox();
                    if (co.ok())
                        std::cout << "BAD ISECTS: " << co << std::endl;

                    const Box& dstBox = isects[i].second;
                    const Box& srcBox = BoxLib::coarsen(dstBox,refRatio[lev-1]);

                    NodeFab dst(dstBox,1);
                    for (IntVect iv(srcBox.smallEnd());
                         iv<=srcBox.bigEnd();
                         srcBox.next(iv))
                    {
                        const IntVect baseIV = refRatio[lev-1] * iv;
                        for (IntVect ivt(rangeBox.smallEnd());ivt<=rangeBox.bigEnd();rangeBox.next(ivt))
                            dst(baseIV + ivt,0) = Node(iv,lev-1,-1,Node::VALID);
                    }
                    const Box ovlp = dstBox & ifab.box();

                    Box mo = ovlp & fai.validbox();
                    if (mo.ok())
                    {
                        std::cout << "BAD OVERLAP: " << mo << std::endl;
                        std::cout << "         vb: " << fai.validbox() << std::endl;
                    }
                    if (ovlp.ok())
                        ifab.copy(dst,ovlp,0,ovlp,0,1);
                }
            }
        }

        // Block out cells covered by finer grid
        if (lev < nLevs-1)
        {
            const BoxArray coarsenedFineBoxes =
                BoxArray(gridArray[lev+1]).coarsen(refRatio[lev]);
                
            for (MFIter fai(nodes[lev]); fai.isValid(); ++fai)
            {
                NodeFab& ifab = nodes[lev][fai];
                const Box& box = ifab.box();
                std::vector< std::pair<int,Box> > isects = coarsenedFineBoxes.intersections(box);
                for (int i = 0; i < isects.size(); i++)
                {
                    const Box& ovlp = isects[i].second;
                    for (IntVect iv=ovlp.smallEnd(); iv<=ovlp.bigEnd(); ovlp.next(iv))
                        ifab(iv,0) = Node(iv,lev,fai.index(),Node::COVERED);
                }
            }
        }

        // Set nodeIds
        nodeIds[lev].setVal(-1);
        for (MFIter fai(nodes[lev]); fai.isValid(); ++fai)
        {
            NodeFab& ifab = nodes[lev][fai];
            IntFab& nfab = nodeIds[lev][fai];
            const Box& box = fai.validbox() & gridArray[lev][fai.index()];
            for (IntVect iv(box.smallEnd()); iv<=box.bigEnd(); box.next(iv))
            {
                if (ifab(iv,0).type == Node::VALID)
                {
                    if (ifab(iv,0).level != lev) 
                        std::cout << "bad level: " << ifab(iv,0) << std::endl;
                    nfab(iv,0) = cnt++;
                }
            }
        }
    }
    nNodes_local = cnt;
    nNodes_global = nNodes_local;

#if BL_USE_MPI
    // Adjust node numbers to be globally unique
    int num_size = ParallelDescriptor::NProcs();
    Array<int> num_nodes_p(num_size);
    BL_MPI_REQUIRE(MPI_Allgather(&nNodes_local, 1, MPI_INT, num_nodes_p.dataPtr(), 1, MPI_INT, 
                                 ParallelDescriptor::Communicator()));

    int offset = 0;
    for (int i=0; i<ParallelDescriptor::MyProc(); ++i) {
        offset += num_nodes_p[i];
    }

    // Adjust node numbers 
    for (int lev=0; lev<nLevs; ++lev)
    {
        for (MFIter mfi(nodeIds[lev]); mfi.isValid(); ++mfi)
        {
            IntFab& nfab = nodeIds[lev][mfi];
            const Box& box = mfi.validbox() & gridArray[lev][mfi.index()];
            for (IntVect iv(box.smallEnd()); iv<=box.bigEnd(); box.next(iv))
            {
                if (nfab(iv,0) >= 0)
                {
                    nfab(iv,0) += offset;
                }
            }
        }
    }

    // Compute total number of nodes
    nNodes_global = 0;
    for (int i=0; i<ParallelDescriptor::NProcs(); ++i) {
        nNodes_global += num_nodes_p[i];
    }

    // Now communicate node numbering to neighbors grow cells
    for (int lev=1; lev<nLevs; ++lev) // This is really concerned with filling c-f grow cells
    {
        BoxArray bndC = BoxArray(bndryCells[lev]).coarsen(refRatio[lev-1]);
        
        MultiIntFab crseIds(bndC,1,0); crseIds.setVal(-1);
        crseIds.copy(nodeIds[lev-1]); // parallel copy

        // "refine" crseIds
        MultiIntFab fineIds(bndryCells[lev],1,0); fineIds.setVal(-1);
        const Box rangeBox = Box(IntVect::TheZeroVector(),
                                 refRatio[lev-1] - IntVect::TheUnitVector());
        for (MFIter mfi(crseIds); mfi.isValid(); ++mfi)
        {
            const Box& cbox = crseIds[mfi].box();
            for (IntVect iv = cbox.smallEnd(), End=cbox.bigEnd(); iv<=End; cbox.next(iv)) {
                int nodeIdx = crseIds[mfi](iv,0);
                const IntVect baseIV = refRatio[lev-1] * iv;
                for (IntVect ivt = rangeBox.smallEnd(), End=rangeBox.bigEnd(); ivt<=End;rangeBox.next(ivt))
                    fineIds[mfi](baseIV + ivt,0) = nodeIdx;
            }
        }

        nodeIds[lev].FillBoundary(0,1);
        BoxLib::FillPeriodicBoundary<IntFab>(geomArray[lev],nodeIds[lev],0,1);

        MultiIntFab ng(BoxArray(nodeIds[lev].boxArray()).grow(nodeIds[lev].nGrow()),1,0);
        for (MFIter mfi(nodeIds[lev]); mfi.isValid(); ++mfi)
        {
            ng[mfi].copy(nodeIds[lev][mfi]); // get valid + f-f (+periodic f-f)
        }
        
        ng.copy(fineIds); // Parallel copy to get c-f from bndryCells

        for (MFIter mfi(nodeIds[lev]); mfi.isValid(); ++mfi)
        {
            nodeIds[lev][mfi].copy(ng[mfi]); // put it all back
        }

    }
#else
    for (int lev=0; lev<nLevs; ++lev)
    {
        nodeIds[lev].FillBoundary(0,1);
        BoxLib::FillPeriodicBoundary<IntFab>(geomArray[lev],nodeIds[lev],0,1);
    }    
#endif

    int n = nNodes_local; // Number of local columns of J
    int m = nNodes_local; // Number of local rows of J
    int N = nNodes_global; // Number of global columns of J 
    int M = nNodes_global; // Number of global rows of J 
    int d_nz = 1 + 2*BL_SPACEDIM; // Number of nonzero local columns of J
    int o_nz = 0; // Number of nonzero nonlocal (off-diagonal) columns of J
    PetscErrorCode ierr = MatCreateMPIAIJ(ParallelDescriptor::Communicator(), m, n, M, N, d_nz, PETSC_NULL, o_nz, PETSC_NULL, &J_mat); CHKPETSC(ierr);
    mats_I_created.push_back(&J_mat);
    ierr = VecCreateMPI(ParallelDescriptor::Communicator(),nNodes_local,nNodes_global,&JRowScale_vec); CHKPETSC(ierr);        
    vecs_I_created.push_back(&JRowScale_vec);
    initialized = true;
}

void
Layout::SetNodeIds(BaseFab<int>& idFab, int lev, int grid) const
{
    BL_ASSERT(initialized);
    BL_ASSERT(lev < nodes.size());
    const MultiIntFab& nodeIds_lev = nodeIds[lev];
    const DistributionMapping dm = nodeIds_lev.DistributionMap();
    BL_ASSERT(grid < nodeIds_lev.size());
    BL_ASSERT(dm[grid] == ParallelDescriptor::MyProc());

    const IntFab& nFab = nodeIds_lev[grid];
    Box ovlp = nFab.box() & idFab.box();
    if (ovlp.ok()) {
        idFab.copy(nFab,ovlp,0,ovlp,0,1);
    }
}

Mat&
Layout::Jacobian()
{
    return J_mat;
}

Vec&
Layout::JRowScale()
{
    return JRowScale_vec;
}

bool
Layout::IsCompatible(const MFTower& mft) const
{
    int myproc = ParallelDescriptor::MyProc();
    bool isio = ParallelDescriptor::IOProcessor();
    bool isok = mft.size()==nLevs;
    if (isio) {
        for (int lev=0; lev<nLevs; ++lev)
        {
            isok &= mft[lev].boxArray()==gridArray[lev];
        }
    }
    return isok;
}

void
Layout::BuildMFTower(MFTower& mft,
                     int      nCompMF,
                     int      nGrowMF) const
{
    mft.define(gridArray,refRatio,nCompMF,nGrowMF);
}


PetscErrorCode
Layout::MFTowerToVec(Vec&           V,
                     const MFTower& mft,
                     int            comp)
{
    BL_ASSERT(initialized);
    BL_ASSERT(IsCompatible(mft));

    int myproc = ParallelDescriptor::MyProc();
    bool isio = ParallelDescriptor::IOProcessor();


    PetscErrorCode ierr;
    ierr = VecCreateMPI(ParallelDescriptor::Communicator(),nNodes_local,nNodes_global,&V); CHKPETSC(ierr);
    vecs_I_created.push_back(&V);

    FArrayBox fab;
    IntFab ifab;
    for (int lev=0; lev<nLevs; ++lev)
    {
        for (MFIter mfi(nodeIds[lev]); mfi.isValid(); ++mfi)
        {
            // NOTE: We built the node numbering so that the following operation works correctly
            const Box& vbox = mfi.validbox();

            BoxArray ba(vbox);
            if (lev<nLevs-1) {
                BoxArray cfba = BoxArray(mft[lev+1].boxArray()).coarsen(refRatio[lev]);
                ba = BoxLib::complementIn(vbox,cfba);
            }

            for (int i=0; i<ba.size(); ++i) {
                const Box& box = ba[i];
                int ni = box.numPts();

                ifab.resize(box,1);
                ifab.copy(nodeIds[lev][mfi]);
                BL_ASSERT(ifab.min()>=0);

                const int* ix = ifab.dataPtr();
                
                fab.resize(box,1);
                fab.copy(mft[lev][mfi]);
                const Real* y = fab.dataPtr();
                
                ierr = VecSetValues(V,ni,ix,y,INSERT_VALUES); CHKPETSC(ierr);                
            }
        }
    }
    ierr = VecAssemblyBegin(V); CHKPETSC(ierr);
    ierr = VecAssemblyEnd(V); CHKPETSC(ierr);
    return ierr;
}

PetscErrorCode
Layout::VecToMFTower(MFTower&   mft,
                     const Vec& V,
                     int        comp) const
{
    BL_ASSERT(initialized);
    BL_ASSERT(IsCompatible(mft));

    int myproc = ParallelDescriptor::MyProc();
    bool isio = ParallelDescriptor::IOProcessor();

    PetscErrorCode ierr;

    FArrayBox fab;
    IntFab ifab;
    ierr = VecAssemblyBegin(V); CHKPETSC(ierr);
    ierr = VecAssemblyEnd(V); CHKPETSC(ierr);
    for (int lev=0; lev<nLevs; ++lev)
    {
        if (mft[lev].nGrow()) {
            mft[lev].setBndry(0);
        }

        for (MFIter mfi(nodeIds[lev]); mfi.isValid(); ++mfi)
        {
            const Box& vbox = mfi.validbox();
            BoxArray ba(vbox);
            if (lev<nLevs-1) {
                BoxArray cfba = BoxArray(mft[lev+1].boxArray()).coarsen(refRatio[lev]);
                ba = BoxLib::complementIn(vbox,cfba);
            }

            for (int i=0; i<ba.size(); ++i) {
                const Box& box = ba[i];
                int ni = box.numPts();

                ifab.resize(box,1);
                ifab.copy(nodeIds[lev][mfi]);
                const int* ix = ifab.dataPtr();
                
                fab.resize(box,1);
                fab.copy(mft[lev][mfi]);
                Real* y = fab.dataPtr();

                ierr = VecGetValues(V,ni,ix,y); CHKPETSC(ierr);
                mft[lev][mfi].copy(fab);
            }
        }
    }
    return ierr;
}

std::ostream& operator<< (std::ostream&  os, const Layout::IntFab& ifab)
{
    os << "IntFab" << '\n';
    const Box& box = ifab.box();
    os << "Box: " << box << '\n';
    for (IntVect iv(box.smallEnd()); iv<=box.bigEnd(); box.next(iv))
    {
        os << iv << " ";
        for (int n=0; n<ifab.nComp(); ++n) {
            os << ifab(iv,n) << " ";
        }
        os << '\n';
    } 
    return os;
}

#endif
