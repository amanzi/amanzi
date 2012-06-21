#include <Layout.H>

std::ostream& operator<< (std::ostream&  os, const Node& node)
{
    os << "iv:" << node.iv << ",l:" << node.level << ",t:"; 
    if (node.type==Node::INIT)
        os << "I";
    else if (node.type==Node::COVERED)
        os << "C";
    else
        os << "V";
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

MFTower::MFTower(const Layout&    _layout,
                 const IndexType& t,
                 int              _nComp,
                 int              _nGrow)
    : layout(_layout), iType(t), nComp(_nComp), nGrow(_nGrow), nLevs(layout.NumLevels()),
      gridArray(layout.GridArray()), geomArray(layout.GeomArray()), refRatio(layout.RefRatio())    
{
    define_alloc();
}

MFTower::MFTower(Layout&           _layout,
                 PArray<MultiFab>& pamf)
    : layout(_layout), nComp(pamf[0].nComp()), nGrow(pamf[0].nGrow()), nLevs(layout.NumLevels()),
    gridArray(layout.GridArray()), geomArray(layout.GeomArray()), refRatio(layout.RefRatio())    
{
    define_noalloc(pamf);
}

static const IndexType CC(IntVect::TheZeroVector());
static const IndexType NC(IntVect::TheUnitVector());
static const IndexType ECI(BoxLib::BASISV(0));
static const IndexType ECJ(BoxLib::BASISV(1));
#if BL_SPACEDIM==3
static const IndexType ECK(BoxLib::BASISV(2));
#endif

void
MFTower::define_alloc()
{
    mft.resize(nLevs,PArrayManage);
    for (int lev=0; lev<nLevs; ++lev)
    {
        BoxArray ba = gridArray[lev];
        if (iType!=CC) {
            if (iType == ECI) {
                ba.surroundingNodes(0);
            }
            else if (iType == ECJ) {
                ba.surroundingNodes(1);
            }
#if BL_SPACEDIM==3
            else if (iType == ECK) {
                ba.surroundingNodes(2);
            }
#endif
            else if (iType == NC) {
                ba.surroundingNodes();
            }
        }
        mft.set(lev,new MultiFab(ba,nComp,nGrow));
    }
}

void
MFTower::define_noalloc(PArray<MultiFab>& pamf)
{
    iType = pamf[0].boxArray()[0].ixType();
    mft.resize(nLevs,PArrayNoManage);
    for (int lev=0; lev<nLevs; ++lev)
    {
        BL_ASSERT(pamf[lev].boxArray() == gridArray[lev]);
        BL_ASSERT(pamf[lev].nComp() == nComp);
        BL_ASSERT(pamf[lev].nGrow() == nGrow);
        mft.set(lev,&(pamf[lev]));
    }
}

void
MFTower::AXPY(const MFTower& rhs,
              Real           p)
{
    BL_ASSERT(IsCompatible(rhs));
    FArrayBox fab;
    for (int lev=0; lev<nLevs; ++lev)
    {
        for (MFIter mfi(mft[lev]); mfi.isValid(); ++mfi)
        {
            Box vbox = mfi.validbox();
            fab.resize(vbox,1);
            fab.copy(rhs[lev][mfi]);
            if (lev<nLevs-1) {
                BoxArray cba = BoxArray(gridArray[lev-1]).refine(refRatio[lev]);
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
                BoxArray cfba = BoxArray(gridArray[lev-1]).coarsen(refRatio[lev-1]);
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
    bool isok = nLevs==rhs.NumLevels();
    for (int lev=0; lev<nLevs && isok; ++lev)
    {
        isok &= gridArray[lev] == rhs.gridArray[lev];
        if (lev < nLevs-1) {
            isok &= rhs.RefRatio()[lev]==refRatio[lev];
        }
    }
    return isok;
}


#ifdef BL_USE_PETSC
std::ostream& operator<< (std::ostream&  os, const Stencil& a)
{
  os << "(size=" << a.size() << ") ";
  for (Stencil::const_iterator it=a.begin(), End=a.end(); it!=End; ++it) {
    os << "(" << it->first << "):(" << it->second << ") ";
  }
  return os;
}


const Stencil operator*(Real val, const Stencil& lhs)
{
  Stencil ret(lhs);
  ret *= val;
  return ret;
}

const Stencil operator*(const Stencil& lhs, Real val)
{
  Stencil ret(lhs);
  ret *= val;
  return ret;
}

const Stencil operator+(const Stencil& lhs, const Stencil& rhs)
{
  Stencil ret(lhs);
  ret += rhs;
  return ret;
}

const Stencil operator-(const Stencil& lhs, const Stencil& rhs)
{
  Stencil ret(lhs);
  ret -= rhs;
  return ret;
}

Stencil&
Stencil::operator+=(const Stencil& rhs)
{
  for (const_iterator it=rhs.begin(), End=rhs.end(); it!=End; ++it) {
    const_iterator mit=find(it->first);
    if (mit==end()) {
      insert(std::pair<Node,Real>(it->first,it->second));
    } else {
      (*this)[it->first] += it->second;
    }
  }
  return *this;
}

Stencil&
Stencil::operator-=(const Stencil& rhs)
{
  for (const_iterator it=rhs.begin(), End=rhs.end(); it!=End; ++it) {
    const_iterator mit=find(it->first);
    if (mit==end()) {
      insert(std::pair<Node,Real>(it->first,-it->second));
    } else {
      (*this)[it->first] -= it->second;
    }
  }
  return *this;
}

Stencil&
Stencil::operator*=(Real val)
{
  for (const_iterator it=begin(), End=end(); it!=End; ++it) {
    (*this)[it->first] *= val;
  }
  return *this;
}

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
    crseNodes.clear();
    nodeIds.clear();
    crseIds.clear();
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

    NodeFab fnodeFab;

    Clear();
    nodes.resize(nLevs,PArrayManage);
    nodeIds.resize(nLevs,PArrayManage);
    bndryCells.resize(nLevs);
    crseNodes.resize(nLevs,PArrayManage);

    int cnt = 0; // local node counter
    for (int lev=0; lev<nLevs; ++lev)
    {
        nodes.set(lev,new MultiNodeFab(gridArray[lev],1,nGrow));
        MultiNodeFab& nodesLev = nodes[lev];
        nodeIds.set(lev,new MultiIntFab(gridArray[lev],1,nGrow));


        for (MFIter fai(nodesLev); fai.isValid(); ++fai)
        {
            NodeFab& ifab = nodesLev[fai];
            const Box box = ifab.box() & gridArray[lev][fai.index()];
            for (IntVect iv=box.smallEnd(); iv<=box.bigEnd(); box.next(iv))
                ifab(iv,0) = Node(iv,lev,Node::VALID);
        }
            
        if (lev > 0)
        {
            BoxArray cba = BoxArray(gridArray[lev]).coarsen(refRatio[lev-1]);
            BoxArray cgba = BoxArray(cba).grow(1);
            crseNodes.set(lev,new MultiNodeFab(cgba,1,0));
            for (MFIter mfi(crseNodes[lev]); mfi.isValid(); ++mfi)
            {
                NodeFab& nfab=crseNodes[lev][mfi];
                const Box box = nfab.box() & geomArray[lev-1].Domain();
                for (IntVect iv=box.smallEnd(); iv<=box.bigEnd(); box.next(iv)) {
                    nfab(iv,0) = Node(iv,lev-1,Node::VALID);
                }
             }

            for (MFIter mfi(crseNodes[lev]); mfi.isValid(); ++mfi)
            {
                NodeFab& nfab=crseNodes[lev][mfi];
                const Box box = nfab.box() & geomArray[lev-1].Domain();
                std::vector< std::pair<int,Box> > isects = cba.intersections(box);
                for (int i=0; i<isects.size(); ++i) {
                    const Box& ibox=isects[i].second;
                    for (IntVect iv=ibox.smallEnd(), End=ibox.bigEnd(); iv<=End; ibox.next(iv)) {
                        nfab(iv,0).type = Node::COVERED;
                    }
                }
            }

            const Box rangeBox = Box(IntVect::TheZeroVector(),
                                     refRatio[lev-1] - IntVect::TheUnitVector());
            bndryCells[lev] = BoxLib::intersect( GetBndryCells(nodesLev.boxArray(),refRatio[lev-1],geomArray[lev]),
                                                 geomArray[lev].Domain() );
            BoxArray bndC = BoxArray(bndryCells[lev]).coarsen(refRatio[lev-1]);

            for (MFIter mfi(nodesLev); mfi.isValid(); ++mfi)
            {
                NodeFab& nodeFab = nodesLev[mfi];
                const Box& box = nodeFab.box();
                Box gbox = Box(mfi.validbox()).grow(refRatio[lev-1]);
                std::vector< std::pair<int,Box> > isects = bndryCells[lev].intersections(gbox);                
                for (int i=0; i<isects.size(); ++i) {
                    const Box& fbox = isects[i].second;
                    Box ovlp = fbox & box;
                    if (ovlp.ok()) {
                        fnodeFab.resize(fbox,1);
                        Box cbox = Box(fbox).coarsen(refRatio[lev-1]);
                        for (IntVect civ = cbox.smallEnd(), End=cbox.bigEnd(); civ<=End; cbox.next(civ)) {
                            const IntVect baseIV = refRatio[lev-1] * civ;
                            for (IntVect ivt = rangeBox.smallEnd(), End=rangeBox.bigEnd(); ivt<=End;rangeBox.next(ivt)) {
                                fnodeFab(baseIV + ivt,0) = crseNodes[lev][mfi](civ,0);
                            }
                        }
                        nodeFab.copy(fnodeFab);
                    }
                }
            }
        }

        // Block out cells covered by finer grid
        if (lev < nLevs-1)
        {
            const BoxArray coarsenedFineBoxes =
                BoxArray(gridArray[lev+1]).coarsen(refRatio[lev]);
                
            for (MFIter fai(nodesLev); fai.isValid(); ++fai)
            {
                NodeFab& ifab = nodesLev[fai];
                const Box& box = ifab.box();
                std::vector< std::pair<int,Box> > isects = coarsenedFineBoxes.intersections(box);
                for (int i = 0; i < isects.size(); i++)
                {
                    const Box& ovlp = isects[i].second;
                    for (IntVect iv=ovlp.smallEnd(); iv<=ovlp.bigEnd(); ovlp.next(iv))
                        ifab(iv,0) = Node(iv,lev,Node::COVERED);
                }
            }
        }

        // Set nodeIds
        nodeIds[lev].setVal(-1);
        for (MFIter fai(nodesLev); fai.isValid(); ++fai)
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
#endif

    // Now communicate node numbering to neighbors grow cells
    crseIds.resize(nLevs-1,PArrayManage);
    for (int lev=0; lev<nLevs; ++lev)
    {
        if (lev>0) 
        {
            BoxArray bndC = BoxArray(bndryCells[lev]).coarsen(refRatio[lev-1]);
        
            crseIds.set(lev-1,new MultiIntFab(bndC,1,0,Fab_allocate)); crseIds[lev-1].setVal(-1);
            crseIds[lev-1].copy(nodeIds[lev-1]); // parallel copy

            // "refine" crseIds
            MultiIntFab fineIds(bndryCells[lev],1,0); fineIds.setVal(-1);
            const Box rangeBox = Box(IntVect::TheZeroVector(),
                                     refRatio[lev-1] - IntVect::TheUnitVector());
            for (MFIter mfi(crseIds[lev-1]); mfi.isValid(); ++mfi)
            {
                const Box& cbox = crseIds[lev-1][mfi].box();
                for (IntVect iv = cbox.smallEnd(), End=cbox.bigEnd(); iv<=End; cbox.next(iv)) {
                    int nodeIdx = crseIds[lev-1][mfi](iv,0);
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
        nodeIds[lev].FillBoundary(0,1);
        BoxLib::FillPeriodicBoundary<IntFab>(geomArray[lev],nodeIds[lev],0,1);
    }

    int n = nNodes_local; // Number of local columns of J
    int m = nNodes_local; // Number of local rows of J
    int N = nNodes_global; // Number of global columns of J 
    int M = nNodes_global; // Number of global rows of J 
    int d_nz = 1 + 2*BL_SPACEDIM; // Number of nonzero local columns of J
    int o_nz = 0; // Number of nonzero nonlocal (off-diagonal) columns of J

    PetscErrorCode ierr = 
        MatCreateMPIAIJ(ParallelDescriptor::Communicator(), m, n, M, N, d_nz, PETSC_NULL, o_nz, PETSC_NULL, &J_mat);
    CHKPETSC(ierr);
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
    bool isok = mft.NumLevels()==nLevs;
    if (isio) {
        for (int lev=0; lev<nLevs; ++lev)
        {
            isok &= mft[lev].boxArray()==gridArray[lev];
        }
    }
    return isok;
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
    //ierr = VecCreateMPI(ParallelDescriptor::Communicator(),nNodes_local,nNodes_global,&V); CHKPETSC(ierr);

    //std::cout << "created PETSc vec: " << &V << std::endl;

    //vecs_I_created.push_back(&V);

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
Layout::VecToMFTower(MFTower& mft,
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

std::ostream& operator<< (std::ostream&  os, const Layout::NodeFab& nfab)
{
    os << "NodeFab" << '\n';
    const Box& box = nfab.box();
    os << "Box: " << box << '\n';
    for (IntVect iv(box.smallEnd()); iv<=box.bigEnd(); box.next(iv))
    {
        os << iv << " ";
        for (int n=0; n<nfab.nComp(); ++n) {
            os << nfab(iv,n) << " ";
        }
        os << '\n';
    } 
    return os;
}

std::ostream& operator<< (std::ostream&  os, const Layout::MultiNodeFab& mnf)
{
    os << "MultiNodeFab" << '\n';
    const BoxArray& ba = mnf.boxArray();
    os << "BoxArray: " << ba << '\n';
    const DistributionMapping& dm = mnf.DistributionMap();
    for (int i=0; i<mnf.size(); ++i) {
        if (dm[i]==ParallelDescriptor::MyProc()) {
            os << mnf[i] << std::endl;
        }
        ParallelDescriptor::Barrier();
    }
    return os;
}

#endif

ABecTower::ABecTower(const Layout& _layout)
    : layout(_layout), gridArray(_layout.GridArray()), 
      geomArray(_layout.GeomArray()), nLevs(_layout.NumLevels()),
      coefs(BL_SPACEDIM,PArrayManage), refRatio(_layout.RefRatio())
{
    Array<IndexType> itype(BL_SPACEDIM);
    D_DECL(itype[0]=ECI,
           itype[1]=ECJ,
           itype[2]=ECK);
    for (int d=0; d<BL_SPACEDIM; ++d) {
        coefs.set(d, new MFTower(layout,itype[d]));
    }
}

/*
     BuildInterpCoefs:
     
     This routine returns the Lagrange interpolating coefficients for a
     polynomial through N points, evaluated at xInt=-1 (see Numerical Recipes,
     v2, p102, e.g.):
     
             (x-x2)(x-x3)...(x-xN)              (x-x1)(x-x2)...(x-x(N-1))
     P(x) = ----------------------- y1  + ... + ------------------------  yN
            (x1-x2)(x1-x3)...(x1-xN)            (x1-x2)(x1-x3)...(x1-xN)

      P(xInt) = sum_(i=1)^(N) y[i]*c[i]
*/

void
BuildInterpCoefs(Real xVal, Array<Real>& coefs)
{
    int N = coefs.size();
    Array<Real> x(N);
    x[0] = -0.5-xVal; // xVal is location of off-grid Dirichlet value
    for (int i=1; i<N; ++i) {
        x[i] = i-1;
    }
    for (int j=0; j<N; ++j) {
        Real num = 1;
        Real den = 1;
        for (int i=0; i<N; ++i) {
            if (i!=j) {
                num *= -1 - x[i];
                den *= x[j] - x[i];
            }
        }
        coefs[j] = num/den;
    }
}

void
ABecTower::BuildCFParallelInterpStencil()
{
    // Some handy intvects
    Array<IntVect> ivp(BL_SPACEDIM), ivpp(BL_SPACEDIM), ivm(BL_SPACEDIM), ivmm(BL_SPACEDIM);
    for (int d=0; d<BL_SPACEDIM; ++d) {
        ivp[d] = BoxLib::BASISV(d);
        ivpp[d] = ivp[d] + BoxLib::BASISV(d);
        ivm[d] = -BoxLib::BASISV(d);
        ivmm[d] = ivm[d] - BoxLib::BASISV(d);
    }

    const PArray<Layout::MultiNodeFab>& nodes = layout.Nodes();
    const PArray<Layout::MultiNodeFab>& crseNodes = layout.CrseNodes();
    const Array<BoxArray>& bndryCells = layout.BndryCells();
    parallelInterpStencil.resize(nLevs);

    for (int lev=1; lev<nLevs; ++lev)
    {
        parallelInterpStencil[lev].resize(BL_SPACEDIM);
        for (MFIter mfi(nodes[lev]); mfi.isValid(); ++mfi)
        {
            const Layout::NodeFab& nodeFab = nodes[lev][mfi];
            const Box& cgbox = crseNodes[lev][mfi].box();
            for (int d=0; d<BL_SPACEDIM; ++d) {

                Array<int> dtan;
                for (int d0=0; d0<BL_SPACEDIM; ++d0) {
                    if (d!=d0) {
                        dtan.push_back(d0);
                    }
                }
                
                Box gdbox = Box(mfi.validbox()).grow(d,1) & geomArray[lev].Domain();
                std::vector< std::pair<int,Box> > isects = bndryCells[lev].intersections(gdbox);
                for (int i=0; i<isects.size(); ++i) {
                    const Box& bndrySect = isects[i].second;
                    for (IntVect iv=bndrySect.smallEnd(), End=bndrySect.bigEnd(); iv<=End; bndrySect.next(iv)) {
                        
                        Node nC = nodeFab(iv,0);
                        BL_ASSERT(cgbox.contains(nC.iv) && nC.type==Node::VALID);
                        std::map<int,Real> x;
                        
                        Stencil c; c[nC] = 1;
                        parallelInterpStencil[lev][d][iv] = c;
                        
                        for (int d0=0; d0<dtan.size(); ++d0) {
                            int itan = dtan[d0];
                            Stencil der, der2;
                            int r = refRatio[lev-1][itan];

                            x[itan] = (iv[itan]%r - 0.5*(r-1))/r;
                            
                            const Node& nR = crseNodes[lev][mfi](nC.iv+ivp[itan],0);
                            const Node& nL = crseNodes[lev][mfi](nC.iv+ivm[itan],0);
                            
                            bool Rvalid = nR.type==Node::VALID;
                            bool Lvalid = nL.type==Node::VALID;
                            
                            if (Rvalid && Lvalid) {
                                // Centered full
                                der[nL]  = -0.5;  der[nR] = +0.5;
                                der2[nL] = +1.0; der2[nC] = -2.0; der2[nR] = +1.0;
                            }
                            else if (Rvalid) {
                                const Node& nRR = crseNodes[lev][mfi](nC.iv+ivpp[itan],0);
                                bool RRvalid = nRR.type==Node::VALID;
                                if (RRvalid) {
                                    // R-shifted full
                                    der[nC]  = -0.5;  der[nRR] = +0.5;
                                    der2[nC] = +1.0; der2[nR]  = -2.0; der2[nRR] = +1.0;
                                } else {
                                    // R-shifted linear
                                    der[nC] = -1.0; der[nR] = +1.0;
                                }
                            } else if (Lvalid) {
                                const Node& nLL = crseNodes[lev][mfi](nC.iv+ivmm[itan],0);
                                bool LLvalid = nLL.type==Node::VALID;
                                if (LLvalid) {
                                    // L-shifted full
                                    der[nLL]  = -0.5;  der[nC] = +0.5;
                                    der2[nLL] = +1.0; der2[nL] = -2.0; der2[nC] = +1.0;
                                } else {
                                    // L-shifted linear
                                    der[nL] = -1.0; der[nC] = +1.0;
                                }
                            } else {
                                // piecewise constant (no derivatives)
                            }
                            parallelInterpStencil[lev][d][iv] += x[itan]*der + 0.5*x[itan]*x[itan]*der2;
                        } // tangential direction

                        if (dtan.size()==2) {
                            const Node& npp = crseNodes[lev][mfi](nC.iv+ivp[dtan[0]]+ivp[dtan[1]],0);
                            const Node& nmm = crseNodes[lev][mfi](nC.iv+ivm[dtan[0]]+ivm[dtan[1]],0);
                            const Node& npm = crseNodes[lev][mfi](nC.iv+ivp[dtan[0]]+ivm[dtan[1]],0);
                            const Node& nmp = crseNodes[lev][mfi](nC.iv+ivm[dtan[0]]+ivp[dtan[1]],0);
                        
                            bool PPvalid = npp.type==Node::VALID;
                            bool MMvalid = nmm.type==Node::VALID;
                            bool PMvalid = npm.type==Node::VALID;
                            bool MPvalid = nmp.type==Node::VALID;
                        
                            if (PPvalid && MMvalid && PMvalid && MPvalid) {
                                Stencil crossterm;
                                crossterm[npp] = +0.25;
                                crossterm[nmm] = +0.25;
                                crossterm[npm] = -0.25;
                                crossterm[nmp] = -0.25;
                                parallelInterpStencil[lev][d][iv] += x[dtan[0]]*x[dtan[1]]*crossterm;
                            }
                        }
                    } // bndry iv
                } // bndry box
            } // bc direction
        } // fine box
    } // lev
}

void
ABecTower::BuildStencil(const BCRec& bc,
                        int maxorder)
{
    int myproc = ParallelDescriptor::MyProc();

    BuildCFParallelInterpStencil();
    const PArray<Layout::MultiNodeFab>& nodes = layout.Nodes();

    // Precompute an often-used interp stencil
    Array<Real> iCoefsZero(maxorder); BuildInterpCoefs(0,iCoefsZero); // value at wall
    Array<Array<Real> > iCoefsCF(BL_SPACEDIM, Array<Real>(maxorder));

    perpInterpStencil.resize(nLevs);
    for (int lev=0; lev<nLevs; ++lev) {
        const Box& dbox = geomArray[lev].Domain();
        Array<Array<Box> > bndry(2, Array<Box>(BL_SPACEDIM));
        for (int d=0; d<BL_SPACEDIM; ++d) {
            bndry[0][d] = BoxLib::adjCellLo(dbox,d);
            bndry[1][d] = BoxLib::adjCellHi(dbox,d);
        }

        // Precompute other often-used interp stencils
        if (lev>0) {
            for (int d=0; d<BL_SPACEDIM; ++d) {
                const IntVect& rat = refRatio[lev-1];
                BuildInterpCoefs(0.5*rat[d],iCoefsCF[d]); // value at wall
            }
        }

        const BoxArray& ba = gridArray[lev];
        const Layout::MultiNodeFab& fmn = nodes[lev];
        const Geometry& gl = geomArray[lev];
        for (MFIter mfi(fmn); mfi.isValid(); ++mfi) {
            const Box& vbox = mfi.validbox();
            const Layout::NodeFab& fn = fmn[mfi];
            for (int d=0; d<BL_SPACEDIM; ++d) {
                Array<Box> myBndry(2);
                myBndry[0] = BoxLib::adjCellLo(vbox,d);
                myBndry[1] = BoxLib::adjCellHi(vbox,d);
                
                for (int hilo=0; hilo<2; ++hilo)
                {
                    Box povlp = myBndry[hilo] & bndry[hilo][d];
                    int bc_flag = (hilo==0 ? bc.lo()[d] : bc.hi()[d]);
                    bool pbc = povlp.ok() && bc_flag==EXT_DIR;
                    if (pbc) {
                        for (IntVect iv=povlp.smallEnd(), End=povlp.bigEnd(); iv<=End; povlp.next(iv)) {
                            Stencil perp;
                            int sgn = (hilo==0 ? +1  : -1); // Direction of interp stencil (inward)
                            for (int k=0; k<iCoefsZero.size(); ++k) {
                                IntVect siv = iv + sgn*k*BoxLib::BASISV(d);
                                perp[fn(siv,0)] = iCoefsZero[k];
                            }
                            std::cout << "For " << iv << " stencil: " << perp << std::endl; 
                        }
                    }
                    else if (lev>0) {
                        // Build c-f stencil
                        BoxArray sba = BoxLib::complementIn(myBndry[hilo],ba);
                        if (gl.isPeriodic(d)) {
                            BoxArray per_ba = BoxLib::intersect(sba,BoxArray(ba).shift(d,dbox.length(d)));
                            for (int j=0; j<per_ba.size(); ++j) {
                                sba = BoxLib::complementIn(per_ba[j],sba);
                            }
                            per_ba = BoxLib::intersect(sba,BoxArray(ba).shift(d,-dbox.length(d)));
                            for (int j=0; j<per_ba.size(); ++j) {
                                sba = BoxLib::complementIn(per_ba[j],sba);
                            }
                        } else {
                            sba = BoxLib::intersect(sba,dbox);
                        }
                        
                        // Now, with coefs create stencil entries
                        const Layout::IVSMap& parStencil = parallelInterpStencil[lev][d];
                        for (int j=0; j<sba.size(); ++j) {
                            const Box& sbox = sba[j];
                            for (IntVect iv=sbox.smallEnd(), End=sbox.bigEnd(); iv<=End; sbox.next(iv)) {
                                
                                Layout::IVScit it = parStencil.find(iv);
                                BL_ASSERT(it!=parStencil.end());
                                const Stencil& parallelStencil = it->second;
                                
                                Stencil totalStencil = iCoefsCF[d][0]*parallelStencil;
                                
                                int sgn = (hilo==0 ? +1  : -1); // Direction of interp stencil (inward)
                                for (int k=1; k<iCoefsCF[d].size(); ++k) {
                                    IntVect siv = iv + sgn*k*BoxLib::BASISV(d);
                                    totalStencil[fn(siv,0)] = iCoefsCF[d][k];
                                }
                                std::cout << "For " << iv << " stencil: " << totalStencil << std::endl; 
                            }
                        }
                    }
                }
            }
        }
    }
}

void
ABecTower::DoCoarseFineParallelInterp(MFTower& mft,
                                      int      sComp,
                                      int      nComp) const
{
    BL_ASSERT(layout.IsCompatible(mft));
    const Array<BoxArray>& bndryCells = layout.BndryCells();

    for (int lev=1; lev<nLevs; ++lev) {

        MultiFab& mf = mft[lev];
        BL_ASSERT(mf.nGrow()>=1);
        BL_ASSERT(sComp+nComp<=mf.nComp());

        BoxArray bacgd = BoxArray(mf.boxArray()).coarsen(refRatio[lev-1]).grow(1);
        MultiFab crseMF(bacgd,nComp,0);
        crseMF.copy(mft[lev-1],sComp,0,nComp); // parallel copy

        BoxArray bnd = bndryCells[lev];
        const Geometry& gl = geomArray[lev];
        const Array<Layout::IVSMap>& parInterpLev = parallelInterpStencil[lev];

        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {

            int gridIdx = mfi.index();

            FArrayBox& crseFab = crseMF[mfi];
            FArrayBox& fineFab = mf[mfi];

            for (int d=0; d<BL_SPACEDIM; ++d) {

                const Layout::IVSMap& parInterp = parInterpLev[d];
                Box boxgd = Box(mfi.validbox()).grow(d,1) & gl.Domain();
                std::vector< std::pair<int,Box> > isects = bnd.intersections(boxgd);
                for (int i=0; i<isects.size(); ++i) {
                    const Box& bndrySect = isects[i].second;
                    for (IntVect iv=bndrySect.smallEnd(), End=bndrySect.bigEnd(); iv<=End; bndrySect.next(iv)) {
                        
                        
                        Layout::IVScit it=parInterp.find(iv);
                        if (it!=parInterp.end()) {
                            const Stencil& s = it->second;

                            for (int n=0; n<nComp; ++n) {
                                Real res = 0;
                                for (Stencil::const_iterator it=s.begin(), End=s.end(); it!=End; ++it) {
                                    const IntVect& ivs=(it->first).iv;
                                    BL_ASSERT(crseFab.box().contains(ivs));
                                    BL_ASSERT((it->first).level==lev-1);
                                    res += crseFab(ivs,sComp+n) * it->second;
                                }
                                fineFab(iv,sComp+n) = res;
                            }
                        }
                    }
                }
            }
        }

        mf.FillBoundary(sComp,nComp);
        gl.FillPeriodicBoundary(mf,sComp,nComp);
    }
}

