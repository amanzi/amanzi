#include <winstd.H>
#include <iostream>
#include <fstream>
#include <unistd.h>

#include <GSLibInt.H>
#include <GSLibInt_F.H>
#include <VisMF.H>
#include <Utility.H>

void  
GSLibInt::parRand(Array <Real> kappaval, 
                  Real         dkappa,
                  Array<int>   n_cell, 
                  int          twoexp, 
                  MultiFab&    mfdata)
{
  int nkpval = kappaval.size();
  int domlo[BL_SPACEDIM];
  int domhi[BL_SPACEDIM];
  
  for (int dir = 0; dir < BL_SPACEDIM; dir++) { 
    domlo[dir] = -3*twoexp;
    domhi[dir] = (n_cell[dir]+3)*twoexp;
  }

  for (MFIter mfi(mfdata); mfi.isValid(); ++mfi) {
    const int  i       = mfi.index();
    const Real* kp_dat = mfdata[mfi].dataPtr();
    const int*  kp_lo  = mfdata[mfi].loVect();
    const int*  kp_hi  = mfdata[mfi].hiVect();

    int iuc = rand()%1000000 + i;

    FORT_PHIRAND(kp_dat,ARLIM(kp_lo),ARLIM(kp_hi),
		 kappaval.dataPtr(),&nkpval,&dkappa, 
		 domlo,domhi,&iuc);
  }
}


void
GSLibInt::seqGaussianSim(const Array <Real>& kappaval, 
                         const Array<int>&   n_cell, 
                         const Array<Real>&  problo,
                         const Array<Real>&  probhi,
                         int                 twoexp,
                         MultiFab&           mfdata,
                         const std::string&  gsfile)
{
  int strl = gsfile.length();

  // Create grids at finest level.
  int nx[3];
  Real dx[3],hdx[3];

  for (int i=0;i<BL_SPACEDIM; i++) {
    nx[i] = (n_cell[i]+6)*twoexp;
    dx[i] = (probhi[i]-problo[i])/(n_cell[i]*twoexp);
    hdx[i] = dx[i]*0.5;
  }

#if (BL_SPACEDIM == 2)
  nx[2] = 1;
  dx[2] = dx[1];
  hdx[2] = hdx[1];
#endif

  int rand_seed;
  FORT_INIT_GSLIB(gsfile.c_str(),&strl,nx,dx,hdx,&rand_seed);
  BoxLib::InitRandom(rand_seed);

  std::cout << "Doing sequential gaussian simulation  ... \n";
  
  Box bx;
#if (BL_SPACEDIM == 3)
  bx = Box(IntVect(0,0,0),
	   IntVect(n_cell[0]-1,n_cell[1]-1,n_cell[2]-1));
#else
  bx = Box(IntVect(0,0),
	   IntVect(n_cell[0]-1,n_cell[1]-1));
#endif
  bx.grow(3);
  bx.refine(twoexp);
  
  BoxArray ba(bx);
  MultiFab mf(ba,1,0);

  if (ParallelDescriptor::IOProcessor())
    std::cout << "starting serial sequential gaussian simulation.\n";

  if (ParallelDescriptor::MyProc() == mf.DistributionMap()[0]) {
    FORT_SGSIM(mf[0].dataPtr(),
               ARLIM(mf[0].loVect()),ARLIM(mf[0].hiVect()),dx);
  }

  ParallelDescriptor::Barrier();
  mfdata.copy(mf);
}

void
GSLibInt::rdpGaussianSim(const Array<Real>& kappaval, 
                         const Array<int>&  n_cell, 
                         const Array<Real>& problo,
                         const Array<Real>& probhi,
                         int                twoexp,
                         MultiFab&          mfdata,
                         int                crse_init_factor,
                         int                max_grid_size_fine_gen,
                         int                ngrow_fine_gen,
                         const std::string& gsfile)
{
  if (ParallelDescriptor::IOProcessor())
    std::cout << "Doing random parallel gaussian simulation.\n";

  std::string kcfile;

  cndGaussianSim(kappaval,n_cell, problo, probhi, twoexp, mfdata,
		 crse_init_factor,max_grid_size_fine_gen,ngrow_fine_gen,kcfile,gsfile);
}

void
GSLibInt::cndGaussianSim(const Array<Real>& kappaval, 
                         const Array<int>&  n_cell, 
                         const Array<Real>& problo,
                         const Array<Real>& probhi,
                         int                twoexp,
                         MultiFab&          mfdata,
                         int                crse_init_factor,
                         int                max_grid_size_fine_gen,
                         int                ngrow_fine_gen,
                         std::string&       kcfile,
                         const std::string& gsfile)
{
  int  nx[BL_SPACEDIM],nxyz;
  Real dx[BL_SPACEDIM],dxc[BL_SPACEDIM];


  int strl = gsfile.length();

  int cond_option = 0;
  int c_sz        = 0;
  int real_sz     = 0;
  int int_sz      = 0;
  int c_idx[10], real_idx[20],int_idx[20];

  //
  // Initializing gslib
  //
  int rand_seed;
  FORT_INIT_GSLIB2(gsfile.c_str(),&strl,&c_sz,c_idx,
		   &real_sz,real_idx,&int_sz,int_idx, &cond_option, &rand_seed);

  BoxLib::InitRandom(rand_seed);
 
 if (cond_option == 1 && c_sz == 0) {
    std::cout << "GSLIB data is missing.  Doing unconditioned simulation.\n";
    cond_option = 0;
  }

  int nGrow = mfdata.nGrow();
  nxyz = 1;
  for (int i=0;i<BL_SPACEDIM; i++) {
    dxc[i] = (probhi[i]-problo[i])/n_cell[i];
    nx[i]  = (n_cell[i]+2*nGrow)*twoexp;
    dx[i]  = dxc[i]/twoexp;
    nxyz  *= nx[i];
  }

  Box bx(IntVect(D_DECL(          0,          0,          0)),
         IntVect(D_DECL(n_cell[0]-1,n_cell[1]-1,n_cell[2]-1)));
  bx.grow(nGrow);

  MultiFab mfc;
  int      cfac = crse_init_factor;
  for (int dir=0; dir<BL_SPACEDIM; dir++) {
    cfac = std::min(cfac,n_cell[dir]);
  }

  Array<int> domloc(BL_SPACEDIM),domhic(BL_SPACEDIM);
  Array <Real> scratch_c(1);
  if (cond_option > 0) { // originally supported other values of cond_option
    scratch_c.resize(c_sz,1.e20);

    if (cond_option == 1) {
      const IntVect ivDum(D_DECL(0,0,0));
      const Real* dDum = scratch_c.dataPtr();
	
      FORT_INTERNAL_DATA(dDum,ARLIM(ivDum),ARLIM(ivDum),
			 scratch_c.dataPtr(),&c_sz,c_idx,
			 &kappaval[0],dxc,problo.dataPtr(),
			 domloc.dataPtr(),domhic.dataPtr());

    }

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMin(scratch_c.dataPtr(),c_sz,IOProc);
    ParallelDescriptor::Bcast(scratch_c.dataPtr(),c_sz,IOProc);
  }

  bx.refine(twoexp);
  BoxArray ba(bx); 
  ba.maxSize(max_grid_size_fine_gen);
  MultiFab mf(ba,1,ngrow_fine_gen);
  mf.setVal(0.);
  
  Array< Array<Real> > scratch_r(mf.size());
  Array< Array<int>  > scratch_i(mf.size());
  Array< Array<int>  > order(mf.size());

  int max_fab_size = 0;
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {    
    const int  i     = mfi.index();
    const int* lo    = mfi.validbox().loVect();
    const int* hi    = mfi.validbox().hiVect();
	  
    const int* k_lo  = mf[mfi].loVect();
    const int* k_hi  = mf[mfi].hiVect();
    const Real* kdat = mf[mfi].dataPtr();
    int nvalid = 1;
    int ntotal = 1;
    for (int d=0; d<BL_SPACEDIM; d++) {
      ntotal *= k_hi[d]-k_lo[d]+1;
      nvalid *= hi[d]-lo[d]+1;
    }
    real_sz = real_sz + ntotal;
    scratch_i[i].resize(int_sz);
    scratch_r[i].resize(real_sz);
    order[i].resize(nvalid);

    srand(int(time(0)));
    int iuc = rand()%1000000 + i;

    max_fab_size = std::max(max_fab_size,nvalid);

    FORT_SGSIM_SETUP(kdat,ARLIM(k_lo),ARLIM(k_hi),
    		     order[i].dataPtr(),&nvalid,
		     scratch_c.dataPtr(),&c_sz,c_idx,
		     scratch_r[i].dataPtr(),&real_sz,real_idx,
		     scratch_i[i].dataPtr(),&int_sz,int_idx,    
		     lo,hi,dx,problo.dataPtr(),&rand_seed);
  }
  ParallelDescriptor::ReduceIntMax(max_fab_size);

  int it = 0;
  while (it < max_fab_size) {

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {    
      const int  i     = mfi.index();
	  
      const int* k_lo  = mf[mfi].loVect();
      const int* k_hi  = mf[mfi].hiVect();
      const Real* kdat = mf[mfi].dataPtr();
      if (it < order[i].size()) 
      {
	int idx_chosen = order[i][it];

	FORT_SGSIM_ITER(kdat,ARLIM(k_lo),ARLIM(k_hi),
			scratch_c.dataPtr(),&c_sz,c_idx,
			scratch_r[i].dataPtr(),&real_sz,real_idx,
			scratch_i[i].dataPtr(),&int_sz,int_idx,
			&idx_chosen);
      }
    }

    mf.FillBoundary();

    it = it + 1;
  }

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const int  i     = mfi.index();
    const int* k_lo  = mf[mfi].loVect();
    const int* k_hi  = mf[mfi].hiVect();
    const Real* kdat = mf[mfi].dataPtr();

    FORT_SGSIM_POST(kdat,ARLIM(k_lo),ARLIM(k_hi),
		    scratch_c.dataPtr(),&c_sz,c_idx,
		    scratch_r[i].dataPtr(),&real_sz,real_idx, 
		    scratch_i[i].dataPtr(),&int_sz,int_idx);

    FORT_LGNORM(kdat,ARLIM(k_lo),ARLIM(k_hi),&kappaval[0]);
  }

  ParallelDescriptor::Barrier();

  FORT_SGSIM_DEALLOC();

  BoxArray gba = BoxArray(mf.boxArray()).grow(nGrow);
  MultiFab mfg(gba,1,0);
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    mfg[mfi].copy(mf[mfi]);
  }
  mf.clear();

  BoxArray gba1 = BoxArray(mfdata.boxArray()).grow(nGrow);
  MultiFab mf1(gba1,1,0);
  mf1.copy(mfg); // No-grow to no-grow parallel copy

  int nComp = mfdata.nComp(); // For now, all components get same data
  for (MFIter mfi(mf1); mfi.isValid(); ++mfi) {
    for (int n=0; n<nComp; ++n) {
      const Box& bx = mf1[mfi].box();
      mfdata[mfi].copy(mf1[mfi],bx,0,bx,n,1);
    }
  }
}
