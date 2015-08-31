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
                         const Box&          domain,
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
    nx[i] = domain.length(i)*twoexp;
    dx[i] = (probhi[i]-problo[i])/nx[i];
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
  
  BoxArray ba(domain);
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
                         const Box&         domain,
                         const Array<Real>& problo,
                         const Array<Real>& probhi,
                         MultiFab&          mfdata,
                         int                max_grid_size_fine_gen,
                         int                ngrow_fine_gen,
                         const std::string& gsfile)
{
  if (ParallelDescriptor::IOProcessor())
    std::cout << "Doing random parallel gaussian simulation.\n";

  std::string kcfile;

  cndGaussianSim(kappaval, domain, problo, probhi, mfdata,
		 max_grid_size_fine_gen,ngrow_fine_gen,kcfile,gsfile);
}

static void
memUsage(const std::string& note)
{
  Real min_alloc_fab_gb = BoxLib::TotalBytesAllocatedInFabs()/(1024.0*1024.0);
  Real max_alloc_fab_gb = min_alloc_fab_gb;
  Real min_fab_gb = BoxLib::TotalBytesAllocatedInFabsHWM()/(1024.0*1024.0);
  Real max_fab_gb = min_fab_gb;

  ParallelDescriptor::ReduceRealMin(min_fab_gb,ParallelDescriptor::IOProcessorNumber());
  ParallelDescriptor::ReduceRealMax(max_fab_gb,ParallelDescriptor::IOProcessorNumber());
  ParallelDescriptor::ReduceRealMin(min_alloc_fab_gb,ParallelDescriptor::IOProcessorNumber());
  ParallelDescriptor::ReduceRealMax(max_alloc_fab_gb,ParallelDescriptor::IOProcessorNumber());
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "\n" << note << " FAB GB spread across MPI nodes: ["
	      << min_fab_gb
	      << " ... "
	      << max_fab_gb
	      << "  alloc: "
	      << min_alloc_fab_gb
	      << " ... "
	      << max_alloc_fab_gb
	      << "]\n";
  }
}

void
GSLibInt::cndGaussianSim(const Array<Real>& kappaval, 
                         const Box&         domain,
                         const Array<Real>& problo,
                         const Array<Real>& probhi,
                         MultiFab&          mfdata,
                         int                max_grid_size_fine_gen,
                         int                ngrow_fine_gen,
                         std::string&       kcfile,
                         const std::string& gsfile)
{
  Real dx[BL_SPACEDIM];
  int strl = gsfile.length();

  int cond_option = 0;
  int c_sz        = 0;
  int real_sz     = 0;
  int int_sz      = 0;
  int c_idx_siz   = 10;
  int r_idx_siz   = 20;
  int i_idx_siz   = 20;
  Array<int> c_idx(c_idx_siz),real_idx(r_idx_siz),int_idx(i_idx_siz);

  //
  // Initializing gslib
  //
  int rand_seed;
  FORT_INIT_GSLIB2(gsfile.c_str(),&strl,&c_sz,c_idx.dataPtr(),&c_idx_siz,
		   &real_sz,real_idx.dataPtr(),&r_idx_siz,&int_sz,int_idx.dataPtr(),&i_idx_siz,
		   &cond_option, &rand_seed);

  BoxLib::InitRandom(rand_seed);
 
  if (cond_option == 1 && c_sz == 0) {
    std::cout << "GSLIB data is missing.  Doing unconditioned simulation.\n";
    cond_option = 0;
  }

  for (int i=0;i<BL_SPACEDIM; i++) {
    dx[i]  = (probhi[i] - problo[i])/domain.length(i);
  }

  Array <Real> scratch_c(1);
  if (cond_option > 0) { // originally supported other values of cond_option
    scratch_c.resize(c_sz,1.e20);
    if (cond_option == 1) {
      FORT_INTERNAL_DATA(scratch_c.dataPtr(),&c_sz,c_idx.dataPtr());
    }
  }

  BoxArray ba(domain);
  ba.maxSize(max_grid_size_fine_gen);
  MultiFab mfg(ba,1,ngrow_fine_gen);
  mfg.setVal(0.);

  Array< Array<Real> > scratch_r(mfg.size());
  Array< Array<int>  > scratch_i(mfg.size());
  Array< Array<int>  > order(mfg.size());

  int max_fab_size = 0;
  for (MFIter mfi(mfg); mfi.isValid(); ++mfi) {    
    const int  i     = mfi.index();
    const int* lo    = mfi.validbox().loVect();
    const int* hi    = mfi.validbox().hiVect();
    const int* dlo   = domain.loVect();
    const int* dhi   = domain.hiVect();
	  
    const int* k_lo  = mfg[mfi].loVect();
    const int* k_hi  = mfg[mfi].hiVect();
    const Real* kdat = mfg[mfi].dataPtr();
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
		     scratch_c.dataPtr(),&c_sz,c_idx.dataPtr(),&c_idx_siz,
		     scratch_r[i].dataPtr(),&real_sz,real_idx.dataPtr(),&r_idx_siz,
		     scratch_i[i].dataPtr(),&int_sz,int_idx.dataPtr(),&i_idx_siz,
		     lo,hi,dx,problo.dataPtr(),dlo,dhi,&rand_seed);
  }
  ParallelDescriptor::ReduceIntMax(max_fab_size);

  for (int i=0; i<order.size(); ++i) {
    for (int j=0; j<order[i].size(); ++j) {
      std::cout << "i, j, o[j]: " << i << ", " << j << ", " << order[i][j] << std::endl;
      if (order[i][j] > 10000) {
	BoxLib::Abort("got it");
      }
    }
  }

  int it = 0;
  while (it < max_fab_size) {

    for (MFIter mfi(mfg); mfi.isValid(); ++mfi) {    
      const int  i     = mfi.index();
	  
      const int* k_lo  = mfg[mfi].loVect();
      const int* k_hi  = mfg[mfi].hiVect();
      const Real* kdat = mfg[mfi].dataPtr();
      if (it < order[i].size()) 
      {
	int idx_chosen = order[i][it];

	//std::cout << "it, idx_chosen: " << it << " " << idx_chosen << " (max = " << mfg[mfi].box().numPts() << ")" << std::endl;
	if (idx_chosen > mfg[mfi].box().numPts()) {
	  std::cout << "idx_chosen: " << idx_chosen << std::endl;
	  BoxLib::Abort("uh oh");
	}

	FORT_SGSIM_ITER(kdat,ARLIM(k_lo),ARLIM(k_hi),
			scratch_c.dataPtr(),&c_sz,c_idx.dataPtr(),&c_idx_siz,
			scratch_r[i].dataPtr(),&real_sz,real_idx.dataPtr(),&r_idx_siz,
			scratch_i[i].dataPtr(),&int_sz,int_idx.dataPtr(),&i_idx_siz,
			&idx_chosen);
      }
    }

    mfg.FillBoundary();

    it = it + 1;
  }

  for (MFIter mfi(mfg); mfi.isValid(); ++mfi) {
    const int  i     = mfi.index();
    const int* k_lo  = mfg[mfi].loVect();
    const int* k_hi  = mfg[mfi].hiVect();
    const Real* kdat = mfg[mfi].dataPtr();

    FORT_SGSIM_POST(kdat,ARLIM(k_lo),ARLIM(k_hi),
		    scratch_c.dataPtr(),&c_sz,c_idx.dataPtr(),&c_idx_siz,
		    scratch_r[i].dataPtr(),&real_sz,real_idx.dataPtr(),&r_idx_siz,
		    scratch_i[i].dataPtr(),&int_sz,int_idx.dataPtr(),&i_idx_siz);
  }

  VisMF::Write(mfg,"JUNK");
  
  for (MFIter mfi(mfg); mfi.isValid(); ++mfi) {
    const int  i     = mfi.index();
    const int* k_lo  = mfg[mfi].loVect();
    const int* k_hi  = mfg[mfi].hiVect();
    const Real* kdat = mfg[mfi].dataPtr();
    const Box& vbox = mfi.validbox();
    FORT_LGNORM(kdat,ARLIM(k_lo),ARLIM(k_hi),&kappaval[0],
		vbox.loVect(),vbox.hiVect());
  }

  mfdata.copy(mfg,0,0,1); // Parallel copy
  mfg.clear();

  FORT_SGSIM_DEALLOC();

  int nComp = mfdata.nComp(); // For now, all components get same data
  for (int n=1; n<nComp; ++n) {
    MultiFab::Copy(mfdata,mfdata,0,n,1,mfdata.nGrow());
  }
}
