#include <winstd.H>
#include "iostream"
#include "fstream"

#include "PGslib.H"
#include "GSLIB_F.H"
#include "POROUS_F.H"
#include "PROB_PM_F.H"
#include <ParallelDescriptor.H>
#include <time.h>

void  
PGslib::seqTurn3D(Array<int> n_cell, int twoexp, MultiFab& mfdata)
{
  srand(int(time(0)));
  int iuc = rand()%1000000;

  // Create grids at finest level.
  int nx,ny,nz;

  nx = (n_cell[0]+6)*twoexp;
  ny = (n_cell[1]+6)*twoexp;
  nz = 1;
#if (BL_SPACEDIM == 3)
  nz = (n_cell[2]+6)*twoexp;
#endif

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
  if (ParallelDescriptor::MyProc() == mf.DistributionMap()[0])
  {
    FORT_TURN(mf[0].dataPtr(),&nx,&ny,&nz,&iuc);
  }

  ParallelDescriptor::Barrier();
  
  mfdata.copy(mf);  // Parallel copy

}


void  
PGslib::parRand(Array <Real> kappaval, 
		Real         dkappa, 
		Array<int>   n_cell, 
		int          twoexp, 
		MultiFab&     mfdata)
{

  int nkpval = kappaval.size();
  int domlo[BL_SPACEDIM];
  int domhi[BL_SPACEDIM];
  
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
  { 
    domlo[dir] = -3*twoexp;
    domhi[dir] = (n_cell[dir]+3)*twoexp;
  }

  for (MFIter mfi(mfdata); mfi.isValid(); ++mfi)
  {
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
PGslib::seqGaussianSim(Array <Real> kappaval, 
		       Real         dkappa, 
		       Array<int>   n_cell, 
		       Array<Real>  problo,
		       Array<Real>  probhi,
		       int          twoexp,
		       MultiFab&    mfdata,
		       std::string  gsfile)
{
  
  int strl = gsfile.length();

  // Create grids at finest level.
  int nx[3];
  Real dx[3],hdx[3];

  for (int i=0;i<BL_SPACEDIM; i++)  
  {
    nx[i] = (n_cell[i]+6)*twoexp;
    dx[i] = (probhi[i]-problo[i])/(n_cell[i]*twoexp);
    hdx[i] = dx[i]*0.5;
  }

#if (BL_SPACEDIM == 2)
  nx[2] = 1;
  dx[2] = dx[1];
  hdx[2] = hdx[1];
#endif

  FORT_INIT_GSLIB(gsfile.c_str(),&strl,nx,dx,hdx);

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

  if (ParallelDescriptor::MyProc() == mf.DistributionMap()[0])
  {
        FORT_SGSIM(mf[0].dataPtr(),
           ARLIM(mf[0].loVect()),ARLIM(mf[0].hiVect()),dx);
        //FORT_LGNORM(mf[0].dataPtr(),
	//  ARLIM(mf[0].loVect()),ARLIM(mf[0].hiVect()),
	// &kappaval[0],&dkappa);
		
  }

  ParallelDescriptor::Barrier();
  
  mfdata.copy(mf);  

}

void  
PGslib::rdpGaussianSim(Array<Real> kappaval, 
		       Real        dkappa,  
		       Array<int>  n_cell, 
		       Array<Real> problo,
		       Array<Real> probhi,
		       int         twoexp,
		       MultiFab&   mfdata,
		       std::string gsfile)
{
  if (ParallelDescriptor::IOProcessor())
    std::cout << "Doing random parallel gaussian simulation.\n";

  std::string kcfile;

  cndGaussianSim(kappaval,dkappa,n_cell, problo, probhi, twoexp, mfdata,
		 kcfile,gsfile);
}

void  
PGslib::cndGaussianSim(Array<Real> kappaval, 
		       Real        dkappa,  
		       Array<int>   n_cell, 
		       Array<Real>  problo,
		       Array<Real>  probhi,
		       int          twoexp,
		       MultiFab&    mfdata,
		       std::string  kcfile,
		       std::string  gsfile)
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
  FORT_INIT_GSLIB2(gsfile.c_str(),&strl,&c_sz,c_idx,
		   &real_sz,real_idx,&int_sz,int_idx, &cond_option);

  if (cond_option == 1 && c_sz == 0) {
    std::cout << "GSLIB data is missing.  Doing unconditioned simulation.\n";
    cond_option = 0;
  }

  int maxgridsize = 32;
  nxyz = 1;
  for (int i=0;i<BL_SPACEDIM; i++)  
  {
    dxc[i] = (probhi[i]-problo[i])/n_cell[i];
    nx[i]  = (n_cell[i]+6)*twoexp;
    dx[i]  = dxc[i]/twoexp;
    nxyz  *= nx[i];
  }

  Box bx;
#if (BL_SPACEDIM == 3)
  bx = Box(IntVect(0,0,0),
	   IntVect(n_cell[0]-1,n_cell[1]-1,n_cell[2]-1));
#else
  bx = Box(IntVect(0,0),
	   IntVect(n_cell[0]-1,n_cell[1]-1));
#endif

  bx.grow(3);

  MultiFab mfc;
  int      cfac = 32;
  for (int dir=0; dir<BL_SPACEDIM; dir++)
    cfac = std::min(32,n_cell[dir]);

  Array<int> domloc(BL_SPACEDIM),domhic(BL_SPACEDIM);
  if (cond_option >= 2) {
    if (kcfile.empty())
    {
      // create a coarse box
	Box bxc;
#if (BL_SPACEDIM == 3)
	bxc = Box(IntVect(0,0,0),
		  IntVect(n_cell[0]/cfac-1,n_cell[1]/cfac-1,n_cell[2]/cfac-1));
#else 
	bxc = Box(IntVect(0,0),
		  IntVect(n_cell[0]/cfac-1,n_cell[1]/cfac-1));
#endif
      // coarse values initialized randomly
      BoxArray bac(bxc);
      mfc.define(bac,1,0,Fab_allocate);
      int nkpval = kappaval.size();
      int domlo[BL_SPACEDIM];
      int domhi[BL_SPACEDIM];

      for (int dir = 0; dir < BL_SPACEDIM; dir++)
      { 
	domlo[dir] = 0;
	domhi[dir] = n_cell[dir]/cfac;
      }
      for (MFIter mfi(mfc); mfi.isValid(); ++mfi)
      {
	const int  i     = mfi.index();
	const int*  kp_lo  = mfc[mfi].loVect();
	const int*  kp_hi  = mfc[mfi].hiVect();
	const Real* kp_dat = mfc[mfi].dataPtr();

	srand(int(time(0)));
	int iuc = rand()%1000000 + i;
	FORT_PHIRAND(kp_dat,ARLIM(kp_lo),ARLIM(kp_hi),
		     kappaval.dataPtr(),&nkpval,&dkappa, 
		     domlo,domhi,&iuc);
      
      }
      
      kcfile = "coarse_kp";
      VisMF::Write(mfc,kcfile);
    }
    else 
    { 
      kcfile += "/kp";
      VisMF::Read(mfc,kcfile);
      if (mfc.size() > 1)  BoxLib::Abort("kcfile's size is > 1.");
      for (int i=0;i<BL_SPACEDIM; i++)  
      {  
	dxc[i] = (probhi[i]-problo[i])/n_cell[i]*8;
      }  
    }

    if (ParallelDescriptor::IOProcessor()) 
      std::cout << "Initialized coarse distribution\n";


    // We assume the multifab containing the conditioning data has
    // no ghost cells.
    for (MFIter mfi(mfc); mfi.isValid(); ++mfi)
    {
      const int* k_lo  = mfc[mfi].loVect();
      const int* k_hi  = mfc[mfi].hiVect();

      c_sz = 1;
      for (int i=0;i<BL_SPACEDIM; i++) {
	c_sz *= k_hi[i]-k_lo[i] + 1;
	domloc[i] = k_lo[i];
	domhic[i] = k_hi[i];

      }
    }
    ParallelDescriptor::ReduceIntSum(c_sz);
    ParallelDescriptor::ReduceIntMin(domloc.dataPtr(),BL_SPACEDIM);
    ParallelDescriptor::ReduceIntMax(domhic.dataPtr(),BL_SPACEDIM);

    c_idx[0] = 1;
    for (int i=1;i<10;i++)
      c_idx[i] = c_idx[i-1] + c_sz;

    c_sz *= 10;
  }
  
  Array <Real> scratch_c(1);
  if (cond_option > 0) {
    scratch_c.resize(c_sz);

    for (int i=0; i<c_sz; i++)
      scratch_c[i] = 1.e20;

    for (MFIter mfi(mfc); mfi.isValid(); ++mfi)
    {
      const int* k_lo  = mfc[mfi].loVect();
      const int* k_hi  = mfc[mfi].hiVect();
      const Real* kdat = mfc[mfi].dataPtr();

      FORT_INTERNAL_DATA(kdat,ARLIM(k_lo),ARLIM(k_hi),
			 scratch_c.dataPtr(),&c_sz,c_idx,
			 &kappaval[0],&dkappa,dxc,
			 domloc.dataPtr(),domhic.dataPtr());
    }

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMin(scratch_c.dataPtr(),c_sz,IOProc);
    ParallelDescriptor::Bcast(scratch_c.dataPtr(),c_sz,IOProc);
  }

  
  bx.refine(twoexp);
  BoxArray ba(bx); 
  ba.maxSize(maxgridsize);
  MultiFab mf(ba,1,10);
  mf.setVal(0.);
  
  Array< Array<Real> > scratch_r(mf.size());
  Array< Array<int>  > scratch_i(mf.size());
  Array< Array<int>  > order(mf.size());

  if (ParallelDescriptor::IOProcessor())
    std::cout << "Start to generate gaussian simulated distribution ...\n";

  int max_fab_size = 0;
  for (MFIter mfi(mf); mfi.isValid(); ++mfi)
  {    
    const int  i     = mfi.index();
    const int* lo    = mfi.validbox().loVect();
    const int* hi    = mfi.validbox().hiVect();
	  
    const int* k_lo  = mf[mfi].loVect();
    const int* k_hi  = mf[mfi].hiVect();
    const Real* kdat = mf[mfi].dataPtr();
    int nvalid = 1;
    int ntotal = 1;
    for (int d=0; d<BL_SPACEDIM; d++) 
    {
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
		     lo,hi,dx,&iuc);
  }
  ParallelDescriptor::ReduceIntMax(max_fab_size);

  int it = 0;
  while (it < max_fab_size)
  {
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {    
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
  
  for (MFIter mfi(mf); mfi.isValid(); ++mfi)
  {    
    const int  i     = mfi.index();
    
    const int* k_lo  = mf[mfi].loVect();
    const int* k_hi  = mf[mfi].hiVect();
    const Real* kdat = mf[mfi].dataPtr();

    FORT_SGSIM_POST(kdat,ARLIM(k_lo),ARLIM(k_hi),
		    scratch_c.dataPtr(),&c_sz,c_idx,
		    scratch_r[i].dataPtr(),&real_sz,real_idx, 
		    scratch_i[i].dataPtr(),&int_sz,int_idx);

    FORT_LGNORM(kdat,ARLIM(k_lo),ARLIM(k_hi),
		&kappaval[0],&dkappa);

  
  }
  ParallelDescriptor::Barrier();

  FORT_SGSIM_DEALLOC();

  // Parallel copy
  mfdata.copy(mf); 

}

// void  
// PGslib::readConductivity(Real density, Real muval, std::string kfile)
// {
//   //
//   // Read in conducitivity values given by Ye Zhang and
//   // convert it to permeability value in mDa.  
//   //

//   const char* infile = "Kf.txt";

//   // Create grids at finest level.
//   int ng;
  
//   ng = 3*twoexp;

//   Box bx;
// #if (BL_SPACEDIM == 3)
//   bx = Box(IntVect(0,0,0),
// 	   IntVect(n_cell[0]-1,n_cell[1]-1,n_cell[2]-1));
// #else
//   bx = Box(IntVect(0,0),
// 	   IntVect(n_cell[0]-1,n_cell[1]-1));
// #endif
//   bx.grow(3);
//   bx.refine(twoexp);

//   BoxArray ba(bx);
//   MultiFab mf(ba,1,0);
//   mf.setVal(5.e-5);

//   // unit in input is mD x 10^-7
//   // unit in file is m/yr
//   // conversion k = kf*(mu*1e-3)/(rho*g*31536000) x 1e15 x 1e-7 mD
//   Real coef = muval/(density*9.81*31536000)*1.e5;
//   if (ParallelDescriptor::MyProc() == mf.DistributionMap()[0])
//   {
//     std::ifstream indata;
//     indata.open(infile);
//     if (!indata.is_open())
//       BoxLib::Error("File cannot be opened");
//     else
//     {
//       std::string line;
//       Real tmp;

//       // skip the first line
//       indata >> line;
//       std::getline(indata, line);

//       FArrayBox& fab = mf[0];
//       const int* lo = fab.loVect();
//       const int* hi = fab.hiVect();

// #if (BL_SPACEDIM == 2)	
// 	for (int iy=lo[1]+ng-1; iy<hi[1]-ng+1; iy++) {
// 	  for (int ix=lo[0]+ng-1; ix<hi[0]-ng+1; ix++) {
// 	    indata >> tmp;
// 	    fab(IntVect(ix,iy),0) = tmp*coef;
// 	  }
// 	}
// #else
// 	for (int iz=lo[2]+ng-1; iz<hi[2]-ng+1; iz++) {
// 	  for (int iy=lo[1]+ng-1; iy<hi[1]-ng+1; iy++) {
// 	    for (int ix=lo[0]+ng-1; ix<hi[0]-ng+1; ix++) {
// 	       indata >> tmp;
// 	       fab(IntVect(ix,iy,iz),0) = tmp*coef;
// 	    }
// 	  }
// 	}

	
// 	for (int iz=lo[2]; iz<hi[2]; iz++) {
// 	  for (int iy=lo[1]; iy<hi[1]; iy++) {
// 	    for (int ix=lo[0]; ix<hi[0]; ix++) {
// 	      if (ix < lo[0]+ng-1)
// 		fab(IntVect(ix,iy,iz),0) = fab(IntVect(lo[0]+ng-1,iy,iz),0);
// 	      if (ix > hi[0]-ng+1)
// 		fab(IntVect(ix,iy,iz),0) = fab(IntVect(hi[0]-ng+1,iy,iz),0);
// 	    }
// 	  }
// 	}
// #endif
//     }
//     indata.close();
//   }
  
    
//   ParallelDescriptor::Barrier();
  
//   ba.maxSize(maxBaseGrid);
//   MultiFab mfsmall(ba,1,0);
//   mfsmall.copy(mf);  
  
//   VisMF::Write(mfsmall,kfile);
// }

// //
// // This constructs the non-uniform layered structure at F-Basin.
// //
// void  
// PGslib::interp_layer(Array <Real>    kappaval, 
// 		     Real            dkappa, 
// 		     std::string     kfile)
// {

//   // Create grids at finest level.
//   Box bx;
// #if (BL_SPACEDIM == 3)
//   bx = Box(IntVect(0,0,0),
// 	   IntVect(n_cell[0]-1,n_cell[1]-1,n_cell[2]-1));
// #else
//   bx = Box(IntVect(0,0),
// 	   IntVect(n_cell[0]-1,n_cell[1]-1));
// #endif
//   bx.grow(3);
//   bx.refine(twoexp);
   
//   BoxArray ba(bx);
//   ba.maxSize(maxBaseGrid);

//   MultiFab mf(ba,1,0);
//   MultiFab zlayer(ba,1,0);

 
//   Real dx[BL_SPACEDIM];
//   for (int i=0;i<BL_SPACEDIM; i++) 
//     dx[i] = (probhi[i]-problo[i])/(n_cell[i]*twoexp);

// #if 0
//   Real problog[BL_SPACEDIM], probhig[BL_SPACEDIM];
 
//   for (int i=0;i<BL_SPACEDIM; i++) 
//   {
//     dx[i] = (probhi[i]-problo[i])/(n_cell[i]*twoexp);
//     problog[i] = problo[i] - 3*twoexp*dx[i];
//     probhig[i] = probhi[i] + 3*twoexp*dx[i];
//     std::cout << problog[i] << ' ' << probhig[i] << '\n';
//   }
  

//   std::cout << ba.size() << '\n';
//   for (int i = 0; i < ba.size(); i++)
//   {
//     RealBox gridloc = RealBox(ba[i],dx,problog);
//     for (int n = 0; n< BL_SPACEDIM; n++)
//       std::cout << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
//   }
// #endif 

//   int nkpval = kappaval.size();
  
//   const int IOProc = ParallelDescriptor::IOProcessorNumber();
//   char infile[10];
//   for (int nl = 0; nl < nkpval+1; nl++)
//   {
//     int ns;
//     Array<Real> x, y, z;

//     if (ParallelDescriptor::IOProcessor())
//     {
//       sprintf(infile,"surf%i.txt",nl);
//       std::cout << "Opening surface file " << infile << std::endl;

//       std::ifstream indata;
//       indata.open(infile);
//       if (!indata.is_open())
// 	BoxLib::Error("File cannot be opened :: check the layer files.");
//       else
//       {
// 	indata >> ns;

// 	x.resize(ns);
// 	y.resize(ns);
// 	z.resize(ns);
// 	for (int nn = 0; nn < ns; nn++) 
// 	{
// 	  indata >> x[nn] >> y[nn] >> z[nn];
// 	}
// 	indata.close();
//       }
      
//     }

//     ParallelDescriptor::Bcast(&ns, 1, IOProc);
//     if (!ParallelDescriptor::IOProcessor())
//     {
//       	x.resize(ns);
// 	y.resize(ns);
// 	z.resize(ns);
//     }

//     ParallelDescriptor::Bcast(x.dataPtr(), ns, IOProc);
//     ParallelDescriptor::Bcast(y.dataPtr(), ns, IOProc);
//     ParallelDescriptor::Bcast(z.dataPtr(), ns, IOProc);

//     for (MFIter mfi(zlayer); mfi.isValid(); ++mfi)
//     {
//       const Real* zl_dat = zlayer[mfi].dataPtr();
//       const int*  zl_lo  = zlayer[mfi].loVect();
//       const int*  zl_hi  = zlayer[mfi].hiVect();


// #if (BL_SPACEDIM == 3)
//       FORT_INTERP_LAYER(zl_dat,ARLIM(zl_lo),ARLIM(zl_hi),
// 			x.dataPtr(),y.dataPtr(),z.dataPtr(),
// 			&ns,&nl,dx,problo.dataPtr(),probhi.dataPtr());
// #else
//       FORT_INTERP_LAYER(zl_dat,ARLIM(zl_lo),ARLIM(zl_hi),
// 			x.dataPtr(),z.dataPtr(),
// 			&ns,&nl,dx,problo.dataPtr(),probhi.dataPtr());
// #endif

//     }
//     std::string zf = "zlayer";
//     VisMF::Write(zlayer,zf);
     
//   }

//   for (MFIter mfi(mf); mfi.isValid(); ++mfi)
//   {
//     const Real* kp_dat = mf[mfi].dataPtr();
//     const int*  kp_lo  = mf[mfi].loVect();
//     const int*  kp_hi  = mf[mfi].hiVect();

//     const Real* zl_dat = zlayer[mfi].dataPtr();
//     const int*  zl_lo  = zlayer[mfi].loVect();
//     const int*  zl_hi  = zlayer[mfi].hiVect();

//     FORT_ASSIGN_LAYER(kp_dat,ARLIM(kp_lo),ARLIM(kp_hi),
// 		      zl_dat,ARLIM(zl_lo),ARLIM(zl_hi),
// 		      kappaval.dataPtr(),&nkpval,&dkappa,
// 		      dx,problo.dataPtr(),probhi.dataPtr());

//   }
  
//   VisMF::Write(mf,kfile);

// }

// //
// // This invokes user-defined permeability distribution defined
// // in PROB_$D.F
// //
// void  
// PGslib::user_def(Array <Real>    kappaval, 
// 		 Real            dkappa, 
// 		 std::string     kfile)
// {

//   // Create grids at finest level.
//   Box bx;
// #if (BL_SPACEDIM == 3)
//   bx = Box(IntVect(0,0,0),
// 	   IntVect(n_cell[0]-1,n_cell[1]-1,n_cell[2]-1));
// #else
//   bx = Box(IntVect(0,0),
// 	   IntVect(n_cell[0]-1,n_cell[1]-1));
// #endif
//   bx.grow(3);
//   bx.refine(twoexp);

//   BoxArray ba(bx);
//   ba.maxSize(maxBaseGrid);

//   MultiFab mf(ba,1,0);
 
//   Real dx[BL_SPACEDIM];
//   for (int i=0;i<BL_SPACEDIM; i++) 
//     dx[i] = (probhi[i]-problo[i])/(n_cell[i]*twoexp);

//   int nkpval = kappaval.size();
//   for (MFIter mfi(mf); mfi.isValid(); ++mfi)
//   {
//     const Real* kp_dat = mf[mfi].dataPtr();
//     const int*  kp_lo  = mf[mfi].loVect();
//     const int*  kp_hi  = mf[mfi].hiVect();

//     FORT_USER_KAPPA(kp_dat,ARLIM(kp_lo),ARLIM(kp_hi),
// 		    kappaval.dataPtr(),&nkpval,&dkappa,
// 		    dx,problo.dataPtr(),probhi.dataPtr());
//   }

//   VisMF::Write(mf,kfile);

// }
