#include <winstd.H>

#ifdef COREREACT

#include "iostream"

#include "PChemDriver.H"
#include "CHEMDRIVER_F.H"
#include <ParallelDescriptor.H>

bool PChemDriver::mInitialized = false;
const Real EQTOL_DEF = 1.e-16;

PChemDriver::PChemDriver (int ncomps, Array<int> num, Array<Real> mw)
  :
  mEQTOL(EQTOL_DEF)
{
   if (!mInitialized)
   {
     FORT_INIT_CHEM_TOUGH();
     
     n_primary = ncomps;
     n_core    = 0;
     for (int i=0; i<4; i++) {
       n_type[i] = num[i];
       n_core += num[i];
       if (i == 2) n_core += num[i];
     }
     mweight.resize(n_primary);
     for (int i=0; i<n_primary; i++) 
       mweight[i] = mw[i];

     mInitialized = true;
   }
}

/* 
VODE solve:
   Y: component densities
   P: pressure
   C: secondary concentrations
   F: function count
*/

void  
PChemDriver::solveTransient(FArrayBox&        Y,
			    FArrayBox&        P,
			    FArrayBox&        C,
			    FArrayBox&        F,
			    const Box&        box,
			    int               ncomps,
			    Real              dt,
			    FArrayBox&        kappa,
			    FArrayBox&        phi) 
{

  FORT_CHEM_TOUGH(Y.dataPtr(),ARLIM(Y.loVect()),ARLIM(Y.hiVect()),
		  P.dataPtr(),ARLIM(P.loVect()),ARLIM(P.hiVect()),
		  C.dataPtr(),ARLIM(C.loVect()),ARLIM(C.hiVect()),
		  F.dataPtr(),ARLIM(F.loVect()),ARLIM(F.hiVect()),
		  &ncomps,box.loVect(),box.hiVect(),&dt,
		  &n_core,n_type,mweight.dataPtr(),
		  kappa.dataPtr(),ARLIM(kappa.loVect()),ARLIM(kappa.hiVect()),
		  phi.dataPtr(),ARLIM(phi.loVect()),ARLIM(phi.hiVect()));
}

void
PChemDriver::solveSaturate(FArrayBox&        Y,
			   FArrayBox&        X,
			   FArrayBox&        P,
			   FArrayBox&        C,
			   const Box&        box,
			   int               ncomps,
			   FArrayBox&        kappa,
			   FArrayBox&        phi) 
{
  FORT_SAT_TOUGH(Y.dataPtr(),ARLIM(Y.loVect()),ARLIM(Y.hiVect()),
		 X.dataPtr(),ARLIM(X.loVect()),ARLIM(X.hiVect()),
		 P.dataPtr(),ARLIM(P.loVect()),ARLIM(P.hiVect()),
		 C.dataPtr(),ARLIM(C.loVect()),ARLIM(C.hiVect()),
		 &ncomps,box.loVect(),box.hiVect(),
		 &n_core,n_type,mweight.dataPtr(),
		 kappa.dataPtr(),ARLIM(kappa.loVect()),ARLIM(kappa.hiVect()),
		 phi.dataPtr(),ARLIM(phi.loVect()),ARLIM(phi.hiVect()));
}

void
PChemDriver::totalToInd(FArrayBox&        Y,
			FArrayBox&        P,
			FArrayBox&        C,
			const Box&        box,
			int               ncomps,
			FArrayBox&        kappa,
			FArrayBox&        phi) 
{

  FORT_INIT_SECONDARY_TOUGH(Y.dataPtr(),ARLIM(Y.loVect()),ARLIM(Y.hiVect()),
			    P.dataPtr(),ARLIM(P.loVect()),ARLIM(P.hiVect()),
			    C.dataPtr(),ARLIM(C.loVect()),ARLIM(C.hiVect()),
			    &ncomps,box.loVect(),box.hiVect(),
			    &n_core,n_type,mweight.dataPtr(),
			    kappa.dataPtr(),ARLIM(kappa.loVect()),ARLIM(kappa.hiVect()),
			    phi.dataPtr(),ARLIM(phi.loVect()),ARLIM(phi.hiVect()));
}
#endif
