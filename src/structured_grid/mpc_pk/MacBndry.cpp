//
// $Id: MacBndry.cpp,v 1.5 2011-01-14 02:45:01 gpau Exp $
//
#include <winstd.H>

#include <LO_BCTYPES.H>
#include <MacBndry.H>

MacBndry::MacBndry ()
    :
    InterpBndryData()
{}

MacBndry::MacBndry (const BoxArray& _grids,
                    int             _ncomp,
                    const Geometry& _geom)
    :
    InterpBndryData(_grids,_ncomp,_geom)
{}

void
MacBndry::setBndryConds (const BCRec& phys_bc,
                         IntVect&     ratio,
			 int          comp)
{
    //
    // ALL BCLOC VALUES ARE NOW DEFINED AS A LENGTH IN PHYSICAL
    // DIMENSIONS *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE
    //
    const BoxArray& grids      = boxes();
    const int ngrds            = grids.size();
    const Real* dx             = geom.CellSize();
    const Box& domain          = geom.Domain();

    for (OrientationIter fi; fi; ++fi)
    {
      Array<Real>& bloc                = bcloc[fi()];
      Array< Array<BoundCond> >& bctag = bcond[fi()];

      int dir    = fi().coordDir();
      Real delta = dx[dir]*ratio[dir];
      int p_bc   = (fi().isLow() ? phys_bc.lo(dir) : phys_bc.hi(dir));

      for (int i = 0; i < ngrds; i++)
      {
	const Box& grd = grids[i];
	
	if (domain[fi()] == grd[fi()] && !geom.isPeriodic(dir))
	{
	  //
	  // All physical bc values are located on face.
	  //
          
	  if (p_bc == EXT_DIR)
	  {
	    bctag[i][comp] = LO_DIRICHLET;
	    bloc[i] = 0.;
	  }
	  else if (p_bc == FOEXTRAP      ||
		   p_bc == HOEXTRAP      || 
		   p_bc == REFLECT_EVEN  ||
		   p_bc == SEEPAGE)
	  {
	    bctag[i][comp] = LO_NEUMANN;
	    bloc[i] = 0.;
	  }
	  else if (p_bc == REFLECT_ODD)
	  {
	    bctag[i][comp] = LO_REFLECT_ODD;
	    bloc[i] = 0.;
	  }
	}
	else
	{
	  //
	  // Internal bndry.
	  //
	  bctag[i][comp] = LO_DIRICHLET;
	  bloc[i] = 0.5*delta;
	}
      }
    }
}

void
MacBndry::setHomogValues (const BCRec& bc,
			   IntVect&     ratio,
                           int          ncomp)
{
  setBndryConds(bc, ratio, ncomp);

    for (OrientationIter fi; fi; ++fi)
    {
      for (FabSetIter fsi(bndry[fi()]); fsi.isValid(); ++fsi)
      {
	bndry[fi()][fsi].setVal(0);
      }
    }
}
