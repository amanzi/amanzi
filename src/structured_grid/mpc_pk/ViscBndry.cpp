#include <winstd.H>

#include <LO_BCTYPES.H>
#include <ViscBndry.H>

void
ViscBndry::setBndryConds (const BCRec&   bc,
                          IntVect& ratio,
			  int comp)
{
    //
    //  NOTE: ALL BCLOC VALUES ARE DEFINED AS A LENGTH IN PHYSICAL
    //        DIMENSIONS *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE
    //
    const Real* dx    = geom.CellSize();
    const Box& domain = geom.Domain();

    for (OrientationIter fi; fi; ++fi)
    {
        Array<Real>& bloc = bcloc[fi()];
        Array< Array<BoundCond> >& bctag = bcond[fi()];

        int dir          = fi().coordDir();
        const Real delta = dx[dir]*ratio[dir];
        int p_bc         = (fi().isLow() ? bc.lo(dir) : bc.hi(dir));
	
        for (int i = 0; i < boxes().size(); i++)
        {
            if (domain[fi()] == boxes()[i][fi()] && !geom.isPeriodic(dir))
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
		    bctag[i][comp] = LO_NEUMANN;//LO_REFLECT_ODD;
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
ViscBndry::setScalarBndryConds (const BCRec&   bc,
				IntVect& ratio,
				int comp)
{
    //
    //  NOTE: ALL BCLOC VALUES ARE DEFINED AS A LENGTH IN PHYSICAL
    //        DIMENSIONS *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE
    //
    const Real* dx    = geom.CellSize();
    const Box& domain = geom.Domain();

    for (OrientationIter fi; fi; ++fi)
    {
        Array<Real>& bloc = bcloc[fi()];
        Array< Array<BoundCond> >& bctag = bcond[fi()];

        int dir          = fi().coordDir();
        const Real delta = dx[dir]*ratio[dir];
        int p_bc         = (fi().isLow() ? bc.lo(dir) : bc.hi(dir));
	
        for (int i = 0; i < boxes().size(); i++)
        {
            if (domain[fi()] == boxes()[i][fi()] && !geom.isPeriodic(dir))
            {
	      //
	      // All physical bc values are located on face.
	      //
	      if (p_bc == EXT_DIR)
                {
                    bctag[i][comp] = LO_DIRICHLET;
                    bloc[i] = 0.;
                }
                else if (p_bc == SEEPAGE       ||
			 p_bc == FOEXTRAP      ||
                         p_bc == HOEXTRAP      || 
                         p_bc == REFLECT_EVEN)
                {
                    bctag[i][comp] = LO_NEUMANN;
                    bloc[i] = 0.;
                }
                else if (p_bc == REFLECT_ODD)
                {
		    bctag[i][comp] = LO_NEUMANN; //LO_REFLECT_ODD;
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
ViscBndry::setHomogValues (const BCRec& bc,
                           IntVect& ratio)
{
    setBndryConds(bc, ratio);

    for (OrientationIter fi; fi; ++fi)
    {
        for (FabSetIter fsi(bndry[fi()]); fsi.isValid(); ++fsi)
        {
            bndry[fi()][fsi].setVal(0);
        }
    }
}

void
ViscBndry::setScalarValues (const BCRec&    bc,
			    IntVect&        ratio,
			    const MultiFab* mf)
{
    setScalarBndryConds(bc, ratio);

    for (OrientationIter fi; fi; ++fi) 
    {      
      for (FabSetIter fsi(bndry[fi()]); fsi.isValid(); ++fsi)
      {	      
	const Box& bx = fsi.validbox();
	int idx = fsi.index();
	bndry[fi()][fsi].copy((*mf)[idx],bx,0,bx,0,1);	    
      }
    }
}

void
ViscBndry::setdeltaSValues (const BCRec& bc,
			    IntVect& ratio)
{
    setScalarBndryConds(bc, ratio);

    for (OrientationIter fi; fi; ++fi) 
    {
      for (FabSetIter fsi(bndry[fi()]); fsi.isValid(); ++fsi)
      {	      
	bndry[fi()][fsi].setVal(0.);
      }      
    }
}
