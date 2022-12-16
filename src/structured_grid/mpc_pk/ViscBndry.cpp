/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <winstd.H>

#include <LO_BCTYPES.H>
#include <ViscBndry.H>

void
ViscBndry::setBndryConds (const BCRec&   bc,
                          const IntVect& ratio,
			  int            comp)
{
    //
    //  NOTE: ALL BCLOC VALUES ARE DEFINED AS A LENGTH IN PHYSICAL
    //        DIMENSIONS *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE
    //
    const Real* dx    = geom.CellSize();
    const Box& domain = geom.Domain();
    //
    // We note that all orientations of the FabSets have the same distribution.
    // We'll use the low 0 side as the model.
    //
    for (FabSetIter fsi(bndry[Orientation(0,Orientation::low)]); fsi.isValid(); ++fsi)
    {
        const int                  i     = fsi.index();
        RealTuple&                 bloc  = bcloc[i];
        Array< Array<BoundCond> >& bctag = bcond[i];

        for (OrientationIter fi; fi; ++fi)
        {
            const Orientation face  = fi();
            const int         dir   = face.coordDir();

            if (domain[face] == boxes()[i][face] && !geom.isPeriodic(dir))
            {
                //
                // All physical bc values are located on face.
                //
                const int p_bc  = (face.isLow() ? bc.lo(dir) : bc.hi(dir));

                if (p_bc == EXT_DIR)
                {
                    bctag[face][comp] = LO_DIRICHLET;
                    bloc[face]        = 0;
                }
                else if (p_bc == FOEXTRAP      ||
                         p_bc == HOEXTRAP      ||
                         p_bc == REFLECT_EVEN  ||
			 p_bc == SEEPAGE)
                {
                    bctag[face][comp] = LO_NEUMANN;
                    bloc[face]        = 0;
                }
                else if (p_bc == REFLECT_ODD)
                {
                    bctag[face][comp] = LO_REFLECT_ODD;
                    bloc[face]        = 0;
                }
            }
            else
            {
                //
                // Internal bndry.
                //
                const Real delta = dx[dir]*ratio[dir];
                if (delta < 0) {
                  std::cout << "hello" << std::endl;
                }
                BL_ASSERT(delta >= 0);

                bctag[face][comp] = LO_DIRICHLET;
                bloc[face]        = 0.5*delta;
            }
        }
    }
}

void
ViscBndry::setHomogValues (const BCRec& bc,
                           const IntVect& ratio)
{
    setBndryConds(bc, ratio);

    for (OrientationIter fi; fi; ++fi)
    {
      const Orientation& face = fi();
      for (FabSetIter fsi(bndry[face]); fsi.isValid(); ++fsi)
        {
	  bndry[face][fsi].setVal(0);
        }
    }
}

void
ViscBndry::setScalarValues (const BCRec&    bc,
			    IntVect&        ratio,
			    const MultiFab* mf)
{
  for (int n=0; n<nComp(); ++n) {
    setBndryConds (bc, ratio, n);
  }

    for (OrientationIter fi; fi; ++fi)
    {
      int face = fi();
      for (FabSetIter fsi(bndry[face]); fsi.isValid(); ++fsi)
      {
	const Box& bx = fsi.validbox();
	int idx = fsi.index();
	bndry[face][fsi].copy((*mf)[idx],bx,0,bx,0,1);
      }
    }
}

void
ViscBndry::setdeltaSValues (const BCRec& bc,
			    IntVect& ratio)
{
  for (int n=0; n<nComp(); ++n) {
    setBndryConds (bc, ratio, n);
  }

    for (OrientationIter fi; fi; ++fi)
    {
      int face = fi();
      for (FabSetIter fsi(bndry[face]); fsi.isValid(); ++fsi)
      {
	bndry[face][fsi].setVal(0.);
      }
    }
}
