#ifndef _PRESSBNDRYDATA_H_
#define _PRESSBNDRYDATA_H_

//
// $Id: PressBndry.H,v 1.1 2011-08-27 18:04:11 gpau Exp $
//

#include <InterpBndryData.H>
#include <BC_TYPES.H>
#include <BC_TYPES_PM.H>

class PressBndry
    :
    public InterpBndryData
{
public:

    PressBndry ();

    PressBndry (const BoxArray& _grids,
		int             _ncomp,
		const Geometry& _geom);

    virtual void setBndryConds (const BCRec& phys_bc,
                                IntVect&     ratio,
                                int          comp=0);

    void setHomogValues (const BCRec& phys_bc,
			 IntVect&     ratio,
			 int          comp=0);

};

#endif
