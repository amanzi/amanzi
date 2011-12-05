#ifndef _PRESSBNDRYDATA_H_
#define _PRESSBNDRYDATA_H_

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
