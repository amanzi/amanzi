
#ifndef _MACBNDRYDATA_H_
#define _MACBNDRYDATA_H_

//
// $Id: MacBndry.H,v 1.3 2011-01-14 02:45:01 gpau Exp $
//


#include <InterpBndryData.H>
#include <BC_TYPES.H>
#include <BC_TYPES_PM.H>

class MacBndry
    :
    public InterpBndryData
{
public:

    MacBndry ();

    MacBndry (const BoxArray& _grids,
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
