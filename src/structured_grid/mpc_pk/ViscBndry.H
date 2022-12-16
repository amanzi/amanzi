/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _VISCBNDRYDATA_H_
#define _VISCBNDRYDATA_H_

#include <InterpBndryData.H>
#include <BC_TYPES.H>
#include <BC_TYPES_PM.H>

class ViscBndry
    :
    public InterpBndryData
{
public:
    ViscBndry()
        :
        InterpBndryData()  {}

    ViscBndry (const BoxArray& _grids,
               int             _ncomp,
               const Geometry& _geom)
        :
        InterpBndryData(_grids,_ncomp,_geom)
    {}

    virtual void setBndryConds (const BCRec&   phys_bc,
                                const IntVect& ratio,
				int            comp=0);

    void setHomogValues (const BCRec& bc,
                         const IntVect& ratio);

    void setScalarValues (const BCRec& bc,
			     IntVect& ratio,
			     const MultiFab* mf);

    void setdeltaSValues (const BCRec& bc,
                          IntVect& ratio);
};

#endif
