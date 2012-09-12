
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

    virtual void setScalarBndryConds (const BCRec&   phys_bc,
				      IntVect& ratio,
				      int comp=0);

    void setHomogValues (const BCRec& bc,
                         IntVect& ratio);

    void setScalarValues (const BCRec& bc,
			     IntVect& ratio,
			     const MultiFab* mf);

    void setdeltaSValues (const BCRec& bc,
                          IntVect& ratio);
};

#endif 
