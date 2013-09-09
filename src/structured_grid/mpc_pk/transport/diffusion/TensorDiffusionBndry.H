
#ifndef _TensorDiffusionBndry_H_
#define _TensorDiffusionBndry_H_

#include <MCInterpBndryData.H>
#include <BC_TYPES.H>
#include <MCLinOp.H>

class TensorDiffusionBndry
    :
    public MCInterpBndryData
{
public:
    TensorDiffusionBndry () {}

    TensorDiffusionBndry (const BoxArray& grids,
                          int             nc,
                          const Geometry& geom)
        : MCInterpBndryData(grids,nc,geom) {}

    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio,
                                int          comp=0);

    void setHomogValues (const Array<BCRec>& phys_bc,
                         int                 ratio);

};

#endif /*_TensorDiffusionBndry_H_*/
