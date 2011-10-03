#ifndef _MACOPERATOR_H_
#define _MACOPERATOR_H_

//
// $Id: MacOperator.H,v 1.9 2011-02-17 19:53:30 gpau Exp $
//

#include <ABecLaplacian.H>
#include <BndryData.H>

class MacOperator
    :
    public ABecLaplacian
{
public:

    MacOperator (const BndryData& mgb,
                 const Real*      _h);

    virtual ~MacOperator ();
    //
    // Define the gradient coeficients used by the multigrid object.
    //
    void setCoefficients (MultiFab*   area,
                          MultiFab*   coef,
                          const Real* dx,
                          const BCRec& phys_bc,
                          const Box& domain);
    //
    // Initialize the rhs for the level solves.
    //
    void defRHS (MultiFab* area,
                 MultiFab& vol,
                 MultiFab& rhs,
                 MultiFab* vel,
		 MultiFab* RhoD,
                 const BCRec&    phys_bc,
                 const Box&      domain,
                 Real      scale);
    //
    // Apply the pressure gradient to the mac velocities.
    //
    void velUpdate (MultiFab*       vel,
                    MultiFab&       phi,
                    const MultiFab* coef,
                    const BCRec&    phys_bc,
                    const Real*     dx,
                    const Box&      domain,
		    int             level,
                    Real            scale);

    //
    // Scale the rhs for the sync solves.
    //
    void syncRhs (const MultiFab& Volume,
                  MultiFab&       Rhs,
		  MultiFab*       u_mac_sync,
		  MultiFab*       rhoD_sync,
		  MultiFab*       area,
                  Real            rhs_scale,
		  const BCRec&    phys_bc,
		  const Box&      domain);
};

#endif
