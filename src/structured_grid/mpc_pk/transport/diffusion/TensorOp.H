#ifndef _TensorOp_H_
#define _TensorOp_H_

#include <Tuple.H>
#include <MCLinOp.H>

class TensorOp
  :
  public MCLinOp
{
public:
  //
  // Constructor for box array, boundary data and scalar mesh spacing.
  //
  TensorOp (const BndryData& bd,
	    Real             _h,
	    int              _nc = 1);
  //
  // Constructor for box array, boundary data and vector mesh spacing.
  //
  TensorOp (const BndryData& bd,
	    const Real*      _h,
            int              _nc = 1);
  //
  // The destructor.
  //
  virtual ~TensorOp ();
  //
  // Set scalar coefficients.
  //
  void setScalars (Real _alpha,
                   Real _beta);
  //
  // Get scalar alpha coefficient.
  //
  Real get_alpha () const;
  //
  // Get scalar beta coefficient.
  //
  Real get_beta () const;
  //
  // Return reference to "a" coefficients for base level.
  //
  const MultiFab& aCoefficients (int level = 0);
  //
  // Return reference to "b" coefficients for base level.
  //
  const MultiFab& bCoefficients (int dir,
                                 int level = 0);
  //
  // Return reference to "b" coefficients for base level.
  //
  const MultiFab& b1Coefficients (int dir,
                                  int level = 0);
  //
  // Copy _a into "a" coeffs for base level.
  //
  void aCoefficients (const MultiFab& _a);
  //
  // Copy _b into "b" coeffs in "dir" coordinate direction for base level.
  //
  void bCoefficients (const MultiFab& _b,
                      int             dir);
  //
  // Copy _b1 into "b1" coeffs in "dir" coordinate direction for base level.
  //
  void b1Coefficients (const MultiFab& _b1,
                       int             dir);
  //
  // Alternative (older) version of above members.
  //
  void setCoefficients (const MultiFab &_a,
                        const MultiFab &_bX,
                        const MultiFab &_b1X,
                        const MultiFab &_bY,
                        const MultiFab &_b1Y
#if BL_SPACEDIM==3
                        ,const MultiFab &_bZ,
                        const MultiFab &_b1Z
#endif
    );
  //
  // Collective version of above members (taking an array of MultiFabs for "b" and "b1").
  //
  void setCoefficients (const MultiFab& _a,
                        const MultiFab* _b,
                        const MultiFab* _b1);
  //
  // Allocate/fill internal data for new level.
  //
  virtual void prepareForLevel (int level);
  //
  // Remove internal data for this level and all levels above.
  //
  virtual void clearToLevel (int level);
  //
  // Set flags so that a coeffs at lev and coarser require recalculation.
  //
  void invalidate_a_to_level (int lev);
  //
  // Set flags so that b coeffs at lev and coarser require recalculation.
  //
  void invalidate_b_to_level (int lev);
  //
  // Set flags so that b1 coeffs at lev and coarser require recalculation.
  //
  void invalidate_b1_to_level (int lev);
  //
  // Computes flux associated with operator.
  //
  virtual void compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
                         MultiFab& in, const MCBC_Mode& bc_mode=MCInhomogeneous_BC,
                         int sComp=0, int dComp=0, int nComp=1, int bndComp=0);
    
protected:
  //
  // Initialize a full set (a,b) of coefficients on the box array.
  //
  void initCoefficients (const BoxArray &_ba);
  //
  // Compute out=L (in) at level=level.
  //
  virtual void Fapply (MultiFab&       out,
                       const MultiFab& in,
                       int             level);
  //
  // Apply smoother to improve residual to L(solnL)=rhsL.
  //
  virtual void Fsmooth (MultiFab&       solnL,
                        const MultiFab& rhsL,
                        int             level,
                        int             phaseflag);
  //
  // Return number of components.  This is virtual since only the derived knows.
  //
  virtual int numberComponents ();
  //
  // Number of relaxation phases.  On derived class knows.
  //
  virtual int numberPhases ();
  //
  // Helper function for object construction.
  //
  void initConstruct (const Real* _h);

protected:
  //
  // Array (on level) of "a" coefficients
  //
  Array< MultiFab* > acoefs;
  //
  // Array (on level) of Tuple (on dimension) of "b", "b1" coefficients
  //
  Array< Tuple< MultiFab*, BL_SPACEDIM> > bcoefs;
  Array< Tuple< MultiFab*, BL_SPACEDIM> > b1coefs;
  //
  // Scalar "alpha" coefficient
  //
  Real alpha;
  //
  // Scalar "beta" coefficient
  //
  Real beta;

private:
  //
  // Flag, can a coeffs be trusted at a level.
  //
  Array<int> a_valid;
  //
  // Flag, can b,b1 coeffs be trusted at a level.
  //
  Array<int> b_valid;
  Array<int> b1_valid;
  //
  // Default value for a (MultiFab) coeficient.
  //
  static Real a_def;
  //
  // Default value for b, b1 (MultiFab) coeficient.
  //
  static Real b_def, b1_def;
  //
  // Default value for alpha (scalar) coeficient.
  //
  static Real alpha_def;
  //
  // Default value for beta (scalar) coeficient.
  //
  static Real beta_def;
  //
  // Disallow copy constructors (for now...to be fixed).
  //
  TensorOp (const TensorOp&);
  TensorOp& operator = (const TensorOp&);
};

inline
Real
TensorOp::get_alpha () const
{
  return alpha;
}

inline
Real
TensorOp::get_beta () const
{
  return beta;
}

inline
void
TensorOp::setCoefficients (const MultiFab &_a,
                           const MultiFab &_bX,
                           const MultiFab &_b1X,
                           const MultiFab &_bY,
                           const MultiFab &_b1Y
#if BL_SPACEDIM==3
                           ,const MultiFab &_bZ
                           ,const MultiFab &_b1Z
#endif
  )
{
  aCoefficients(_a);
  bCoefficients(_bX, 0);
  b1Coefficients(_b1X, 0);
  bCoefficients(_bY, 1);
  b1Coefficients(_b1Y, 1);
#if BL_SPACEDIM==3
  bCoefficients(_bZ, 2);
  b1Coefficients(_b1Z, 2);
#endif
}

inline
void
TensorOp::setCoefficients (const MultiFab& _a,
                           const MultiFab* _b,
                           const MultiFab* _b1)
{
  aCoefficients(_a);
  for (int n = 0; n < BL_SPACEDIM; ++n) {
    bCoefficients(_b[n], n);
    b1Coefficients(_b1[n], n);
  }
}

#endif /*_TensorOp_H_*/
