
How Amanzi-S works!

This file will be considered somewhat of a living document, serving as
a location where some of the more esoteric aspects of the structured
grid implementation are outlined.  This is not to be considered
official documentation, but rather a set of notes to pass around
between developers and "expert users" (a group that is likely the same
as the aforementioned).  Anyway, let's get to it...


The "Multi-level Richard Solve" (aka, the Newton solver for the
backward-Euler Richards equation steps):

The goal of the RichardSolver class is to orchestrate the solution of
the following equation:

      d(rho.phi.s)/dt + Div(q) = 0,                      (1)

subject to a combination of Dirichlet and Neumann type boundary
conditions.  Here rho is the fluid (water) density and s is the water
saturation of the material.  The Darcy flux, q, is

      q = -(rho.k.kr/mu)( Grad(p) + rho.g)

p is the water pressure, k and kr=kr(s) are the intrinsic and relative
permeabilty, mu is the viscosity and g is the gravity vector.  For
unsaturated conditions, p = -pcap, where pcap=pcap(s) is the capiilary
pressure.


We solve (1) in finite-volume form.  Integrating over each discrete cell
(using the divergence theorem, taking the state at the center of
each cell to represent the cell-average value), and using a backward-
Euler (in time) discretization over time interval, dt:

(rho.phi.s)^{new}-(rho.phi.s)^{old} + (dt/vol).Sum(q^{new}.Area) = Res

As a result, Res=0 is a nonlinear equation for the state at "new" time.
We solve this equation using damped Newton iterations, specifically 
using p as the primary variable.  There are many details....


The big picture:

In Amanzi-speak, the "MPC" (multi-process coordinator) and flow PK
(process kernal) work together to drive the solution of this system,
the sense that they provide the initial data (p^{old}) and the time
interval, dt.  The solver attempts to solve (1) for dt, and reports
back on its success (or failure) to do so.  The time-step is adjusted,
and the solver is called to integrate the next time interval (or retry
over a shorter interval, as needed).

The solution process itself is managed in the RichardSolver class.



Calling the solver

There are a number of options for the solver piece of this task,
including the ability to select different solver approaches.  By far,
the most developed approach is to let PETSc's SNES nonlinear solver
software drive.  In this case, the code in RichardSolver.cpp provides
all the functions that are expected by SNES to do its work.  Well, 
actually, that is where the code exists to manage this construction.
The gory details of evaluating material properties, computing 
saturations from pressures, and that sort of thing, are done by the
original worker functions buried in the PorousMedia class.  That
being said, RichardSolver is really devoted as an interface between
PorousMedia and PETSc.  As such, it does need quite a few support
class and routines.  We'll mention a few as we go, but this is not
an exaustive manual by any stretch.  Here's a brief rundown of what
the RichardSolver does when asked to solve Richards equation:

1. Copy the initial state from the data structures native to Amanzi-S
into those expected by vanilla PETSc (note that we might have chosen
to extend PETSc's structures to BoxLib data directly, but that would
have been a larger effort not yet warranted....ie, premature
optimization!).

2. Driving PETSc:
- Set a pointer identifying the residual function
- Set a "post check" function (which provides an alternative update)
- Optionally scale the unknowns 
- Call SNES->Solve
- Optionally unscale the unknowns
- Return to caller, with some clues about success and difficulty


Building the residual:

There is one residual function, RichardRes_DtDt, called by SNES.  It
takes x, the pressure and returns f, the residual.  The pressure
enters as a PETSc "Vec", it is optionally scaled,and then converted to
a BoxLib MFTower (more on those later).  DpDtResidual then does the
work on BoxLib structures, passing the the hardest part to the routine
DivRhoU.  Once DivRhoU is computed, the final residual is assembled in
a very simple fortran routine cleverly named, FORT_RS_PDOTRES.

Unless told otherwise, PETSc will use RichardRes_DtDt also to
construct a Jacobian matrix, using finite differences, ie:
J_ij=dR_i/dp_j.  It does so by twiddling each cell's pressure
indepenently and recomputing R+ = R(p + p_eps).  Then J_ij=(R+ -
R)/p_eps.  Well, it's a little smarter, in that we've set it up to
divide the complete set of unknowns into independent sets (colors) so
that the number of actual calls to the residual function for this
purpose is optimally minimal.  As a variation on this strategy, if
semi_analytic_J is true, the J is computed using the routine
"RichardR2" instead.  RichardR2 only computes the DivRhoU piece of R.
The time derivative part is computed analytically, and added to the
matrix elements inside the FD driver, paradoxically (but maybe less so
now?) named SemiAnalyticMatFDColoringApply.


Options, options, options....

As usual, PETSc brings in a whole raft of options, from how to solve
the linear systems, form the J (or not), using a sparse data structure
for J, whether and how to precondition the linear solve, etc... In
vanilla PETSc, you can add these options to the command line.  For
Amanzi, I added a parameter "PETSc Options File" to the XML input,
which is formatted in the usual PETSc way.


One particularly important option is related to "damping" the Newton
step, that is, modifying the dp that is computed by Newton prior to
actually updating the solution for that iteration.  PETSc provides a
few approaches that are all based on scaling the entire update to
guarantee a reduction in the in the residual.  I have found that this
is too restrictive.  The code here includes a couple of other
approaches.  The first is activated by choosing "-snes_ls_basic".  By
default this would do nothing and would simply take the entire dp.
However, the SNES solver points to the function PostCheck in
RichardSolver.cpp in this case, where I've implemented a simple
geometrical strategy to determine the scaling factor such that the new
residual norm is less than some factor, ls_acceptance_factor, times
the original residual norm.  If the factor is > 1, for example, you're
ok with the residual to increase a little before dropping.

A second approach we tried is a bit more complicated.  First, using
the chain rule, we compute the change in saturation associated with
the Newton dp (using the same code was discussed above in the
semi-analytic J formulation).  If the new saturation is larger than a
threshold value (defined at the top of RichardSolver.cpp, and via
runtime options), a new dp is computed such that s remains bounded
between 0 and 1.  Note that this is essentially just a clever way to
implement "variable switching", and tries to keep the updates
reasonable.  It is activated if the saturation threshold is > 0 by
setting the PostCheck function for SNES to PostCheckAlt


Another biggie is the linear solver to use for the Newton system.  It
turns out that as the time step gets arbitrarily large, the Jacobian
matrix becomes really badly conditioned, horribly asymmetric and
dominated by off-diagonal terms.  Because of this, standard
Krylov-based linear solvers just don't do well. One option is direct
inversion via SuperLU and SuperLU-dist.  These are selected in the
PETSc options file by setting "-pc_type lu",
"-pc_factor_mat_solver_package superlu_dist" and "-ksp_type preonly"

It should be said that at this time, there's no clear winner for all
this stuff.  The set of problems I've been looking at are difficult
enough that getting a solution requires horsing with all of these
until one works.


More grungy detail:

Temporary data

Without being nazi-like about wasted space, I tried to be economical
about memory usage.  The solution p over the AMR levels comes in,
effectively, as an array of MultiFabs, one per level.  These are
loaded directly (without copy) into an "MFTower" data structure.  The
MFTower structure basically formalizes the multilevel relationship, by
tying them together with specified refinement ratios, and geometry
objects (that specify the domain size, grid-spacing, and periodicity,
thus providing things like face areas and cell volumes).  Much of this
information is contained in (and actually pulled directly from) the
Amr/AmrLevel objects in the base library, but MFTower also implements
the notion of a mapping of all "uncovered" cells to a unique index
scheme.  This is necessary so that we can build and use the "matrix"
data required for SNES to represent the Jacobian explicitly.

Similarly, the SNES framework interacts with Vec structures.  Once we
have defined a mapping from the (IntVect,level) concept accessed by
MultiFabs to a single unique integer index, the transfer of data
between the two is straightforward; we even reuse the original BoxLib
data distribution, so the translation to Vec structures is simply a
local copy.  Another way to thik of this is that we are taking a 
very structured data format (MultiFab array) and generating a map to 
a completely unstructured format....ick!  But, in doing so we get
access to a nice set of solvers to play with.  The other thing we do 
pay however is the space to hold all these temporary Vecs used in
translation.

So, as I said, without going overboard, I tried to make MFTower
objects (by reference) for as much of the work space as possible,
minimizing data allocations and copies.  Look in the RichardSolver
constructor to see how well (or badly) I did.  In particular, there's
MFTowers for 3 edge-based quantities (which could probably be cut down
to one, or could maybe even borrow memory allocated but unused for the
PorousMedia class.  Not going there just yet!)  Also, there is a fair
amount of meta-data associated with the MFTowers and such, and
associated with storing the PDE operators (particularly as coarse-fine
interfaces) as a matrix operation, but that is shared between all the
structures as much as possible.


Two other comments regarding RichardSolver....

As you wade through the solver, you'll probably find that there are
really only two aspects of this that drove the design of everything
that followed: the needs of the ComputeDarcyVelocity routine in
RichardSolver.cpp and the interface to the PETSc SNES and linear
solver software.  You might start looking there to understand how
things are laid out and where/when things are computed.

Very little optimization/profiling has been done with this suite of
software.  It is not yet known how much energy is wasted with all the
data copying, and how much memory is available with smarter
organization of the work.  This means there's probably lots of
opportunity for the bright young mind to improve all this!

