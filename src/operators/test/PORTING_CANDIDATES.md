# operators/test port-feasibility survey

Context: while debugging `MPCCoupledWaterNonInjected` (a new tightly-coupled
surface/subsurface water MPC using `TreeOperator` + `Op_MeshInjection` +
`Operator_Schema` for off-diagonal cross-mesh coupling blocks), two real,
previously-undetected bugs were found and fixed in `src/operators/TreeOperator.cc`:

1. `TreeOperator::AssembleMatrix()` was missing a `Amat_->resumeFill()` call
   before `Amat_->zero()` (present in `Operator::AssembleMatrix()`'s
   equivalent no-arg overload, but missing here) -- calling `AssembleMatrix()`
   a second time on an already-`fillComplete()`'d matrix failed a
   `MatrixFE::sumIntoLocalValues` assertion.
2. `TreeOperator::apply()` unconditionally called `Y.putScalar(0.0)` on
   *every* call, including recursive per-block calls used to accumulate
   multiple blocks' contributions (e.g. a diagonal block AND an off-diagonal
   coupling block both writing into the same output rows) into the same `Y`.
   Each subsequent block's call wiped out the previous block's contribution.
   Fixed by respecting the `scalar` argument's documented zero-vs-accumulate
   convention (already followed correctly by `Operator::apply`/
   `applyUnassembled`, just not by `TreeOperator::apply`).

Both bugs went undetected because **no unit tests currently exercise
`TreeOperator` in this build** -- `grep -rl "TreeOperator" src/operators/test/`
finds 3 files, and none of the 3 are referenced in
`src/operators/CMakeLists.txt`. The entire legacy (Epetra-era)
`src/operators/test/` unit test suite was effectively disabled wholesale
during the Epetra->Tpetra migration (either fully commented out, or in a
couple of cases gated behind `if(0)`), and most of those files were never
actually ported to the current Tpetra/Kokkos API, so they would not compile
as-is today.

This file surveys `src/operators/test/*.cc` for port feasibility, prioritized
by relevance to the bugs above (TreeOperator / Op_MeshInjection /
Operator_Schema coupling-block machinery first).

No code has been changed as part of this survey -- recommendations only.


## Priority 1: `src/operators/test/operators_tree.cc`

**Not referenced anywhere in `src/operators/CMakeLists.txt` -- dead/unbuilt.**
Contains 4 `UnitTest++` tests in `SUITE(SURFACE_SUBSURFACE)`, written against
the pre-port Epetra API (`Epetra_Vector`, `Epetra_Map`, `->SubVector()`,
`->Data()`, etc.). This is the single highest-value file to port: it directly
exercises a 2x2 `TreeOperator` with cross-mesh `Op_MeshInjection`/
`Operator_Schema` off-diagonal coupling blocks -- structurally the *exact*
shape of what `MPCCoupledWaterNonInjected`
(`src/physics/ats/src/pks/mpc/mpc_coupled_water_noninjected.cc`) now does in
production. That file is a working, already-ported, in-repo template for the
injection-map-construction and `Operator_Schema`/`Op_MeshInjection`/
`TreeOperator` wiring idiom.

### Per-test verdict

| Test | Location | Verdict |
|---|---|---|
| `SURFACE_SUBSURFACE_SUPERMAP` | lines 33-99 | **Mechanical rename, low effort.** |
| `SURFACE_SUBSURFACE_DIFFUSION` | lines 102-230 | **Moderate.** Every construct maps 1:1 to current API; `mpc_coupled_water_noninjected.cc` is a direct template to copy the injection-map-build and `Operator_Schema`/`Op_MeshInjection`/`TreeOperator` wiring from. No `CHECK`/`CHECK_EQUAL` assertions exist in the original (it only prints the matrix) -- a faithful port should add real assertions. |
| `SURFACE_SUBSURFACE_LOTS_OF_DIAGONALS` | lines 233-512 | **Moderate.** Same fixes as above, at larger scale (4-leaf tree). **This test already asserts that the assembled-matrix apply and the matrix-free apply agree** (lines 506-511) -- i.e. it already encodes the exact apply-vs-assembled consistency check that would have caught bug #2 above directly. Porting this one test is the single highest-value item in this whole survey. |
| `THREE_LEVEL_HIERARCHY` | lines 515-692 | **Blocked.** Depends on `PDE_CouplingFlux`, which does not exist anywhere in this checkout (`grep -rl "PDE_CouplingFlux" src/` finds only this test file). Not a rename target for `Op_MeshInjection`/`Operator_Schema` -- those are deliberately restricted (per their own doc comments) to single-scalar-per-row, one-to-one cross-mesh cell injection, whereas `PDE_CouplingFlux` supports arbitrary same-mesh index-list-based coupling (`std::vector<std::vector<int>>` row/col index lists), a materially more general capability. Porting this test requires writing new production code first, not just porting test code -- out of scope for a test-only port. |

### Construct mapping table (old -> current), for whoever does the port

| Old (Epetra-era) | Current equivalent | Evidence |
|---|---|---|
| `Epetra_Vector`, `Epetra_Map` | `Vector_type` / `Map_type` (via `Map_ptr_type`) | `src/utils/AmanziTypes.hh:169,177,182,193` |
| `->SubVector(i)` (TreeVector/TreeVectorSpace) | `->getSubVector(i)` | `TreeVector.hh:87,90`, `TreeVectorSpace.hh:72` |
| `->Data()` (TreeVector) | `->getData()` | `TreeVector.hh:96,99` |
| `->PutScalar(x)` | `->putScalar(x)` | `TreeVector.hh:110` |
| `Operator::get_row_map()`/`get_col_map()` | **Doesn't exist on `Operator`** -- use `getRowMap()` (row) / `getDomainMap()` (col); note asymmetry, there is no `getColMap()` | `Operator.hh:228,232,233` |
| `TreeOperator::get_row_map()`/`get_col_map()` | **Unchanged**, these do exist on `TreeOperator` | `TreeOperator.hh:79-80` |
| `op.DomainMap()` / `.RangeMap()` (capitalized) | `getDomainMap()` / `getRangeMap()` | `Operator.hh:228,232`, `TreeOperator.hh:70,74` |
| `TreeOperator(tvs)` / `TreeOperator(row_map, col_map)` | unchanged | `TreeOperator.hh:56-59` |
| `TreeVectorSpace(RCP<const CompositeVectorSpace>)` | ctor takes `RCP<const CompositeSpace>`; **no ctor from a raw Map exists (old or new)** -- `get_row_map()`/`getDomainMap()` on `Operator` already return `CompositeSpace`-shaped objects, so this is a rename, not a structural gap | `TreeVectorSpace.hh:37`; confirmed working in `mpc_coupled_water_noninjected.cc:114-115` |
| `Operator_Schema(row_map, col_map, plist, schema_row, schema_col)` | same shape, args are `RCP<const CompositeSpace>` | `Operator_Schema.hh:36-43`; live call in `mpc_coupled_water_noninjected.cc:155-157` |
| `Op_MeshInjection(schema_row, mesh_row, schema_col, mesh_col, injection, transpose)` | **unchanged constructor shape** | `Op_MeshInjection.hh:58-66`; live call in `mpc_coupled_water_noninjected.cc:158-160,174-176` |
| `Epetra_Map(-1, parents.size(), parents.data(), 0, *comm)` | `Map_type(-1, parent_gids_host /* Kokkos::View<GO*,HostSpace> */, 0, comm)` | `mpc_coupled_water_noninjected.cc:131-138` |
| `cell_map.GID(lid)` (Epetra) | `cell_map.getGlobalElement(lid)` (Tpetra), or prefer `mesh->getEntityGID(CELL, lid)` directly | `Operator_Schema.cc:50-51,91-92`; `mpc_coupled_water_noninjected.cc:136` |
| `global_operator()->Init()` (reset before refill) | **`Init()` no longer exists on `Operator`** -- use `->Zero()` | confirmed absent via full-file grep of `Operator.hh`; live replacement used in `operators_diffusion_discretizations.cc:200,257` (`fix.global_op->Zero()`) |
| `BCs(mesh, Entity_kind::FACE, WhetStone::DOF_Type::SCALAR)` | unchanged | `src/data_structures/BCs.hh:137-139` |
| `PDE_DiffusionFV(plist, mesh)` / `.Setup(K, k, dkdp)` / `.SetBCs(bc, bc)` | unchanged | `PDE_DiffusionFV.hh:64`, `PDE_Diffusion.hh:92-94` |
| `PDE_Accumulation(Entity_kind::CELL, mesh)` / `.AddAccumulationTerm(cv, "cell")` | unchanged | `PDE_Accumulation.hh:54,111` |
| `A()->Print(std::cout)` (`Epetra_CrsMatrix::Print`) | no direct equivalent on `Tpetra::CrsMatrix` -- use `Tpetra::MatrixMarket::Writer<Matrix_type>::writeSparseFile(...)`, or drop (it's diagnostic-only) | `mpc_coupled_water_noninjected.cc:363-364` |
| `A()->Apply(ones, result)` (capital `Apply`) | `Tpetra::CrsMatrix::apply(...)` (lowercase) | standard Tpetra rename |
| `supermap->getMap()` | unchanged | `SuperMap.hh:46` |
| `createSuperMap(tvs)`, `copyToSuperVector(map, tv, sv)` | unchanged, used identically in currently-building `TreeOperator.cc` | `SuperMap.hh:113-139`, `TreeOperator.cc:304-305` |
| `getNumTreeVectorLeaves(tvs)`, `Impl::collectTreeOperatorLeaves(...)` | unchanged, used identically in currently-building `TreeOperator.cc` | `DataStructuresHelpers.hh:181`; `TreeOperator.hh`'s `Impl` namespace |
| `TreeVector_Utils.hh` header | doesn't exist as a separate file -- `getNumTreeVectorLeaves` lives in `DataStructuresHelpers.hh`; `collectTreeOperatorLeaves` is declared in `TreeOperator.hh`'s `Impl` namespace | confirmed via grep/find |
| `TreeOperator::get_supermap()` | `TreeOperator` has no bare `get_supermap()`; use `get_row_supermap()` / `get_col_supermap()`, or `getSuperMap()` (alias for `get_row_supermap()`, valid when square) | `TreeOperator.hh:85` |
| `op.SetOperatorBlock(...)` / `.PutScalar(...)` (Stokes-style capitalization, N/A to this file but noted for `operators_stokes.cc`) | `set_operator_block(...)` / `putScalar(...)` | mechanical rename |

Meshing note: the surface-extraction call used throughout this test,
`meshfactory.create(mesh_subsurf, setnames, AmanziMesh::Entity_kind::FACE, true)`,
is **still fully supported** in the current `MeshFactory` API
(`src/mesh/mesh_factory/MeshFactory.hh:60-63`) and is actively used by other
(currently-disabled-but-otherwise-intact) tests
(`operators_nonlinear_surface.cc:125`, `operators_advdiff_surface.cc:75`,
`operators_diffusion_dfn.cc:77`, `operators_fractures.cc:76`,
`operators_parabolic_surface.cc:73`) and by production code
(`ats_mesh_factory.cc:244`). This is a different (and still valid) mesh-pairing
idiom than the two-top-level-State-mesh-entries pattern
`mpc_coupled_water_noninjected.cc`'s own test XML configs use -- both coexist
in the API; neither is stale.


## Priority 2: other zero-or-near-zero-Epetra-hit candidates (not TreeOperator-specific)

These were found by grepping `Epetra_` hit-counts across every file in
`src/operators/test/`. Not individually investigated as deeply as
`operators_tree.cc`, but flagged as cheap, high-signal candidates.

- **`operators_marshak_assembly.cc`, `operators_marshak.cc`** -- zero Epetra
  hits. Unlike almost everything else in `CMakeLists.txt`, these are gated
  behind `if(0) ... endif()` (`src/operators/CMakeLists.txt:157-171`) rather
  than fully commented out, strongly suggesting they build/ran recently and
  were just toggled off for some unrelated reason, not abandoned mid-port.
  **Cheapest possible win in this whole survey: just flip `if(0)` to `if(1)`
  and see what happens.** (Do include `Mesh_MSTK.hh` directly, an older
  pattern vs. generic `MeshFactory`, but that header still exists and is
  built in this checkout.)
- **`operators_evaluators_op.cc`, `operators_evaluators_op_inverse.cc`** --
  effectively Epetra-free; their only `Epetra_` references are two dead/unused
  `#include` lines each (`Epetra_MpiComm.h`, `Epetra_Vector.h` -- confirmed
  not used anywhere else in either file). Not referenced in `CMakeLists.txt`
  at all. Not `TreeOperator`-related (suite name is `EVALUATOR_ON_OP`).
- **`operators_upwind_manifolds.cc`** (`TEST(UPWIND_FLUX_MANIFOLDS)`) -- zero
  Epetra hits, modern-looking includes (`MeshExtractedManifold.hh`,
  `MeshAlgorithms.hh`). Directly relevant to the upwinding/mobility-derivative
  evaluator chain (`EvaluatorUpwindTotalFlux` et al.) investigated at length
  this session. Not in `CMakeLists.txt`.
- **`operators_diffusion_consistent_face.cc`** (`TEST(OPERATOR_DIFFUSION_MIXED)`)
  -- zero Epetra hits, uses `PDE_DiffusionMFD`. Not in `CMakeLists.txt`.
- **`operators_state_container.cc`** -- already fully modern/idiomatic
  (`Op_Cell_Cell`, `Op_Factory`, `State::Require<...>`, `Comm_ptr_type`,
  positional `MeshFactory::create(...)`). Good reference for idiomatic current
  test style generally. Not in `CMakeLists.txt`; narrow in scope (not
  `TreeOperator`-related).
- **`operators_diffusion_low_order.cc`** -- only Epetra-adjacent residue is a
  dead `#include "EpetraExt_RowMatrixOut.h"` line (an old debug-dump tool);
  otherwise clean.


## Not recommended / out of scope

- **`operators_stokes.cc`** -- only 1 raw `Epetra_` hit, but still heavily
  old-API in substance (`SetOperatorBlock`, `DomainMap()`/`RangeMap()`
  capitalized, `PutScalar`, `SubVector`/`Data()`, `InitBlockDiagonalPreconditioner`,
  a standalone `LinearOperatorGMRES<...>` template rather than the port's
  `AmanziSolvers::createInverse` factory pattern). Depends on
  `PDE_Elasticity`/`PDE_Abstract`, whose port status wasn't verified. Uses
  `TreeOperator` but for a same-mesh multi-PDE (Stokes velocity/pressure)
  block structure, not cross-mesh coupling -- lower relevance to the bugs
  found this session, and a heavier, cross-cutting rewrite. Deprioritized.
- **`operators_coupled_diffusion.cc`** -- 26 Epetra hits, plus uses the old
  `CompositeVector` call-operator indexing style (`u("face", 0, f)`), a
  fundamentally different access pattern than the port's `viewComponent(...)`.
  Depends on `AnalyticNonlinearCoupledBase.hh` (4 Epetra hits of its own).
  Heavy rewrite; deprioritized relative to `operators_tree.cc`.
- **`operators_fct.cc`** -- has a real (non-dead) `Epetra_MultiVector`
  allocation, not just stale includes. Flux-corrected-transport, unrelated to
  `TreeOperator`/preconditioning. Deprioritized.
- All other files in `src/operators/test/` were not individually surveyed in
  this pass; the Epetra-hit-count table below (full list, sorted ascending)
  is a reasonable starting point for extending this triage if more coverage
  is wanted later.


## Full Epetra-hit-count table (all `src/operators/test/operators_*.cc` and `operator_*.cc`)

```
0	operators_boundary_flux_solver.cc
0	operators_diffusion_consistent_face.cc
0	operators_diffusion_discretizations.cc   (already active/building)
0	operators_diffusion_low_order.cc
0	operators_marshak_assembly.cc            (if(0)-gated)
0	operators_marshak.cc                     (if(0)-gated)
0	operators_mini1D_nonlinear.cc
0	operators_mini1D.cc
0	operators_solvers.cc                     (already active/building)
0	operators_state_container.cc
0	operators_upwind_manifolds.cc
1	operator_elasticity_weak_symmetry.cc
1	operators_advdiff_surface.cc
1	operators_advection_dfn.cc
1	operators_elasticity.cc
1	operators_fractures.cc
1	operators_gravity.cc
1	operators_laplace_beltrami.cc
1	operators_parabolic_surface.cc
1	operators_stokes.cc
2	operator_diffusion_curved_face.cc
2	operator_diffusion_dfn_nlfv.cc
2	operators_advection_dg.cc
2	operators_commute.cc
2	operators_diffusion_curved.cc
2	operators_diffusion_dg.cc
2	operators_diffusion_edges.cc
2	operators_evaluators_op_inverse.cc       (both hits are dead includes)
2	operators_evaluators_op.cc               (both hits are dead includes)
3	operators_advdiff.cc
3	operators_diffusion_dfn_fv.cc
3	operators_diffusion_dfn.cc
3	operators_diffusion_fractured_matrix.cc
3	operators_diffusion_nodal.cc
3	operators_diffusion_nonsymmetric.cc
3	operators_diffusion_tpfa.cc
3	operators_electromagnetics.cc
3	operators_fct.cc
3	operators_limiters_dg.cc
3	operators_utils.cc
5	operators_marshak_logical.cc
5	operators_nonlinear_surface.cc
5	operators_remap_curved.cc
5	operators_stability.cc
6	operators_diffusion_divk.cc
6	operators_remap_gcl.cc
6	operators_remap.cc
6	operators_upwind.cc
7	operators_remap_gcl_reinit.cc
8	operators_diffusion_convergence.cc
8	operators_diffusion_high_order.cc
8	operators_magnetic_diffusion.cc
8	operators_reconstruction.cc
9	operators_diffusion_nlfv.cc
10	operator_reconstruction.cc
10	operators_tree.cc
11	operators_advection_dg_transient.cc
26	operators_coupled_diffusion.cc
34	operators_limiters.cc
```

(`operators_diffusion_discretizations.cc` and `operators_solvers.cc` are the
only two test source files currently registered and building in
`src/operators/CMakeLists.txt`, confirming the zero-hit heuristic above is a
real signal of "already ported," not coincidence.)
