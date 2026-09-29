---
name: amanzi-style-guide
description: Amanzi/ATS C++ coding style conventions (naming, formatting, file layout). Load before writing or editing C++ code in this repository.
---

# Amanzi/ATS Style Guide

Source of truth: https://github.com/amanzi/ats/wiki/StyleGuide (work in progress; only
covers things NOT handled by clang-format).

## Formatting

- The codebase is `clang-format`-formatted. Run the repo's formatting script on any file
  you add or edit before considering the change done:

  ```
  cd $AMANZI_SRC_DIR/src
  . ../tools/formatting/clang-format.sh
  ```

- Everything not covered by clang-format (naming, file layout) is below. These are
  guidelines, not hard rules — the existing codebase does not universally follow them.
  Prefer to follow them for new code; when editing existing code, match the surrounding
  convention unless you are doing a broad refactor of that code, in which case bring it
  into line.

## Naming

- **Classes, structs, enums, typedefs, template parameters, namespaces**: `CamelCase`
  (e.g. `State`, `MeshFunction`, `SolverNKA`). Should be nouns. When inheriting, prefer to
  prefix with the base class name (`SolverNKA` : `Solver`). Underscores may separate a
  modifier, e.g. `Mesh_Helpers`.
- **Enums**: `enum class`, named with a `_kind` suffix (e.g. `Entity_kind`,
  `Parallel_kind`). Provide a `std::string to_string(Entity_kind)` and an
  `Entity_kind createEntityKind(const std::string&)` in the same namespace.
- **Functions** (member and free): `lowerCamelCase()`, always a verb (`getMyThing()`,
  `doSomethingImportant()`). This convention is recent and most existing code doesn't
  follow it. Which style to use depends on the size of the change:
  - **Small, localized edit** (tweaking one function, adding a small helper to an
    existing old-style class): match the existing surrounding style. Don't rename one
    function to `lowerCamelCase` in the middle of a class that's otherwise
    `UpperCamelCase` — that's pure churn without fixing the whole unit.
  - **Significant refactor of an existing unit** (rewriting a whole class, splitting a
    file, materially changing a class's structure): bring that unit's functions into
    `lowerCamelCase()` as part of the refactor.
  - **Brand new code** (a new class, new file, new namespace, new free functions with no
    prior existing-style anchor): DEFINITELY use `lowerCamelCase()` from the start, even
    if it lives alongside/replaces old-style code in the same header or is introduced as
    part of a branch that also touches old-style files. "New from the perspective of
    this file/branch" counts as new, even if the surrounding project still has plenty of
    old-style code elsewhere.
- **Variables**: lower case with underscores (`my_mesh`, `surface_mesh`,
  `my_member_variable`).
- **Private/protected members and functions**: trailing underscore, always — e.g.
  `my_mesh_`, `getMyPrivateThing_()`. This applies to private/protected free-standing
  helper methods on a class just as much as to data members.

## Files

- Name the file after the class it implements: `Mesh.hh` declares `Mesh`, `Mesh.cc`
  implements it.
- Templated/header-only classes: split into `Class_decl.hh` (declaration),
  `Class_impl.hh` (implementation), and `Class.hh` (includes both).
- Header guards: prefer `#pragma once` in new code (`#ifdef`/`#define` guards are the
  older style, still present in old files).

## Tests

- Test files are named lower_case_with_underscores (breaking the CamelCase class-naming
  rule on purpose), based on the library/subcomponent under test — e.g. `state_dag.cc`
  tests the DAG in the `state` library, `mesh_geometry.cc` tests geometry in `mesh`.
