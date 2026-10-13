---
name: amanzi-style-guide
description: Amanzi/ATS C++ coding style conventions (naming, formatting, file layout). Load before writing or editing C++ code in this repository.
---

# Amanzi/ATS Style Guide

Source of truth: https://github.com/amanzi/ats/wiki/StyleGuide (work in progress; only
covers things NOT handled by clang-format).

## Formatting

- The codebase is `clang-format`-formatted, but **only a specific major version produces
  correct output against this repo's `src/.clang-format`**. `tools/formatting/clang-format.sh`
  records the required version in `EXPECTED_CLANG_FORMAT_VERSION` (as of this writing: 19).
  Before running `clang-format` in this repo, for any reason, check the installed version
  against that script's expectation:

  ```
  clang-format --version
  grep EXPECTED_CLANG_FORMAT_VERSION tools/formatting/clang-format.sh
  ```

  **If the major version doesn't match exactly, do NOT run `clang-format` at all** — not
  even on a single file. A mismatched version (confirmed: v23 against a v19-authored style
  file) does not just produce slightly different whitespace — it actively collapses things
  the style file is written to keep multi-line, e.g. turning a multi-line template
  parameter list into one line, or collapsing an out-of-line short function body onto a
  single malformed line. These are real behavior changes in newer clang-format releases,
  not configuration bugs to "fix" by editing `.clang-format`. If the version doesn't match,
  hand-write the correct multi-line form for just the lines you touched instead, and tell
  the user the installed clang-format is the wrong version rather than installing a
  different one yourself or leaving mis-formatted output in place.

- When the version does match, format only the specific files you added or edited, by
  running `clang-format -i` directly on each one, followed by the repo's short-block
  post-pass (which exists specifically to correct the out-of-line short-function-body
  collapse described above):

  ```
  clang-format -i path/to/FileYouEdited.cc path/to/FileYouEdited.hh
  python tools/formatting/format_short_blocks.py path/to/FileYouEdited.cc
  python tools/formatting/format_short_blocks.py path/to/FileYouEdited.hh
  ```

- **Do NOT** source/run `tools/formatting/clang-format.sh` itself (or any variant that
  walks a whole directory tree, e.g. `cd $AMANZI_SRC_DIR/src && . ../tools/formatting/clang-format.sh`)
  to format a change, even when the version matches. That script reformats every file
  under `src/` — in a repo with pre-existing uncommitted work, this silently rewrites
  hundreds of unrelated files (including other in-progress, uncommitted edits) in the
  same pass as your intended change, with no way to cleanly separate "the formatter's
  reflow" from "someone's real edit" after the fact if there's no commit to diff against.
  If a repo-wide reformat is ever genuinely wanted, that is the user's call to make and
  run themselves, not something to do as a byproduct of finishing an edit. Use the script
  only to read its `EXPECTED_CLANG_FORMAT_VERSION` value — never execute it.
- Before formatting, confirm the target file list is exactly the files you touched —
  nothing broader — and prefer passing explicit paths over any wildcard/recursive
  invocation.

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

## Comments

Default to no comment. Only write one when the code alone would leave a future reader
missing something they need — a hidden constraint, a subtle invariant, a non-obvious
reason, a pre/post-condition that isn't visible at the call site. If removing the comment
would not confuse that reader, remove it.

- **Don't restate the API.** If a signature, type, or name already says it, a comment
  repeating it in prose adds nothing. (`// Returns the number of DoFs` above
  `int getNumDofs() const` is noise.)
- **Don't document implementation details of a dependency the caller doesn't need to
  know.** A class's own `.hh` should describe its contract, not the library underneath it
  (e.g. `FileHDF5`'s public comments should not explain ASCEMIO's internal limitations —
  only the consequences that actually constrain the caller, like "views must be
  contiguous and LayoutRight").
- **Don't narrate development history.** No "originally this was X, refactored to do Y",
  no references to a design plan, a report, or a prior version of the code. A comment
  describes the code as it is now; git history already owns the as-it-was.
- **Don't describe paths not taken.** No "an alternative would be...", "a real fix would
  add...", "this doesn't do X". Comments document what the code does, not what it
  doesn't do or could have done instead.
- **Keep long comments for subtle, high-stakes details only.** A multi-line comment
  earns its size by explaining something whose consequences are easy to get wrong
  (why a broadcast must follow a rank-0-only read to avoid a deadlock; why a 1e-6 nudge
  is needed to avoid landing on the wrong side of a roundoff boundary). A multi-line
  comment that just walks through obvious code is a sign the comment should be cut, not
  trimmed.
- **Don't say what's self-evident from the code itself.** `// this code is SIMD` above
  code that obviously runs identically on every rank is empty. State the actual contract
  instead, e.g. "every rank must pass an empty view, or every rank must pass one already
  sized to match."
- **Say it once.** If the same design point applies to several overloads or several
  files, state it fully in one place and have the others reference it ("see the unpaged
  overload's comment above") rather than repeating the full explanation each time.
- **A short comment on a short, non-obvious helper is fine** — e.g. naming what a oneline
  predicate decides when the name alone doesn't carry it — but don't pad it with the
  reasoning for why the helper's behavior was chosen that way; that belongs in the
  class-level doc if anywhere.

## Documenting input specs

- Any class whose constructor (or a factory function building it) takes a
  `Teuchos::ParameterList` and reads more from it than just passing it along must
  document every parameter it reads, in that class's own `.hh` file (the
  `.. admonition:: ...-spec` block). This includes parameters with defaults — document
  the default too.
- When adding a new parameter to an existing object, or wiring an existing parameter
  through to a new object, update that object's spec doc in the same change. Don't leave
  a parameter the code accepts undocumented, and don't document a parameter the code
  doesn't actually read (check the constructor actually consumes it — a parameter that's
  accepted by a helper class but never read from this object's own plist isn't really
  supported here, and should either be wired through for real or left out of the doc).

## Tests

- Test files are named lower_case_with_underscores (breaking the CamelCase class-naming
  rule on purpose), based on the library/subcomponent under test — e.g. `state_dag.cc`
  tests the DAG in the `state` library, `mesh_geometry.cc` tests geometry in `mesh`.
- **Always prefix the test file (and its registered test name) with the library it
  belongs to**: `<library>_<subcomponent>.cc`, e.g. `utils_file_hdf5.cc`/
  `output_inputoutput_hdf5.cc`, not `file_hdf5.cc`. This is what lets you tell which
  library a test exercises from its name alone, and keeps it correct if a file ever moves
  libraries (e.g. `src/output/test/output_file_hdf5.cc` -> `src/utils/test/utils_file_hdf5.cc`
  when `FileHDF5` itself moves from `output` to `utils` — the file moves with its code,
  and is renamed to match its new library, not left with its old prefix).
