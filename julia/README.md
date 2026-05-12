# SPM.jl — Julia port of SPM (pilot)

This subdirectory contains an in-progress Julia translation of selected
functions from the surrounding MATLAB SPM (Statistical Parametric Mapping)
toolbox. The intent is to remove the MATLAB runtime dependency while
preserving the SPM API and file organization.

This is a **pilot**. The original SPM is ~830,000 lines of MATLAB across
~4,850 `.m` files (plus 700+ C/MEX kernels in `src/` and large vendored
toolboxes in `external/`). A full translation is a multi-year, multi-expert
effort. This commit demonstrates the translation pipeline end-to-end on three
small, self-contained functions, and documents the conventions and the
empirical comparison method that any future translation should follow.

## Layout

```
julia/
├── Project.toml             Julia package manifest (deps: MAT, SpecialFunctions)
├── src/
│   ├── SPM.jl               umbrella module: `using SPM` brings exports
│   ├── spm_Npdf.jl          one function per file, filename = function name
│   ├── spm_Ncdf.jl          same as the parent SPM/MATLAB convention
│   └── spm_Tcdf.jl
├── test/
│   ├── runtests.jl          in-Julia property tests (fast, no Octave needed)
│   ├── properties/          mathematical-property tests
│   └── parity/              empirical parity vs the original MATLAB SPM
│       ├── gen_inputs.jl    Julia: generate adversarial + random input grids
│       ├── run_reference.m  Octave: run the original spm_*.m on those inputs
│       ├── check_parity.jl  Julia: diff Julia translation vs Octave reference
│       └── data/            generated .mat files (gitignored)
└── scripts/
    └── parity.sh            one-button driver: gen -> reference -> check
```

## Translation conventions

- **One function per file** matching the MATLAB filename
  (`spm_Npdf.m` → `spm_Npdf.jl`).
- **Docstring is a verbatim port** of the MATLAB header `%` block, edited
  only to mark the Julia API (defaults, types) and to record the source.
- **Default arguments** mirror MATLAB `if nargin < N` blocks via Julia
  default values (`spm_Npdf(x, u=0, v=1)`).
- **NaN-for-out-of-domain**: SPM returns NaN (with a warning) for inputs
  like `v <= 0`. The Julia port returns NaN silently. (A future change can
  add `@warn` but the pilot keeps the kernel pure.)
- **Size compatibility**: SPM allows scalar-or-matching-shape arguments
  but not NumPy-style singleton broadcasting between non-scalars. The
  port enforces the same contract via `_check_compatible_sizes`.

## Running the parity harness

```
bash julia/scripts/parity.sh
```

This:

1. Generates input grids in Julia (random + adversarial: 0, ±Inf, NaN,
   domain boundaries, algorithm-switch regions for `betainc`) and writes
   `julia/test/parity/data/inputs_*.mat`.
2. In **Octave**, loads those `.mat` files, calls the **unmodified**
   `spm_Npdf.m`, `spm_Ncdf.m`, `spm_Tcdf.m` from the SPM repo root, and
   writes outputs to `reference_*.mat`.
3. In Julia, runs the translations on the same inputs and reports the
   maximum absolute and relative error per case, exiting non-zero if any
   case exceeds the per-function tolerance.

Tolerances are per-function (see `check_parity.jl`):

| Function   | rtol  | atol   | Rationale                                          |
| ---------- | ----- | ------ | -------------------------------------------------- |
| `spm_Npdf` | 1e-14 | 0      | closed form via `exp`; bit-identical in practice   |
| `spm_Ncdf` | 1e-14 | 1e-15  | closed form via `erf`; 1-ulp slack in deep tail    |
| `spm_Tcdf` | 1e-9  | 1e-12  | iterative `betainc`; see "Known parity gaps" below |

## Property tests

`julia/test/properties/test_distributions.jl` checks invariants
**independent of the MATLAB reference**: symmetry of the Gaussian PDF,
`F(0) = 0.5` for symmetric CDFs, monotonicity, Cauchy/Normal limits of
Student-t, NaN propagation for invalid variance/df. Run with
`julia --project=julia julia/test/runtests.jl`.

## Known parity gaps

- **`spm_Tcdf` vs MATLAB**: Julia's `SpecialFunctions.beta_inc` and
  MATLAB's `betainc` produce values that agree to roughly 9–10
  significant figures (worst observed relative error ≈ 4×10⁻¹⁰ on a
  ~2,400-point grid). This is below the threshold of any neuroimaging
  inferential procedure (typical significance thresholds operate on
  log-scale `p` values), but it is **not** bit-identical. To tighten,
  call MATLAB's `betainc` algorithm directly (the SPM `src/` directory
  does not currently include a C implementation; a port of DiDonato &
  Morris's algorithm 708 would be the canonical source). The pilot does
  not pursue this.

- **`spm_Ncdf` deep tail**: For `|z| > ~9`, MATLAB's `erf` saturates to
  `±1` exactly; Julia's libm leaves 1 ulp of slack. The CDF then reads
  `5.55×10⁻¹⁷` instead of `0.0`. Sub-machine-epsilon; no consequence.

## What is *not* translated

The pilot deliberately stays inside pure-numerics territory. The
following are deferred or out of scope:

| Area                                    | Status        | Notes                                                                       |
| --------------------------------------- | ------------- | --------------------------------------------------------------------------- |
| `src/` (C/MEX kernels)                  | keep as C     | Call from Julia via `ccall`; do not rewrite                                 |
| `external/` (FieldTrip, etc.)           | leave alone   | Separate upstreams, license/provenance issues                               |
| `matlabbatch/`, `config/`               | re-design     | MATLAB UI/serialization framework; Julia needs a different design           |
| `spm.m` GUI, `spm_orthviews/`           | re-design     | MATLAB graphics; Julia would use Makie/Gtk/Web                              |
| `@nifti/`, `@gifti/`, `@meeg/`, `@xmltree/`, `@file_array/`, `@slover/` | future        | MATLAB OOP classes; Julia uses multiple dispatch — design pass needed       |
| `compat/`                               | drop          | MATLAB-version shims                                                        |
| `toolbox/` (DCM, DARTEL, MEEGtools)     | each separate | Each is effectively its own project                                         |
| Remaining `spm_*.m` numerics            | follow pilot  | Add file under `julia/src/`, register in parity harness, run                |

## Why Julia, not R or Python?

Per the project request, Julia is preferred. Specifically:

- **Numerics + dispatch**: SPM uses element-wise broadcasting and many
  helper functions dispatch on dimensionality / argument count. Julia
  multiple dispatch fits this idiom directly.
- **Calling existing C/MEX**: Julia's `ccall` is zero-overhead, useful
  for keeping `src/` kernels as-is.
- **Speed without rewriting in C**: relevant for the iterative
  estimators (EM, VB) in DCM/DEM.

R would have been chosen if Julia were unavailable; Python with
NumPy/SciPy is the third-choice fallback. None of those choices is made
in this pilot.

## Trustworthy MATLAB references consulted

When translating, the following references are used as the source of
truth — **not** memory:

1. The function's own `%` header comment in the original `.m` file
   (treated as a specification).
2. MathWorks function reference (`https://www.mathworks.com/help/matlab/ref/`)
   for built-ins (`betainc`, `erf`, `gammaln`, etc.).
3. The Octave manual for documented MATLAB-compat edge cases.
4. The SPM book (Friston et al., 2007) for statistical formulae.
5. **Empirical**: running the unmodified `.m` file in Octave on a wide
   adversarial input grid and comparing the output element-by-element.
   This is the final arbiter when docs disagree.

## Multi-agent plan summary (rationale)

The plan was developed under three perspectives in adversarial dialogue:

- **Statistician** (skeptical about numerics): demanded per-function
  tolerances, NaN-aware comparison, adversarial inputs at algorithm
  boundaries, and **property tests independent of the MATLAB reference**.
- **Beginner neuroimager** (cares about discoverability): demanded that
  function names, filenames, and the `spm_*` prefix remain unchanged.
- **Expert neuroimager** (cares about realism of scope): vetoed
  translating `matlabbatch/`, `external/`, the GUI, and the OOP classes
  in this pass.

The translation strategy implemented here reflects the intersection of
those concerns.
