# actflower Vision

## What We Are Building
`actflower` is the R implementation of activity-flow modeling that should feel inevitable in hindsight: fast, exacting, modular, and pleasant to use.

We are not doing a sentimental translation of Python into R syntax. We are building the version that should have existed first.

## Standard
If the original ActflowToolbox authors ever see this package, the reaction target is:
"This is annoyingly better."

That means:
- Faster on real workloads.
- More correct on edge cases.
- Better tested than the source implementation.
- Cleaner interfaces and sharper failure messages.
- Internals that are understandable without archaeology.

## Non-Negotiables
- `actflower` uses `hdf5r` for HDF5 IO.
- We do not use `rhdf5`.
- We do not preserve inefficient loops just because they existed before.
- We will apply linear algebra shortcuts aggressively when mathematically equivalent.
- We will reuse proven code from `~/code/ariadne` when it is better than rewriting, and we will harden it under actflower contracts.

## Engineering Posture
- Correctness first, then speed, then ergonomics. No compromises on the first two.
- Every exported function gets strong validation and tests.
- Every heavy path gets benchmarked.
- Every stochastic path is seed-controlled and reproducible.
- Every module has one job.

## Where We Intend To Beat The Baseline
- Nodewise multiregression FC estimated through precision-matrix identities rather than repeated model fitting.
- Batched matrix/tensor prediction paths replacing per-target loops.
- Deterministic and inspectable CV logic for lasso/glasso/PC methods.
- Better handling of mask-based/non-circular workflows and missing external tooling.

## What "Done" Looks Like
A researcher can move from raw HDF5 arrays to actflow results with a compact, predictable API; CI proves numerical parity where intended, documents intentional deviations where improvements exist, and benchmark gates show we are materially faster.

If we find obvious opportunities to make the method more efficient or more correct, we implement them without apology.
