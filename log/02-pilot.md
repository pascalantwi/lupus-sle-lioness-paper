# Log, stage 2: project spec and pilot

Workflow 0.2.4-md. Started 2026-11-13.

## Setup

- The user's request (2026-11-13): redo the pilot for idea 1 in this repository and report when it passes. Read as approval to run stage 1 and stage 2 through to the verdict, and to download the GSE235508 series matrix (16 KB, GEO) for the data description. The counts (25 MB) are not downloaded; their approval is an open question (proposal section 9).
- `tools/init.R` was run before stage 1 (package `lionesspreg`); `renv::restore()` ran clean; no Python.

## What was read

- `workflow/02-pilot.md`, `workflow/rules.md`, `templates/pilot-report/README.md`.
- The first project's pilot-c code, pre-registration, verdict script, and report, as the model for this pilot; the LIONESS edge functions and their tests were copied, the mean-shift helper left out.
- The GSE235508 series matrix: the codes of `samplegroup`, `timepoint`, and `pregnancyid`, as in the series record opened at stage 1.

## Steps and timing

| Step | What | Commit |
|---|---|---|
| 1 | `project-spec.md` | 69c7beb |
| 2 | `pilot/preregistration.md`, committed alone | df3e6bb |
| 3 | code, tests (37 expectations), `devtools::check()` 0 errors, 0 warnings, 0 notes; data description | 65893a8 |
| 3 | timing run: 10 replicates per cell, three workers, 84 s; amendment to 200 replicates | 2963242 |
| 4 | the run: `Rscript pilot/pilot.R --reps 200`, three workers | after the run |

## Decisions

The numbered decisions are in `pilot/decisions.md`.

## Process feedback

- The roxygen markdown parser read an inline formula starting with `` `r `` as R code; the formula now starts with `e =`.
