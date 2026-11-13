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

## The run and the verdict

- The run: 8 cells x 200 replicates in 1,622 seconds on three workers (`pilot/results/run_info.rds`), exit 0, no restart.
- `Rscript pilot/verdict.R`: `piloted-close-call` (`pilot/results/verdict.rds`); readings in `pilot/decisions.md`, D1.
- The report rendered; `check-report.py`: no errors or warnings. A Markdown copy, `pilot/pilot-report.md`, was rendered at the user's request ("let the file be md"), with its three figures committed beside it.

## After the first pilot's stop

- D-S1 (2026-11-13). The user said the plan should have held three aims: pregnancy stage, disease activity at each visit, and whether edges change differently in SLE than in RA across pregnancy. The agent's own description of idea 1 had named all three; `idea.md` and the proposal kept only the first. The agent's error, recorded in proposal section 12.1.
- D-S2 (2026-11-13). The user chose a re-pilot on the three aims ("go"), and approved the counts download (25 MB). Defaults taken: RA as one group; all pregnancy samples against all non-pregnant ones.
- The counts were downloaded (`analysis/00_download.R`, not committed) and the factor change measured (`analysis/01_factor_shift.R`, `results/factor_shift.rds`); no edge was computed on real data.
- The revision (proposal section 12) has its own stop before pilot-b, as the agent promised the user.
- D-S3 (2026-11-13, during the pilot-b run). The user: one paper, not two. The pregnancy work is to become a second application in the first project's paper (`work-flow-for-lioness-paper`), each dataset analysed on its own, never merged. No change to the first repository until the user says so; how the pregnancy part enters that paper (a protocol amendment there) is proposed after pilot-b's verdict.

## Pilot-b

- Pre-registration 927c000; code 3924b81; timing 52 s for 10 replicates per cell; amendment to 200 replicates 807b6cb, pushed before the run.
- The run: 10 cells x 200 replicates in 965 seconds on three workers, exit 0. The machine had restarted once in the session before the run; the run itself completed in one pass.
- `Rscript pilot-b/verdict.R`: `piloted-close-call` (`pilot-b/results/verdict.rds`); readings in `pilot-b/decisions.md`, D1. Every deciding clause passed (1(b), 3(a), 3(b)); 2(b) and 3(c) failed; 2(a) and 2(c) close.
- The report rendered with no checker errors; a Markdown copy `pilot-b/pilot-report.md` with its figures.
- One search for other pregnancy datasets with SLE and RA (web search, 2026-11-13): none beyond GSE235508; Hong 2019 (SLE and healthy only) and a 2025 monocyte study (small) noted.

## Pilot-c

- The user: "go to C" (2026-11-13), straight on without a separate stop (D-S4). Revision 2 (proposal section 13): every visit, marker-gene cell mix, activity within woman; pass rule unchanged. The "fewer edges" change the agent had proposed was dropped before pre-registration: each share counts pairs tested one by one, so the set tested does not change any pair's result.
- Marker-gene shift measured (`analysis/02_factor_shift_markers.R`); pre-registration 430b0d4; code and tests; timing 60 s for 10 replicates per cell; amendment to 400 replicates b23e94d, pushed before the run.
- The run: 10 cells x 400 replicates in 2,233 seconds, one pass. Verdict `piloted-close-call` (`pilot-c/results/verdict.rds`); readings in `pilot-c/decisions.md`, D1; the unequal-groups account of the small artifact in D2.
