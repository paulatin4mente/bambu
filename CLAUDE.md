# Bambu — Codebase Guide for AI Agents

## What bambu does

Bambu is an R/Bioconductor package for **multi-sample transcript discovery and quantification** from long-read RNA-Seq data (Oxford Nanopore, PacBio). It takes aligned BAM files plus a reference genome and annotation (GTF/TxDb), and outputs a `SummarizedExperiment` with expression estimates for both known and novel transcripts and genes.

The core challenge bambu solves: long reads are noisy (sequencing errors, alignment artefacts), so naive read-to-transcript assignment produces many false novel isoforms. Bambu addresses this through junction error correction, an XGBoost-based transcript scoring model, NDR (Novel Discovery Rate) control, and an EM algorithm for multi-mapping read assignment.

---

## Pipeline overview

```
BAM files
    │
    ▼
[1] prepareAnnotations      — convert GTF/TxDb into internal GRangesList format
    │
    ▼
[2] bambu.processReads      — per-sample: extract reads, correct junctions,
    │                         build read classes, score with XGBoost
    ▼
[3] bambu.extendAnnotations — cross-sample: combine read class candidates,
    │                         filter, assign NDR score, extend reference annotations
    ▼
[4] assignReadClasstoTranscripts — build equivalence classes, map read classes
    │                              to transcripts (distance-based)
    ▼
[5] bambu.quantify          — EM algorithm to resolve multi-mapping reads,
                              compute counts / CPM / fullLengthCounts
    │
    ▼
SummarizedExperiment output
```

The top-level orchestrator is `bambu()` in [R/bambu.R](R/bambu.R). All stages are called from there. Stages 2–5 can be run per-sample in parallel via `BiocParallel`.

---

## Module map

### Entry points (user-facing exported functions)

| File | Function | Purpose |
|------|----------|---------|
| [R/bambu.R](R/bambu.R) | `bambu()` | Main entry point; orchestrates all pipeline stages |
| [R/prepareAnnotations.R](R/prepareAnnotations.R) | `prepareAnnotations()` | Convert GTF/TxDb to `GRangesList` annotation object; must be called before `bambu()` |
| [R/readWrite.R](R/readWrite.R) | `writeBambuOutput()`, `importBambuResults()`, `writeToGTF()`, `readFromGTF()` | Save/load bambu results to disk (GTF + count matrices) |
| [R/plotBambu.R](R/plotBambu.R) | `plotBambu()` | Visualise expression and isoform structure |
| [R/compareTranscripts.R](R/compareTranscripts.R) | `compareTranscripts()` | Compare alternatively-spliced transcripts between query and subject |
| [R/transcriptToGeneExpression.R](R/transcriptToGeneExpression.R) | `transcriptToGeneExpression()` | Aggregate transcript-level SE to gene-level counts |

### Module 1 — Annotation preparation

| File | Key functions | Role |
|------|--------------|------|
| [R/prepareAnnotations.R](R/prepareAnnotations.R) | `prepareAnnotations()` | Loads GTF or TxDb, extracts exon-by-transcript GRangesList, adds metadata columns (`TXNAME`, `GENEID`, etc.) |
| [R/prepareAnnotations_utilityFunctions.R](R/prepareAnnotations_utilityFunctions.R) | internal helpers | Intron extraction, annotation metadata formatting |

### Module 2 — Read processing (per sample)

| File | Key functions | Role |
|------|--------------|------|
| [R/bambu-processReads.R](R/bambu-processReads.R) | `bambu.processReads()`, `bambu.processReadsByFile()`, `bambu.readsByFile()`, `constructReadClasses()` | Top-level per-sample dispatcher; iterates over BAM files; calls junction table construction, error correction, read class construction |
| [R/bambu-processReads_utilityCreateJunctionTables.R](R/bambu-processReads_utilityCreateJunctionTables.R) | `isore.constructJunctionTables()`, `createJunctionTable()`, `junctionStrandCorrection()` | Extract splice junctions from alignments; compute strand scores from splice motifs |
| [R/bambu-processReads_utilityJunctionErrorCorrection.R](R/bambu-processReads_utilityJunctionErrorCorrection.R) | `junctionErrorCorrection()`, `fitXGBoostModel()`, `findHighConfidenceJunctions()` | Correct noisy splice junctions using an XGBoost model trained on annotated junctions; core noise-reduction step |
| [R/bambu-processReads_utilityConstructReadClasses.R](R/bambu-processReads_utilityConstructReadClasses.R) | `isore.constructReadClasses()`, `constructSplicedReadClasses()`, `constructUnsplicedReadClasses()`, `assignGeneIds()` | Group corrected reads into read classes (unique exon-chain signatures); assign gene IDs by overlap with annotation |
| [R/bambu-processReads_scoreReadClasses.R](R/bambu-processReads_scoreReadClasses.R) | `trainBambu()` + scoring functions | XGBoost-based transcript scoring; `trainBambu()` is exported for training a custom model on a different species/dataset |
| [R/prepareDataFromBam.R](R/prepareDataFromBam.R) | `prepareDataFromBam()` | Low-level BAM → GRangesList reader; handles CIGAR parsing, clipping, chunked reading via `yieldSize` |

### Module 3 — Annotation extension (cross-sample)

| File | Key functions | Role |
|------|--------------|------|
| [R/bambu-extendAnnotations.R](R/bambu-extendAnnotations.R) | `bambu.extendAnnotations()` | Orchestrates cross-sample combination then extension; calls combine then extend |
| [R/bambu-extendAnnotations-utilityCombine.R](R/bambu-extendAnnotations-utilityCombine.R) | `isore.combineTranscriptCandidates()`, `combineSplicedTranscriptModels()`, `extractFeaturesFromReadClassSE()` | Merge read class candidates across samples into a unified set; compute per-position start/end weighted medians |
| [R/bambu-extendAnnotations-utilityExtend.R](R/bambu-extendAnnotations-utilityExtend.R) | `isore.extendAnnotations()`, `recommendNDR()`, `setNDR()`, `filterTranscriptsByAnnotation()`, `calculateDistToAnnotation()` | Filter candidates by NDR threshold; merge with reference; assign novel transcript IDs; `setNDR()` is exported to re-threshold post-hoc |

### Module 4 — Read class to transcript assignment

| File | Key functions | Role |
|------|--------------|------|
| [R/bambu_utilityFunctions.R](R/bambu_utilityFunctions.R) | `assignReadClasstoTranscripts()` (via `calculateDistTable()`), `combineCountSes()`, `generateColData()`, `checkInputs()`, `setBiocParallelParameters()` | Compute distance table mapping read classes to transcripts; build equivalence classes; validate inputs; manage parallelism |
| [R/bambu-assignDist.R](R/bambu-assignDist.R) | `assignReadClasstoTranscripts()`, `generateUniqueCounts()`, `generateNonUniqueCounts()`, `generateIncompatibleCounts()` | Build count matrices (unique, non-unique, incompatible) from read class equivalence classes |
| [R/bambu-quantify_utilityFunctions.R](R/bambu-quantify_utilityFunctions.R) | `genEquiRCs()`, `modifyIncompatibleAssignment()`, `processIncompatibleCounts()` | Generate equivalence read classes; handle reads incompatible with all transcripts in a gene |

### Module 5 — Quantification

| File | Key functions | Role |
|------|--------------|------|
| [R/bambu-quantify.R](R/bambu-quantify.R) | `bambu.quantify()`, `bambu.quantDT()` | Run EM per sample; aggregate counts; compute CPM; assemble sparse output vectors |
| [src/em.cpp](src/em.cpp) | `em_theta()` (C++/Rcpp) | Expectation-Maximisation algorithm implemented in C++ with RcppArmadillo; resolves multi-mapping reads to transcripts using a probability matrix |

### Module 6 — Visualization, comparison & output

| File | Key functions | Role |
|------|--------------|------|
| [R/readWrite.R](R/readWrite.R) | `writeBambuOutput()`, `importBambuResults()`, `writeToGTF()`, `readFromGTF()` | Serialise/load bambu results (GTF + count matrices) |
| [R/plotBambu.R](R/plotBambu.R) | `plotBambu()` | Visualise expression and isoform structure |
| [R/plotBambu_utilityFunctions.R](R/plotBambu_utilityFunctions.R) | internal helpers | Plot construction helpers for `plotBambu()` |
| [R/compareTranscripts.R](R/compareTranscripts.R) | `compareTranscripts()` | Compare alternatively-spliced transcripts between query and subject |
| [R/compareTranscripts_utilityFunctions.R](R/compareTranscripts_utilityFunctions.R) | internal helpers | GRanges manipulation utilities for transcript comparison |
| [R/transcriptToGeneExpression.R](R/transcriptToGeneExpression.R) | `transcriptToGeneExpression()` | Aggregate transcript-level SE to gene-level counts |
| [R/transcriptToGeneExpression_utilityFunctions.R](R/transcriptToGeneExpression_utilityFunctions.R) | internal helpers | Helper functions for gene-level aggregation |

### Shared utilities

| File | Purpose |
|------|---------|
| [R/utility_spliceHelper_functions.R](R/utility_spliceHelper_functions.R) | Distance-based splice overlap finding (`findSpliceOverlapsByDist()`); used in annotation extension and assignment |
| [R/globals.R](R/globals.R) | `globalVariables()` declarations to suppress R CMD check notes for data.table/dplyr column names |
| [R/RcppExports.R](R/RcppExports.R) | Auto-generated Rcpp bindings — do not edit manually |

---

## Key data structures

- **`bambuAnnotation` / `GRangesList`** — exons grouped by transcript (`names` = transcript IDs); `mcols` carries `TXNAME`, `GENEID`, `NDR`, `newTx` flags. Created by `prepareAnnotations()`.
- **Read class `SummarizedExperiment`** — one row per read class, one column per sample; `rowRanges` holds the exon-chain GRanges; `assays` include `counts` (read counts per class); `metadata` carries `countMatrix`, `incompatibleCountMatrix`, `readClassDist`, `eqClassById`.
- **`readClassDt` (data.table)** — equivalence-class table with columns `txid`, `eqClassId`, `eqClassById`, `gene_sid`, `nobs`, `multi_align`, `aval`; fed directly to `em_theta()`.
- **`distTable`** — mapping of read class to transcript(s) with distance scores; computed by `calculateDistTable()` and stored in `metadata(readClassList)$readClassDist`.
- **Output `SummarizedExperiment`** — `rowRanges` = extended annotations; `assays` = `counts`, `CPM`, `fullLengthCounts`, `uniqueCounts`; `metadata` = `incompatibleCounts`, `nonuniqueCounts`, optionally `readToTranscriptMap` and `distTable`.

---

## Important parameters

- **`NDR`** (Novel Discovery Rate) — analogous to FDR for novel transcripts; default is auto-recommended based on `baselineFDR = 0.1`; can be re-applied post-hoc with `setNDR()`.
- **`opt.discovery`** — list controlling isoform reconstruction (min read counts, min sample number, distance thresholds); see `setIsoreParameters()` in [R/bambu_utilityFunctions.R](R/bambu_utilityFunctions.R).
- **`opt.em`** — list controlling EM convergence (`maxiter`, `conv`, `sig.digit`); see `setEmParameters()`.
- **`rcFiles`** — read class files (`.rds`) saved mid-run; allow skipping `processReads` on re-runs.

---

## Code review guidelines

Standard review concerns apply (correctness, tests, clarity, etc.). In
addition, reviewers should check the bambu-specific items below. Extend this
list as shared conventions emerge.

### Bambu-specific review checklist
- **Is `CLAUDE.md` still accurate given this change?** If the PR adds,
  removes, renames, or moves a file under `R/`, or materially changes what a
  module does, flag the `CLAUDE.md` sections that would need updating (module
  map tables, pipeline overview, key data structures, important parameters).
  Suggest the wording in a review comment; do not edit `CLAUDE.md` as part of
  the review.

---

## On-request annotation passes

The two passes below are **opt-in**: only run them when the user explicitly
asks ("tag code issues", "add function headers", etc.). They produce noisy
diffs and are not part of normal development.

### Pass 1 — Tag code issues with typed `TODO:` comments
Surface possible bugs, dead code, and poor naming inline. Do **not** fix the
issues — just flag them. Place each comment on the line **above** the flagged
code. Use one of these typed prefixes so findings are greppable:

- `# TODO: [BUG] ...` — possible bugs (off-by-one, NULL deref, wrong operator
  precedence, vacuous `all()` on empty input, etc.)
- `# TODO: [UNUSED CODE] ...` — unreachable code after `return()`, unused
  variables, functions defined but never called, large commented-out blocks
- `# TODO: [POOR NAMING] ...` — generic names (`temp`, `tmp`, `x`, `ov`,
  `length_tmp`); the same single letter reused for different concepts
- `# TODO: [OTHER] ...` — magic numbers, debug `print()`s, typos, misleading
  comments

### Pass 2 — Add function header blocks with module & call counts
Insert a header immediately above every named top-level function definition
(above the `#'` roxygen block if one exists). Format:

```r
# --- functionName ---
# Module: [module name] | [filename.R]
# Called by: caller1.R, caller2.R
# Call count: N calls, M files
```

For functions with no internal callers, use
`# Call count: 0 internal calls (exported or not called internally)` and
`# Called by: (not called anywhere — user-facing entry point)` for exported
ones.

**Call-graph accuracy** — when counting callers, **exclude matches inside
roxygen `#'` blocks, regular comments, and string literals**. A naïve
`grep 'funcName('` over-counts every `@seealso`, `@examples`, and
`message("… funcName(…) …")`. Filter lines whose first non-whitespace
character is `#` before counting, and sanity-check matches preceded by an
unclosed `"` on the same line. Exported entry points (`bambu`, `plotBambu`,
`writeBambuOutput`, `readFromGTF`, `importBambuResults`, `compareTranscripts`,
`prepareAnnotations`, `trainBambu`, `setNDR`, `writeToGTF`,
`transcriptToGeneExpression`) are user-facing — their roxygen examples and
`@seealso` references are not call sites.

---

## Development notes

- **R package conventions**: uses roxygen2 for documentation; internal functions are tagged `@noRd`. Run `devtools::document()` after changing roxygen headers.
- **C++ code**: `src/em.cpp` uses RcppArmadillo; rebuild with `devtools::compileAttributes()` then `devtools::build()`. `R/RcppExports.R` and `src/RcppExports.cpp` are auto-generated.
- **Tests**: located in [tests/testthat/](tests/testthat/); run with `devtools::test()`. Test data (small chr9 region) is in `inst/extdata/`.
- **Bioconductor branch**: `devel_pre_v4` is the active development branch; `devel` is the main integration branch used for PRs.
- **Parallelism**: `BiocParallel` is used for multi-sample processing; `bpParameters` is configured in `setBiocParallelParameters()` based on `ncore` and platform.
