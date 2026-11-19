# Handoff: BGC datasets, frozen features and baseline clustering

Last updated: 2026-11-19 (Europe/Paris).

This file preserves the discussion and investigation for a new Codex conversation.
The user communicates in Chinese. Explain decisions and results in Chinese; use
clear, reproducible documentation and explicit provenance for experiments.

## Objective and current authorization

The user is developing candidate-independent BGC embeddings for similarity
retrieval and auxiliary clustering. Their biological priorities are biosynthetic
function, product relationships, and module composition, rather than organism or
coarse biosynthetic class alone.

The research objective is to evaluate protein and DNA embeddings using the nine
manually curated GCF datasets distributed with BiG-SCAPE 2.0, then test whether
DNA adds useful information to protein representations. The current project
stage includes completed B/E BASELINE CLUSTERING and G/H coordinate-repaired
input preparation, plus G/H encoder-interface verification and server orchestration.
H's existing Geneious aSDomain annotations have now been recovered in a new
processed snapshot: 2363 restored, two UOC915 annotations remain needs_review.
The previous H queue passed CPU parsing while silently omitting these annotations;
it must not be treated as domain-complete input.
The current authorized follow-up is the server-side bridge from recovered H (and
optionally repaired G) to the existing extraction workflow. All 90 H samples,
including UOC915, remain included with explicit known-incomplete annotation
status. Local work implements and CPU-tests this bridge; the user will prepare
the formal queue and run the models on the server.
B/E server feature extraction is complete and the six exports
have been returned locally. Preparation of A/C/D/F/I remains pending; G/H model
extraction and clustering have not been run.

- On 2026-11-09 the user authorized the first dataset-preparation step: create
  `data/raw/` and `data/processed/`, move original data into raw, and establish a
  traceable sample inventory. This supersedes the earlier documentation-only
  scope. The user subsequently requested continuing this same task after a pause.
- Raw relocation and the initial inventory are now implemented. No original
  sequences/labels were repaired, filtered, or deduplicated; no final evaluation
  cohort spanning all nine datasets has been selected.
- On 2026-11-10 the user selected B and E and authorized server-call scripts for
  frozen ESM2, BGC-MAC and DNA embeddings. This supersedes the earlier restriction
  to inventory work for these two datasets. Implemented entrypoint:
  `scripts/pilot_features.py` (`prepare`, `extract`, `validate`), documented in
  `docs/pilot_features.md`. This step did not launch a server/GPU run; the user
  subsequently completed extraction and returned the outputs, as recorded below.
- On 2026-11-11 the user reported completing server extraction, returned the B/E
  outputs, and authorized two clustering scripts: original cosine + average
  linkage and original cosine + kNN + Leiden. Live inspection found the files
  in `data/processed/features_be_v1/` (not `data/feature/`). All six native exports
  passed provenance/coverage validation, including the extraction run's artifact
  hashes and shared ESM/MAC token-cache provenance.
- The user subsequently requested B/E V-measure comparison across these routes
  and BiG-SCAPE. `scripts/summarize_vmeasure.py` now rechecks existing memberships
  and summarizes them alongside author-published BiG-SCAPE outputs; no local
  BiG-SCAPE run or new clustering sweep is claimed.
- On 2026-11-12 the user authorized a PCA preprocessing experiment with
  average-linkage on the existing B/E features, explicitly excluding further
  Leiden tests for this step. This supersedes the design-only PCA scope below
  for these six feature exports. Entry point: `scripts/cluster_pca_average.py`;
  instructions: `docs/pca_average.md`.
- On 2026-11-19 the user authorized repairing G/H coordinate annotations in
  processed copies, explicitly preserving raw files as backups. This is an
  input-preparation task, not authorization to launch G/H encoding or clustering.
  Entry point: `scripts/repair_gh.py`; evidence policy: `docs/repair_gh.md`.
- The user subsequently requested checking encoder compatibility to measure
  ESM2, BGC-MAC and DNA on G/H using the same workflow as B/E. CPU parsing passes
  all 101 repaired inputs. `pilot_features.py prepare-repaired` now bridges the
  repair snapshot to the existing extraction/validation workflow; instructions
  are in `docs/gh_features.md`. No G/H model inference or clustering has run.
- The user then authorized restoring H's existing aSDomain annotations and
  explicitly deferred Pfam-coordinate splitting. Entry point:
  `scripts/reconstruct_h.py`; evidence policy: `docs/reconstruct_h.md`.
  Preserve raw and all previous processed snapshots; do not deduplicate recovered
  domains or add new Pfam-derived cuts as part of this task. The old H queue must
  not be used for the planned domain-based benchmark without this recovery and
  newly validated prepared inputs.
- The user subsequently authorized server preparation and extraction support via
  `pilot_features.py prepare-reconstructed`, preserving the remaining UOC915
  review state rather than deleting that sample. This supersedes the previous
  lack of a recovery-schema adapter: prepare new server inputs and keep annotation
  completeness separate from preparation/inference success. Do not run local GPU
  inference, publish a formal local queue or alter the GBKs again as part of
  this follow-up; none of those actions was requested.
- The user subsequently requested pushing these changes to GitHub. Commit and
  push the H recovery scripts, server adapter, tests and documentation to the
  existing BGC_retrieve repository. Keep data/ and checkpoints excluded under
  the existing .gitignore; the processed snapshots require separate transfer.
- Frozen encoders remain the initial comparison baseline. The user now supports
  comparing multiple clustering methods and considering trained aggregation and
  ESM2 fine-tuning, especially if frozen representations are insufficient.
  Training is a planned future research branch, not restricted to a final fallback.
- Beyond the two selected baselines and the authorized PCA + average experiment,
  the research directions below remain future
  design options. Specific training choices, losses and splits are not final.
- B/E extraction orchestration and baseline clustering scripts are implemented.
  The user ran actual feature extraction on the server. G/H coordinate repairs
  are published in processed copies; retrieval evaluation and model training
  have not been run, and the encoder parsers themselves were not modified.
- The implementation choices below are recommendations, not a finalized protocol
  approved by the user. Follow the new conversation's instructions to proceed.
  This file does not itself authorize launching server jobs or modifying models.
- Do not repeatedly request permission for work the user subsequently authorizes.
  Preserve the distinction between design-only requests and implementation work.

## Workspace map

- Benchmark workspace: `/home/qiju/projects/BGC_retrieve/`, renamed from
  `/home/qiju/projects/BIG-SCAPE/` before the 2026-11-11 clustering work.
  It contains data, inventory scripts/tests and documentation, not a checkout or
  installation of the BiG-SCAPE clustering software.
- Raw benchmark data: `data/raw/curated_datasets/` (moved from
  `data/curated_datasets/` on 2026-11-09; do not use the obsolete path).
- Dataset source references: `data/raw/curated_datasets/README.md`.
- All derived data: `data/processed/`; current inventory: `data/processed/inventory_v1/`.
- G/H repair output: `data/processed/repaired_gh_v1/`, preserving original dataset,
  family and filename layout, with separate raw/processed hashes and repair evidence.
  Verify publication status in its metadata before using the copies.
- G/H encoder queue: `data/processed/pilot_gh_v1/` is published with metadata
  status complete, all 101 stable-ID GBK copies and embedded repair provenance;
  see `summary.md` and `input_checks.tsv`. This is the historical coordinate-only
  queue: H's hidden Geneious aSDomain annotations were not used. Its complete
  status describes that preparation run, not complete domain recovery; preserve
  it as evidence and do not use its H inputs for a new benchmark. Suggested server feature
  output for the new workflow is `data/processed/features_gh_reconstructed_v1/`,
  which has not been generated locally.
- H aSDomain recovery output: `data/processed/reconstructed_h_v1/`, sourced from
  `repaired_gh_v1/` and preserving all 90 H entries with original names and labels.
  G is not changed. Published with status `needs_review`: 2363/2365 annotations
  restored, 89/90 samples pass recovery checks, UOC915 retains two unresolved
  annotations. This snapshot is not itself a stable-ID encoder queue.
- New server queue location: `data/processed/pilot_gh_reconstructed_v1/`, created
  by `pilot_features.py prepare-reconstructed`. This is a suggested output path,
  not a locally published artifact; the user runs preparation on the server.
  New features go to `data/processed/features_gh_reconstructed_v1/`.
- Migration evidence: `data/processed/raw_relocation/`. All 2941 source files,
  including GBK/TSV/README and supplied HTML/JSON, were hashed before and after
  the directory move. The original file contents and internal layout are preserved.
- Reproduction and schema: `README.md`, `scripts/prepare_raw_layout.py`,
  `scripts/build_inventory.py`, `scripts/validate_inventory.py`, and `tests/`.
- Existing encoder project: `/home/qiju/projects/BGC_annotation/`.
- Original encoder upstream: <https://github.com/EvoCatalysis/BGC_annotation>.
- Server hardware reported by the user: NVIDIA L40S, 48 GB VRAM.
  DNA checkpoint weights and expensive extraction are on the server; config and
  previously exported embeddings are available locally. Verify paths before use.
- The local interpreter `/home/qiju/miniforge3/envs/BGC/bin/python` successfully
  ran the read-only Biopython audits. Do not assume it has the server's GPU stack.

Useful existing files in `../BGC_annotation/`:

| Purpose | Files |
| --- | --- |
| Frozen protein extraction | `experiment/extract_features.py`, `model/features.py`, `model/encoders.py` |
| GBK protein/domain parsing | `data_preparation/feature_inputs.py` |
| ESM inference and token cache | `data_preparation/frozen_esm.py`, `data_preparation/feature_cache.py` |
| Frozen DNA extraction | `experiment/extract_dna_features.py`, `data_preparation/frozen_dna.py`, `data_preparation/dna_inputs.py` |
| DNA full/window validation | `experiment/validate_dna_encoding.py` |
| Existing documentation | `docs/frozen-features.md`, `docs/dna-features.md`, `docs/umap-features.md` |
| Visualization | `visualization/umap_features.py` |
| Previous controlled DNA comparison | `outputs/analysis_dna_context_control/analysis.md` |
| Previous five-way comparison | `outputs/analysis_dna_full_comparison/analysis.md` |

Inspect the current checkout before changing anything. Files under `outputs/`
may be local/ignored rather than tracked. Do not assume benchmark implementation
or reports exist merely because proposed artifact names appear below.

## Existing representations and lessons from prior experiments

The existing feature format is `embeddings.npy`, `bgc_ids.tsv`, `metadata.json`,
and `failures.tsv`. Align rows using explicit row numbers and IDs, not directory
order. Validate dimensions, finite values, normalization, and provenance.

Protein routes:

- ESM2: raw domain/protein vectors -> masked mean -> L2 normalization.
  The usual exported dimension is 1280.
- MAC: the same token vectors -> frozen trained MAC `BGCEncoder` -> masked mean
  -> L2 normalization (usually 256 dimensions).
- MAP: the same token vectors -> frozen trained MAP `BGCEncoder` -> masked mean
  -> L2 normalization (usually 1024 dimensions).
- Dimensions must be read from actual checkpoints/metadata. The original class
  queries and prediction heads are not needed for this independent BGC vector.
  Reusing architecture is different from loading trained weights; do not create
  a randomly initialized encoder and call it a pretrained feature extractor.
- The parser orders CDSs by genomic position and domains within a CDS by protein
  position. A CDS with `aSDomain` annotations contributes its domains; a CDS without
  them contributes its whole translation. Gaps outside annotated domains are
  omitted for domain-annotated proteins. `PFAM_domain` is not interchangeable
  with `aSDomain`; the legacy Pfam script does not supply this domain slicing.
- Default ESM pooling is the final-layer BOS representation, not residue mean.
  The current frozen extractor covers long domain/protein inputs with disjoint
  windows of at most 1021 residues and combines window vectors by residue count
  into one vector per original token. This is not full-context protein encoding.

DNA routes:

- A trained `ModernBertForMaskedLM` with BPE tokenization is used as a frozen
  backbone; existing extraction pools final content-token states and exports a
  normalized 768-dimensional BGC vector. The MLM head is not the BGC embedding.
- The user confirmed with the model owner that training fed entire sequences
  through one forward without truncation/windowing. Record this as owner-reported
  provenance; it does not establish unseen training data or length generalization.
- Config `max_position_embeddings=2500` sets the default windowed baseline size.
  Windowed mode uses 2498 content tokens plus special tokens, with full coverage.
  This config value is not the full-mode input limit. Full mode explicitly sends
  the entire tokenized BGC through one forward.
- On 2026-11-04 the user explicitly confirmed that the DNA encoder supports
  arbitrary-length sequence input and that DNA length filtering is no longer
  needed. Keep the full labeled region; do not filter BGCs by DNA bp/token length,
  truncate long inputs, or inherit the historical 40 kb preparation step.
- Full mode is the primary DNA route. Windowed mode remains a same-input context
  control, not a length-admission policy or an automatic fallback from full mode.
- Compare full and windowed inference with the same checkpoint, tokenizer,
  FlashAttention-2 backend, BF16 precision, and input cohort. The historical
  windowed run used SDPA/FP32 and is a separate control.
- Saved FP32 embedding arrays do not imply FP32 model inference. The CLI still
  defaults to windowed mode; choose `--sequence-mode full` explicitly when needed.
- In the locally inspected CLI, `--max-bp` is optional and full-mode `--max-tokens`
  defaults to unlimited. Leave both length caps unset for this benchmark; skip
  `filter_gbk_by_length.py`. Record length for resource/performance stratification.
  Resource failures remain explicit execution failures; arbitrary-length input
  support does not imply unlimited GPU memory or resolve multi-record semantics.
- Ambiguous DNA bases/UNK tokens are retained; report their counts/proportions.
  Do not silently delete them or concatenate separate GenBank records into an
  invented contiguous chromosome.

The previous report, rechecked on 2026-11-04, evaluated 2635 shared labeled MIBiG
BGCs with coarse multi-label class Jaccard@15 in original cosine space:

| Representation | Class Jaccard@15 |
| --- | ---: |
| MAC | 0.750639 |
| MAP | 0.688476 |
| ESM2 | 0.662507 |
| DNA full, FA2/BF16 | 0.480240 |
| DNA windowed, FA2/BF16 | 0.477076 |
| Historical DNA windowed, SDPA/FP32 | 0.477316 |

These are not GCF, product-chemistry, or independent generalization scores. In the
controlled comparison, the 2041 multi-window queries did not improve on average
with full context (0.512853 full vs 0.514935 windowed). The 594 single-window
vectors were bitwise identical between the two routes; changes in their retrieval
with the full candidate pool came from changed candidate representations. DNA
distances were strongly associated with GC/source composition. These findings
motivate better relevance labels; they do not prove DNA lacks useful complementary
information. Read the report for limits and exact cohorts before quoting results.
Do not confuse the old DNA experiment letters A/B/C with benchmark datasets A/B/C.

## Papers and benchmark semantics

Primary sources already consulted:

1. BiG-SCAPE 1.0 / CORASON: *A computational framework to explore large-scale
   biosynthetic diversity*, DOI `10.1038/s41589-019-0400-9`.
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC6917865/>
2. BiG-SCAPE 2.0 / BiG-SLiCE 2.0: *BiG-SCAPE 2.0 and BiG-SLiCE 2.0: scalable,
   accurate and interactive sequence clustering of metabolic gene clusters*,
   Nature Communications (2026) 17:2000, DOI `10.1038/s41467-026-68733-5`.
   <https://doi.org/10.1038/s41467-026-68733-5>
   Accessible full-text PDF:
   <https://escholarship.org/content/qt2cj4p0m0/qt2cj4p0m0_noSplash_591ebe4fa5c4126d2496c8f75a999be7.pdf>
3. Published dataset/output archive, versioned record:
   <https://zenodo.org/records/16838820>.
   It includes `curated_datasets.tar.gz`, `curated_outputs.tar.gz`, and
   `curated_benchmarked.tar.gz`. The benchmark archive is now downloaded under
   `data/processed/vmeasure_comparison_v1/sources/`; the larger
   `curated_outputs.tar.gz` was not downloaded. Effective IDs remain unverified;
   published scores are external references, not a local reproduction.
4. Official benchmark implementation inspected during the investigation:
   <https://github.com/medema-group/BiG-SCAPE/blob/master/big_scape/benchmarking/benchmark_data_loader.py>
   <https://github.com/medema-group/BiG-SCAPE/blob/master/big_scape/benchmarking/benchmark_metrics.py>

Verified interpretation:

- The nine-dataset suite is from the 2.0 paper, including a curated MIBiG 1.3
  dataset inherited from earlier work. Do not say both original papers used all
  nine datasets.
- The main clustering metric is V-measure, the harmonic mean of homogeneity and
  completeness. It is not chance-adjusted. Homogeneity diagnoses mixed families;
  completeness diagnoses splitting a curated family across predicted clusters.
- For each tool/parameter configuration, the paper selects the scanned cutoff
  whose predicted GCF count is closest to the curated GCF count. It does NOT
  simply choose maximum V-measure. This uses ground-truth family count and is a
  matched-resolution diagnostic, not an independently selected deployment cutoff
  or a mathematical performance upper bound.
- The benchmark includes singletons and mixes biosynthetic classes without
  pre-partitioning by class. Reproduction used BiG-SCAPE 1.1.9 global/glocal and
  BiG-SCAPE 2.0.0-beta.1 combinations of global/glocal/local and extension modes.
  Pin versions/settings before claiming paper reproduction.
- BiG-SCAPE uses domain content, order, and sequence similarities, with pairwise
  alignment choices. Native clustering uses thresholded connected components
  followed by affinity propagation; 2.0 adjusts AP behavior for dense components.
  A global fixed vector and a pair-dependent local alignment retain different
  information. A benchmark can compare them without pretending they are equivalent.
- Dataset B targets similar architecture/metabolites across distant taxa plus
  singleton negatives. H targets subtypes within glycopeptide NRPSs. F has
  product-related curation. These are useful diagnostic cases, not interchangeable
  definitions of chemical similarity. B and especially G are small datasets.
- BiG-SCAPE 1 weights were tuned using MIBiG BGC/compound distances. Dataset E
  therefore is not entirely external to BiG-SCAPE algorithm development either.

Implementation pitfalls observed in the inspected `master` benchmark code:

- Short TSV labels are stored as `data[bgc] = family`: a later duplicate ID
  overwrites an earlier label.
- Metrics are computed on the intersection of curated and computed IDs, so failed
  inputs can disappear from scores unless coverage is checked explicitly.
- Do not infer that published scores necessarily used the same code behavior:
  verify the publication version and its effective inputs before reproducing it.

## Local dataset audit

Layout: each `Dataset_*` has `curated_GCF_assignments.tsv` (`BGC`, `GCF` columns)
and `GCF_assignments/<family>/*.gbk`. Counts below were rechecked on 2026-11-04.
They describe raw files/labels, before resolving duplicates or label ambiguity.

| Dataset directory | GBK / TSV rows | Distinct TSV BGC IDs | GCF labels | Singleton labels |
| --- | ---: | ---: | ---: | ---: |
| Dataset_A_Amycolatopsis | 1301 | 1200 | 439 | 294 |
| Dataset_B_Divergent | 39 | 39 | 23 | 16 |
| Dataset_C_Fusarium | 392 | 392 | 123 | 41 |
| Dataset_D_JK1 | 279 | 279 | 35 | 0 |
| Dataset_E_MIBiG_BS1 | 281 | 281 | 90 | 12 |
| Dataset_F_Metcalf | 110 | 110 | 10 | 0 |
| Dataset_G_Planomonospora | 11 | 11 | 3 | 0 |
| Dataset_H_Reclassification | 90 | 90 | 10 | 1 |
| Dataset_I_Rhodococcus | 408 | 393 | 33 | 1 |

Totals: 2911 file entries; 766 dataset-local family labels; 365 singleton labels.
The totals are not a global family vocabulary or 2911 independent BGCs.

Detailed read-only audit findings from 2026-11-02 (not rerun in full for this
handoff; revalidate when implementing the persistent inventory):

- A: 90 repeated-basename groups / 101 extra ID rows. In 74 groups the ordered DNA
  records are identical; 16 groups contain different sequences under the same
  basename. Some identical sequences have different family assignments.
- A: 18 multi-record files, containing 23 extra records. All 2911 GBKs could be
  read with multi-record parsing, yielding 2934 sequence records overall.
- I: 15 repeated-basename groups, all with identical DNA sequences across family
  assignments. Example: `CP008947.1.cluster021` occurs with
  `heterobactin_rhequibactin` and `nrps_6` (TSV lines 284 and 397 in the audited
  copy). Seven files are in folder `other_11.2` while their TSV label is `other_11`.
- The record-level forward-DNA audit found 2757 unique sequences and 21 sequence
  groups shared across datasets. This was NOT reverse-complement or homology
  deduplication, and record counts are not file-level sample counts.
- B has 16 explicit singleton `neg_*` labels. Preserve their negative-candidate
  role. Some labels (e.g. holomycin/desferrioxamine) have trailing nonbreaking
  whitespace; preserve raw labels and log any canonicalization.
- C has no `aSDomain` annotations at all. Current protein parsing uses whole-CDS
  tokens there, including a successfully parsed token of 12074 aa. Annotation
  differences are a potential input-representation confounder.
- H sequences range from 43475 to 99572 bp. The historical 40 kb filter would
  remove the entire dataset and is no longer part of this benchmark's policy.
  The audited per-record maximum across the suite was 281670 bp in A.

Existing-parser checks, one file at a time (not GPU encoding):

| Dataset | DNA parse passes | Protein parse passes |
| --- | ---: | ---: |
| A | 1283 / 1301 | 1112 / 1301 |
| B | 39 / 39 | 39 / 39 |
| C | 392 / 392 | 391 / 392 |
| D | 279 / 279 | 274 / 279 |
| E | 281 / 281 | 281 / 281 |
| F | 110 / 110 | 110 / 110 |
| G | 11 / 11 | 8 / 11 |
| H | 90 / 90 | 75 / 90 |
| I | 408 / 408 | 404 / 408 |
| Total | 2893 / 2911 | 2694 / 2911 |

Bulk discovery rejects duplicate basenames before parsing. The table does not
mean the current extraction CLI can process A/I or all nine directories directly.

Protein failures include `*` in translations, missing translation, domain-parent
assignment failures, inconsistent coordinates, repeated-domain ambiguity, and
multi-record inputs. Do not globally suppress validation or silently drop CDSs.
Distinguish terminal stops, internal stops, pseudogenes, and recoverable metadata.

A targeted check of G/H domains with translations, protein coordinates, and a
uniquely containing CDS found 32/618 mismatching intervals, respectively;
all checked mismatches matched the parent translation exactly when
the annotated end was increased by one. Other domains already matched half-open
coordinates. This supports an evidence-based legacy coordinate adapter, NOT a
blanket `end += 1` transformation. Require sequence agreement, preserve raw and
resolved coordinates, and verify parent identity. Other failure types remain open.

Do not flatten family directories, use basename alone as identity, select the
first/last conflicting label, or assume family folders settle every conflict.
Determine whether a conflict represents multiple biological units in one region,
overlapping family membership, a naming collision, or a curation inconsistency.

## Current stage: B/E baseline clustering and remaining dataset preparation

The first inventory step, B/E server frozen feature extraction and G/H processed
coordinate repairs are complete. G/H encoder input compatibility and extraction
orchestration are now checked; actual G/H inference remains for the server.
Label adjudication, compatibility handling and
cohort selection for A/C/D/F/I remain open.
The B/E workflow preserves all 39/281 entries, including singletons and cross-set
duplicates, without length filtering or annotation repair. It is not a train/test split.
Confirm consequential unresolved biological policies through source evidence and
the user's instructions; routine implementation choices within an authorized
stage do not need repeated approval.

### Persistent inventory and input policy

Maintain the reversible inventory before encoding. Keep original data unchanged
in `data/raw/` and write derived inputs/results under `data/processed/`. Include dataset, opaque unique sample ID,
original path/ID/labels, record index/count, sequence and source hashes, length,
organism/GC, CDS/domain availability, parsing status, repairs, and exclusions.
Encoding must not read family labels or use them to crop/reorder/select features.

Resolve A/I identities and label conflicts. Same sequence with validated multiple
labels belongs in a multi-label retrieval analysis; it cannot be forced into the
single-partition ground truth assumed by V-measure/ARI. Preserve unresolved cases
in an explicit inventory and report exclusions; do not silently manufacture a
new single family from a label combination.

Validate GBK compatibility while preserving the labeled region's boundaries.
Do not regenerate antiSMASH boundaries or concatenate multiple records without a
documented biological policy. Any inferred translations must respect strand,
joined CDS locations, translation table, partial features, and annotation evidence.
Prefer validated annotated translations where available.

Audit duplicate sequences (including reverse complements), parent/reference IDs,
protein duplicates, and feasible homology overlap. Cross-dataset duplicates must
not leak across development/evaluation partitions. Keep gene-family homology
checks distinct from simple organism/genus grouping.

Distinguish source-file identity, sequence duplicate groups, and final evaluation
samples. Identical DNA with different annotations need not yield identical protein
inputs. Do not automatically turn conflicting labels into validated multi-label
truth, or split a multi-record file into samples that all inherit its file label.
Keep unresolved cases in the inventory with explicit coverage consequences.

Dataset-preparation deliverables are the versioned input/protocol specification,
persistent manifest, duplicate and label-conflict inventories, evidence-based
compatibility policy and repair/exclusion records, and a proposed evaluation
cohort with coverage counts. B/E additionally has a prepared cohort, returned
server feature exports, two completed clustering sweeps and a V-measure summary.
Retrieval scores, model training and DNA/protein fusion results remain pending.

## Future research roadmap (design only; not current tasks)

The following stages record discussed options for later work. B/E server feature
scripts and the original-cosine average/Leiden baselines are authorized and implemented.
PCA + average-linkage is additionally authorized as recorded above; the other directions
do not authorize starting experiments. A candidate
sequence is frozen baselines with multiple clusterers, then learned BGC aggregation
or metric learning, then ESM2 fine-tuning comparisons, followed by renewed DNA
complementarity evaluation. Frozen fusion can also be tested as an initial control.

### 1. Five frozen baselines on a fixed cohort

Encode ESM2, MAC, MAP, DNA full FA2/BF16, and DNA windowed FA2/BF16. Preserve the
current pooling first so that the new relevance definition is the main change.
Reuse one verified ESM token cache across protein routes. Existing MIBiG vectors
can only be reused after checking sequence, annotation, preprocessing, and model
provenance; matching an accession alone is insufficient.

Fix a shared comparison cohort per dataset after the input audit, and require
all methods to account for it. Report route-specific success/failure coverage
alongside any common-success analysis; do not score different hidden intersections.
Report domain/protein tokenization and sequence coverage differences. A uniform
annotation or whole-CDS control may later help separate model and input effects.

Run small server validation before bulk extraction, covering long inputs and
legacy annotation cases. Record checkpoint/tokenizer hashes, pooling, backend,
precision, software versions, GPU, elapsed time, peak memory, and failures. OOM
is an execution failure, not a biological reason to quietly remove a sample.

### 2. Family retrieval as the primary embedding evaluation

Within each dataset, query each BGC against the other members using original
normalized-vector cosine similarity, excluding self and handling duplicate groups
according to the fixed input policy. Compute these retrieval metrics in original
embedding space, not on visualization coordinates. This does not prohibit a
separately declared dimensionality-reduction-plus-clustering experiment below.

Recommended primary metric: mAP; complementary metrics: Recall@5 and Recall@10
where the candidate pool permits it. Define tie handling and small-pool behavior
before running comparisons. Report query-macro and family-macro results, followed
by an equal-dataset summary; retain individual dataset results and sample counts.

Same curated GCF defines binary relevance for validated single-label cohorts.
For verified multi-label cases, predefine shared-family relevance and report a
separate panel. Singletons stay in the candidate pool and clustering task, but
have no positive retrieval targets: omit those queries from AP/recall averages
and explicitly report eligible-query coverage, rather than scoring them as zero.

### 3. Multiple clustering methods and BiG-SCAPE baselines

The user explicitly wants multiple clustering methods in the longer-term plan.
The current PCA experiment is limited to average-linkage; do not extend its grid
to Leiden. Average linkage remains
one reference method, not the sole permitted clustering algorithm. Candidate
routes with different assumptions are:

| Route | Question / parameters to control |
| --- | --- |
| Cosine-distance hierarchical clustering, average and complete linkage | Distance hierarchy; linkage and cut threshold |
| Spherical k-means | Directional prototype structure; cluster count and initialization |
| Similarity-weighted kNN graph + Leiden | Local-neighborhood communities; neighbors, edge weights, resolution, seed |
| HDBSCAN | Density structure and unassigned samples; minimum cluster size, density settings, noise rate |
| PCA + GMM, exploratory | Soft membership; dimensionality, covariance constraints, component count |

Compare a representation-by-clustering-method matrix on the same audited cohort.
Predeclare tuning budgets and use development data to select parameters; do not
report each representation's best test-label-selected algorithm as independent
test performance. Different representations may need different numeric settings.
Ordinary k-means on normalized inputs is not automatically spherical k-means.

Report V-measure, homogeneity, completeness, ARI, predicted cluster count,
singleton count, coverage, and stochastic stability where applicable. Do not
silently discard HDBSCAN noise or treat all noise points as one biological family.
A proposed full-cohort scoring convention assigns each unassigned point its own
singleton label, alongside noise rate and explicitly conditional non-noise scores;
finalize this convention before experiments.

Treat original-space, PCA, and UMAP preprocessing as explicit experimental
conditions. UMAP + HDBSCAN may be studied as a complete pipeline, with declared
dimensions, neighborhood settings, seeds, and sensitivity analysis. A cleaner 2D
plot is not selection evidence. In inductive evaluation, fit preprocessing on
training/development inputs and apply its fixed transform to held-out inputs;
independent train/test UMAP fits do not share a coordinate system. Any transductive
clustering of an entire held-out collection must be identified as such. Supervised
UMAP uses labels and belongs to a separately evaluated supervised pipeline.

Keep two distinct analyses:

- Paper-style matched family count: select the scanned cut/resolution setting
  yielding the closest curated family count, with a predetermined tie policy.
  For methods that take a cluster count directly, setting it from the truth also
  uses that information. This is a matched-resolution diagnostic, not a deployed
  parameter-selection score; do not select by maximum test V-measure.
- Transferable settings: choose thresholds, cluster-count selection rules, graph
  settings, and other parameters on development datasets, then evaluate on held-out
  datasets. A leave-one-dataset-out scheme is a candidate, provided cross-dataset
  duplicates and related groups are controlled.

Use native BiG-SCAPE 1/2 results as external pipeline baselines. Check the author's
published outputs first; download only what is needed once that work is requested.
Recompute metrics on matching audited IDs or rerun when input policies differ.
Do not attribute an end-to-end difference entirely to embeddings when clustering
algorithms/input processing also differ. A common-clusterer comparison on distance
matrices is a possible later control.

### 4. Pooling and learned BGC representations

After baselines, compare ESM/MAC pooling choices on the fixed protocol. Candidate
directions discussed include gene-balanced hierarchical pooling and functional
blocks. Final BGC vectors cannot be repooled: cache ordered token states, masks,
parent/domain mappings, input hashes, and checkpoint identity if needed. Never
align a sidecar solely by matching its length.

Treat representation learning as a distinct experiment from replacing the
clusterer. First consider frozen ESM2 with a trainable BGC aggregator/projection
head, supervised by BGC-level family or independently supported functional
relationships. Contrastive/metric learning is a candidate objective aligned with
retrieval; family classification can be an auxiliary objective. Neither loss nor
architecture has been finalized. Distinguish training a head from adapting ESM2.

### 5. ESM2 adaptation / fine-tuning

The user is open to fine-tuning ESM2 if frozen representations are insufficient.
Plan comparisons rather than assuming frozen weights are a permanent constraint:

1. Frozen ESM2 + trainable BGC aggregation/projection.
2. The same aggregation/projection and objective + ESM2 LoRA or partial unfreezing.
3. Broader or full fine-tuning if validation gains and resource feasibility justify it.

Keep training data, sampling, loss, output dimension, and evaluation comparable
when measuring the incremental effect of adapting the backbone. Compute the loss
after BGC aggregation and backpropagate to the protein encoder when adapting it.
Do not copy a BGC family/product label onto every constituent protein as if it
were that protein's own functional label. A BGC contains multiple complementary
enzymes and potentially shared transport/regulatory proteins.

True ESM2 adaptation requires raw amino-acid inputs and fresh forward passes;
fixed token caches only support frozen-backbone experiments. Preserve provenance
and re-export features for each adapted checkpoint. Do not assume DNA's arbitrary
length support changes protein input/window handling or fine-tuning memory costs.

Separate low-homology new BGCs from seen families and retrieval/clustering of
families absent from training. Keep all derived proteins/domains of a BGC and
related duplicate/homology groups on the same side of relevant split boundaries.
Audit training overlap and selection leakage before claiming generalization.
Improved training-family classification alone does not establish a better space
for unseen-family retrieval or independent product-chemistry relationships.

### 6. DNA/protein complementarity

Evaluate frozen fusion first, starting with separately normalized representations:
`s_fused = (1 - alpha) * s_protein + alpha * s_DNA`.
This equals cosine similarity of the candidate-independent concatenated vector
`z = [sqrt(1 - alpha) * p; sqrt(alpha) * d]` for unit-normalized `p` and `d`.
Choose a small, declared alpha grid or fixed values using development data. The
grid and exact pooling choices have NOT been finalized. Select protein route,
pooling, cutoff, and fusion settings without tuning on the evaluation labels.

Compare against pure protein and an ESM+MAC protein-fusion control. Include simple
GC/length/k-mer composition controls where useful. Examine both correct family
neighbors rescued by DNA and new errors introduced by fusion, at equal retrieval
budget. Different nearest neighbors alone do not establish useful complementarity.
Report cross-taxon, length, annotation-coverage, and fine-grained-family results
with coverage. Use paired family/dataset-aware uncertainty summaries; do not treat
correlated BGC queries as independent experimental replicates.

After protein adaptation, reassess DNA's incremental value against the adapted
protein-only baseline too. For any learned fusion, include a comparable-capacity
protein-only control; extra trainable parameters are not evidence of DNA value.

Attention pooling over DNA windows was considered earlier, but adding a random
attention module is not a justified training-free improvement. Keep it outside the
initial fixed baseline; any learned aggregation needs a separate explicit design.

### Research references informing these future options

These sources were discussed on 2026-11-04; they motivate experiments, not an
assertion that a particular method will win on this BGC benchmark.

- HieVi, *Phage evolutionary relationships emerge from protein language model-based
  proteome representation* (2025): frozen ESM2 protein vectors aggregated into a
  proteome vector, then HDBSCAN in original high-dimensional space. Its taxonomy
  and evolution results do not establish BGC product/GCF performance.
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC12541379/>
- *Protein Set Transformer: a protein-based genome language model to power
  high-diversity viromics* (Nature Communications, 2025): frozen ESM2 inputs,
  trained contextual aggregation, and weighted kNN + Leiden. Gains include the
  learned aggregator; its singleton exclusions are not our benchmark policy.
  <https://www.nature.com/articles/s41467-025-66049-4>
- *Fine-tuning protein language models boosts predictions across diverse tasks*
  (Nature Communications, 2024): direct frozen/full-tuning/LoRA comparisons,
  including ESM2. Supports testing adaptation, without guaranteeing BGC transfer.
  <https://www.nature.com/articles/s41467-024-51844-2>
- *Exploring sequence landscape of biosynthetic gene clusters with protein language
  models* (ICML 2024 workshop): ESM2-8M partial fine-tuning with contrastive learning,
  using coarse biosynthetic categories and leave-one-class-out evaluation. This
  is relevant prior work, not an independent curated-GCF benchmark result.
  <https://openreview.net/forum?id=UguakxTKG5>
- User-provided PLMCCL-TP (Cluster Computing, 2025): therapeutic-peptide multilabel
  classification with PLM, contrastive learning, clustering features, and a CNN.
  Only the publisher abstract/figure previews and public code were inspected;
  the subscription-only full text and Table 1 were not verified. The README
  points to ProtBERT and code uses `BertModel`, not ESM2; the precise local
  checkpoint is not established. Public code intends backbone + head training,
  then UMAP/GMM-derived classification features. Its independently fitted
  train/test UMAP spaces require an implementation audit before reuse; no
  reproduction was run and this observation does not invalidate published results.
  <https://link.springer.com/article/10.1007/s10586-025-05580-x>
  <https://github.com/jiangkun-lin/PLMCCL-TP>
- DeepT2, *A deep learning model for type II polyketide natural product prediction
  without sequence alignment* (Digital Discovery, 2023): ESM2 KS-beta embeddings
  with supervised UMAP + HDBSCAN. Product labels inform the space; do not cite it
  as evidence of wholly unsupervised recovery of independent product labels.
  <https://pubs.rsc.org/en/content/articlehtml/2023/dd/d3dd00107e>

## Leakage and interpretation limits

- Freezing weights does not remove training overlap. MAC/MAP were trained with
  MIBiG-derived data. Audit exact training IDs/checkpoint member and homologs where
  available; dataset E is especially relevant, and other sets contain MIBiG refs.
- DNA pretraining used genomic sequences; the actual training corpus overlap is
  unknown. Label it unknown unless verified. Pretraining exposure and supervised
  target leakage are different issues; report both accurately.
- A family benchmark is closer to the user's goal than six-class retrieval, but
  curated GCF membership is not a universal continuous product-chemistry metric.
  Claims about chemical complementarity need independent chemical relevance.
- B/F/H are informative biological case studies; G has only 11 raw entries. Do
  not generalize from a few cases or hide heterogeneous performance in one mean.
- Dataset, family, and homology grouping have different meanings. Cross-genus
  evaluation is not proof of homology independence.

## Concrete next step and expected artifacts

The full-suite inventory is present
at `data/processed/inventory_v1/`; read its `summary.md` and the project README
before repeating work. It contains 2911 source-file entries, 2934 sequence records,
and all 2911 TSV rows. There are 2904 exact file/label matches and 7 explicit unique
basename matches with folder/TSV label disagreement. The 18 multi-record files and
217 protein-parser failures remain in the inventory; DNA parsing passes 2893 files
and protein parsing passes 2694. This reproduces the earlier single-file audit.

Stable sample IDs identify archive-relative file entries, not confirmed independent
biological BGCs. `label_rows.tsv` preserves source line numbers and raw strings;
duplicate groups retain sample/record membership. No entries were removed or
repaired and there is no DNA length filtering. A preliminary generated inventory
is retained under `data/processed/_history/`; `inventory_v1/` is the current output.

For A/C/D/F/I, subsequent preparation should inspect conflicts and failure
records to define evidence-based label/compatibility policies, then produce an
audited evaluation cohort. Do not jump to future clustering, fine-tuning, or fusion
because those branches are documented above.
B and E are the user-selected pilots (39 and 281 files). F remains a candidate
requiring an explicit duplicate policy. G/H repaired inputs now pass CPU parsers
and have an extraction preparation path. Extend to C/D after compatibility
handling and to A/I after identity/label resolution. If pilot labels influence model selection, do not later call those
same datasets untouched test sets.

### Implemented B/E server workflow

- `data/processed/pilot_be_v1/` contains 320 byte-identical GBK copies renamed to
  inventory sample IDs, `cohort.tsv` with original IDs/paths/raw labels, and per-set
  label TSVs outside the encoder input directories. Preparation verifies GBK and
  label hashes against the inventory and preserves Unicode label whitespace.
- `scripts/pilot_features.py extract` calls existing BGC_annotation CLIs in
  subprocesses with local server checkpoints. ESM2 uses BOS, 1021-aa windows and
  masked-mean/L2 BGC pooling; MAC consumes that exact new ESM cache, using one
  trained ensemble member (default 0). DNA explicitly uses full mode, defaults
  to FA2/BF16, and has no bp/token caps or automatic window fallback.
- Outputs are separate per dataset/route under a user-chosen processed directory.
  Each has native arrays, IDs, metadata, failures and an execution log. `run.json`
  records commands, checkpoint/code hashes, elapsed times and validation results.
- `--resume` checks unchanged configuration and completed output hashes, skips
  validated stages, archives failed attempts and restarts those stages. There is
  no per-sample resume. Changes to models/settings/code require a new run directory.
- Output validation requires every expected sample, finite FP32 unit vectors,
  matching source hashes/cache provenance, no exclusions, and full-DNA coverage.
- Real local B/E inputs are prepared. Unit tests use synthetic outputs/mocked
  child processes; they do not establish model accuracy or server GPU compatibility.

Created: source-file and sample manifests, record inventory, raw TSV-row inventory,
duplicate groups/members, issue inventory, per-dataset summary, provenance metadata,
and reproducible scripts/tests. Policies for this inventory are documented in the
README and metadata. Still pending: biological conflict resolution, actual repair
or exclusion decisions, and the final input/evaluation protocol and cohort.

Verified local outputs: B/E per-route feature exports, both baseline clustering
sweeps and their V-measure comparison. Future outputs beyond these baselines
include per-query retrieval metrics, additional clusterers and adaptation/fusion
ablations. These remain outside the current comparison task.

### Implemented original-cosine clustering baselines

- `scripts/clustering_common.py` reads native `embeddings.npy`, `bgc_ids.tsv`,
  `metadata.json`, `failures.tsv` and joins with `pilot_be_v1/cohort.tsv`. Explicit
  row indices locate vectors; stable IDs determine common ordering. Server paths
  in provenance need not exist locally. No hidden intersections or label stripping.
- `scripts/cluster_average.py`: SciPy average linkage on condensed original cosine
  distances, default cuts at zero and all distinct merge heights. Saves full tree,
  row/leaf mappings, memberships, metrics and diagnostic figures. No truth-based
  cut selection; ties need not yield every possible cluster count.
- `scripts/cluster_leiden.py`: exact cosine kNN, default union graph, positive
  cosine weights counted once per undirected edge, RBConfiguration Leiden.
  Defaults k=5/10/15, resolution=0.25/0.5/1/2/4/8, seeds=0..4. Isolates remain
  separate singletons. Reports every setting and seed, plus pairwise seed ARI.
- Both operate in original embedding space with no PCA/UMAP. Metrics include
  homogeneity, completeness, V-measure, ARI, AMI, cosine silhouette, cluster and
  singleton counts, and full coverage. Raw GCF strings are evaluation-only.
- On 2026-11-11 both default scripts completed on all six real B/E exports:
  `data/processed/clustering_average_v1/` contains 960 cut results, and
  `data/processed/clustering_leiden_v1/` contains 540 parameter/seed results.
  Every setting retains B=39/E=281. Sixty-six tests passed; independent output
  checks confirmed memberships/IDs/labels/cluster counts, artifact/source hashes,
  and bitwise-identical cosine matrices between the two baselines. These are
  diagnostic sweeps; no best-performing configuration or model has been selected.
- `docs/clustering_baselines.md` contains commands and interpretation limits;
  `requirements-clustering.txt` pins the tested CPU dependencies. No test-label
  optimum is automatically selected by the clustering scripts. The separate
  summary below uses true family counts for a diagnostic comparison.

### Implemented B/E V-measure comparison

- `scripts/summarize_vmeasure.py` verifies clustering artifact hashes and recomputes
  homogeneity, completeness, V-measure and ARI for all 1500 membership partitions
  against unchanged cohort labels, agreeing within `1e-12`. Coverage is complete:
  B has 39 BGCs / 23 GCFs; E has 281 BGCs / 90 GCFs.
- Average selects the nearest true GCF count, breaking ties by smaller threshold;
  all six results match the target exactly. For each fixed Leiden k, selection
  minimizes mean `abs(K - true_GCF_count)` across all five seeds, breaking ties
  by smaller resolution. The short table displays k=10, the original grid's middle
  value; k=5/10/15 are all retained. Neither rule maximizes V-measure. This is
  label-informed diagnosis, not independent parameter selection or generalization.
- Average BGC-MAC V-measure is 0.870008 on B and 0.813323 on E, higher than ESM2
  and DNA for this selection rule. E's existing Leiden grid reaches at most 50
  clusters, below its 90 true families; it is not yet a matched-resolution test.
  All-singleton controls and ARI are retained because V-measure is not adjusted
  for chance. Four additional summary selection tests passed.
- `data/processed/vmeasure_comparison_v1/summary.md` and `comparison.png` contain
  the report/figure; TSVs retain all fixed-k results and exploratory max-V results
  separately. Reproduction commands are in README. Choose a new `--output` path
  because the summarizer refuses to overwrite an existing report.
- Eighteen author-published B/E scores are imported from Zenodo record 16838820:
  BiG-SCAPE 1 global/glocal, and BiG-SCAPE 2 global plus local/glocal with three
  extension modes. This is an external reference, not a local BiG-SCAPE rerun.
  Logs have matching nominal sample counts, but effective sample IDs are not
  verified and some exported confusion matrices have incomplete counts. Consult
  `sources/provenance.json` before making any strict same-cohort claim.

### B/E PCA + average-linkage experiment

- `scripts/cluster_pca_average.py` uses the same strict native feature loader and
  unchanged B/E cohort as the raw-space baselines. It fits mean-centered full SVD
  PCA separately per dataset/route, without per-feature scaling or whitening,
  then L2-normalizes projected rows before cosine + existing average linkage.
- Default fixed dimensions are 2/4/8/16/32/64/128/256 and variance targets are
  0.90/0.95/0.99. Raw, centered-only and full numerical-rank PCA controls separate
  centering from truncation effects. Variance ratios describe centered input
  variance before post-projection L2, not family information or performance.
- Rank is bounded by `min(N-1,D)` and numerical SVD tolerance. Oversized dimension
  requests are recorded as skipped, never silently clamped. Duplicate selectors
  reaching the same dimension share one condition. Zero/near-zero projected rows
  fail the entire condition explicitly; no sample is dropped to salvage a score.
- Every condition scans all distinct merge heights plus zero, then summarizes
  the closest true family count with smaller-threshold ties. Labels do not enter
  PCA or tree building. All dimensions remain reported; 95% is a predefined
  short-table display, not the best V-measure dimension.
- Fitting PCA on the full unlabeled cohort is transductive. This diagnostic uses
  family counts in cut selection and is not an independent generalization score.
  The current step runs average-linkage only, with no Leiden or model inference.
- Output defaults to `data/processed/clustering_pca_average_v1/`; retains fitted
  transforms, aligned sample mappings, condition vectors/distances/trees,
  memberships, all-cut metrics, matched metrics, plots and provenance. See
  `docs/pca_average.md`.
- On 2026-11-12 the default real B/E run completed: 72 conditions and 12,488
  all-height cut partitions; 9 oversized B dimension requests skipped, with no
  failed conditions. Every partition retains all B=39/E=281 samples. Six raw
  reference checks passed against the prior average run, including bitwise-equal
  distances and metrics within `1e-12`. Independent checks recomputed H/C/V/ARI
  for all 72 selected memberships, verified nearest-family-count selection and
  output hashes; full-rank versus centered-only cosine errors were below `3e-15`.
  Twenty-seven related tests passed (PCA core/pipeline, shared inputs and average
  hierarchy); this step ran no Leiden analysis or Leiden tests.
- The fixed 95% variance display does not improve every route: for B/E, ESM2
  V-measure is 0.8327/0.8057, MAC 0.8674/0.8011 and DNA 0.7754/0.7588. E/ESM2's
  increase over raw (0.7904) is already mostly present in centered-only (0.8060).
  Consult the complete dimension curves and centered controls before attributing
  a gain to dimensionality truncation; no best-performing dimension was adopted.

### G/H coordinate repair in processed copies

- The implemented entry point is `scripts/repair_gh.py`, using
  `scripts/domain_coordinate_repair.py`. The real run was published on 2026-11-19
  to `data/processed/repaired_gh_v1/`, with status complete and all 101 G/H GBKs.
- G has 266 repaired domains in 5/11 files (50 aSDomain, 216 PFAM_domain); H has
  618 repaired aSDomains in 15/90 files. Across 884 genomic endpoint corrections,
  866 also increment explicit protein_end; 18 G domains receive only a genomic
  correction and retain absent protein-coordinate qualifiers. The other 81 GBKs
  and both label TSVs are byte-identical to raw.
- Protein parsing improves from G=8/11 and H=75/90 to G=11/11 and H=90/90;
  DNA parsing passes all 101. All 2941 raw files match inventory and before/after
  SHA256 snapshots. Checks cover 5039 domain annotations; 2 remain needs_review.
  A second repair plan on each output proposes no changes. Original inventory
  failures still describe raw inputs; read repaired samples.tsv for current status.
  Thirteen focused core/workflow tests passed, covering strand-aware byte edits,
  unsupported evidence, preservation checks, output path protection and failure
  without publishing a partial result.
- Repair only sequence-proven `aSDomain` / `PFAM_domain` C-terminal endpoint errors: require a
  unique same-strand containing CDS, agreement with its protein interval, raw
  domain DNA translating to the annotated sequence minus its final residue, and
  a one-codon extension translating exactly to the complete annotated sequence.
  The proposed location must match parent CDS coordinate/frame mapping.
- For explicit inclusive protein ends, jointly increment `protein_end` by one
  and extend genomic end by three bases on the plus strand or genomic start by
  minus three on the minus strand. Domains lacking protein coordinates can have
  genomic-only fixes when unique translation matching proves the interval;
  do not add absent protein-coordinate qualifiers or apply a blanket increment.
- Preserve feature types: repaired `PFAM_domain` entries remain Pfam annotations,
  not new ESM aSDomain tokens. Two orphan Pfam annotations in H `ApSA.gbk` have no
  supported parent CDS and are retained unchanged as `needs_review`; protein-parser
  compatibility does not establish that all biological annotation issues are solved.
- Preserve raw files and the original inventory; make only necessary numeric
  byte edits in GBK copies. Preserve DNA, CDS/domain translations, other
  annotations, original filenames, family labels and complete G/H membership.
  Retain byte-identical copies for unaffected files and labels.
- Publish to a new directory strictly under `data/processed/` only after all
  repaired GBKs pass existing protein and DNA parsers, semantic preservation
  checks and raw before/after integrity verification. Record `samples.tsv`,
  `domain_checks.tsv`, `repairs.tsv`, `raw_integrity.tsv`, per-dataset counts and
  provenance hashes. No encoder parser/model modifications or sample exclusions.
- The repair snapshot itself is not a prepared encoder directory. Use the new
  `pilot_features.py prepare-repaired` bridge described below for the historical
  coordinate-only queue. H requires the later aSDomain recovery before a new
  domain-based benchmark. No G/H extraction or clustering is claimed by this step.

### G/H encoder compatibility and server workflow

- The real preparation has completed at `data/processed/pilot_gh_v1/`, metadata
  status complete. Loading the published cohort validates all G=11/H=90 rows;
  a temporary relocated copy was also checked without its original source paths.
  See [input summary](data/processed/pilot_gh_v1/summary.md). Eleven G/H preparation
  tests and eighteen existing extraction/validator tests pass. Existing B/E six
  exports also pass the updated validator, and the nine local encoder source
  hashes match their recorded B/E run versions.
- Existing BGC_annotation ESM2/MAC/DNA interfaces accept all repaired inputs in
  CPU parser checks: G=11 BGCs/3 GCFs, H=90 BGCs/10 GCFs. This establishes parsing
  compatibility only, not correct recovery of every domain or checkpoint/GPU
  inference success. The H Geneious omission discovered afterward invalidates
  the earlier inference that no further input repair was needed.
- `scripts/prepare_repaired_features.py`, called by
  `pilot_features.py prepare-repaired`, validates the complete repair artifact
  hashes, sample/repair/label consistency and GBK membership. It copies processed
  bytes under stable sample IDs, checks protein/DNA inputs with current parsers,
  and prepares `bgc_repaired_pilot_inputs_v1` metadata, `cohort.tsv`,
  `input_checks.tsv`, a summary and embedded `repair_provenance/` manifests/labels.
- In that historical cohort, `source_sha256` is the actual processed GBK hash used by all
  output validators; `raw_source_sha256` is separate original-backup provenance.
  Preserve all 101 entries and raw Unicode labels; no filtering or exclusions.
  The complete prepared directory is portable to the server: validation does not
  require original local raw/repaired absolute paths. Existing B/E format and
  commands remain supported; raw `prepare` still accepts only B/E.
- G DNA spans 15,484–76,182 bp with 15–71 protein/domain tokens per BGC; H spans
  43,475–99,572 bp with 13–92 tokens. G has no token above 1021 aa. H has 156 such
  tokens across 42 BGCs, maximum 9508 aa. ESM2 covers them with the existing
  disjoint 1021-aa window policy, then residue-weighted pooling restores one
  vector per input token; MAC still receives at most 92 tokens per H BGC.
  These are pre-recovery counts. All 156 H long tokens were unsplit CDSs with
  hidden Geneious aSDomain annotations, not domains that remain long after
  applying their existing annotations. The actual recovery reduces >1021-aa
  tokens to zero; this does not resolve the two missing UOC915 annotations.
- H retains 246 N bases; G has no ambiguous bases in the audited inputs. DNA
  remains full/FA2/BF16 with no bp/token cap. Prior B/E full-DNA exports include
  successful inputs up to 185,250 bp; G/H actual BPE lengths/VRAM remain to be
  measured on the server. Do not silently switch to windowed after OOM.
- Reuse the B/E ESM2 checkpoint, MAC member 0, DNA checkpoint/tokenizer and server
  environment; preserve code/model hashes in the new run. MAC consumes this run's
  ESM cache. `extract` and `validate` require all G=11/H=90 rows per route and
  the same native arrays/ID/metadata/failure format. Commands: `docs/gh_features.md`.
- No G/H embeddings, V-measure, BiG-SCAPE execution or model training are claimed.
  The clustering/PCA/score scripts still register only B/E; extend their cohort
  handling after G/H features return before analysis. Continue with average/PCA
  under the user's current preference; do not add further Leiden tests.

### H Geneious aSDomain recovery

- The user authorized recovering existing H aSDomain annotations, while explicitly
  deferring Pfam-coordinate splitting. Continue using the existing rule:
  recognized aSDomain annotations replace a CDS with its domain tokens; a CDS
  without recognized domains contributes its whole translation. Do not introduce
  a new segmentation policy or make recovered domains from PFAM_domain features.
- The actual run is published at `data/processed/reconstructed_h_v1/` with schema
  `bgc_h_asdomain_reconstruction_v1`, status `needs_review`, snapshot_complete true
  and ready_for_full_h_encoding false. It retains all 90 H samples/10 GCFs,
  modifies 45 GBKs and restores 2363/2365 annotations. There are 91 parent-ID
  corrections, 187 duplicate groups/190 extra entries retained, and two unresolved
  annotations in UOC915. Protein/DNA parsing passes 90/90; 89 samples pass the
  aSDomain recovery check. Long tokens (>1021 aa) decrease from 156 to zero, but
  UOC915's known missing annotations still prevent claiming a complete H baseline.
  Total tokens increase from 4836 to 6791 (4718 domains + 2073 whole proteins);
  the maximum token is now 1008 aa. The other 45 GBKs are byte-identical to source.
- Sixteen core and six workflow tests pass. The full run validates rereading all
  GBKs, idempotence, non-target feature preservation and the expected token changes.
  All 2941 raw files and the complete source snapshot retain their before/after
  hashes. Independent verification also checks all 106 artifacts and current
  script hashes. No model inference or new prepared encoder queue was created.
- H contains 2365 `misc_feature` annotations in 45/90 files with the exact note
  `Geneious type: aSDomain`. They lack translation and protein-coordinate qualifiers
  and are ignored by the current parser. Missing aSDomain family names are not
  permission to invent names, HMM scores or Pfam identifiers.
- Read from `data/processed/repaired_gh_v1/`, preserving prior sequence-proven
  endpoint repairs. Write a new H-only snapshot under
  `data/processed/reconstructed_h_v1/`; retain all 90 H entries, family paths,
  filenames and labels. Do not modify G, raw, inventory or earlier processed data.
- Restore only with an unambiguous same-strand containing parent CDS and exact
  agreement among nucleotide interval, parent CDS reading frame/codon_start,
  translated domain DNA and the corresponding CDS protein slice. Add zero-based,
  half-open protein coordinates and the proven translation. Evidence establishes
  interval consistency, not independent biological validation of old HMM boundaries.
- The completed recovery restores 2272 annotations directly plus 91 annotations whose
  old locus_tag must be mapped to their sequence/coordinate-proven parent. Record
  the old ID, new ID and evidence explicitly; never silently match by old tag alone.
  Two UOC915 annotations retain the same 1154-bp interval and an out-of-frame
  start. The original COMMENT records a GG insertion at 17584–17585 within their
  interval; the current CDS DNA and translation agree, but the old domain bounds
  do not fit that reading frame. No exact parent-protein match with transferable
  trusted boundaries was found across the nine raw datasets. Retain these two
  unresolved annotations; do not round coordinates or apply a blanket extension.
- Preserve duplicate multiplicity: the audit found 187 repeated parent/location
  groups with 190 extra annotations in 40 H BGCs. Report groups but do not merge,
  drop or downweight entries in this recovery. Any future deduplication would be
  a separately evaluated input-policy change.
- Preserve DNA, CDS translations, record IDs, GCF labels, existing native aSDomain
  and PFAM annotations and all unrelated features. Use minimal text edits, retain
  old qualifiers in the audit and record source/output hashes. Never apply the
  earlier endpoint-repair rule without its original translation evidence.
- Unresolved annotations stay unchanged as review items. If any remain unresolved,
  a complete, inspectable copy may be published with `needs_review`, but this is
  not an assertion of benchmark readiness. In schema
  `bgc_h_asdomain_reconstruction_v1`, `snapshot_complete: true` refers to validated
  output files only; `ready_for_full_h_encoding: false` explicitly blocks a claim
  of complete H input readiness. Do not silently exclude their BGCs.
- The old `pilot_gh_v1/` stays unchanged. Its CPU parser success did not test whether
  Geneious aSDomain annotations had been recovered. Recovered H needs a new
  stable-ID cohort with its own provenance; do not overwrite GBKs inside the old
  cohort or pass the recovery snapshot to a command expecting the coordinate-only
  repair schema. Existing `prepare-repaired` retains the coordinate-only schema;
  the newly authorized `prepare-reconstructed` handles the recovery snapshot.
  The recovery run itself created no new encoder queue. No new H model inference,
  clustering or V-measure is claimed.

### Server bridge for recovered H and repaired G

- `pilot_features.py prepare-reconstructed` reads H from
  `--reconstructed-root data/processed/reconstructed_h_v1`; for `--datasets G H`,
  it also reads G from `--repaired-root data/processed/repaired_gh_v1`.
  `--encoder-project` selects the CPU parser checkout, and `--output` must be a
  new directory, suggested `data/processed/pilot_gh_reconstructed_v1`.
  H-only preparation does not require G inputs. Server commands and checkpoints
  are documented in `docs/gh_features.md`.
- Prepared schema is `bgc_reconstructed_pilot_inputs_v1`; extraction run schema is
  `bgc_reconstructed_feature_run_v1`. The full H snapshot can be copied to the
  server without rerunning recovery; server-side recovery instead requires its
  raw/inventory/repaired inputs. H-only preparation may omit `--repaired-root`.
- Preserve the complete selected cohorts: G=11/3 GCFs and H=90/10 GCFs. Keep
  UOC915, all labels/singletons and duplicate aSDomain entries. The adapter copies
  GBK bytes, verifies provenance and runs CPU input checks; it performs no further
  domain repair, Pfam segmentation, length filtering or sample exclusions.
- A prepared metadata status of `complete` means mechanical preparation succeeded,
  not that H domain annotations are complete. Carry `needs_review` and
  the evidence for the two unresolved UOC915 annotations through
  prepared metadata, cohort members and extraction `run.json`. Do not mutate
  the reconstruction snapshot's status or ready_for_full_h_encoding flag.
  Prepared/config/run top-level fields are annotation_status, annotation_complete,
  unresolved_asdomain_count, unresolved_sample_ids and review_policy. Current H
  uses needs_review, false, 2, UOC915's stable ID, and retain, respectively.
  The annotation fields cover this Geneious aSDomain reconstruction only;
  earlier PFAM etc. review items remain separate in inherited_review_count.
- Preserve new actual-source GBK hashes, raw hashes and embedded audit provenance.
  H records are embedded in reconstruction_provenance/, including source_provenance/
  for the earlier repair chain; optional G records use repair_provenance/.
  Do not replace files in the old `pilot_gh_v1/` or reuse its H ESM cache. Existing
  `extract` / `validate` accept the new prepared format and require every selected
  sample. Use a new `features_gh_reconstructed_v1/` output directory.
- Keep B/E settings: ESM2 BOS/1021-aa windows, MAC member 0 consuming the same
  run's ESM cache, DNA full/FA2/BF16 without a length cap or silent windowed fallback.
  Current restored H has 48–123 tokens per BGC, total 6791, maximum 1008 aa; G has
  372 total tokens, maximum 945 aa. Neither currently needs ESM token windows.
- Local work is code and temporary CPU tests only; no formal new queue has been
  published locally and no G/H model inference has run. The user runs preparation,
  extraction and validation on the server. Keep all known annotation limitations
  in future performance reports, and adapt average/PCA/scoring only after features
  return. No additional Leiden experiments are part of this task. GitHub
  publication was subsequently authorized explicitly, as recorded above.
  The six existing B/E real feature exports have passed the updated validation.
- Actual CPU preparation and loading were verified in /tmp using all 101 real
  G/H samples and separately all 90 H-only samples. Tokens are G=372/H=6791,
  no token exceeds 1021 aa, and both unresolved annotations remain recorded.
  No torch import or model inference occurred; real input hashes are unchanged.
  The formal data/processed/pilot_gh_reconstructed_v1 directory does not exist.
- A portability check changes a temporary copy of G's recorded raw_root to a simulated server path,
  relocates the prepared queue and renames its original input copies. Intercepted
  read_json/read_tsv/sha256 operations then permit only paths inside the relocated
  queue, and loading still succeeds. Source identity relies on content hashes,
  while historical source paths are retained only as provenance.
- The initial bridge's combined run passed 47 tests: 13 new preparation, 5 new run
  orchestration, 11 existing G/H preparation and 18 existing extraction/native
  validator tests. Model subprocesses use stubs in orchestration tests; real
  input checks use CPU parsers only. These are interface/provenance/workflow
  checks, not GPU model inference tests. Log: /tmp/h-bridge-final-tests.log.
- The user subsequently reported a server-side failure after successfully rerunning
  H recovery: `Reconstruction raw provenance root differs`. The adapter incorrectly
  equated the recovery's current raw_root with its embedded coordinate-repair
  ancestor's historical raw_root. These paths may differ after moving the source
  snapshot between machines. Keep both recorded paths without editing existing
  snapshots; validate the metadata hash, sample identities and raw/processed/DNA/
  protein hashes instead. Preparation and extraction protect both raw locations.
  The already generated H snapshot remains valid and can be prepared after updating
  the scripts; the failure occurs before any prepared output is created. Previous
  portability tests covered moving G or a completed queue, but missed rerunning
  H recovery on a new machine from an older coordinate-repair snapshot.
- The regression now passes 52 tests, including five new server-rerun and
  integrity/protection cases. A temporary full copy of raw/inventory/coordinate
  inputs was used to rerun actual H recovery with a different current raw_root,
  then prepare and relocate all G=11/H=90 inputs successfully. All 2941 original
  raw files and the original source snapshots retain their hashes; no torch or
  model inference ran. Evidence: /tmp/h-raw-root-regression-tests.log and
  /tmp/h-raw-root-real-regression.json. Existing server recovery snapshots do
  not need regeneration or manual metadata edits for this fix.

The 2026-11-02 audit was originally read-only inline Python. Its major counts and
single-file parser outcomes have now been revalidated by the persistent inventory
script on 2026-11-09. Record-level RC duplicate groups are also reported; this is
not near-homology analysis and does not settle biological labels.
