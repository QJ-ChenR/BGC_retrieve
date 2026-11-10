# Handoff: BGC dataset preparation and future embedding evaluation roadmap

Last updated: 2026-11-10 (Europe/Paris).

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
stage is B/E PILOT FEATURE EXTRACTION PREPARATION, with dataset preparation still
pending for the other seven sets. The user will run real model extraction on the server.

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
  `docs/pilot_features.md`. No server/GPU run is authorized or claimed here.
- Frozen encoders remain the initial comparison baseline. The user now supports
  comparing multiple clustering methods and considering trained aggregation and
  ESM2 fine-tuning, especially if frozen representations are insufficient.
  Training is a planned future research branch, not restricted to a final fallback.
- The research directions below are future design options, not current execution
  tasks. Specific model choices, losses, splits, and hyperparameters are not final.
- B/E extraction orchestration is implemented; actual benchmark model extraction,
  clustering, retrieval evaluation and parser repairs have not been run.
- The implementation choices below are recommendations, not a finalized protocol
  approved by the user. Follow the new conversation's instructions to proceed.
  This file does not itself authorize launching server jobs or modifying models.
- Do not repeatedly request permission for work the user subsequently authorizes.
  Preserve the distinction between design-only requests and implementation work.

## Workspace map

- Benchmark workspace: `/home/qiju/projects/BIG-SCAPE/`.
  It contains data, inventory scripts/tests and documentation, not a checkout or
  installation of the BiG-SCAPE clustering software.
- Raw benchmark data: `data/raw/curated_datasets/` (moved from
  `data/curated_datasets/` on 2026-11-09; do not use the obsolete path).
- Dataset source references: `data/raw/curated_datasets/README.md`.
- All derived data: `data/processed/`; current inventory: `data/processed/inventory_v1/`.
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
   `curated_benchmarked.tar.gz`. The latter two were not present locally in the
   audited directory. Verify versions and effective IDs before reusing outputs.
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

## Current stage: B/E pilot and remaining dataset preparation

The first inventory step is complete. B/E are selected for a frozen feature pilot;
label adjudication, parser repairs and cohort selection for other sets remain open.
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
cohort with coverage counts. B/E additionally has server extraction scripts and a
prepared cohort. No GPU encoding results, clustering, retrieval scores, model
training or DNA/protein fusion results have been produced.

## Future research roadmap (design only; not current tasks)

The following stages record discussed options for later work. Only the B/E server
feature-script subset is currently authorized and implemented; the other directions
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

The user explicitly wants multiple clustering methods. Average linkage remains
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

For datasets other than B/E, subsequent preparation should inspect conflicts and failure
records to define evidence-based label/compatibility policies, then produce an
audited evaluation cohort. Do not jump to future clustering, fine-tuning, or fusion
because those branches are documented above.
B and E are the user-selected pilots (39 and 281 files). F remains a candidate
requiring an explicit duplicate policy. Extend to C/D/G/H after compatibility
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

Pending server outputs: B/E per-route feature exports. Future outputs beyond this pilot: per-query retrieval metrics,
per-family/per-dataset summaries, multi-method clustering comparisons and parameter
curves, adaptation/fusion ablations, and an analysis report with reproduction
commands. These later analyses are outside the current B/E feature-script task.

The 2026-11-02 audit was originally read-only inline Python. Its major counts and
single-file parser outcomes have now been revalidated by the persistent inventory
script on 2026-11-09. Record-level RC duplicate groups are also reported; this is
not near-homology analysis and does not settle biological labels.
