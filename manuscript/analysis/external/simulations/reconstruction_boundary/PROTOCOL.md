# Reconstruction-boundary exploration

This is an exploratory mechanism study, not a new speed benchmark or proof of a universal DNA error-tolerance threshold. Production DBGPS2 source is copied unchanged into `src/`. The diagnostic profiler reuses its counter and mirrors its bounded search. Accepted sequences are independently compared with the production search on smoke inputs and deterministic query samples. References are used only to score outcomes and compute oracle diagnostics, never to prune or repair search candidates.

## Locked experiments

The original matrix (`design_lock.json`) has 1,080 counted graphs, each evaluated at k=21 or 31, c=0/1/3/5 and r=1.25/1.5/2/2.55/3/4/8/effectively-off. The latter uses the native finite value 1e6; absent seed coverage still prevents entry unless initial ratios are skipped. It is an adapted-rule ablation, not the published DBGPS binary or an independently optimized legacy baseline.

Main pools: three independent payload seeds, 128 designs each, with 200-nt strands: 18-nt primer, 16-nt known index, 140-nt payload, 8-nt CRC16 and 18-nt terminal primer. Uniform random payloads and standard CRC16/XMODEM are used; these are not experimentally synthesized libraries or a biochemical-constraint design study.

For the main fragmentation model, each original molecule is partitioned at a random phase with fixed internal spacing 25/45/80 nt; terminal fragments may be shorter. The 200 setting is an intact molecule without cutting. Internal fragments are retained independently with probability q=1/0.25/0.05; fragments touching either original end are always retained. All retained pieces enter k-mer counting, including pieces shorter than k (which cannot contribute a k-mer). This isolates a recovered-library bias; it is not a calibrated chemical-decay or library-preparation model.

Sampling continues per strand until its retained-base target reaches 10/30/100 times 200 bases. The final molecule is retained in full, so actual mean base depth can slightly exceed the target; actual depth, physical molecule count, fragment-length distribution and center/end depth ratio are recorded. q comparisons match the target retained sequencing budget, not initial physical molecule number. Nested per-strand random streams pair depths and recovery conditions. Sequencing substitutions are added afterwards at probability 0/0.01/0.03/0.10 per source base, using paired random draws across error levels.

Independent confirmation uses three separate payload seeds and 512 designs. Predefined stress cohorts use 256 designs with independent mixed edit events (half substitutions, quarter deletions, quarter insertions), or lognormal strand abundance sigma=0.9. Abundance shifts match expected mean depth; finite-sample actual means are recorded and are not forced equal. These validation pools were specified before outcome analysis and do not establish generalization to independent sequencing chemistries.

Two additional geometry arms were locked before their outcomes were inspected: Bernoulli cutting with probability 1/L at each internal bond, and matched regular partitions, both with three identical 256-design pools and the same target budget, q, error levels and k. Each arm has 144 graphs. This compares length-distribution shapes at matched cut density, rather than assuming that mean fragment length alone determines recoverability. The model arms run sequentially after the raw HT176 validation.

Raw HT176 validation uses both original paired FASTQ files, without merging or trimming, k=31, the legacy mapping, c=1/5 and all eight r settings. Prior full-raw profile outcome counts are checked where operating points overlap. It is separate from the merged-molecule reconstruction speed comparison. The local 176-day sample's public accession remains unverified.

## Exact diagnostic decomposition

For a given graph, k, c and r, let N be the number of designed strands. All coverage counts are exact canonical counts subject to the native 14-bit saturation.

* Ucov: all 148 output-defining extension k-mers (payload plus CRC) have count at least c+1. The known index/primer seed and the unreported look-ahead base are treated separately.
* Uend: Ucov holds and at least one of the four terminal successors has count at least c+1. The successor need not be the reference primer base: the production validator does not report or check that base.
* Useed: Uend holds and the finite native ratio calculation can start from its seed (positive seed coverage when skip=0).
* Uadm: the true output path passes every actual integer bound and at least one terminal successor passes the final bound. At a parent count a, child count b must satisfy max(floor(a/r),c+1) <= b <= min(floor(a*r),16383), with native floating arithmetic and initial skip policy.
* S: the native search uniquely accepts exactly the designed codeword.

Thus S <= Uadm <= Useed <= Uend <= Ucov. Losses decompose into N-Ucov (output coverage), Ucov-Uend (terminal support), Uend-Useed (seed support), Useed-Uadm (ratio pruning), and Uadm-S (search/uniqueness). Oracle values are diagnostic bounds for this specified reconstruction interface, not information-theoretic limits for all decoders or outer-code recovery.

The continuous symmetric ratio condition is recorded only as a conservative comparison. Integer truncation can retain a child rejected by that continuous approximation (for example parent 6, child 2, r=2.55). An initial diagnostic implementation incorrectly required the reference look-ahead base; its aborted development runs are segregated under `excluded_initial_fixed_lookahead_oracle/` and excluded from all analysis. The locked primary parameter matrix was retained, and every condition was restarted with the corrected oracle.

## Interpretation safeguards

Near-ceiling examples must have a non-trivial ceiling; zero recovery with a zero ceiling is not evidence of practical success. A two-percentage-point gap is an exploratory reporting convention, not a learned universal threshold. Best-r or best-c choices made using truth are descriptive oracle choices, not deployable tuning rules. The simple clean-k-mer survival factor (1-e)^k applies to independent substitutions only; overlapping window counts are correlated. The synthetic pools are small, so conclusions about cross-strand collision pressure require larger-scale validation. Native timing fields are diagnostics and are not used for new speed claims.
