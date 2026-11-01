# Deep validation: prospective protocol

This phase follows the first exploration, whose outcomes were already examined. The new design, policy definitions and measurement sources are locked before measuring these independent pools. It is not a retrospective blind test of the first-phase results. No production algorithm is changed.

## Questions and locked tests

1. Can a reference-free rule select a useful ratio window? Evaluate five fixed policies: c=1/r=2.55; maximize the count of uniquely CRC-accepted codewords over r at c=1; maximize the same objective over the c/r grid; per-index union of all c=1 successes only if their complete accepted codewords agree; and the latter with at least two distinct r settings supporting a candidate. Ties select the smaller finite r, with off last; joint c/r ties first select the higher c. The repeated settings are correlated, so two-setting support is not two independent confirmations or a calibrated error probability. Sweeps require more computation than a single fixed setting.
2. Does error structure matter beyond average error rate? Compare independent substitutions, a stationary two-state burst model (mean error run 6 nt, stationary error fraction 3%), and independent template-position substitutions copied identically into every molecular descendant. Positions include the primers and index. These are controlled extremes, not fitted chemical or sequencing models. Matched mean error probability does not imply exactly identical realized error fractions. Evaluate mixed insertions/deletions separately with the first-phase definition.
3. Does scale change the window? Three independent random-payload families contain nested 1,024/8,192/65,536-design pools. Conditions: regular cut spacing 45 nt, target retained depth 100×, q=0.05, independent substitution probability 1%. Nested sizes share the initial design sequences and per-strand random streams; sizes are paired, not independent replicates. The three families are the replicates. This remains below million-design scale.
4. Can model-specific coverage bounds predict failure? Use a Poisson number of physical molecules per strand, with mean chosen to give expected retained-base depth 10/30/100×. Cross regular/Bernoulli cuts, nominal spacing 25/45/80 nt, q=1/0.05, and independent substitution probability 1%. Record origin-correct k-window counts separately from all native graph counts. Analyze the latter separately because collisions or error-derived matching words may add support.

The exact job list has 165 synthetic read inputs, each counted at k=21 and 31: 330 graphs, each with c=1/3/5 and r=1.25/1.5/2/2.55/3/4/8/off (native ratio 10^6). The five policies are frozen in `blind_policy.py`, which receives only native acceptance statuses, accepted full codewords, and the c/r settings. It does not receive reference sequences, correctness labels, coverage ceilings, true-path features, or expected error parameters.

Holdout, error-structure, mixed-error and abundance cohorts use three new seeds (7701/8802/9903), 1,024 designs each. Scale families use 17701/18802/19903. Poisson theory uses 11011/12012/13013, 1,024 designs each. The abundance stress uses lognormal sigma=1.5 with expected abundance one and Poisson molecules; its finite realized mean is recorded. Homogeneous Poisson theory must not be applied to that mixture without conditioning on abundance.

## Theory to test

For a window beginning at j (zero-based) in a 200-nt strand, under independent Bernoulli cutting p=1/L, its retention probability is

`w_j = (1-p)^(k-1) [q + (1-q)(A_j+B_j-A_j B_j)]`,

where `A_j=(1-p)^j`, `B_j=(1-p)^(200-j-k)`. A fragment touching either original end is always retained; other fragments are retained with probability q. For regular random-phase partitions, w_j is the exact average over the L possible offsets. Intact molecules have w_j=1.

Let `E_retained = sum_j w_j(k=1)` and `mu=200 D/E_retained`. Poisson thinning gives origin-correct count marginal `C_j ~ Poisson(lambda_j)`, `lambda_j=mu w_j (1-e)^k`, for independent substitutions. This does not assert independence between neighboring windows.

For output-defining windows J and cutoff c, `f_j=P(Poisson(lambda_j)<=c)`. Without assuming window independence,

`max(0,1-sum_J f_j) <= P(all C_j>=c+1) <= 1-max_J f_j`.

These are model-specific probability bounds on origin-correct output-path coverage, not per-library deterministic bounds and not universal decoder/file-recovery limits. Seed and terminal availability are separate. Compare empirical frequencies with uncertainty; a finite random library may fall outside an expectation bound by sampling variation. A product of marginal probabilities is retained only as a naive dependence-ignoring comparator.

## Validation and interpretation

The new profiler only adds pre-scoring native codeword/status recording. Smoke tests compare outputs with the frozen parent profiler and the production search. A second reference containing different valid payloads but the same index range must leave all blind statuses and candidate bytes unchanged. Independent Python scoring must match native correct/wrong totals. Preserve the first phase's exact integer and terminal-look-ahead definitions.

Evaluate every locked policy on every condition, including zero-ceiling failures and erroneous accepted outputs. The default and oracle-best single setting are comparators; only truth-free policies are deployable candidates. Do not tune policies after reading these held-out outcomes and call the same data independent validation. Report false-positive counts, yield, regret to the truth-selected best grid point, and conflicts; maximizing CRC acceptance alone does not guarantee correctness. Confidence calculations across strands do not substitute for independent experimental libraries.

Native diagnostic times are not a new speed benchmark. No speed advantage is claimed for a multi-setting sweep. The original 176-day raw library can validate reference-free policy behavior externally, but it was already inspected in phase one and is not an unseen test library.
