"""Preserve prior manuscript OOXML while integrating validated boundary work."""
from copy import deepcopy
from pathlib import Path
from zipfile import ZipFile
import hashlib,importlib.util,json
from lxml import etree as E
from PIL import Image
import pandas as pd

R=Path(__file__).resolve().parent;O=R.parents[2]
SOURCE=O/'docx/DBGPS2_manuscript_speed_comparison_serial_3workers.docx'
DEST=O/'docx/DBGPS2_manuscript_reconstruction_boundaries.docx'
spec=importlib.util.spec_from_file_location('h',O/'benchmarks/2026-10-19/speed_extension/docx_ooxml_helpers.py')
h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h);NS=h.NS;W=h.W


def text(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))


M='{http://schemas.openxmlformats.org/officeDocument/2006/math}'
def mr(s):
    n=E.Element(M+'r');t=E.SubElement(n,M+'t');t.text=str(s);return n
def sub(base,index):
    n=E.Element(M+'sSub');E.SubElement(n,M+'e').append(mr(base));E.SubElement(n,M+'sub').append(mr(index));return n
def power(base,index):
    n=E.Element(M+'sSup');e=E.SubElement(n,M+'e')
    for a in base:e.append(a)
    E.SubElement(n,M+'sup').append(mr(index));return n
def frac(a,b):
    n=E.Element(M+'f');E.SubElement(n,M+'num').append(mr(a));E.SubElement(n,M+'den').append(mr(b));return n
def equation(*nodes):
    p=E.Element(W+'p');pr=E.SubElement(p,W+'pPr');pr.extend([h.element('jc',val='center'),h.element('spacing',after=100)])
    para=E.SubElement(p,M+'oMathPara');math=E.SubElement(para,M+'oMath')
    for n in nodes:math.append(n)
    return p


def main():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()=='0756cd06875b268aeaaf8ccbc733505feb6bf7f3702fd9eaf499332cce187d6e'
    v=json.loads((R/'validation_summary.json').read_text());cv=json.loads((R/'crossing_validation_summary.json').read_text())
    assert v['independent_python_candidate_scores_match_native'] and cv['all_new_candidate_scores_match_native']
    assert json.loads((R/'figure_validation.json').read_text())['all_passed']
    p=pd.concat([pd.read_csv(R/'policy_results.csv'),pd.read_csv(R/'crossing_policy_results.csv')],ignore_index=True)
    sim=p[p.cohort!='real_HT176_raw'];real=p[p.cohort=='real_HT176_raw'];NGRAPH=366
    assert len(sim)==NGRAPH*5
    scale=sim[(sim.cohort=='scale')&(sim.n==65536)&(sim.k==21)]
    getscale=lambda name:scale[scale.policy==name].correct_pct
    default=getscale('default');blind=getscale('global_crc_r')
    false=lambda name:int(scale[scale.policy==name].false_accepted.sum())
    nr=sim[sim.policy=='global_crc_r'];matches=int((nr.regret_r_pp==0).sum())
    error=sim[(sim.cohort=='error_structure')&(sim.policy=='default')&(sim.k==31)&(sim.depth==30)]
    er=lambda channel:error[error.channel==channel].correct_pct
    b=pd.read_csv(R/'crossing_bounds_validation.csv');example=b[(b.cohort=='poisson_theory')&(b.geometry=='bernoulli')&(b.length==80)&(b.depth==10)&(b.recovery==.05)&(b.k==21)&(b.c==1)].iloc[0]
    with ZipFile(SOURCE) as z:parts={n:z.read(n) for n in z.namelist()};infos=z.infolist()
    root=E.fromstring(parts['word/document.xml']);body=root.find('w:body',NS);ps=body.findall('w:p',NS)
    find=lambda s:next(p for p in ps if text(p).startswith(s))
    template=find('We compared both graph pruning strategies with the prefix-clustering')
    heading=find('Controlled speed comparison with clustering and multiple alignment')
    caption_template=find('Table 7 ');legend_template=find('Fig. 13 |')
    para=lambda s:h.paragraph(s,template)
    insert_before=lambda anchor,nodes:[body.insert(body.index(anchor),node) for node in nodes]

    title=ps[0];h.replace_text(title,'Fast DNA storage reconstruction with local coverage constraints')
    pp=title.find('w:pPr',NS)
    if pp is None:pp=E.Element(W+'pPr');title.insert(0,pp)
    st=pp.find('w:pStyle',NS)
    if st is None:st=h.element('pStyle',val='Title');pp.insert(0,st)
    else:st.set(W+'val','Title')
    rp=E.Element(W+'rPr');title.find('w:r',NS).insert(0,rp);rp.append(h.element('color',val='000000'))
    abstract=(
        'DNA fragmentation and sequence errors limit archival data retrieval by disrupting sequence continuity and introducing competing reconstruction paths. '
        'We present DBGPS2, which combines an absolute k-mer count threshold with local coverage constraints and optional query and two-pass counting Bloom filters. '
        'In the most fragmented sample of a 210,000-strand thermal-aging series, exact recovery increased from 48.6% with the adapted DBGPS rule to 95.4% at matched k, reaching 97.1% with shorter k-mers, without incorrect outputs at these settings. '
        'Three high-error photolithographic libraries yielded 70.02–98.43% exact recovery and two completely recovered files. '
        'On DNA Fountain, Organick File 16 and a 176-day aging sample, fully serial or three-primary-worker controls were 8.9–31.2-fold faster than the tested prefix-clustering, alignment and voting implementation; these reconstruction timings exclude shared preprocessing and outer decoding. '
        'Controlled simulations show that fragment-length tails and error correlations change reconstruction boundaries beyond their mean length or error rate. '
        'We derive model-specific coverage bounds using adjacent-window dependence and check them in new simulated libraries. '
        'In larger simulations, a separate reference-free CRC-count sweep raised correct yield while retaining erroneous outputs, exposing a precision limit at 65,536-design scale. '
        'The results separate missing sequence support, ratio-induced pruning and bounded-search failure, while identifying both the recovery gains and precision limits of local coverage constraints.')
    h.replace_text(find('DNA degradation and sequence errors obstruct'),abstract)
    intro=(
        'To test which aspects of damage determine reconstruction boundaries, we additionally vary fragment geometry, recovered positional coverage, substitution structure, mixed edits and strand abundance in controlled simulations. '
        'An exploratory matrix is followed by independently seeded validation, larger nested libraries and a prospective check of model-specific probability bounds. '
        'The analysis separates the existence of an output path from its acceptance by the actual integer ratio rule and its recovery by a bounded search. '
        'A reference-free parameter sweep is evaluated separately from the production default and from the speed benchmarks; parameter-independent CRC validity is an observable objective, not proof of exact reconstruction.')
    insert_before(find('Results'),[para(intro)])

    # Correct the interpretation of legacy continuous/reference-interval metrics
    # without replacing their measured values or changing their original media.
    h.replace_text(find('We defined the coverage-limited ceiling'),
        'As a descriptive coverage benchmark for the original thermal series, we measured the fraction of reference intervals whose k-mers all had coverage ≥2, the minimum child count at c = 1. '
        'This DStoreQC statistic ranged from 97.3% (HT5) to 99.7% (control) (Supplementary Table 4). DBGPS2 was within 0.12–1.8 percentage points of it, whereas the gap for the adapted DBGPS rule increased to 25.6 and 48.6 points in HT4 and HT5 (Fig. 5d). '
        'The original profiles classified the remaining losses using low reference-path counts and a continuous symmetric coverage-ratio statistic (Fig. 5e,f). '
        'That reference interval includes the seed and a specified look-ahead base; the continuous ratio statistic is conservative under integer rounding. '
        'Neither is an exact upper bound or minimum ratio for the production interface, which can use any admissible unreported terminal extension and does not impose the child cutoff on its seed. '
        'The exact nested bounds below treat these conventions separately. The twofold coverage step created by overlapping mates remains a practical reason for a default above two, while the original HT4/HT5 sweeps favored a slightly looser bound near three.')
    h.replace_text(find('Recovery approaches the limit set by the k-mer coverage'),'Recovery follows reference path coverage')
    h.append_text(find('To diagnose decoding outcomes we wrote dbgps2_profile'),
        ' The quantity labelled rneeded in the original profiles is the continuous symmetric ratio burden on the specified reference path, not the exact smallest native ratio under integer truncation; the diagnostic also includes the specified reference look-ahead. The boundary-validation profiler introduced below evaluates native admissibility directly and permits any eligible terminal successor.')
    h.replace_run_text(find('To diagnose decoding outcomes we wrote dbgps2_profile'),
        'the smallest ratio bound that admits the entire true path',
        'the continuous symmetric coverage-ratio burden on the specified true path')

    nodes=[h.paragraph('Coverage availability and search limits across controlled damage models',heading)]
    nodes.append(para(
        'We first explored 1,368 counted synthetic graphs and then evaluated 330 independently seeded validation graphs at k = 21 or 31, c = 1, 3 or 5 and eight ratio settings. '
        'An additional 36 graphs from new 4,096-design pools tested a downcrossing bound derived after the initial simulations. The new phases therefore comprise 366 graphs and 8,784 parameter points. '
        'Conditions have three independently seeded design families; paired conditions and nested library sizes are not independent replicates. The 200-nt designs contain known indices, random payloads and CRC16. '
        'The models isolate fragmentation, fragment recovery and error structure, rather than fitting a biochemical digital twin.'))
    nodes.append(para(
        'For each counted graph we distinguish output support Ucov, support for any terminal extension Uend, usable seed support Useed, admissibility under the native integer ratio rule Uadm, and unique correct output S. '
        'They satisfy S ≤ Uadm ≤ Useed ≤ Uend ≤ Ucov. The output-support condition covers the 148 payload-plus-CRC extension words; the seed and unreported terminal look-ahead are separate. '
        'The successive differences identify coverage, terminal, seed, ratio-pruning and search/uniqueness losses. References score outputs and calculate these diagnostics but never guide traversal. '
        'These are interface-specific, reference-assisted bounds, not information-theoretic or complete-file recovery limits.'))
    nodes.append(para(
        'Fragment-length distributions affected support beyond their mean. In the matched exploratory comparison, regular partitions and Bernoulli cuts had nearly identical mean retained lengths of 22.32 and 22.31 nt, but their k = 31 output-support ceilings and correct recovery were 0% and 100%, respectively, at target depth 100×, q = 1 and substitution probability 1%. '
        'Regular spacing 25 nt provides no 31-mer; the Bernoulli length tail does. The length distributions and expected positional window support explain why mean fragment length below k is not itself a universal reconstruction barrier (Fig. 14a,b).'))
    nodes.append(para(
        'Under the Poisson physical-molecule model with independent substitutions, origin-correct window counts have Poisson marginals determined by whole-window retention and error survival. '
        'We derive a lower probability bound using the initial coverage failure and adjacent transitions below the cutoff, which accounts for shared-window dependence without assuming independent windows (Methods). '
        f'Across {cv["homogeneous_conditions"]} model conditions, including {cv["fresh_conditions"]} conditions from the new confirmation pools, '
        f'{cv["compatible_conditions"]} Bonferroni-adjusted binomial intervals intersected the derived probability ranges (Fig. 14c). '
        f'For Bernoulli spacing 80 nt, depth 10×, q = 0.05, e = 1%, k = 21 and c = 1, the original three-pool complete-support frequency was {100*example.empirical:.2f}%; '
        f'the downcrossing range was {100*example.model_lower:.2f}–{100*example.model_upper:.2f}%, whereas a product of marginals predicted {100*example.naive_product:.2f}% (Fig. 14d). '
        'Some probability ranges remain broad; compatibility does not establish a calibrated predictor for real DNA damage.'))
    nodes.append(para(
        'Error correlation changed recovery at similar mean error fractions (Fig. 15). With regular spacing 80 nt, q = 1 and target depth 30×, the default k = 31 reconstruction recovered '
        f'{er("substitution").mean():.2f} ± {er("substitution").std():.2f}% under independent 3% substitutions and '
        f'{er("burst").mean():.2f} ± {er("burst").std():.2f}% under the stationary burst model, but only '
        f'{er("persistent").mean():.2f} ± {er("persistent").std():.2f}% when the same template-position mutations were present in all copies. '
        'Here ± denotes sample s.d. across three pools. At 100×, independent-error recovery reached 100%, while persistent-error recovery remained unchanged. '
        'The persistent model is an extreme mechanistic control, not a fitted representation of the aged samples.'))
    nodes.append(para(
        'Larger libraries exposed a precision limit missed by the small exploratory pools (Fig. 16). At 65,536 designs, k = 21, regular spacing 45 nt, target depth 100×, q = 0.05 and substitution probability 1%, '
        f'default recovery was {default.mean():.3f} ± {default.std():.3f}%, with {false("default")} incorrect outputs across three independent families. '
        f'The fixed reference-free rule maximizing uniquely CRC-accepted outputs over r raised recovery to {blind.mean():.3f} ± {blind.std():.3f}%, '
        f'but still accepted {false("global_crc_r")} incorrect outputs. The consensus union accepted {false("union_consensus")} incorrect outputs at this scale; requiring two ratio settings retained {false("two_support_consensus")}. '
        'Agreement across correlated settings therefore does not certify correctness. Post hoc donor tracing identifies long exact segments from other designs in false candidates, consistent with graph-path mixing.'))
    nodes.append(para(
        f'The CRC-count choice attained the truth-selected best single-r correct count in {matches}/{NGRAPH} simulated grids, including {int((nr.coverage_ceiling_pct>=50).sum())} grids with output-support ceilings of at least 50%; zero-ceiling cases alone do not demonstrate practical success. '
        'This is compatible with residual false positives: maximizing total accepted output guarantees only that correct-yield regret is no larger than the number of false accepted outputs, which is unknown in deployment. '
        'For raw HT176, the same reference-free objective selected c = 1, r = 3 and recovered 208,688 designs, compared with 208,593 at the default, without incorrect selected outputs (Table 8). '
        'The full query range, including unused indices, entered the objective. This real library had already been inspected and is an external behavior check, not an unseen experimental test. '
        'The sweep is a separate diagnostic wrapper; it is not an automatic feature of the production default or a new equal-budget speed comparison.'))
    nodes.append(h.paragraph('Table 8 Reference free decisions on the original paired HT176 reads',caption_template))
    labels={'default':'Fixed r = 2.55','global_crc_r':'CRC count over r','global_crc_cr':'CRC count over c and r','union_consensus':'Consensus union','two_support_consensus':'Two setting support'}
    rows=[['Decision','Correct','False outputs','c','r','Gap to Ucov pp']]
    for name,label in labels.items():
        row=real[real.policy==name].iloc[0]
        rows.append([label,f'{int(row.correct):,}',str(int(row.false_accepted)),
                     'varies' if pd.isna(row.selected_c) else str(int(row.selected_c)),
                     'varies' if pd.isna(row.selected_r) else f'{row.selected_r:g}',f'{row.gap_to_cov_pp:.4f}'])
    table=h.table(rows);widths=[2900,1400,1300,700,700,1260]
    table.find('w:tblPr/w:tblW',NS).set(W+'w',str(sum(widths)))
    for col,width in zip(table.findall('w:tblGrid/w:gridCol',NS),widths):col.set(W+'w',str(width))
    for row in table.findall('w:tr',NS):
        for cell,width in zip(row.findall('w:tc',NS),widths):cell.find('w:tcPr/w:tcW',NS).set(W+'w',str(width))
    nodes.extend([table,para('Correct recovery is scored against 210,000 designs after all decisions. False outputs include incorrect sequences and accepted unused indices. The output-support ceiling at c = 1 is 208,855 designs. These raw paired reads were not merged or trimmed; this table is separate from the successful-merge speed benchmark. The consensus policies combine complete candidates only when all successful c = 1 settings agree; two-setting support additionally requires success at two distinct r values, which are not independent observations.')])
    insert_before(find('Discussion'),nodes)

    discussion=(
        'The controlled models identify a more general reconstruction boundary than a single mean fragment length or error rate. Whole-window retention, surviving strand abundance and error correlation jointly determine the availability of output k-mers. '
        'Adjacent windows share most of their bases and often the same molecules, so multiplying marginal support probabilities can underestimate intact-path availability severely. '
        'The Poisson thinning and downcrossing analysis gives explicit model-specific bounds without this independence assumption. Its assumptions are restrictive: molecules are sampled independently, fragmentation and recovery follow the specified channel, and substitutions are independent for the analytical result. '
        'The persistent and mixed-edit controls illustrate why these bounds must not be transferred to systematic errors or arbitrary chemistry unchanged.')
    insert_before(find('The Bloom filters address implementation cost independently'),[para(discussion),para(
        'Reference-free CRC-count selection provides an observable tuning objective, but the scale experiment shows why correct yield and precision must remain separate endpoints. '
        'A CRC-valid false path can survive several correlated ratio settings, especially if the true path is pruned or the larger search is inconclusive. '
        'The selected grid may attain the largest correct count while retaining erroneous outputs, and consensus unions can increase both yield and false acceptance. '
        'A deployable adaptive method therefore requires stronger error control and evaluation at a common total reconstruction budget; a multi-setting sweep cannot inherit the speed ratio of a single operating point.')])
    h.replace_text(find('The two evaluation regimes do not establish a universal tolerance'),
        'The experimental libraries and simulations do not establish a universal tolerance to DNA damage or a single maximum error rate. The thermal-aging analysis covers one strand design and a damage series that was not strictly monotonic. '
        'Public photolithographic data broaden the sequence-error regime but do not test the same libraries after aging. Joint severe fragmentation and high errors are examined here only through controlled models, not a newly prepared experimental library. '
        'Synthetic payloads are random, the largest new simulated pools contain 65,536 designs, and the models are not fitted to synthesis or chemical-decay measurements. '
        'The analytical comparison tests compatibility with specified models, not a universal phase transition. Prospective validation on independent experimental libraries, higher-scale specificity tests and byte-identical file recovery across sequencing budgets remain necessary.')

    methods=[h.paragraph('Controlled simulations and reference free validation',heading)]
    methods.append(para(
        'The exploratory phase comprised 1,368 counted graphs, followed by a source-locked 330-graph validation phase on new seeds and a separately locked 36-graph confirmation of the downcrossing bound. '
        'Exploratory designs used 128, 256 or 512 strands per pool. The main validation used 1,024 strands, nested scale families used 1,024, 8,192 and 65,536 strands, and the fresh theory confirmation used 4,096 strands. '
        'All synthetic designs used the same 200-nt layout, direct two-bit ACGT mapping, random payloads and CRC16/XMODEM. The index range began at 1,000. Three random design families provided independent replicates within each condition; matched conditions and nested sizes shared their underlying designs.'))
    methods.append(para(
        'Regular fragmentation partitions each original molecule with uniformly random phase and fixed internal spacing L; fragments touching either original end are retained and internal fragments are independently retained with probability q. '
        'Bernoulli fragmentation instead cuts independently at each internal bond with probability 1/L. Uncut 200-nt controls are separate from these distributions. '
        'Budget-sampled arms retain complete molecular passes until the per-design target retained-base budget is reached, recording overshoot. Poisson arms draw each physical molecule count independently with mean μ chosen for the expected retained-base budget. '
        'Abundance stress uses lognormal strand multipliers with mean one and sigma 1.5; finite realized depths, physical molecule counts, fragment histograms and positional depth are retained.'))
    methods.append(para(
        'Independent substitutions change each source base with probability e to one of the other three bases uniformly. The burst model is a stationary two-state chain with mean error run six bases and stationary error fraction e: recovery probability per error base is 1/6 and entry probability per correct base is e/[6(1−e)]. '
        'Persistent errors independently change template positions once and copy those same changes into every descendant, including primer/index positions. Mixed edits allocate half the event probability to substitutions and one quarter each to deletions and insertions. '
        'All retained reads enter native counting; no reference-based read repair, partitioning or path selection is used. Realized error fractions are measured from the generated reads.'))
    methods.append(para(
        'The diagnostic wrapper records full accepted codeword bytes and truth-independent search statuses before comparing with any reference. '
        'Five policies were frozen before the main validation: fixed c = 1/r = 2.55; maximum unique CRC-accepted count over r at c = 1; maximum count over c/r; per-index consensus union across c = 1 settings; and that union requiring at least two supporting r values. '
        'Counts include all queried indices. Ties favor a larger c for the joint policy, then a smaller finite r; effectively-off r sorts last. The cutoffs are c = 1, 3 or 5; r takes 1.25, 1.5, 2, 2.55, 3, 4, 8 or 1,000,000. '
        'The module receives no truth labels, oracle counts or channel parameters. Independent Python scoring subsequently checks correct, wrong and unused-index outputs against the native classifications.'))
    methods.append(para(
        'For the exact nested decomposition, Ucov requires every output-defining extension k-mer to have count ≥c+1; Uend additionally requires any terminal successor with this count; Useed additionally requires an available seed under skip = 0; and Uadm applies the actual integer bounds at every step and accepts any eligible unreported successor. '
        'The difference Uadm−S separates bounded-search and uniqueness failure from true-path pruning. The counter remains the unchanged 14-bit saturating native counter. '
        'Smoke tests, deterministic per-grid samples and reference-permutation checks test the diagnostic search against production search and confirm that reference changes leave blind candidate/status bytes unchanged. Seven representative generated inputs were replayed byte-identically.'))
    methods.append(para(
        'For a length-k window at zero-based position j in a 200-nt template under Bernoulli fragmentation, its retention probability is given below, with A and B describing connection to either template end. '
        'Regular-partition retention is the exact average over the L possible phases. Expected retained bases per molecule are the sum of the length-one retention probabilities; physical-molecule intensity is then set from target expected depth D as shown below. '
        'With independent substitutions and Poisson physical molecules, the origin-correct window-count marginal follows the Poisson law below. Output windows are numbered j = 0 to 147 in search order; their genomic start positions are 35 − k + j, which supplies the retention-function subscript in their count means. These counts are recorded by the generator independently of the canonical graph, whose counts can additionally include matching words from other origins.'))
    methods.append(equation(sub('w','j'),mr('(k) = '),power([mr('(1 − '),frac('1','L'),mr(')')],'k − 1'),mr('[q + (1 − q)('),sub('A','j'),mr(' + '),sub('B','j'),mr(' − '),sub('A','j'),sub('B','j'),mr(')]')))
    methods.append(equation(sub('A','j'),mr(' = '),power([mr('(1 − '),frac('1','L'),mr(')')],'j'),mr(',    '),sub('B','j'),mr(' = '),power([mr('(1 − '),frac('1','L'),mr(')')],'200 − j − k')))
    methods.append(equation(mr('μ = '),frac('200D','Eretained'),mr(',    '),sub('C','j'),mr(' ∼ Poisson('),sub('λ','j'),mr('),    '),sub('λ','j'),mr(' = μ '),sub('w','35−k+j'),mr('(k)'),power([mr('(1 − e)')],'k')))
    methods.append(para(
        'For the probability bounds, number the 148 output windows from zero. Let f be each window’s lower-tail probability. A basic union bound gives the range below. '
        'For adjacent windows, shared clean support has the stated mean. Poisson splitting decomposes the pair into a shared count and two independent private counts. '
        'The downcrossing probability is obtained exactly from the bivariate Poisson lower-tail sum. Every failed output path either starts below threshold or has a first downcrossing, giving the final lower bound below, without assuming a Markov or independence model across the full path. '
        'The bound was derived after inspecting the original window-product discrepancy, then checked on separately locked new seeds; this chronology is retained in the source record.'))
    methods.append(equation(sub('f','j'),mr(' = P('),sub('C','j'),mr(' ≤ c),    max(0, 1 − '),sub('∑','j'),sub('f','j'),mr(') ≤ P(all '),sub('C','j'),mr(' > c) ≤ 1 − '),sub('max','j'),sub('f','j')))
    methods.append(equation(sub('λ','shared,j'),mr(' = μ '),sub('w','35−k+j'),mr('(k + 1)'),power([mr('(1 − e)')],'k + 1'),mr(',    '),sub('d','j'),mr(' = P('),sub('C','j'),mr(' > c, '),sub('C','j+1'),mr(' ≤ c)')))
    methods.append(equation(mr('max(0, 1 − '),sub('f','0'),mr(' − '),sub('∑','j=0,…,146'),sub('d','j'),mr(') ≤ P(all '),sub('C','j'),mr(' > c)')))
    methods.append(para(
        'Homogeneous-model frequencies combine independent origin-support indicators across the three pools and use exact Clopper–Pearson intervals with Bonferroni adjustment over all compared homogeneous conditions (family alpha 0.05). '
        'A separate post hoc abundance extension conditions on each simulated strand multiplier and uses simultaneous Hoeffding intervals for independent, non-identically distributed indicators; it is not a homogeneous Poisson fit. '
        'Interval intersection is a compatibility check, not proof of calibration. Native graph output errors are not assigned binomial confidence intervals, because pooled graph paths couple different strands. '
        'Fresh HT176 runs use both original raw FASTQ files without trimming or merging, the legacy map, k = 31 and the same 24-point grid. The previously analyzed real library is reported separately from unseen synthetic validation and from merged-input timings.'))
    insert_before(find('Optional Bloom filters for lookup and counting'),methods)
    h.append_text(find('DBGPS2 is available at https://github.com/'),
        ' The reconstruction_boundary and boundary_deep_validation archives contain locked simulation sources, exact diagnostic extensions, independent validation settings, blind candidate bytes, origin-support arrays, per-index false-output traces, probability calculations and source tables and vector artwork for Figures 14–16 and Table 8. The production core remains unchanged; reference-free sweeps are separate analysis wrappers.')

    # Clarify the existing Fig. 5 diagnostic without altering its artwork.
    fig5=find('Fig. 5 |')
    h.replace_text(fig5,
        'Fig. 5 | Reference path support and reconstruction outcomes in the original thermal series. a, Mean graph nodes per index; an unbranched path requires 150. b, Path-overflow frequency without ratio constraint. c, Minimum true-path coverage among HT5 strands recovered by both rules, DBGPS2 only or neither; the dashed line is the best adapted-DBGPS cutoff. d, Legacy reference-interval support (all k-mers ≥2) compared with observed correct recovery. This statistic includes seed and specified reference look-ahead support and is a descriptive benchmark, rather than an exact upper bound for every production output. e, Original reference-path loss classification by low counts or steep continuous ratio transitions. f, The statistic labelled rneeded is the continuous symmetric ratio burden on the specified path, not the exact smallest native integer-window ratio. The new exact decomposition separates seed, output, any terminal successor and native admissibility (Methods). Original measured values are retained.')

    rels=E.fromstring(parts['word/_rels/document.xml.rels']);new_images={}
    pictemplate=next(p for p in reversed(ps) if p.find('.//w:drawing',NS) is not None)
    def figure(number,name,title,legend):
        path=R/f'{name}.png';pic=deepcopy(pictemplate);ppr=pic.find('w:pPr',NS)
        if ppr is None:ppr=E.Element(W+'pPr');pic.insert(0,ppr)
        for tag in ['pageBreakBefore','keepNext']:
            if ppr.find('w:'+tag,NS) is None:ppr.append(h.element(tag))
        im=Image.open(path);width=5943600;height=int(width*im.height/im.width)
        pic.find('.//wp:extent',NS).set('cx',str(width));pic.find('.//wp:extent',NS).set('cy',str(height))
        for ext in pic.findall('.//a:ext',NS):ext.set('cx',str(width));ext.set('cy',str(height))
        dp=pic.find('.//wp:docPr',NS);dp.set('id',str(max(int(x.get('id')) for x in root.findall('.//wp:docPr',NS))+1));dp.set('name',f'Figure {number}');dp.set('descr',title)
        rid='rId'+str(max(int(x.get('Id')[3:]) for x in rels if x.get('Id').startswith('rId'))+1)
        filename=f'{name}.png';E.SubElement(rels,'{http://schemas.openxmlformats.org/package/2006/relationships}Relationship',Id=rid,Type=NS['r']+'/image',Target='media/'+filename)
        pic.find('.//a:blip',NS).set('{'+NS['r']+'}embed',rid);pic.find('.//pic:cNvPr',NS).set('name',filename)
        caption=h.paragraph(f'Fig. {number} | {title}',legend_template)
        first=caption.find('w:r',NS);rp=first.find('w:rPr',NS)
        if rp is None:rp=E.Element(W+'rPr');first.insert(0,rp)
        if rp.find('w:b',NS) is None:rp.append(h.element('b'))
        h.append_text(caption,' '+legend);run=caption.findall('w:r',NS)[-1];rp=E.Element(W+'rPr');rp.extend([h.element('b',val=0),h.element('sz',val=22)]);run.insert(0,rp)
        cp=caption.find('w:pPr',NS)
        for kn in cp.findall('w:keepNext',NS):cp.remove(kn)
        cp.append(h.element('keepNext',val=0))
        ix=len(body)-1 if body[-1].tag==W+'sectPr' else len(body);body.insert(ix,pic);body.insert(ix+1,caption)
        new_images['word/media/'+filename]=path.read_bytes()
    figure(14,'Fig14_effective_coverage_model','Fragment geometry and dependent window support constrain coverage availability.',
        'a, Retained fragment-length distributions at matched nominal cut spacing 25 nt, expected depth 100× and q = 1; counts pool three Poisson-sampled libraries per model. The vertical line is k = 31. b, Expected origin-correct window counts under the two cut models at L = 45 nt, D = 100×, q = 0.05, e = 1% and k = 21. c, Empirical complete output-window support versus the model probability range from adjacent downcrossings and marginal upper bounds, at c = 1. Horizontal ranges are probability bounds, not confidence intervals. d, The same empirical frequencies versus a product of marginal probabilities. Panels c,d pool three independent libraries per condition and include the fresh bound-confirmation pools. The diagonal marks equality. Bounds concern origin-correct coverage under the specified models; seed, terminal and search decisions are separate.')
    figure(15,'Fig15_error_structure','Equal mean error fractions can produce different reconstruction boundaries.',
        'a, Realized substituted-base fractions for independent, stationary burst and persistent-template substitution models, all with expected rate 3%. b, Correct recovery with default c = 1 and r = 2.55; circles indicate k = 21 and open squares k = 31. Regular partition spacing is 80 nt and q = 1. Bars/points show means with sample s.d. across three independent 1,024-design families at target retained depths 30× and 100×. Persistent positions are copied identically into all molecules of a design, including its primers and index, and are an extreme mechanistic control rather than a calibrated decay channel. Shared mean error rate does not imply shared error correlations or clean-window availability.')
    figure(16,'Fig16_blind_selection_scale','Scale exposes false acceptance in default and reference free decisions.',
        'a,c, Correct recovery at k = 21 and 31. b,d, Erroneous accepted outputs, including any unused-index outputs. Conditions are regular spacing 45 nt, target depth 100×, q = 0.05 and independent substitution probability 1%. Three independently seeded families contain nested pools of 1,024, 8,192 and 65,536 designs. Recovery points are means with sample s.d.; error panels show means and small points for each family. CRC-count choice maximizes uniquely valid output count over the fixed r grid at c = 1. Consensus union requires agreement of all successful c = 1 settings; two-setting support also requires success at two distinct r values, which are correlated. The dotted line in recovery panels is the reference-assisted output-support ceiling. The joint c/r count policy and all individual points are in the source tables. No new speed or complete-file recovery comparison is implied.')
    updates={'word/document.xml':E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True),
             'word/_rels/document.xml.rels':E.tostring(rels,xml_declaration=True,encoding='UTF-8',standalone=True)}
    with ZipFile(DEST,'w') as z:
        for info in infos:z.writestr(info,updates.get(info.filename,parts[info.filename]))
        for name,data in new_images.items():z.writestr(name,data)
    with ZipFile(DEST) as z:
        assert z.testzip() is None
        for name,data in parts.items():
            if name.startswith('word/media/'):assert z.read(name)==data,name
    mn={'m':'http://schemas.openxmlformats.org/officeDocument/2006/math'};old=E.fromstring(parts['word/document.xml'])
    original_math=[E.tostring(n) for n in old.findall('.//m:oMath',mn)]
    final_math=[E.tostring(n) for n in root.findall('.//m:oMath',mn)]
    assert all(n in final_math for n in original_math)
    result=dict(source=str(SOURCE),output=str(DEST),sha256=hashlib.sha256(DEST.read_bytes()).hexdigest(),
                preserved_prior_media=True,preserved_prior_equations=True,figures=len(root.findall('.//w:drawing',NS)),
                tables=len(root.findall('.//w:tbl',NS)),new_simulated_graphs=NGRAPH,new_parameter_points=8784,
                changes=['title and abstract','Introduction scope','exact versus legacy coverage diagnostics',
                         'controlled-boundary Results and Table 8','error and selection limitations in Discussion',
                         'reproducible models and downcrossing Methods','Figures 14 to 16','Code availability'])
    assert result['figures']==16 and result['tables']==8
    (R/'manuscript_structural_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':main()
