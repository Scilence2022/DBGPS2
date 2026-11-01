"""Stream paired libraries; add Poisson sampling and structured substitution models."""
import math
import random
from pathlib import Path
import numpy as np
from simulate import designs, seed_of, corrupt


def window_retention(geometry, length, q, k):
    n = 200
    j = np.arange(n-k+1)
    if length >= n:
        return np.ones(len(j))
    if geometry == 'bernoulli':
        z = 1-1/length
        left, right = z**j, z**(n-j-k)
        return z**(k-1) * (q+(1-q)*(left+right-left*right))
    total = np.zeros(len(j))
    for offset in range(1, length+1):
        cuts = [0, *range(offset, n, length), n]
        for a, b in zip(cuts, cuts[1:]):
            if b-a >= k:
                total[a:b-k+1] += 1 if a == 0 or b == n else q
    return total/length


def fragment_one(seq, index, job, abundance):
    L, q, D = job['length'], job['recovery'], job['depth']
    rng = random.Random(seed_of(job['seed'], L, index, job['geometry'], 'deep-fragments'))
    expected_retained = float(window_retention(job['geometry'], L, q, 1).sum())
    molmean = D*200*abundance/expected_retained
    nrng = np.random.default_rng(seed_of(job['seed'], L, index, job['geometry'], 'poisson-M'))
    molecules = int(nrng.poisson(molmean)) if job['sampling'] == 'poisson' else None
    rows = []; coords = []; lengths = []; retained = 0; m = 0
    while (m < molecules if molecules is not None else retained < D*200*abundance):
        m += 1
        if L >= 200:
            cuts = [0, 200]
        elif job['geometry'] == 'regular':
            cuts = [0, *range(rng.randrange(1, L+1), 200, L), 200]
        else:
            cuts = [0, *[p for p in range(1, 200) if rng.random() < 1/L], 200]
        for a, b in zip(cuts, cuts[1:]):
            if a and b < 200 and rng.random() > q:
                continue
            rows.append((f'{index}_{m}_{a}_{b}', seq[a:b])); coords.append(a)
            lengths.append(b-a); retained += b-a
    return rows, coords, lengths, m, molmean


def structured_errors(rows, coords, seq, index, job):
    channel, e = job['channel'], job['error']
    # Seeds intentionally exclude library size and depth, pairing nested pools
    # and the common initial molecular/read sequence at different budgets.
    es = seed_of(job['seed'], job['geometry'], job['length'], job['recovery'], index, channel, 'deep-errors')
    if channel in ['substitution', 'mixed']:
        return corrupt(rows, es, e, channel)
    rng = np.random.default_rng(es)
    extra = rng.integers(1, 4, 200, dtype=np.uint8)
    lut = {b:i for i,b in enumerate('ACGT')}
    if channel == 'persistent':
        mask = rng.random(200) < e
        altered = ''.join('ACGT'[(lut[b]+int(extra[p]))%4] if mask[p] else b for p,b in enumerate(seq))
        out = [(name, altered[a:a+len(s)]) for (name,s),a in zip(rows, coords)]
    elif channel == 'burst':
        # Stationary two-state chain: error fraction e; mean error run 6 nt.
        beta = 1/6; alpha = e*beta/(1-e); out = []
        for name,s in rows:
            state = bool(rng.random() < e); text = []
            u = rng.random(len(s)); delta = rng.integers(1, 4, len(s))
            for p,b in enumerate(s):
                text.append('ACGT'[(lut[b]+int(delta[p]))%4] if state else b)
                state = u[p] >= beta if state else u[p] < alpha
            out.append((name, ''.join(text)))
    else:
        raise ValueError(channel)
    nbase = sum(len(s) for _,s in rows)
    nsub = sum(sum(a != b for a,b in zip(s,t)) for (_,s),(_,t) in zip(rows,out))
    return out, dict(substitutions=nsub, deletions=0, insertions=0,
                     source_bases=nbase, output_bases=nbase, empty_reads=0)


def generate(job, readpath, refpath, clean_prefix):
    ref = designs(job['seed'], job['n'])
    arng = np.random.default_rng(seed_of(job['seed'], 'deep-abundance'))
    sigma = job['abundance_sigma']
    abundance = arng.lognormal(-sigma*sigma/2, sigma, job['n']) if sigma else np.ones(job['n'])
    clean = {k:np.zeros((job['n'], 201-k), dtype=np.uint16) for k in [21,31]} if job['sampling']=='poisson' else {}
    totals = dict(reads=0, retained_bases=0, output_bases=0, substitutions=0, deletions=0,
                  insertions=0, physical_molecules=0, zero_molecule_strands=0)
    frag_hist = np.zeros(201, dtype=np.int64); base_sum=np.zeros(200, dtype=np.int64)
    strand_depth=[]; expected_m=[]
    with Path(readpath).open('w') as fo, Path(refpath).open('w') as fr:
        for row,(index,seq) in enumerate(ref.items()):
            fr.write(f'>{index}\n{seq}\n')
            pieces, coords, lengths, molecules, molmean = fragment_one(seq, index, job, abundance[row])
            reads, errors = structured_errors(pieces, coords, seq, index, job)
            totals['reads'] += len(reads); totals['retained_bases'] += sum(lengths)
            totals['physical_molecules'] += molecules; totals['zero_molecule_strands'] += molecules==0
            for name in ['output_bases','substitutions','deletions','insertions']:
                totals[name] += errors[name]
            frag_hist += np.bincount(lengths, minlength=201)
            diff=np.zeros(201,dtype=np.int64)
            for (name,text),a,L,(_,original) in zip(reads,coords,lengths,pieces):
                if text:fo.write(f'>{name}\n{text}\n')
                diff[a]+=1; diff[a+L]-=1
                if clean:
                    assert len(text)==len(original)
                    mismatch=np.r_[0,np.cumsum([a!=b for a,b in zip(text,original)])]
                    for k in [21,31]:
                        if L>=k:clean[k][row,a:a+L-k+1] += (mismatch[k:]-mismatch[:-k]==0).astype(np.uint16)
            depth=np.cumsum(diff)[:-1];base_sum+=depth;strand_depth.append(float(depth.mean()));expected_m.append(molmean)
    for k, arr in clean.items():
        arr.astype('<u2').tofile(str(clean_prefix)+f'.clean_k{k}.u16')
    observed = totals['retained_bases']/(job['n']*200)
    totals.update(observed_base_depth=observed, mean_fragment_length=totals['retained_bases']/max(1,totals['reads']),
                  realized_substitution_rate=totals['substitutions']/max(1,totals['retained_bases']),
                  realized_edit_rate=(totals['substitutions']+totals['insertions']+totals['deletions'])/max(1,totals['retained_bases']),
                  between_strand_depth_cv=float(np.std(strand_depth)/max(1e-12,np.mean(strand_depth))),
                  mean_expected_molecules=float(np.mean(expected_m)), fragment_length_histogram=frag_hist.tolist(),
                  mean_base_depth_by_position=(base_sum/job['n']).tolist())
    return totals
