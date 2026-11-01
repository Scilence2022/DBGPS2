"""Score locked blind outputs, then evaluate reference-assisted mechanisms."""
from pathlib import Path
import csv
import binascii
import hashlib
import json
import math
import numpy as np
import pandas as pd
from blind_policy import select
from deep_simulate import window_retention
from simulate import seed_of
from run_deep import CS,RS,digest

R=Path(__file__).resolve().parent


def packed_reference(path, minid, n):
    expected=np.zeros((n,41),dtype=np.uint8);present=np.zeros(n,dtype=bool)
    ident=None;seq=[]
    def add(ident,text):
        i=int(ident)-minid
        if not 0<=i<n:return
        word=text[18:182];assert len(word)==164
        val=np.array(['ACGT'.index(a) for a in word],dtype=np.uint8).reshape(41,4)
        expected[i]=(val[:,0]<<6)|(val[:,1]<<4)|(val[:,2]<<2)|val[:,3];present[i]=True
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                if ident is not None:add(ident,''.join(seq))
                ident=line[1:].split()[0];seq=[]
            else:seq.append(line.strip())
    if ident is not None:add(ident,''.join(seq))
    return expected,present


def score_record(path):
    rec=json.loads(path.read_text());prefix=str(path.with_suffix(''))
    meta=json.loads(Path(prefix+'.blind_meta.json').read_text());n=meta['n_id']
    public=[dict(c=c,r=r) for r in RS for c in CS]
    status=np.fromfile(prefix+'.blind_status.u8',dtype=np.uint8).reshape(24,n)
    candidates=np.fromfile(prefix+'.accepted.u8',dtype=np.uint8).reshape(24,n,41)
    # Complete all policy decisions before reference loading or scoring.
    choices=select(status,candidates,public)
    reference=rec['command'][rec['command'].index('-x')+1]
    local=R/'references'/Path(reference).name
    if local.exists():reference=local
    elif rec['settings']['cohort']=='real_HT176_raw' and (R/'reference_evaluation_only/HT176_ref.fa').exists():
        reference=R/'reference_evaluation_only/HT176_ref.fa'
    expected,present=packed_reference(reference,meta['min_id'],n)
    N=int(present.sum());assert N==rec['settings']['n']
    correct_matrix=(status==0)&present[None,:]&np.all(candidates==expected[None,:,:],axis=2)
    grid=[]
    for j,row in enumerate(rec['grid']):
        row={k:float(v) if k not in ['c','skip'] else int(v) for k,v in row.items()}
        assert row['c']==public[j]['c'] and row['r']==public[j]['r']
        correct=int(correct_matrix[j].sum());wrong=int(((status[j]==0)&present&~correct_matrix[j]).sum())
        spur=int(((status[j]==0)&~present).sum())
        assert (correct,wrong,spur)==(int(row['correct']),int(row['wrong']),int(row['spurious_nonref']))
        assert int((status[j]==0).sum())==correct+wrong+spur
        row.update(rec['settings'],k=rec['k'],case_id=path.stem)
        for field in ['correct','wrong','ceil_cov','ceil_adm','ratio_losses','search_gap','overflow_adm','ambiguous_adm']:
            row[field+'_pct']=100*row[field]/N
        grid.append(row)
    c1=[j for j,s in enumerate(public) if s['c']==1]
    oracle_r=max(int(correct_matrix[j].sum()) for j in c1)
    oracle_cr=max(int(correct_matrix[j].sum()) for j in range(24))
    cov1=int(rec['grid'][next(j for j in c1 if public[j]['r']==2.55)]['ceil_cov'])
    rows=[]
    for name,choice in choices.items():
        accepted=choice['accepted'];matches=np.all(choice['candidates']==expected,axis=1)&present
        S=int((accepted&matches).sum());W=int((accepted&present&~matches).sum());spur=int((accepted&~present).sum())
        j=choice['selected_combo']
        rows.append(dict(**rec['settings'],k=rec['k'],case_id=path.stem,policy=name,
                         selected_c=public[j]['c'] if j is not None else None,
                         selected_r=public[j]['r'] if j is not None else None,
                         correct=S,wrong=W,spurious_nonref=spur,false_accepted=W+spur,
                         accepted=int(accepted.sum()),conflicts=choice.get('conflicts',0),
                         correct_pct=100*S/N,false_accepted_pct=100*(W+spur)/N,
                         coverage_ceiling=cov1,coverage_ceiling_pct=100*cov1/N,
                         gap_to_cov_pp=100*(cov1-S)/N,
                         oracle_r_correct=oracle_r,oracle_cr_correct=oracle_cr,
                         regret_r_pp=100*(oracle_r-S)/N,regret_cr_pp=100*(oracle_cr-S)/N))
        assert S<=cov1
    forensic=[]
    truth=np.memmap(prefix+'.truth.u16',dtype='<u2',mode='r',shape=(n,150))
    alts=np.memmap(prefix+'.alts.u16',dtype='<u2',mode='r',shape=(n,149,3))
    for j,setting in enumerate(public):
        false=(status[j]==0)&~correct_matrix[j]
        for i in np.flatnonzero(false):
            word=candidates[j,i];v=truth[i];c=setting['c'];r=setting['r'] or 1e6
            coverage=bool(np.all(v[1:149]>=c+1))
            def eligible(a,b):return max(int(int(a)/r),c+1)<=int(b)<=min(int(int(a)*r),16383)
            admissible=all(eligible(a,b) for a,b in zip(v[:148],v[1:149])) and any(eligible(v[148],b) for b in [v[149],*alts[i,148]])
            assert not admissible or not present[i]
            differing=None
            if present[i]:
                delta=word^expected[i]
                differing=int(sum(np.count_nonzero((delta>>shift)&3) for shift in [0,2,4,6]))
                assert np.array_equal(word[:4],expected[i,:4])
            if rec['settings']['cohort']!='real_HT176_raw':
                raw=word.tobytes();assert binascii.crc_hqx(raw[:39],0)==int.from_bytes(raw[39:41],'big')
            policies=[name for name,ch in choices.items() if ch['accepted'][i] and np.array_equal(ch['candidates'][i],word)]
            forensic.append(dict(case_id=path.stem,cohort=rec['settings']['cohort'],seed=rec['settings']['seed'],
                n=N,k=rec['k'],c=c,r=setting['r'],index=int(meta['min_id']+i),
                reference_present=bool(present[i]),nucleotide_hamming=differing,
                true_output_min_count=int(v[1:149].min()),coverage_eligible=coverage,
                true_rule_admissible=bool(admissible),candidate_hex=word.tobytes().hex(),
                policies_accepting_this_candidate=';'.join(policies)))
    return rows,grid,rec,forensic


def poisson_failure(lam,c):
    term=np.ones_like(lam);series=term.copy()
    for j in range(1,c+1):term=term*lam/j;series+=term
    return np.exp(-lam)*series


def theory_record(rec,path):
    job=rec['settings'];k=rec['k'];key=path.stem.split('_k')[0]
    clean=np.fromfile(R/'clean_counts'/f'{key}.clean_k{k}.u16',dtype='<u2').reshape(job['n'],201-k)
    output=clean[:,35-k:35-k+148]
    native=np.fromfile(str(path.with_suffix(''))+'.truth.u16',dtype='<u2').reshape(job['n'],150)[:,1:149]
    assert np.all(native>=np.minimum(output,16383))
    ret=window_retention(job['geometry'],job['length'],job['recovery'],k)[35-k:35-k+148]
    expected_bases=float(window_retention(job['geometry'],job['length'],job['recovery'],1).sum())
    mu=200*job['depth']/expected_bases;lam=mu*ret*(1-job['error'])**k
    rows=[]
    for c in CS:
        f=np.clip(poisson_failure(lam,c),0,1)
        clean_ok=(output>=c+1).all(axis=1);native_ok=(native>=c+1).all(axis=1)
        assert np.all(~clean_ok|native_ok)
        rows.append(dict(**job,k=k,c=c,case_id=path.stem,
                         clean_complete=int(clean_ok.sum()),native_complete=int(native_ok.sum()),
                         native_extra_complete=int(native_ok.sum()-clean_ok.sum()),
                         clean_complete_pct=100*float(clean_ok.mean()),
                         native_complete_pct=100*float(native_ok.mean()),
                         lower=max(0,1-f.sum()),upper=1-f.max(),naive_product=float(np.prod(1-f)),
                         min_lambda=float(lam.min()),max_lambda=float(lam.max()),
                         marginal_mean_rmse=float(np.sqrt(np.mean((output.mean(axis=0)-lam)**2))),
                         observed_mean_count=float(output.mean()),predicted_mean_count=float(lam.mean())))
    return rows


def exact_binomial_interval(count,n,alpha):
    """Clopper-Pearson interval by stable log-CDF inversion, no SciPy needed."""
    j=np.arange(n+1,dtype=float)
    logchoose=np.array([math.lgamma(n+1)-math.lgamma(a+1)-math.lgamma(n-a+1) for a in range(n+1)])
    def cdf(x,p):
        if p<=0:return 1.
        if p>=1:return float(x>=n)
        terms=logchoose[:x+1]+j[:x+1]*math.log(p)+(n-j[:x+1])*math.log1p(-p)
        m=terms.max();return min(1.,float(math.exp(m)*np.exp(terms-m).sum()))
    def invert(x,target,lo,hi):
        for _ in range(60):
            mid=(lo+hi)/2
            if cdf(x,mid)>target:lo=mid
            else:hi=mid
        return (lo+hi)/2
    low=0. if count==0 else invert(count-1,1-alpha/2,0.,count/n)
    high=1. if count==n else invert(count,alpha/2,count/n,1.)
    return low,high


def abundance_theory_record(rec,path):
    """Post hoc conditional extension; no homogeneous-binomial assumption."""
    job=rec['settings'];k=rec['k'];key=path.stem.split('_k')[0];n=job['n']
    clean=np.fromfile(R/'clean_counts'/f'{key}.clean_k{k}.u16',dtype='<u2').reshape(n,201-k)[:,35-k:35-k+148]
    native=np.fromfile(str(path.with_suffix(''))+'.truth.u16',dtype='<u2').reshape(n,150)[:,1:149]
    assert np.all(native>=np.minimum(clean,16383))
    sigma=job['abundance_sigma'];rng=np.random.default_rng(seed_of(job['seed'],'deep-abundance'))
    abundance=rng.lognormal(-sigma*sigma/2,sigma,n)
    w=window_retention(job['geometry'],job['length'],job['recovery'],k)[35-k:35-k+148]
    mu=200*job['depth']/window_retention(job['geometry'],job['length'],job['recovery'],1).sum()
    lam=abundance[:,None]*mu*w[None,:]*(1-job['error'])**k
    out=[]
    for c in CS:
        f=np.clip(poisson_failure(lam,c),0,1)
        low=np.maximum(0,1-f.sum(axis=1));high=1-f.max(axis=1)
        out.append(dict(**job,k=k,c=c,case_id=path.stem,
            clean_complete=int((clean>=c+1).all(axis=1).sum()),
            native_complete=int((native>=c+1).all(axis=1).sum()),
            model_lower=float(low.mean()),model_upper=float(high.mean()),
            sampled_mean_abundance=float(abundance.mean())))
    return out


def main():
    lock=json.loads((R/'design_and_source_lock.json').read_text())
    for file,sha in lock['measurement_files'].items():assert digest(R/file)==sha,file
    paths=sorted(p for p in (R/'runs').glob('*_k*.json') if not p.name.endswith(('.oracle.json','.blind_meta.json')))
    assert len(paths)==330,len(paths)
    policy=[];grid=[];theory=[];forensic=[];abundance=[]
    for i,path in enumerate(paths):
        p,g,rec,fp=score_record(path);policy.extend(p);grid.extend(g);forensic.extend(fp)
        if rec['settings']['cohort']=='poisson_theory':theory.extend(theory_record(rec,path))
        if rec['settings']['cohort']=='abundance_holdout':abundance.extend(abundance_theory_record(rec,path))
        if i%40==0:print('ANALYZED',i+1,'/330',flush=True)
    # The real sample includes queried IDs not present in the known design pool.
    p,g,_,fp=score_record(R/'real_HT176_raw.json');policy.extend(p);grid.extend(g);forensic.extend(fp)
    pd.DataFrame(policy).to_csv(R/'policy_results.csv',index=False)
    pd.DataFrame(grid).to_csv(R/'native_grid_results.csv',index=False)
    pd.DataFrame(forensic).to_csv(R/'false_accept_forensics.csv',index=False)
    t=pd.DataFrame(theory);t.to_csv(R/'theory_per_pool.csv',index=False)
    keys=['geometry','length','depth','recovery','error','k','c']
    rows=[];groups=list(t.groupby(keys))
    for group,part in groups:
        N=int(part.n.sum());count=int(part.clean_complete.sum());lo,hi=exact_binomial_interval(count,N,.05/len(groups))
        lower=float(part.lower.iloc[0]);upper=float(part.upper.iloc[0])
        assert np.allclose(part.lower,lower) and np.allclose(part.upper,upper)
        rows.append(dict(zip(keys,group),N=N,clean_complete=count,empirical=count/N,
                         simultaneous_cp_low=lo,simultaneous_cp_high=hi,
                         model_lower=lower,model_upper=upper,
                         naive_product=float(part.naive_product.iloc[0]),
                         confidence_interval_intersects_bounds=hi>=lower and lo<=upper,
                         native_complete_pct=float(part.native_complete.sum()/N*100),
                         native_extra_complete=int(part.native_extra_complete.sum()),
                         min_lambda=float(part.min_lambda.iloc[0])))
    bounds=pd.DataFrame(rows);bounds.to_csv(R/'theory_bounds_validation.csv',index=False)
    at=pd.DataFrame(abundance);at.to_csv(R/'abundance_conditional_theory_per_pool.csv',index=False)
    groups_a=list(at.groupby(['length','depth','recovery','error','k','c','abundance_sigma']))
    abundance_rows=[]
    for key,part in groups_a:
        N=int(part.n.sum());observed=part.clean_complete.sum()/N
        eps=math.sqrt(math.log(2*len(groups_a)/.05)/(2*N))
        low=float(np.average(part.model_lower,weights=part.n));high=float(np.average(part.model_upper,weights=part.n))
        abundance_rows.append(dict(zip(['length','depth','recovery','error','k','c','abundance_sigma'],key),
            N=N,empirical=float(observed),model_lower=low,model_upper=high,
            hoeffding_epsilon=eps,simultaneous_low=max(0,observed-eps),simultaneous_high=min(1,observed+eps),
            interval_intersects_bounds=observed+eps>=low and observed-eps<=high))
    abounds=pd.DataFrame(abundance_rows);abounds.to_csv(R/'abundance_conditional_bounds.csv',index=False)
    p=pd.DataFrame(policy);g=pd.DataFrame(grid);sim=p[p.cohort!='real_HT176_raw']
    for name,oracle in [('global_crc_r','oracle_r_correct'),('global_crc_cr','oracle_cr_correct')]:
        selected=p[p.policy==name]
        regret=selected[oracle]-selected.correct
        assert ((regret>=0)&(regret<=selected.false_accepted)).all()
    rlock=json.loads((R/'real_design_and_source_lock.json').read_text())
    assert digest(R/'bin/blind_profile')==rlock['binary_sha256']
    assert digest(R/'blind_policy.py')==rlock['policies_sha256']
    summary=dict(graphs=330,synthetic_operating_points=len(g[g.cohort!='real_HT176_raw']),
                 real_operating_points=24,policy_evaluations=len(p),
                 new_native_search_comparisons=int(g.native_checks.sum())+json.loads((R/'smoke_validation.json').read_text())['native_search_comparisons'],
                 independent_python_candidate_scores_match_native=True,all_source_locks_match=True,
                 crc_count_objective_regret_bound_verified=True,
                 truth_permutation_blind_output_invariant=True,
                 native_counts_dominate_origin_clean_counts=True,
                 theoretical_conditions=len(bounds),theory_intersections=int(bounds.confidence_interval_intersects_bounds.sum()),
                 conditional_abundance_conditions=len(abounds),conditional_abundance_intersections=int(abounds.interval_intersects_bounds.sum()),
                 simultaneous_family_error=.05,interval_method='Bonferroni Clopper-Pearson exact binomial',
                 policy_false_accepts_over_distinct_runs=sim.groupby('policy').false_accepted.sum().astype(int).to_dict(),
                 maximum_library_size=int(sim.n.max()),cohorts=sim[sim.policy=='default'].groupby('cohort').size().astype(int).to_dict())
    (R/'validation_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)
    print('DEEP_ANALYSIS_COMPLETE',flush=True)


if __name__=='__main__':main()
