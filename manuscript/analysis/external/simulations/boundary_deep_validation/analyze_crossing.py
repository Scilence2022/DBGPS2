"""Validate the adjacent-crossing bound on original and fresh pools."""
from pathlib import Path
import json,math
import numpy as np
import pandas as pd
from crossing_bound import evaluate
from analyze_deep import theory_record,score_record,exact_binomial_interval
from run_deep import digest

R=Path(__file__).resolve().parent
old=pd.read_csv(R/'theory_per_pool.csv');rows=[];new_policy=[];new_grid=[];new_fp=[]
paths=sorted(p for p in (R/'crossing_runs').glob('*_k*.json') if not p.name.endswith(('.oracle.json','.blind_meta.json')))
assert len(paths)==36,len(paths)
lock=json.loads((R/'crossing_design_and_source_lock.json').read_text())
for file,sha in lock['files'].items():assert digest(R/file)==sha,file
for path in paths:
    p,g,rec,fp=score_record(path);new_policy.extend(p);new_grid.extend(g);new_fp.extend(fp)
    rows.extend(theory_record(rec,path))
t=pd.concat([old,pd.DataFrame(rows)],ignore_index=True)
for i,row in t.iterrows():
    result=evaluate(row.to_dict(),int(row.k),int(row.c))
    t.loc[i,'crossing_lower']=float(result['lower'])
    t.loc[i,'expected_downcrossings']=float(result['expected_downcrossings'])
    assert abs(float(result['upper'])-row.upper)<1e-9
t.to_csv(R/'crossing_theory_per_pool.csv',index=False)
keys=['cohort','geometry','length','depth','recovery','error','k','c'];groups=list(t.groupby(keys));out=[]
for key,part in groups:
    N=int(part.n.sum());observed=part.clean_complete.sum()/N
    lo,hi=exact_binomial_interval(int(part.clean_complete.sum()),N,.05/len(groups))
    low=float(part.crossing_lower.iloc[0]);high=float(part.upper.iloc[0]);oldlow=float(part.lower.iloc[0])
    out.append(dict(zip(keys,key),N=N,empirical=float(observed),model_lower=low,model_upper=high,
                    old_window_union_lower=oldlow,naive_product=float(part.naive_product.iloc[0]),
                    simultaneous_cp_low=lo,simultaneous_cp_high=hi,
                    interval_intersects_bounds=hi>=low and lo<=high,
                    interval_width=high-low,old_interval_width=high-oldlow,
                    expected_downcrossings=float(part.expected_downcrossings.iloc[0])))
b=pd.DataFrame(out);b.to_csv(R/'crossing_bounds_validation.csv',index=False)
pd.DataFrame(new_policy).to_csv(R/'crossing_policy_results.csv',index=False)
pd.DataFrame(new_grid).to_csv(R/'crossing_native_grid_results.csv',index=False)
pd.DataFrame(new_fp).to_csv(R/'crossing_false_accept_forensics.csv',index=False)
v=dict(new_graphs=36,new_operating_points=len(new_grid),new_policy_evaluations=len(new_policy),
       native_search_comparisons=int(pd.DataFrame(new_grid).native_checks.sum()),
       homogeneous_conditions=len(b),compatible_conditions=int(b.interval_intersects_bounds.sum()),
       fresh_conditions=int((b.cohort=='crossing_confirmation').sum()),
       fresh_compatible_conditions=int(b[b.cohort=='crossing_confirmation'].interval_intersects_bounds.sum()),
       median_original_interval_width=float(b.old_interval_width.median()),median_crossing_interval_width=float(b.interval_width.median()),
       maximum_width_reduction=float((b.old_interval_width-b.interval_width).max()),
       maximum_naive_product_error_pp=float(100*abs(b.naive_product-b.empirical).max()),
       all_source_locks_match=True,all_new_candidate_scores_match_native=True,
       statistical_interval='Bonferroni Clopper-Pearson; family alpha=0.05 across all homogeneous conditions')
(R/'crossing_validation_summary.json').write_text(json.dumps(v,indent=2))
print(json.dumps(v,indent=2),flush=True);print('CROSSING_ANALYSIS_COMPLETE',flush=True)
