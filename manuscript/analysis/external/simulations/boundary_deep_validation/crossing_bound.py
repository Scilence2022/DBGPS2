"""Model-specific rigorous coverage bound using adjacent downcrossings."""
import numpy as np
from deep_simulate import window_retention


def poisson_cdf(lam,c):
    term=np.ones_like(lam);series=term.copy()
    for i in range(1,c+1):term=term*lam/i;series+=term
    return np.exp(-lam)*series


def evaluate(job,k,c,abundance=None):
    start=35-k
    mu=200*job['depth']/window_retention(job['geometry'],job['length'],job['recovery'],1).sum()
    lam=mu*window_retention(job['geometry'],job['length'],job['recovery'],k)[start:start+148]*(1-job['error'])**k
    shared=mu*window_retention(job['geometry'],job['length'],job['recovery'],k+1)[start:start+147]*(1-job['error'])**(k+1)
    if abundance is not None:
        lam=abundance[:,None]*lam[None,:];shared=abundance[:,None]*shared[None,:]
    assert np.all(shared<=np.minimum(lam[...,:-1],lam[...,1:])+1e-9)
    left=np.maximum(0,lam[...,:-1]-shared);right=np.maximum(0,lam[...,1:]-shared)
    f=poisson_cdf(lam,c);joint=np.zeros_like(shared);pmf=np.exp(-shared)
    for z in range(c+1):
        if z:pmf=pmf*shared/z
        joint+=pmf*poisson_cdf(left,c-z)*poisson_cdf(right,c-z)
    drops=np.maximum(0,f[...,1:]-joint)
    lower=np.maximum(0,1-f[...,0]-drops.sum(axis=-1))
    old=np.maximum(0,1-f.sum(axis=-1));upper=1-f.max(axis=-1)
    assert np.all(lower>=old-1e-9) and np.all(lower<=upper+1e-9)
    return dict(lower=lower,upper=upper,window_union_lower=old,
                first_failure=f[...,0],expected_downcrossings=drops.sum(axis=-1),
                naive_product=np.prod(1-f,axis=-1))
