"""Dimension-aware recovery of a prespecified generating factor subspace."""
import numpy as np


def basis(x):
    x=np.asarray(x,dtype=float)-np.mean(x,axis=0)
    u,s,_=np.linalg.svd(x,full_matrices=False)
    if not len(s) or s[0]==0:return u[:,:0]
    tol=max(x.shape)*np.finfo(float).eps*s[0]
    return u[:,s>tol]


def full_subspace_recovery(scores,truth):
    """Sum squared principal cosines / true dimension, padding missing axes by 0."""
    qa=basis(scores);qb=basis(truth)
    if not qa.shape[1]:return 0.0
    return float(np.sum(np.linalg.svd(qa.T@qb,compute_uv=False)**2)/truth.shape[1])
