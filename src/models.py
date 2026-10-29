"""Official MOFA+ and a Gaussian copula factor sampler with extended rank likelihood.

The copula implementation is a benchmark implementation, not the bfa package.
Proper loading and residual-variance priors regularise latent rank scales.
"""
import contextlib
import io
import time
import numpy as np
from scipy.special import ndtr, ndtri
from scipy.stats import rankdata


def correlation(loadings, residual=None):
    c = loadings @ loadings.T
    c += np.diag(np.ones(len(c)) if residual is None else residual)
    sd = np.sqrt(np.diag(c))
    return c / sd[:, None] / sd[None, :]


def rank_metadata(x):
    result = []
    for j in range(x.shape[1]):
        obs = np.flatnonzero(np.isfinite(x[:, j]))
        order = obs[np.argsort(x[obs, j], kind='stable')]
        _, starts, groups = np.unique(x[order, j], return_index=True, return_inverse=True)
        result.append((order, starts, groups))
    return result


def latent_update(z, mean, metadata, rng):
    """Gibbs updates of alternating value groups; ties share bounds, not values.

    Nonadjacent groups have conditionally independent updates. Two passes keep
    every strict observed ordering valid, unlike simultaneous adjacent updates.
    Missing observations have unconstrained normal full conditionals.
    """
    n, p = z.shape
    for j, (order, starts, groups) in enumerate(metadata):
        for parity in (0, 1):
            values = z[order, j]
            lows = np.r_[-np.inf, np.maximum.reduceat(values, starts)[:-1]]
            highs = np.r_[np.minimum.reduceat(values, starts)[1:], np.inf]
            sel = groups % 2 == parity
            rows = order[sel]
            a = ndtr(lows[groups[sel]] - mean[rows, j])
            b = ndtr(highs[groups[sel]] - mean[rows, j])
            u = a + (b-a)*rng.uniform(1e-10, 1-1e-10, len(rows))
            z[rows, j] = mean[rows, j] + ndtri(np.clip(u, 1e-14, 1-1e-14))
        missing = np.ones(n, dtype=bool)
        missing[order] = False
        z[missing, j] = mean[missing, j] + rng.normal(size=missing.sum())
    return z


def pack_metadata(metadata,n,p):
    orders=[];starts=[];ids=[];first=[];last=[];parities=[]
    offset=0;group_offset=0
    for j,(order,localstarts,groups) in enumerate(metadata):
        orders.extend(order*p+j);starts.extend(localstarts+offset)
        ids.extend(groups+group_offset)
        first.append(group_offset);last.append(group_offset+len(localstarts)-1)
        parities.extend(groups%2)
        offset+=len(order);group_offset+=len(localstarts)
    orders=np.asarray(orders);mask=np.ones(n*p,bool);mask[orders]=False
    return orders,np.asarray(starts),np.asarray(ids),first,last,np.asarray(parities),np.flatnonzero(mask)


def latent_update_packed(z,mean,packed,rng,sd=None):
    """Same valid checkerboard Gibbs update, vectorised across all features."""
    order,starts,groups,first,last,parities,missing=packed
    flat=z.ravel();mu=mean.ravel()
    sigma=np.ones_like(flat) if sd is None else np.broadcast_to(sd,z.shape).ravel()
    for parity in (0,1):
        values=flat[order]
        lo=np.r_[-np.inf,np.maximum.reduceat(values,starts)[:-1]]
        hi=np.r_[np.minimum.reduceat(values,starts)[1:],np.inf]
        lo[first]=-np.inf;hi[last]=np.inf
        sel=parities==parity;rows=order[sel]
        a=ndtr((lo[groups[sel]]-mu[rows])/sigma[rows]);b=ndtr((hi[groups[sel]]-mu[rows])/sigma[rows])
        u=a+(b-a)*rng.uniform(1e-10,1-1e-10,len(rows))
        flat[rows]=mu[rows]+sigma[rows]*ndtri(np.clip(u,1e-14,1-1e-14))
    flat[missing]=mu[missing]+sigma[missing]*rng.normal(size=len(missing))
    return z


def assert_order(z, metadata, tol=1e-7):
    for j, (order, starts, groups) in enumerate(metadata):
        lo = np.maximum.reduceat(z[order, j], starts)
        hi = np.minimum.reduceat(z[order, j], starts)
        if np.any(lo[:-1] > hi[1:] + tol):
            raise AssertionError('Extended rank constraints violated')


def split_rhat(chains):
    # Rotation/sign-invariant correlation entries, split each independent chain.
    a = np.asarray(chains)
    half = a.shape[1]//2
    a = np.concatenate([a[:, :half], a[:, -half:]], axis=0)
    within = np.mean(np.var(a, axis=1, ddof=1), axis=0)
    between = half*np.var(a.mean(1), axis=0, ddof=1)
    return np.sqrt(((half-1)/half*within+between/half)/np.maximum(within, 1e-12))


def fit_copula(x, k, seed, burn=400, draws=600, thin=3, chains=2):
    start = time.perf_counter()
    n, p = x.shape
    meta = rank_metadata(x)
    packed=pack_metadata(meta,n,p)
    monitored = [(j, (j+17) % p) for j in range(0, p, 3)]
    means, traces = [], []
    all_moments = []
    for chain in range(chains):
        rng = np.random.default_rng(seed + chain*100003)
        z = np.zeros_like(x)
        for j, (order, _, _) in enumerate(meta):
            z[order, j] = ndtri((rankdata(x[order, j])-0.5)/len(order))
        h = rng.normal(size=(n, k))
        loadings = rng.normal(0, .2, (p, k))
        psi=np.ones(p)*.5
        stored, moments = [], []
        accumulator = np.zeros((p, p))
        for iteration in range(burn+draws):
            v = np.linalg.inv(np.eye(k)[None,:,:]*2 + (h.T@h)[None,:,:]/psi[:,None,None])
            means_l=np.einsum('pk,pkj->pj',z.T@h/psi[:,None],v)
            loadings=means_l+np.einsum('pij,pj->pi',np.linalg.cholesky(v),rng.normal(size=(p,k)))
            v = np.linalg.inv(np.eye(k)+(loadings.T/psi)@loadings)
            h = (z/psi)@loadings@v + rng.normal(size=(n, k))@np.linalg.cholesky(v).T
            residual=z-h@loadings.T
            psi=(.5+.5*np.sum(residual*residual,axis=0))/rng.gamma(2+n/2,size=p)
            z = latent_update_packed(z, h@loadings.T, packed, rng,np.sqrt(psi))
            if iteration >= burn and (iteration-burn) % thin == 0:
                c = correlation(loadings,psi)
                accumulator += c
                stored.append([c[a,b] for a,b in monitored])
                moments.append(np.diag(loadings.T@loadings).sum())
        assert_order(z, meta)
        traces.append(stored)
        all_moments.append(moments)
        means.append(accumulator/len(stored))
    rh = split_rhat(traces)
    return {'correlation': np.mean(means, axis=0), 'seconds': time.perf_counter()-start,
            'rhat_max': float(np.max(rh)), 'rhat_median': float(np.median(rh)),
            'trace': np.asarray(traces), 'iterations': burn+draws,
            'converged': bool(np.max(rh)<1.1)}


def fit_mofa(views, k, seed, iterations=500, gaussian=False):
    from mofapy2.run.entry_point import entry_point
    start = time.perf_counter()
    log = io.StringIO()
    # Continuous RNA is log2(CPM+1); protein is log2 abundance. Genomic dosage
    # is treated as continuous Gaussian; mutation indicators use Bernoulli.
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        ent = entry_point()
        ent.set_data_options(scale_views=False, center_groups=False)
        ent.set_data_matrix([[v.copy()] for v in views],
                            likelihoods=(['gaussian']*4 if gaussian else ['gaussian','bernoulli','gaussian','gaussian']),
                            views_names=['germline','somatic','RNA','protein'])
        ent.set_model_options(factors=k, spikeslab_weights=True, ard_weights=True,
                              spikeslab_factors=False, ard_factors=False)
        ent.set_train_options(iter=iterations, convergence_mode='medium',
                              seed=seed, verbose=False, quiet=True, dropR2=None)
        ent.build()
        ent.run()
        w = ent.model.nodes['W'].getExpectation()
        z = ent.model.nodes['Z'].getExpectation()
        tau = [np.mean(t,axis=0) if np.ndim(t)==2 else np.asarray(t)
               for t in ent.model.nodes['Tau'].getExpectation()]
        stats = ent.model.getTrainingStats()
    valid_elbo = np.asarray(stats.get('elbo', [])); valid_elbo=valid_elbo[np.isfinite(valid_elbo)]
    return {'scores': z, 'weights': w, 'precision':tau,
            'seconds': time.perf_counter()-start,
            'rhat_max': np.nan, 'rhat_median': np.nan, 'converged': 'Converged!' in log.getvalue(),
            'elbo': valid_elbo, 'log': log.getvalue()}


def normal_scores_fit(x):
    sorted_columns = [np.sort(col[np.isfinite(col)]) for col in x.T]
    return normal_scores_apply(x, sorted_columns), sorted_columns


def normal_scores_apply(x, sorted_columns):
    out = np.zeros_like(x)
    for j, vals in enumerate(sorted_columns):
        good = np.isfinite(x[:,j])
        left = np.searchsorted(vals,x[good,j],side='left')
        right = np.searchsorted(vals,x[good,j],side='right')
        u = (left+right+1)/(2*(len(vals)+1))
        out[good,j] = ndtri(np.clip(u, .5/(len(vals)+1),1-.5/(len(vals)+1)))
    return out
