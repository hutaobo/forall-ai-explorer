"""Group-adjusted association; no extra model fitting dependencies."""
import numpy as np
import pandas as pd
from scipy import stats

def group_adjusted_cohorts(clean):
    """Keep the full cohort and label the restricted two-lineage sensitivity."""
    primary=clean.loc[~clean.replicate].copy() if 'replicate' in clean else clean.copy()
    summary,residuals=group_adjusted(primary)
    cohorts={'all':{'statistics':summary,'residuals':residuals}}
    restricted=primary.loc[primary.group.isin(['BCP-ALL','T-ALL'])]
    if set(restricted.group)=={'BCP-ALL','T-ALL'} and len(restricted)<len(primary):
        summary,residuals=group_adjusted(restricted)
        summary['cohort_selection']='BCP-ALL and T-ALL only; other groups excluded, not reassigned'
        summary['excluded_from_full_n']=len(primary)-len(restricted)
        cohorts['bcp_t']={'statistics':summary,'residuals':residuals}
    return cohorts

def group_adjusted(clean, permutations=9999, bootstraps=2000, seed=20261022):
    data=clean.dropna(subset=['x','y','group']).copy()
    x=data.x.to_numpy(float); y=data.y.to_numpy(float)
    groups=[np.flatnonzero(data.group.to_numpy()==g) for g in sorted(data.group.unique())]
    n=len(x); k=len(groups); df=n-k-1
    if df<3: raise ValueError('Too few residual degrees of freedom for group adjustment.')
    if any(len(ix)<2 for ix in groups): raise ValueError('Each included group needs at least two complete primary pairs.')
    def center(a):
        out=np.array(a,dtype=float,copy=True)
        for ix in groups: out[...,ix]-=out[...,ix].mean(axis=-1,keepdims=True)
        return out
    def corr(a,b):
        den=np.sqrt(np.sum(a*a,axis=-1)*np.sum(b*b,axis=-1))
        return np.divide(np.sum(a*b,axis=-1),den,out=np.full(np.shape(den),np.nan),where=den>0)
    xr,yr=center(x),center(y)
    rx,ry=center(stats.rankdata(x)),center(stats.rankdata(y))
    r=float(corr(xr,yr)); rho=float(corr(rx,ry))
    if not np.isfinite(rho): raise ValueError('No within-group variation remains.')
    rng=np.random.default_rng(seed)
    perm=np.empty((permutations,n))
    for ix in groups:
        order=np.argsort(rng.random((permutations,len(ix))),axis=1)
        perm[:,ix]=ry[ix][order]
    null=corr(np.broadcast_to(rx,perm.shape),perm)
    p_perm=float((1+np.sum(np.abs(null)>=abs(rho)-1e-12))/(permutations+1))
    # Stratified bootstrap preserves observed group sizes; rerank each resample.
    ids=np.concatenate([rng.choice(ix,size=(bootstraps,len(ix)),replace=True) for ix in groups],axis=1)
    bx=stats.rankdata(x[ids],axis=1); by=stats.rankdata(y[ids],axis=1)
    offset=0
    for ix in groups:
        sl=slice(offset,offset+len(ix)); offset+=len(ix)
        bx[:,sl]-=bx[:,sl].mean(axis=1,keepdims=True)
        by[:,sl]-=by[:,sl].mean(axis=1,keepdims=True)
    boot=corr(bx,by); ci=np.quantile(boot[np.isfinite(boot)],[0.025,0.975])
    z=pd.get_dummies(data.group,drop_first=True,dtype=float).to_numpy()
    design=np.column_stack([np.ones(n),x,z])
    if np.linalg.matrix_rank(design)!=design.shape[1]: raise ValueError('Singular adjusted design.')
    inv=np.linalg.inv(design.T@design); beta=inv@design.T@y
    error=y-design@beta; leverage=np.sum((design@inv)*design,axis=1)
    hc3=inv@(design.T@((error/(1-leverage))[:,None]**2*design))@inv
    se=float(np.sqrt(hc3[1,1])); critical=stats.t.ppf(0.975,df)
    slope=float(beta[1]); p_slope=float(2*stats.t.sf(abs(slope/se),df))
    p_pearson=float(2*stats.t.sf(abs(r)*np.sqrt(df/max(1-r*r,1e-15)),df))
    total=float(np.sum((y-y.mean())**2)); group_sse=float(yr@yr); full_sse=float(error@error)
    summary={'n':n,'group_counts':data.groupby('group').size().to_dict(),'residual_df':df,
        'partial_Spearman':{'coefficient':rho,'permutation_p':p_perm,'bootstrap_95_ci':ci.tolist(),
            'definition':'Pearson correlation of global average-rank residuals after categorical group adjustment',
            'permutations':permutations,'bootstrap_resamples':bootstraps,'seed':seed},
        'partial_Pearson':{'coefficient':r,'nominal_t_p':p_pearson},
        'linear_model':{'formula':'y ~ x + categorical(group); common within-group slope',
            'slope':slope,'slope_HC3_se':se,'slope_HC3_95_ci':[slope-critical*se,slope+critical*se],
            'slope_HC3_t_approx_p':p_slope,'group_only_R2':1-group_sse/total,
            'full_R2':1-full_sse/total,'incremental_R2':(group_sse-full_sse)/total},
        'scope':'Exploratory adjustment for supplied group only; independent cell-line observations and within-group exchangeability assumed. Not causal independence. Spearman primary; Pearson/regression sensitivity; p values not multiplicity adjusted.'}
    data['x_group_residual']=xr; data['y_group_residual']=yr
    data['x_rank_group_residual']=rx; data['y_rank_group_residual']=ry
    return summary,data
