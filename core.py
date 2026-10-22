"""Independent FORALL processed-table analysis; no upstream R code copied."""
import hashlib
import json
import os
import re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

FILES = {'protein': '41467_2022_29224_MOESM4_ESM.xlsx', 'drug': '41467_2022_29224_MOESM14_ESM.xlsx',
         'cells': '41467_2022_29224_MOESM2_ESM.xlsx', 'drugs': '41467_2022_29224_MOESM7_ESM.xlsx'}
UNITS = {'protein': 'log2 relative protein abundance', 'drug': 'sDSS'}
DOI = 'https://doi.org/10.1038/s41467-022-29224-5'
# Dataset-specific label correspondence from Supplementary Data 1, 3 and 13.
# Not a general rule: h/non-h COG derivatives may be distinct cell lines elsewhere.
FORALL_PROTEIN_LABELS = {'COG-LL-319':'COG-LL-319h','COG-LL-332':'COG-LL-332h',
                       'COG-LL-394':'COG-LL-394h','COG-LL-402':'COG-LL-402h','TMD-5':'TMD5'}
MODULES=['Profile','Pairwise','Correlation scan','Cell-line profile','Group comparison','PCA / UMAP','Data tables']

def correlation_scan(frames,metadata,kind,feature,selected,second_kind,method='Spearman',harmonize=True):
    template=pairwise_table(frames,metadata,kind,feature,selected,second_kind,frames[second_kind].index[0],harmonize)
    y=frames[second_kind].reindex(columns=template.second_sample.tolist()).to_numpy(float)
    x=np.broadcast_to(frames[kind].loc[feature,selected].to_numpy(float),y.shape).copy()
    valid=np.isfinite(x)&np.isfinite(y); n=valid.sum(axis=1)
    x[~valid]=np.nan; y[~valid]=np.nan
    if method=='Spearman': x=stats.rankdata(x,axis=1,nan_policy='omit'); y=stats.rankdata(y,axis=1,nan_policy='omit')
    x-=np.divide(np.nansum(x,axis=1),n,out=np.zeros(len(n)),where=n>0)[:,None]
    y-=np.divide(np.nansum(y,axis=1),n,out=np.zeros(len(n)),where=n>0)[:,None]
    den=np.sqrt(np.nansum(x*x,axis=1)*np.nansum(y*y,axis=1))
    r=np.divide(np.nansum(x*y,axis=1),den,out=np.full(len(n),np.nan),where=(den>0)&(n>=3))
    r=np.clip(r,-1,1); p=2*stats.t.sf(np.abs(r)*np.sqrt(np.maximum(n-2,0)/np.maximum(1-r*r,1e-15)),np.maximum(n-2,1))
    out=pd.DataFrame({'feature':frames[second_kind].index,'n':n,'coefficient':r,'nominal_p':p,'q_BH':bh(p)})
    return out.sort_values('coefficient',ascending=False,na_position='last').reset_index(drop=True)

def pairwise_table(frames, metadata, kind, feature, selected, second_kind, second, harmonize=True):
    def canonical(name, layer):
        return FORALL_PROTEIN_LABELS.get(name,name) if harmonize and layer=='protein' else name
    other=frames[second_kind].loc[second]
    lookup={canonical(s,second_kind):s for s in other.index if '_BR' not in s}
    if len(lookup)!=sum('_BR' not in s for s in other.index):
        raise ValueError('Canonical sample collision; refusing ambiguous pairing.')
    out=metadata[kind].loc[selected].copy()
    out['x']=frames[kind].loc[feature,selected]
    # Preserve same-layer replicate behavior; never match a replicate to a primary.
    out['second_sample']=[s if s in other.index else lookup.get(canonical(s,kind)) for s in selected]
    out['y']=[other.get(s,float('nan')) for s in out.second_sample]
    out['label_harmonized']=[isinstance(t,str) and s!=t for s,t in zip(selected,out.second_sample)]
    return out

def within_group_correlations(clean, method):
    rows=[]
    for group, block in clean.groupby('group'):
        if len(block)<3 or block.x.nunique()<2 or block.y.nunique()<2: continue
        s=(stats.spearmanr if method=='Spearman' else stats.pearsonr)(block.x,block.y)
        rows.append({'group':group,'n':len(block),'coefficient':float(s.statistic),'nominal_p':float(s.pvalue)})
    result=pd.DataFrame(rows)
    if len(result): result['q_BH_within_groups']=bh(result.nominal_p)
    return result
GLOSSARY = '''FORALL data dictionary: SEM and REH are named leukemia CELL LINES, never error measures.
sDSS means SELECTIVE drug sensitivity score, normalized to normal bone marrow response using a modified AUC metric.
Do not call sDSS scaled/standardized DSS, invent a universal numerical range, interpret negative values as stimulated growth,
equate near-zero with statistical nonsignificance, or convert sDSS to IC50. Protein values are log2 abundance relative to
a pooled control, not RNA expression. A table of cell lines supplies no error bars, uncertainty or significance unless
explicitly computed and supplied. A higher sDSS is a higher published score; it does not establish clinical efficacy.'''

def result_context(result, parameters, sources, statistics=None):
    """Full selected observations for profiles/pairs; full aggregates for large tables."""
    full = parameters.get('module') in ('Profile', 'Pairwise', 'Cell-line profile')
    context = {'parameters': parameters, 'rows': len(result), 'sources': sources,
               'statistics': statistics or {}, 'doi': DOI}
    if full:
        context.update(data_scope='All selected rows, including missing/unmatched observations; no preview truncation.',
                       records=json.loads(result.to_json(orient='records')))
    else:
        context.update(data_scope='Aggregates computed across all selected rows; individual rows are not supplied.',
                       numeric_summary=json.loads(result.describe().to_json()) if len(result.select_dtypes('number').columns) else {},
                       missing_counts=result.isna().sum().to_dict())
    return context

def explain_result(summary, question, api_key=None):
    if summary.get('statistics',{}).get('group_adjusted'):
        return explain_adjusted_result(summary,question,api_key)
    system=GLOSSARY+' Return JSON only: {"sentences":["...","..."]}. Use concise factual sentences in the user language, enough to address the requested analyses. '
    system+='Only describe supplied data, parameters and computed statistics. Report supplied correlation method, coefficient, paired n and nominal p when available; never invent tests. No extra background, drug mechanisms, cell-line genotypes, growth claims or treatment advice. '
    system+='A numerically higher profile value alone does not establish significance. Use inferential language only when computed tests are supplied. '
    system+='When group_adjusted is supplied, compare the pooled and adjusted coefficients, report its n, permutation p and bootstrap CI. If BCP_T_sensitivity is supplied, report that separate cohort, its n and exclusions, rather than calling it the full cohort. '
    system+='These are cell-line lineage groups, not molecular subtypes. Distinguish adjusted association from individual within-group tests; use supplied BH q for those tests. A large p/q does not prove no association, and a residual association does not prove causal independence. '
    system+='Respect data_scope: full records are not a preview. Distinguish total selected, paired and excluded counts. Correlation is exploratory association, not causality; nominal p is unadjusted. Cite the supplied DOI.'
    messages=[{'role':'system','content':system},{'role':'user','content':question+'\n'+json.dumps(summary,ensure_ascii=False)}]
    forbidden=('scaled drug sensitivity','standard error of the mean','均值的标准误','缩放后的药物敏感性',
               '刺激细胞生长','促进生长','stimulates growth')
    for attempt in range(2):
        raw=chat(messages,True,api_key)
        try:
            sentences=json.loads(raw)['sentences']
            if not isinstance(sentences,list) or not sentences or any(not isinstance(x,str) for x in sentences): raise ValueError()
            text='\n\n'.join(sentences)
            if any(t in text.casefold() for t in forbidden): raise ValueError()
            return text
        except (ValueError,KeyError,TypeError):
            messages.append({'role':'user','content':'Previous output violated the evidence rules. Return JSON sentences using only supplied numeric results, cohort counts, computed uncertainty and association limitations. Do not invent missing statistics or biology.'})
    raise ValueError('Model explanation failed evidence checks; computed values unchanged.')

def explain_adjusted_result(summary, question, api_key=None):
    """Render cohort facts from computation; model supplies qualitative interpretation only."""
    s=summary['statistics']; a=s['group_adjusted']; parts=[]
    parts.append(f"Pooled {s['method']}: n={s['paired_n']}, coefficient={s['coefficient']:.3f}, nominal p={s['nominal_p']:.4g}.")
    for label,c in [('Full selected primary cohort',a),('BCP-ALL / T-ALL sensitivity',s.get('BCP_T_sensitivity'))]:
        if not c:continue
        r=c['partial_Spearman'];lo,hi=r['bootstrap_95_ci']
        groups=', '.join(f'{g}: {n}' for g,n in c['group_counts'].items())
        parts.append(f"{label}: n={c['n']} ({groups}); group-adjusted partial Spearman={r['coefficient']:.3f}; stratified permutation p={r['permutation_p']:.4g}; bootstrap 95% CI [{lo:.3f}, {hi:.3f}].")
        if 'excluded_from_full_n' in c:
            parts.append(f"Sensitivity-only exclusions: {c['excluded_from_full_n']} other-group pairs; retained in the full cohort and unadjusted scatter. No group reassignment.")
    for g in s.get('within_group',[]):
        parts.append(f"Within {g['group']}: n={g['n']}, {s['method']}={g['coefficient']:.3f}, nominal p={g['nominal_p']:.4g}, BH q={g['q_BH_within_groups']:.4g}.")
    system=GLOSSARY+' Return JSON only: {"sentences":["..."]}. Explain qualitatively in the user language. The application already renders all numeric and cohort facts verbatim. Do not repeat any numbers, group names, cohort membership or exclusions. Compare pooled and group-adjusted association, uncertainty and within-group evidence. Never infer excluded group counts from missing records. A residual association is not causal independence, does not control all confounders, and weak within-group evidence does not prove absence. These are categorical cell-line lineage groups, never molecular subtypes. Do not call the association independent; say conditional on the supplied grouping. A similar pooled adjusted estimate in the sensitivity cohort does not demonstrate a robust association inside each group. No mechanisms or treatment advice.'
    language='Chinese' if re.search(r'[\u4e00-\u9fff]',question) else 'English'
    messages=[{'role':'system','content':system},{'role':'user','content':'Computed evidence: '+json.dumps(s,ensure_ascii=False)+'\nInterpret only the pooled versus group-adjusted association and within-group uncertainty in '+language+'. Numeric facts and cohort inventory are already displayed. Omit all numbers, names, subtypes and independence/robustness claims.'}]
    for attempt in range(2):
        try:
            sentences=json.loads(chat(messages,True,api_key))['sentences']
            if not isinstance(sentences,list) or not sentences or any(not isinstance(x,str) for x in sentences):raise ValueError()
            text='\n\n'.join(sentences)
            if re.search(r'\d|BCP|T-ALL|B-ALL|EBV|亚型|subtype|独立关联|independent association|稳健|robust',text.replace('BCL2',''),re.I):raise ValueError()
            return '\n\n'.join(parts+[text,'Source: '+DOI])
        except (ValueError,KeyError,TypeError):
            messages.append({'role':'user','content':'Return qualitative sentences only. No numbers or group names: the application renders those computed facts separately.'})
    fallback='AI wording did not pass the cohort checks. The verified numeric results above remain available; retry for an AI interpretation.'
    return '\n\n'.join(parts+[fallback,'Source: '+DOI])

def load(root):
    root = Path(root)
    expected=json.loads((Path(__file__).parent/'sources.json').read_text())
    frames, metadata, sources = {}, {}, {}
    for kind, filename in FILES.items():
        path = root / filename
        spec=next(s for s in expected['datasets'] if s['filename']==filename)
        if path.stat().st_size != spec['bytes'] or hashlib.md5(path.read_bytes()).hexdigest()!=spec['md5']:
            raise ValueError('Published workbook integrity check failed: '+filename)
        raw = pd.read_excel(path, header=None, keep_default_na=False)
        sources[kind] = {'filename': filename, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                         'doi': DOI, 'license': 'CC BY 4.0'}
        if kind in ('protein', 'drug'):
            first = 8 if kind == 'protein' else 4
            labels = raw.iloc[2, 1:].astype(str).tolist()
            if len(set(labels)) != len(labels):
                raise ValueError('Duplicate published sample headers require explicit review.')
            df = raw.iloc[first:, 1:].copy()
            df.columns = labels
            df.index = raw.iloc[first:, 0].astype(str)
            if not df.index.is_unique:
                raise ValueError('Duplicate identifiers require explicit review.')
            frames[kind] = df.apply(pd.to_numeric, errors='coerce')
            metadata[kind] = pd.DataFrame({'sample': labels, 'group': raw.iloc[3, 1:].astype(str).tolist(),
                'replicate': ['_BR' in s for s in labels], 'excel_column': [excel_column(i+2) for i in range(len(labels))]}).set_index('sample',drop=False)
        else:
            frames[kind] = raw.iloc[3:].copy()
            frames[kind].columns = raw.iloc[2].astype(str)
    return frames, metadata, sources

def excel_column(n):
    label = ''
    while n:
        n, rem = divmod(n-1, 26)
        label = chr(65+rem)+label
    return label

def samples(meta, groups, include_replicates):
    mask = meta.group.isin(groups)
    if not include_replicates:
        mask &= ~meta.replicate
    return meta.index[mask].tolist()

def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    result = np.full(len(p), np.nan)
    good = np.flatnonzero(np.isfinite(p)); order = good[np.argsort(p[good])]
    if len(order):
        result[order] = np.minimum(1, np.minimum.accumulate((p[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1])
    return result

def compare(frame, meta, left, right, method):
    a, b = frame.loc[:, meta.index[meta.group == left]], frame.loc[:, meta.index[meta.group == right]]
    rows = []
    for name in frame.index:
        x, y = a.loc[name].dropna().to_numpy(), b.loc[name].dropna().to_numpy()
        if min(len(x),len(y)) < 3: continue
        if method == 'Welch t-test':
            p = stats.ttest_ind(x,y,equal_var=False).pvalue
        else:
            p = stats.mannwhitneyu(x,y,alternative='two-sided',method='auto').pvalue
        rows.append({'feature':name,'mean_difference':float(x.mean()-y.mean()),'p':p,'n_left':len(x),'n_right':len(y)})
    out = pd.DataFrame(rows)
    if len(out): out['q_BH'] = bh(out.p); out = out.sort_values('q_BH')
    return out

def embedding(frame, method, features, neighbors, distance, seed):
    matrix = frame.T
    # Complete-case features only; never silently impute expression.
    matrix = matrix.loc[:, matrix.notna().all(axis=0)]
    matrix = matrix.loc[:,matrix.var().sort_values(ascending=False).head(features).index]
    if matrix.shape[0]<4 or matrix.shape[1]<2: raise ValueError('Too few complete samples/features.')
    if method == 'PCA':
        from sklearn.decomposition import PCA
        model = PCA(n_components=2,random_state=seed)
        xy = model.fit_transform(matrix.to_numpy())
        extra = {'explained_variance_ratio':model.explained_variance_ratio_.tolist()}
    else:
        import umap
        model = umap.UMAP(n_components=2,n_neighbors=min(neighbors,len(matrix)-1),min_dist=distance,random_state=seed,n_jobs=1)
        xy = model.fit_transform(matrix.to_numpy())
        extra = {'neighbors':min(neighbors,len(matrix)-1),'min_dist':distance}
    return pd.DataFrame(xy,index=matrix.index,columns=['axis1','axis2']), {'method':method,'features':matrix.shape[1],
        'seed':seed,'scaling':'published values, centered for PCA; no unit variance scaling','missing':'drop incomplete features',**extra}

def validate_controls(raw, frames, metadata):
    """Return only executable UI fields; never display model free-form extras."""
    plan = json.loads(raw) if isinstance(raw, str) else raw
    if not isinstance(plan,dict):raise ValueError('AI plan must be a JSON object')
    def layer(value, feature=None):
        aliases={'protein':'protein','proteins':'protein','proteomics':'protein','protein expression':'protein','protein abundance':'protein','蛋白':'protein','蛋白质':'protein','蛋白表达':'protein','蛋白质表达':'protein',
                 'drug':'drug','drugs':'drug','drug sensitivity':'drug','drug response':'drug','sdss':'drug','药物':'drug','药敏':'drug','药物敏感性':'drug'}
        requested=str(value).strip().casefold()
        normalized=aliases.get(requested)
        if normalized:return normalized
        if value is not None and requested not in ('','protein/drug'):return None
        matches=[k for k in metadata if feature is not None and str(feature).strip().casefold() in {str(f).casefold() for f in frames[k].index}]
        return matches[0] if len(matches)==1 else None
    def feature_name(value, kind):
        if value is None:return None
        value=str(value).strip()
        if value in frames[kind].index:return value
        matches=[f for f in frames[kind].index if str(f).casefold()==value.casefold()]
        return matches[0] if len(matches)==1 else None
    kind=layer(plan.get('kind'),plan.get('feature'))
    module=next((m for m in MODULES if m.casefold()==str(plan.get('module')).strip().casefold()),None)
    if kind not in frames or kind not in metadata: raise ValueError('Unknown layer')
    if module not in MODULES: raise ValueError('Unknown module')
    feature = feature_name(plan.get('feature'),kind)
    if plan.get('feature') is not None and feature is None: raise ValueError('Unknown feature')
    if module in ('Profile','Pairwise','Correlation scan') and feature is None:raise ValueError('A feature is required')
    available = metadata[kind].group.unique().tolist()
    groups = plan.get('groups') or available
    if not isinstance(groups,list) or any(not isinstance(g,str) or g not in available for g in groups): raise ValueError('Unknown groups')
    result={'kind':kind,'module':module,'groups':list(dict.fromkeys(groups)),**({'feature':feature} if feature else {})}
    if module=='Pairwise':
        sk=layer(plan.get('second_kind'),plan.get('second_feature'))
        if sk not in metadata:raise ValueError('Unknown second layer')
        sf=feature_name(plan.get('second_feature'),sk)
        if sf is None: raise ValueError('Unknown second feature')
        result['second_kind']=sk
        if sf:result['second_feature']=sf
        method=str(plan.get('method',plan.get('correlation_method',''))).strip().casefold()
        if method not in ('spearman','pearson'):raise ValueError('Unknown correlation method')
        result['correlation_method']=method.title()
    if module=='Cell-line profile' and plan.get('cell_line') in metadata[kind].index:result['cell_line']=plan['cell_line']
    return result

def plan_controls(question, frames, metadata, api_key=None):
    """Normalize unambiguous names and repair an invalid model plan once."""
    options={k:metadata[k].group.unique().tolist() for k in metadata}
    instruction=GLOSSARY+' Translate the question into a single executable JSON object, without explanation. '
    instruction+='kind is exactly protein or drug; module is one of '+json.dumps(MODULES)+'. '
    instruction+='Use Profile for one protein/drug expression plot, Pairwise for a two-feature correlation, Correlation scan for one versus many, and Group comparison only for an explicit test between at least two groups. '
    instruction+='Set kind, module, feature, groups (array). For Pairwise require second_kind, second_feature, method (Spearman/Pearson). Keep the first mentioned variable on x and the second on y. '
    instruction+='A request for all paired cell lines must use all available groups in the first layer; mentioning BCP-ALL/T-ALL adjustment does not exclude other groups from the full cohort. The app separately computes the BCP-ALL/T-ALL sensitivity. '
    instruction+='Example: {"kind":"protein","module":"Pairwise","feature":"BCL2","groups":[],"second_kind":"drug","second_feature":"Venetoclax","method":"Spearman"}. Empty groups means all available groups. '
    instruction+='Available groups: '+json.dumps(options)+'. For Cell-line profile optionally specify cell_line. Never invent data or statistics.'
    messages=[{'role':'system','content':instruction},{'role':'user','content':question}]
    for attempt in range(2):
        raw=chat(messages,True,api_key)
        try:
            plan=validate_controls(raw,frames,metadata)
            requested=re.search(r'(?:全部|所有|all)\s*(\d+)\s*(?:个)?\s*(?:细胞系|(?:paired\s+)?cell[\s-]*lines)',question,re.I)
            if plan['module']=='Pairwise' and requested:
                selected=samples(metadata[plan['kind']],plan['groups'],False)
                paired=pairwise_table(frames,metadata,plan['kind'],plan['feature'],selected,plan['second_kind'],plan['second_feature']).dropna(subset=['x','y'])
                if len(paired)!=int(requested.group(1)):
                    raise ValueError(f"The question requests all {requested.group(1)} paired cell lines, but these groups yield {len(paired)}. Use all available groups for the full cohort; BCP/T adjustment is a separate sensitivity.")
                plan['reps']=False
            return plan
        except (ValueError,KeyError,TypeError) as error:
            if attempt:raise ValueError('AI settings could not be validated. Your analysis was not changed.') from None
            messages += [{'role':'assistant','content':raw if isinstance(raw,str) else json.dumps(raw)},
                         {'role':'user','content':f'Correct the JSON settings: {error}. Use exactly protein/drug, an exact module name, both feature names for Pairwise, and valid groups. No explanation.'}]

def chat(messages, json_mode=False, api_key=None):
    import urllib.request, urllib.error
    key = api_key or os.environ.get('SCILIFELAB_OPENLLM_API_KEY','')
    if not key and os.environ.get('OPENLLM_KEY_FILE'):
        key = Path(os.environ['OPENLLM_KEY_FILE']).read_text().strip()
    if not key: raise ValueError('OpenLLM credential not configured on the server.')
    body = {'model':os.environ.get('OPENLLM_MODEL','Mistral-Large-3-675B-Instruct-2512-NVFP4'),
            'messages':messages,'temperature':0,'max_tokens':1600,'stream':False}
    if json_mode: body['response_format']={'type':'json_object'}
    req = urllib.request.Request('https://openllm.scilifelab.se/api/chat/completions',
        data=json.dumps(body,ensure_ascii=False).encode(), headers={'Authorization':'Bearer '+key,'Content-Type':'application/json'})
    try:
        with urllib.request.urlopen(req,timeout=90) as r: result=json.loads(r.read(2000000))
    except urllib.error.HTTPError as e: raise RuntimeError(f'OpenLLM HTTP {e.code}; credential/response not logged.') from None
    return result['choices'][0]['message']['content']
