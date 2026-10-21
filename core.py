"""Independent FORALL processed-table analysis; no upstream R code copied."""
import hashlib
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

FILES = {'protein': '41467_2022_29224_MOESM4_ESM.xlsx', 'drug': '41467_2022_29224_MOESM14_ESM.xlsx',
         'cells': '41467_2022_29224_MOESM2_ESM.xlsx', 'drugs': '41467_2022_29224_MOESM7_ESM.xlsx'}
UNITS = {'protein': 'log2 relative protein abundance', 'drug': 'sDSS'}
DOI = 'https://doi.org/10.1038/s41467-022-29224-5'
GLOSSARY = '''FORALL data dictionary: SEM and REH are named leukemia CELL LINES, never error measures.
sDSS means SELECTIVE drug sensitivity score, normalized to normal bone marrow response using a modified AUC metric.
Do not call sDSS scaled/standardized DSS, invent a universal numerical range, interpret negative values as stimulated growth,
equate near-zero with statistical nonsignificance, or convert sDSS to IC50. Protein values are log2 abundance relative to
a pooled control, not RNA expression. A table of cell lines supplies no error bars, uncertainty or significance unless
explicitly computed and supplied. A higher sDSS is a higher published score; it does not establish clinical efficacy.'''

def explain_result(summary, question, api_key=None):
    system=GLOSSARY+' Return JSON only: {"sentences":["...","..."]}. At most three short factual sentences in the user language. '
    system+='Only describe the supplied table values and parameters. No extra background, drug mechanisms, cell-line genotypes, growth claims, treatment advice, or significance claims. '
    system+='A numerically higher score is simply higher, never "significantly higher". State that significance cannot be established from profile values alone. '
    system+='When the preview is partial, say it covers only the supplied preview. Cite the supplied DOI.'
    messages=[{'role':'system','content':system},{'role':'user','content':question+'\n'+json.dumps(summary,ensure_ascii=False)}]
    forbidden=('scaled drug sensitivity','standard error of the mean','均值的标准误','缩放后的药物敏感性',
               '刺激细胞生长','促进生长','显著高于','无统计学意义','无显著效应','不显著','significantly higher','stimulates growth')
    for attempt in range(2):
        raw=chat(messages,True,api_key)
        try:
            sentences=json.loads(raw)['sentences']
            if not isinstance(sentences,list) or not 1<=len(sentences)<=3 or any(not isinstance(x,str) for x in sentences): raise ValueError()
            text='\n\n'.join(sentences)
            if any(t in text.casefold() for t in forbidden): raise ValueError()
            return text
        except (ValueError,KeyError,TypeError):
            messages.append({'role':'user','content':'Previous output violated the evidence rules. Return only 2 short factual sentences about actual numeric values and lack of uncertainty information. No biological interpretation.'})
    raise ValueError('Model explanation failed evidence checks; computed values unchanged.')

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
    kind, module = plan.get('kind'), plan.get('module')
    if kind not in frames or kind not in metadata: raise ValueError('Unknown layer')
    if module not in ['Profile','Pairwise','Group comparison','PCA / UMAP','Data tables']: raise ValueError('Unknown module')
    feature = plan.get('feature')
    if feature is not None and feature not in frames[kind].index: raise ValueError('Unknown feature')
    available = metadata[kind].group.unique().tolist()
    groups = plan.get('groups') or available
    if not isinstance(groups,list) or any(not isinstance(g,str) or g not in available for g in groups): raise ValueError('Unknown groups')
    return {'kind':kind,'module':module,'groups':list(dict.fromkeys(groups)),**({'feature':feature} if feature else {})}

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
