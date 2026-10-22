"""Full-data payload and exact-name BCL2/Venetoclax scientific acceptance."""
import json
import os
from pathlib import Path
from scipy import stats
from core import load, samples, result_context

frames, metadata, sources = load(os.environ['FORALL_DATA_ROOT'])
drug = frames['drug'].loc['Venetoclax']
profile = metadata['drug'].copy()
profile['value'] = drug
context = result_context(profile, {'module':'Profile'}, sources)
assert len(context['records']) == 43
assert {r['sample'] for r in context['records']} == set(drug.index)
assert 'preview' not in context
primary = samples(metadata['protein'], metadata['protein'].group.unique(), False)
paired = metadata['protein'].loc[primary].copy()
paired['x'] = frames['protein'].loc['BCL2',primary]
paired['y'] = drug.reindex(primary)
clean = paired.dropna(subset=['x','y'])
assert clean.index.is_unique
out = {'drug_total':len(drug),'protein_primary':len(primary),'paired_n':len(clean),
       'drug_unmatched_names':sorted(set(drug.index)-set(primary)),
       'protein_unmatched_names':sorted(set(primary)-set(drug.index)),
       'pairing':'exact original names; primary columns only; no alias merging',
       'by_group':clean.groupby('group').size().to_dict()}
for method, fn in [('Spearman',stats.spearmanr),('Pearson',stats.pearsonr)]:
    s = fn(clean.x,clean.y)
    out[method] = {'coefficient':float(s.statistic),'nominal_p':float(s.pvalue)}
payload = result_context(paired,{'module':'Pairwise'},sources,out)
assert len(payload['records'])==len(primary)
assert sum(r['x'] is not None and r['y'] is not None for r in payload['records'])==len(clean)
target=Path(os.environ.get('FORALL_REPORT_ROOT','.'))
target.mkdir(parents=True,exist_ok=True)
clean.to_csv(target/'bcl2-venetoclax-pairs.csv',index=False)
(target/'bcl2-venetoclax-statistics.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out))
