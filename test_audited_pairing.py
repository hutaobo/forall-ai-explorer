"""Publisher-inventory-supported, dataset-specific pairing and stratified analysis."""
import os,json
from pathlib import Path
import pandas as pd
import plotly.express as px
from scipy import stats
from core import load,samples,pairwise_table,within_group_correlations,FORALL_PROTEIN_LABELS
f,m,s=load(os.environ['FORALL_DATA_ROOT'])
cells=f['cells'].copy(); cells.columns=[' '.join(str(c).split()) for c in cells.columns]
evidence=[]
for protein,drug in FORALL_PROTEIN_LABELS.items():
    row=cells[cells['Cell Line Name']==drug]
    assert len(row)==1 and not (cells['Cell Line Name']==protein).any()
    assert protein in f['protein'].columns and drug in f['drug'].columns
    assert m['protein'].loc[protein,'group']==m['drug'].loc[drug,'group']==row.iloc[0]['Lineage']
    evidence.append({'protein_label':protein,'drug_label':drug,'canonical_inventory_label':drug,
                     'TMT10_channel':row.iloc[0]['TMT10 channel'],'lineage':row.iloc[0]['Lineage'],
                     'rna_sample':row.iloc[0]['RNA-seq sample ID'],
                     'scope':'FORALL tables only, not a universal biological synonym'})
selected=samples(m['drug'],m['drug'].group.unique(),False)
exact=pairwise_table(f,m,'drug','Venetoclax',selected,'protein','BCL2',False)
paired=pairwise_table(f,m,'drug','Venetoclax',selected,'protein','BCL2',True)
assert len(exact.dropna(subset=['x','y']))==38
assert len(paired.dropna(subset=['x','y']))==43 and paired.second_sample.is_unique
assert paired.label_harmonized.sum()==5
# Plot with protein on x, drug on y, retain both published labels.
paired=paired.rename(columns={'x':'sDSS','y':'BCL2','sample':'drug_sample','second_sample':'protein_sample'})
paired['x']=paired.BCL2; paired['y']=paired.sDSS
out={'pairing_n':43,'added_pairs':evidence,'methods':{}}
for method,fn in [('Spearman',stats.spearmanr),('Pearson',stats.pearsonr)]:
    stat=fn(paired.x,paired.y)
    out['methods'][method]={'overall':{'n':43,'coefficient':float(stat.statistic),'nominal_p':float(stat.pvalue)},
                            'within_groups':json.loads(within_group_correlations(paired,method).to_json(orient='records'))}
root=Path(os.environ.get('FORALL_REPORT_ROOT','.')); root.mkdir(parents=True,exist_ok=True)
(root/'audited-pairing-statistics.json').write_text(json.dumps(out,indent=2))
paired.drop(columns=['x','y']).to_csv(root/'audited-43-pairs.csv',index=False)
colors={'BCP-ALL':'#0072B2','T-ALL':'#E69F00','B-ALL':'#CC79A7'}
for group in ['All','BCP-ALL','T-ALL']:
    block=paired if group=='All' else paired[paired.group==group]
    stat=stats.spearmanr(block.BCL2,block.sDSS)
    fig=px.scatter(block,x='BCL2',y='sDSS',color='group',hover_name='drug_sample',
                   hover_data=['protein_sample','label_harmonized'],color_discrete_map=colors,
                   labels={'BCL2':'BCL2 · log2 relative protein abundance','sDSS':'Venetoclax · sDSS'})
    fig.update_layout(template='plotly_white',font_family='Arial',height=620,
                      title=f'{group} · n={len(block)} · Spearman ρ={stat.statistic:.3f} · nominal p={stat.pvalue:.3g}')
    fig.write_html(root/f'audited-{group}-interactive.html',include_plotlyjs=True)
print(json.dumps(out))
