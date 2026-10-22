"""Focused real-data acceptance for reviewer trial coverage and AI settings."""
import os
import numpy as np
from scipy import stats
from core import load,samples,correlation_scan,pairwise_table,validate_controls
f,m,s=load(os.environ['FORALL_DATA_ROOT'])
selected=samples(m['protein'],m['protein'].group.unique(),False)
scan=correlation_scan(f,m,'protein','BCL2',selected,'drug')
assert len(scan)==528
v=scan.set_index('feature').loc['Venetoclax']
assert v['n']==43 and np.isclose(v.coefficient,0.7383883582013337)
for drug in ['Venetoclax','Dexamethasone','Idasanutlin']:
    pairs=pairwise_table(f,m,'protein','BCL2',selected,'drug',drug).dropna(subset=['x','y'])
    ref=stats.spearmanr(pairs.x,pairs.y)
    row=scan.set_index('feature').loc[drug]
    assert np.isclose(row.coefficient,ref.statistic) and np.isclose(row.nominal_p,ref.pvalue)
assert scan.q_BH.dropna().between(0,1).all()
plan=validate_controls({'kind':'protein','module':'Pairwise','feature':'BCL2','groups':['BCP-ALL','T-ALL','B-ALL'],
                        'second_kind':'drug','second_feature':'Venetoclax','method':'Spearman'},f,m)
assert plan['second_feature']=='Venetoclax' and plan['correlation_method']=='Spearman'
from streamlit.testing.v1 import AppTest
os.environ['FORALL_SEEKYOO_AUTH']='0' # Test process only, no listening server.
at=AppTest.from_file('app.py',default_timeout=60).run()
assert not at.exception and at.text_area[0].label=='Ask AI'
at.selectbox(key='module').set_value('Pairwise').run()
at.selectbox(key='second_kind').set_value('drug').run()
at.selectbox(key='second_feature').set_value('Venetoclax').run()
assert not at.exception and any('0.382' in str(x.value) for x in at.metric)
at.selectbox(key='module').set_value('Correlation scan').run()
next(b for b in at.button if b.label=='Compute correlation scan').click().run()
assert not at.exception and any(len(d.value)==528 for d in at.dataframe)
at.selectbox(key='kind').set_value('drug').run()
at.selectbox(key='module').set_value('Cell-line profile').run()
assert not at.exception and any(len(d.value)==528 for d in at.dataframe)
print('Passed: 528-drug scan/independent scipy checks, AI pair settings, English label, adjusted UI, scan UI, cell-line UI.')
