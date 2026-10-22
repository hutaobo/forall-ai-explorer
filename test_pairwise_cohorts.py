"""Real-data regression of the full and restricted cohort, plus UI contract."""
import json,os
from pathlib import Path
from unittest.mock import patch
import numpy as np
from scipy import stats
from core import load,samples,pairwise_table,result_context
from adjusted import group_adjusted_cohorts

frames,metadata,sources=load(os.environ['FORALL_DATA_ROOT'])
selected=samples(metadata['protein'],metadata['protein'].group.unique(),False)
clean=pairwise_table(frames,metadata,'protein','BCL2',selected,'drug','Venetoclax').dropna(subset=['x','y'])
assert len(clean)==43 and clean.index.is_unique and clean.second_sample.is_unique
assert not clean.replicate.any() and np.isfinite(clean[['x','y']]).all().all()
cohorts=group_adjusted_cohorts(clean)
full=cohorts['all']['statistics'];restricted=cohorts['bcp_t']['statistics']
assert full['group_counts']=={'BCP-ALL':25,'T-ALL':16,'B-ALL':2}
assert restricted['n']==41 and restricted['excluded_from_full_n']==2
assert np.isclose(stats.spearmanr(clean.x,clean.y).statistic,0.7383883582013337)
assert np.isclose(full['partial_Spearman']['coefficient'],0.38193791142489736)
assert np.isclose(restricted['partial_Spearman']['coefficient'],0.38044992885852663)
assert full['partial_Spearman']['permutation_p']==0.0152
assert restricted['partial_Spearman']['permutation_p']==0.017
for entry in cohorts.values():
    residuals=entry['residuals']
    assert np.allclose(residuals.groupby('group')[['x_rank_group_residual','y_rank_group_residual']].mean(),0,atol=1e-12)
    assert all(block.x.nunique()>1 and block.y.nunique()>1 for g,block in residuals.groupby('group') if g in ('BCP-ALL','T-ALL'))
context=result_context(clean,{'module':'Pairwise'},sources,{'group_adjusted':full,'BCP_T_sensitivity':restricted})
assert len(context['records'])==43 and context['statistics']['BCP_T_sensitivity']['n']==41
out=Path(os.environ['FORALL_REPORT_ROOT']);out.mkdir(parents=True,exist_ok=True)
variation={g:{'n':len(block),'unique_x':block.x.nunique(),'unique_y':block.y.nunique()} for g,block in clean.groupby('group')}
(out/'cohort-statistics.json').write_text(json.dumps({'all43':full,'BCP_T41':restricted,'within_group_variation':variation},indent=2))
clean.to_csv(out/'complete-43-pairs.csv',index=False)
from streamlit.testing.v1 import AppTest
os.environ['FORALL_SEEKYOO_AUTH']='0'
app=AppTest.from_file('app.py',default_timeout=60).run()
app.selectbox(key='module').set_value('Pairwise').run()
app.selectbox(key='second_kind').set_value('drug').run()
app.selectbox(key='second_feature').set_value('Venetoclax').run()
assert not app.exception
assert any(x.label=='Adjusted scatter cohort' for x in app.selectbox)
table=next(x.value for x in app.dataframe if 'partial_Spearman' in x.value.columns)
assert table.n.tolist()==[43,41]
assert app.get('plotly_chart') and len(app.get('plotly_chart'))==2
next(x for x in app.selectbox if x.label=='Adjusted scatter cohort').set_value('bcp_t').run()
assert not app.exception and any('n=41' in x.value for x in app.caption)
print('Passed: 43/41 independent reference values, unique finite primary pairs, residual centering, full AI context and both interactive chart/cohort controls.')
