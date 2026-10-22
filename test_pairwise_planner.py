"""Synthetic schema-repair tests; no provider calls or real data."""
import pandas as pd
from unittest.mock import patch
from core import validate_controls, plan_controls, explain_result
import json

frames={'protein':pd.DataFrame(index=['BCL2']),'drug':pd.DataFrame(index=['Venetoclax'])}
metadata={k:pd.DataFrame({'group':['BCP-ALL','T-ALL','B-ALL'],'replicate':[False]*3}) for k in frames}
plan={'kind':'蛋白','module':'pairwise','feature':'bcl2','groups':[],
      'second_kind':'药敏','second_feature':'venetoclax','method':'spearman'}
fixed=validate_controls(plan,frames,metadata)
assert fixed['kind']=='protein' and fixed['feature']=='BCL2'
assert fixed['second_kind']=='drug' and fixed['second_feature']=='Venetoclax'
assert fixed['groups']==['BCP-ALL','T-ALL','B-ALL'] and fixed['correlation_method']=='Spearman'
for broken in [{**plan,'kind':'RNA'},{**plan,'second_feature':None},{**plan,'method':'invented'}]:
    try:validate_controls(broken,frames,metadata)
    except ValueError:pass
    else:raise AssertionError('Unsupported or incomplete settings were accepted')
with patch('core.chat',side_effect=[{'kind':'RNA'},plan]) as model:
    assert plan_controls('全部细胞系的BCL2与Venetoclax相关性',frames,metadata)==fixed
    assert model.call_count==2
subset={**plan,'groups':['BCP-ALL','T-ALL']}
def fake_pairs(*args):
    n=43 if len(args[4])==3 else 41
    return pd.DataFrame({'x':range(n),'y':range(n)})
with patch('core.chat',side_effect=[subset,plan]) as model,patch('core.pairwise_table',side_effect=fake_pairs):
    assert plan_controls('全部43个细胞系并校正BCP-ALL/T-ALL',frames,metadata)=={**fixed,'reps':False}
    assert model.call_count==2
with patch('core.chat',return_value={'kind':'RNA'}) as model:
    try:plan_controls('synthetic request',frames,metadata)
    except ValueError:pass
    else:raise AssertionError('Invalid model response did not fail closed')
    assert model.call_count==2
print('Passed: bilingual exact-name normalization, required pair/method, unsupported RNA rejection and one repair attempt.')
sentences=['已提供43个完整配对。','校正结果与全体相关性分开报告。','组内检验不显著不代表没有关联。','这些是细胞系分组，不是分子亚型。']
with patch('core.chat',return_value=json.dumps({'sentences':sentences},ensure_ascii=False)):
    assert explain_result({'statistics':{'paired_n':43}},'synthetic formatting check').count('\n\n')==3
print('Passed: requested multi-cohort explanations are not truncated to three sentences or rejected solely for calibrated Chinese test terminology.')
fixture={'method':'Spearman','paired_n':43,'coefficient':.738,'nominal_p':1.59e-8,'group_adjusted':{'n':43,'group_counts':{'BCP-ALL':25,'T-ALL':16,'B-ALL':2},'partial_Spearman':{'coefficient':.382,'permutation_p':.0152,'bootstrap_95_ci':[.055,.638]}}}
with patch('core.chat',side_effect=[json.dumps({'sentences':['EBV has 2 excluded cases.']}),json.dumps({'sentences':['Adjustment weakens the association; residual association does not establish causality.']})]) as model:
    answer=explain_result({'statistics':fixture},'synthetic cohort guard')
    assert model.call_count==2 and 'B-ALL: 2' in answer and 'EBV' not in answer
    assert 'n=43' in answer and 'p=0.0152' in answer
print('Passed: computed cohort facts rendered verbatim; invented model counts rejected and repaired.')
