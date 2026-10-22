"""Independent residualization checks and real-data group-adjusted analysis."""
import os,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import plotly.express as px
from adjusted import group_adjusted
from core import load,samples,pairwise_table

rng=np.random.default_rng(71)
fixture=pd.DataFrame({'x':rng.normal(size=40),'group':np.repeat(['A','B'],20)})
fixture['y']=2*fixture.x+fixture.group.map({'A':-10,'B':10})+rng.normal(scale=.1,size=40)
res,table=group_adjusted(fixture,permutations=199,bootstraps=100)
assert abs(res['linear_model']['slope']-2)<.1
assert np.allclose(table.groupby('group')[['x_group_residual','y_group_residual']].mean(),0,atol=1e-12)
ranked=fixture[['x','y']].rank(method='average')
ranked-=ranked.groupby(fixture.group).transform('mean')
assert np.isclose(stats.pearsonr(ranked.x,ranked.y).statistic,res['partial_Spearman']['coefficient'])
shifted=fixture.copy(); shifted.x+=shifted.group.map({'A':100,'B':-200}); shifted.y+=shifted.group.map({'A':-7,'B':31})
other,_=group_adjusted(shifted,199,100)
assert np.isclose(other['partial_Pearson']['coefficient'],res['partial_Pearson']['coefficient'])
assert np.isclose(other['linear_model']['slope'],res['linear_model']['slope'])
assert res['partial_Spearman']['permutation_p']==group_adjusted(fixture,199,100)[0]['partial_Spearman']['permutation_p']

frames,meta,sources=load(os.environ['FORALL_DATA_ROOT'])
selected=samples(meta['protein'],meta['protein'].group.unique(),False)
out={}; root=Path(os.environ.get('FORALL_REPORT_ROOT','.'));root.mkdir(parents=True,exist_ok=True)
for label,harmonize,groups in [('all43',True,['BCP-ALL','T-ALL','B-ALL']),('BCP_T41',True,['BCP-ALL','T-ALL']),('strict38',False,['BCP-ALL','T-ALL','B-ALL'])]:
    data=pairwise_table(frames,meta,'protein','BCL2',selected,'drug','Venetoclax',harmonize)
    data=data[data.group.isin(groups)].dropna(subset=['x','y'])
    summary,residuals=group_adjusted(data)
    assert summary['n']=={'all43':43,'BCP_T41':41,'strict38':38}[label]
    out[label]=summary
    residuals.to_csv(root/f'adjusted-{label}-data.csv',index=False)
    for tag,x,y in [('ranks','x_rank_group_residual','y_rank_group_residual'),('values','x_group_residual','y_group_residual')]:
        fig=px.scatter(residuals,x=x,y=y,color='group',hover_name='sample',hover_data=['second_sample'],
            labels={x:'BCL2 · '+('rank residual' if tag=='ranks' else 'group-demeaned log2 abundance'),
                    y:'Venetoclax · '+('rank residual' if tag=='ranks' else 'group-demeaned sDSS')},
            color_discrete_map={'BCP-ALL':'#0072B2','T-ALL':'#E69F00','B-ALL':'#CC79A7'})
        fig.add_hline(y=0,line_dash='dot',line_color='gray');fig.add_vline(x=0,line_dash='dot',line_color='gray')
        fig.update_layout(template='plotly_white',font_family='Arial',height=620,title=f'Group-adjusted association · {label}')
        fig.write_html(root/f'adjusted-{label}-{tag}.html',include_plotlyjs=True)
(root/'group-adjusted-statistics.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out))
