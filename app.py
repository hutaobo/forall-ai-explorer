import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from scipy import stats
from core import load, samples, compare, embedding, chat, explain_result, UNITS, DOI

st.set_page_config(page_title='FORALL AI Explorer',page_icon='🧬',layout='wide')
st.title('FORALL · AI Explorer')
st.caption('Interactive published-data research · SciLifeLab OpenLLM · Independent implementation, not endorsed by the original authors')
@st.cache_data(show_spinner='Loading verified published workbooks…')
def dataset(root): return load(root)
root = os.environ.get('FORALL_DATA_ROOT','/data/forall')
try: frames, metadata, sources = dataset(root)
except Exception as e: st.error(str(e)); st.stop()
if 'pending_controls' in st.session_state:
    for key,value in st.session_state.pop('pending_controls').items(): st.session_state[key]=value

with st.sidebar:
    st.header('Analysis controls')
    kind = st.selectbox('Data layer',['protein','drug'],key='kind')
    module = st.selectbox('Analysis',['Profile','Pairwise','Group comparison','PCA / UMAP','Data tables'],key='module')
    group_names = list(dict.fromkeys(metadata[kind].group))
    groups = st.multiselect('Cell-line groups',group_names,default=group_names,key='groups')
    reps = st.checkbox('Include published replicate columns',value=False,key='reps')
    st.caption('Replicates retain their original columns. Group tests use only primary columns.')
    st.markdown(f'[Original paper]({DOI}) · [Original FORALL](https://proteomics.se/forall/)')
    st.caption('Public tables: CC BY 4.0. Protein is log2 relative abundance; sDSS is not IC50.')

meta = metadata[kind]; selected=samples(meta,groups,reps)
if not selected: st.info('Select at least one group.'); st.stop()
frame = frames[kind].loc[:,selected]
result = None; params={'layer':kind,'module':module,'groups':groups,'replicates':reps}
left, right=st.columns([3,1])
with left:
    st.subheader(module)
    if module in ['Profile','Pairwise']:
        term = st.selectbox('Protein / drug',frame.index.tolist(),index=frame.index.get_loc('BCL2' if kind=='protein' else 'Venetoclax'),key='feature')
        params['feature']=term
    if module=='Profile':
        chart=st.selectbox('Chart',['Points','Bars','Boxplot'])
        sort=st.checkbox('Sort by value',value=True)
        result=meta.loc[selected].copy();result['value']=frame.loc[term];result['missing']=result.value.isna()
        if sort:result=result.sort_values('value',ascending=False,na_position='last')
        clean=result.dropna(subset=['value'])
        if chart=='Boxplot':fig=px.box(clean,x='group',y='value',color='group',points='all',hover_data=['sample','excel_column'])
        elif chart=='Bars':fig=px.bar(clean,x='sample',y='value',color='group',hover_data=['excel_column'])
        else:fig=px.scatter(clean,x='sample',y='value',color='group',symbol='replicate',hover_data=['excel_column'])
        fig.update_layout(yaxis_title=UNITS[kind],xaxis_title='Published sample / group')
        st.caption(f'{len(clean)} plotted values; {result.missing.sum()} missing values retained in export. No averaging or imputation.')
    elif module=='Pairwise':
        second_kind=st.selectbox('Second data layer',['protein','drug'],index=0 if kind=='protein' else 1)
        other=frames[second_kind]
        second=st.selectbox('Second protein / drug',other.index.tolist(),index=other.index.get_loc('NR3C1' if second_kind=='protein' else 'Dexamethasone'))
        method=st.selectbox('Correlation',['Spearman','Pearson'])
        result=meta.loc[selected].copy();result['x']=frame.loc[term];result['y']=other.loc[second].reindex(selected)
        clean=result.dropna(subset=['x','y'])
        if len(clean)<3 or clean.x.nunique()<2 or clean.y.nunique()<2:st.warning('At least three paired nonconstant measurements required.');st.stop()
        stat=(stats.spearmanr if method=='Spearman' else stats.pearsonr)(clean.x,clean.y)
        st.metric(f'{method} r',f'{stat.statistic:.3f}');st.caption(f'n={len(clean)} paired sample columns; nominal p={stat.pvalue:.4g}. Exploratory, not causality; replicate inclusion affects independence.')
        fig=px.scatter(clean,x='x',y='y',color='group',symbol='replicate',hover_data=['sample','excel_column'])
        fig.update_layout(xaxis_title=f'{term} · {UNITS[kind]}',yaxis_title=f'{second} · {UNITS[second_kind]}')
        st.caption(f'Cross-layer pairing uses exact published sample names only; {len(result)-len(clean)} unmatched/missing pairs excluded. Aliases are not silently merged.')
        params.update(second=second,second_layer=second_kind,method=method)
    elif module=='Group comparison':
        if len(groups)<2:st.warning('Select two groups in the sidebar.');st.stop()
        a=st.selectbox('Group A',groups);b=st.selectbox('Group B',groups,index=1)
        test=st.selectbox('Test',['Welch t-test','Mann–Whitney U'])
        if a==b:st.warning('Choose different groups.');st.stop()
        primary=samples(meta,groups,False)
        result=compare(frames[kind].loc[:,primary],meta.loc[primary],a,b,test)
        if result.empty:st.warning('Each group needs at least three primary measurements per feature.');st.stop()
        fig=px.scatter(result,x='mean_difference',y=-np.log10(result.q_BH.clip(lower=1e-300)),hover_name='feature')
        fig.update_layout(xaxis_title=f'Mean {a} − {b} · {UNITS[kind]}',yaxis_title='−log10(BH q)')
        st.caption('Exploratory across cell lines; replicates excluded. BH adjustment over all tested features. This is not DEqMS/edgeR reproduction.')
        st.caption(f'{result.p.isna().sum()} untestable/constant features have undefined p/q values, not evidence of significance.')
        params.update(group_a=a,group_b=b,test=test)
    elif module=='PCA / UMAP':
        method=st.selectbox('Embedding',['PCA','UMAP'])
        features=st.number_input('Most-variable complete features',min_value=2,max_value=len(frame),value=min(1000,len(frame)))
        neighbors=st.slider('Neighbors',2,50,15);distance=st.slider('Minimum distance',0.0,0.99,0.1);seed=st.number_input('Random seed',0,999999,42)
        signature=(kind,tuple(selected),method,features,neighbors,distance,seed)
        if st.session_state.get('embedding_signature') != signature: st.session_state.pop('embedding_result',None)
        if st.button('Compute embedding'):
            with st.spinner('Computing…'):
                result,details=embedding(frame,method,int(features),neighbors,distance,int(seed));params.update(details)
                result=result.join(meta[['group','replicate']]);result['sample']=result.index
                fig=px.scatter(result,x='axis1',y='axis2',color='group',symbol='replicate',hover_name='sample')
                st.session_state['embedding_result']=(result,params,fig)
                st.session_state['embedding_signature']=signature
        if 'embedding_result' not in st.session_state:st.info('Set parameters and compute.');st.stop()
        result,params,fig=st.session_state['embedding_result'];st.json(params)
    else:
        table=st.selectbox('Table',['cells','drugs']);result=frames[table]
        st.dataframe(result,hide_index=True);fig=None
    if fig is not None:
        fig.update_layout(template='plotly_white',font_family='Arial',height=530,legend_title_text='Group')
        st.plotly_chart(fig,use_container_width=True,config={'toImageButtonOptions':{'format':'svg'}})
        st.download_button('Download interactive HTML',fig.to_html(include_plotlyjs=True),'forall-plot.html','text/html')
    st.dataframe(result,width='stretch')
    st.download_button('Download data CSV',result.to_csv(index=True,index_label='row_id'),'forall-data.csv','text/csv')
    st.download_button('Download parameters & provenance',json.dumps({'parameters':params,'sources':sources},indent=2),'forall-provenance.json')
with right:
    st.subheader('OpenLLM assistant')
    session_key=st.text_input('OpenLLM API key (optional if configured by host)',type='password',help='Used only in this browser session; not saved to disk.')
    st.caption('Your question and the displayed public-data summary are sent to SciLifeLab OpenLLM only when you click a button. Keys stay server-side.')
    prompt=st.text_area('Ask / 操作或解释',placeholder='绘制 BCL2 蛋白表达图，选择 T-ALL；或解释当前图中的趋势。')
    if st.button('AI: suggest controls'):
        try:
            suggestion=chat([{'role':'system','content':'Return JSON with only kind (protein/drug), module (Profile/Pairwise/Group comparison/PCA / UMAP/Data tables), feature, groups (array), explanation. Suggest controls for this public FORALL explorer. Do not invent results or executable code. Available groups: '+json.dumps(group_names)}, {'role':'user','content':prompt}],True,session_key)
            st.session_state['suggestion']=json.loads(suggestion)
        except Exception:st.error('OpenLLM did not return valid controls. No analysis changed.')
    if 'suggestion' in st.session_state:
        plan=st.session_state['suggestion'];st.json(plan)
        if st.button('Apply suggested controls'):
            newkind=plan.get('kind',kind);newmodule=plan.get('module',module);newfeature=plan.get('feature')
            if newkind not in frames or newkind not in metadata:st.error('Unknown layer.');st.stop()
            if newmodule not in ['Profile','Pairwise','Group comparison','PCA / UMAP','Data tables']:st.error('Unknown module.');st.stop()
            validgroups=plan.get('groups') or metadata[newkind].group.unique().tolist()
            if any(g not in metadata[newkind].group.values for g in validgroups):st.error('Unknown groups.');st.stop()
            if newfeature and newfeature not in frames[newkind].index:st.error('Feature not found in original table.');st.stop()
            st.session_state['pending_controls']={'kind':newkind,'module':newmodule,'groups':validgroups,**({'feature':newfeature} if newfeature else {})}
            st.rerun()
    if st.button('AI: explain current result'):
        try:
            summary={'parameters':params,'rows':len(result),'preview':json.loads(result.head(12).to_json(orient='records')),'source':sources[kind], 'scope':'Preview only, not all rows. No clinical advice or causal conclusions.'}
            answer=explain_result(summary,prompt,api_key=session_key)
            st.write(answer)
        except Exception:st.error('OpenLLM unavailable; your computed result remains available.')
    with st.expander('Data provenance'):st.json(sources)
