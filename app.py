import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
import streamlit.components.v1 as components
from scipy import stats
from core import load, samples, compare, embedding, chat, explain_result, result_context, pairwise_table, within_group_correlations, correlation_scan, MODULES, validate_controls, GLOSSARY, UNITS, DOI
from seekyoo_access import gate, hosted_request
from adjusted import group_adjusted

st.set_page_config(page_title='FORALL AI Explorer',page_icon='🧬',layout='wide')
identity=gate(st)
st.title('FORALL · AI Explorer')
st.caption('Interactive published-data research · SciLifeLab OpenLLM · Independent implementation, not endorsed by the original authors')
@st.cache_data(show_spinner='Loading verified published workbooks…')
def dataset(root): return load(root)
root = os.environ.get('FORALL_DATA_ROOT','/data/forall')
try: frames, metadata, sources = dataset(root)
except Exception as e: st.error(str(e)); st.stop()
if 'pending_controls' in st.session_state:
    for key,value in st.session_state.pop('pending_controls').items(): st.session_state[key]=value
run_requested=st.session_state.pop('run_requested',False)

with st.sidebar:
    st.header('Analysis controls')
    kind = st.selectbox('Data layer',['protein','drug'],key='kind')
    module = st.selectbox('Analysis',MODULES,key='module')
    group_names = list(dict.fromkeys(metadata[kind].group))
    if 'groups' not in st.session_state: st.session_state['groups']=group_names
    else: st.session_state['groups']=[g for g in st.session_state['groups'] if g in group_names]
    groups = st.multiselect('Cell-line groups',group_names,key='groups')
    reps = st.checkbox('Include published replicate columns',value=False,key='reps')
    st.caption('Replicates retain their original columns. Group tests use only primary columns.')
    st.markdown(f'[Original paper]({DOI}) · [Original FORALL](https://proteomics.se/forall/)')
    st.caption('Public tables: CC BY 4.0. Protein is log2 relative abundance; sDSS is not IC50.')

meta = metadata[kind]; selected=samples(meta,groups,reps)
if not selected: st.info('Select at least one group.'); st.stop()
frame = frames[kind].loc[:,selected]
result = None; computed_statistics={}; params={'layer':kind,'module':module,'groups':groups,'replicates':reps}
left, right=st.columns([3,1])
with left:
    st.subheader(module)
    if module in ['Profile','Pairwise','Correlation scan']:
        if st.session_state.get('feature') not in frame.index: st.session_state['feature']='BCL2' if kind=='protein' else 'Venetoclax'
        term = st.selectbox('Protein / drug',frame.index.tolist(),key='feature')
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
        if 'second_kind' not in st.session_state:st.session_state['second_kind']=kind
        second_kind=st.selectbox('Second data layer',['protein','drug'],key='second_kind')
        other=frames[second_kind]
        if st.session_state.get('second_feature') not in other.index:st.session_state['second_feature']='NR3C1' if second_kind=='protein' else 'Dexamethasone'
        second=st.selectbox('Second protein / drug',other.index.tolist(),key='second_feature')
        method=st.selectbox('Correlation',['Spearman','Pearson'],key='correlation_method')
        harmonize=st.checkbox('Use audited FORALL label harmonization',value=True)
        result=pairwise_table(frames,metadata,kind,term,selected,second_kind,second,harmonize)
        clean=result.dropna(subset=['x','y'])
        if len(clean)<3 or clean.x.nunique()<2 or clean.y.nunique()<2:st.warning('At least three paired nonconstant measurements required.');st.stop()
        stat=(stats.spearmanr if method=='Spearman' else stats.pearsonr)(clean.x,clean.y)
        computed_statistics={'method':method,'coefficient':float(stat.statistic),'nominal_p':float(stat.pvalue),'paired_n':len(clean),'selected_n':len(result),'excluded_n':len(result)-len(clean),'multiple_testing_adjustment':'none; exploratory single selected pair'}
        grouped=within_group_correlations(clean,method)
        computed_statistics['within_group']=json.loads(grouped.to_json(orient='records'))
        computed_statistics['within_group_adjustment']='BH across eligible displayed groups, separately from overall test; groups with n<3 omitted'
        st.metric(f'{method} r',f'{stat.statistic:.3f}');st.caption(f'n={len(clean)} paired sample columns; nominal p={stat.pvalue:.4g}. Exploratory, not causality; replicate inclusion affects independence.')
        fig=px.scatter(clean,x='x',y='y',color='group',symbol='replicate',hover_data=['sample','second_sample','label_harmonized','excel_column'])
        fig.update_layout(xaxis_title=f'{term} · {UNITS[kind]}',yaxis_title=f'{second} · {UNITS[second_kind]}')
        st.caption(f'{len(result)-len(clean)} unmatched/missing pairs excluded; {int(clean.label_harmonized.sum())} pairs use audited dataset-specific labels. Original names retained in hover and CSV. No general h-suffix removal.')
        st.caption('Correspondence: published Supplementary Data 1 sample inventory vs Data 3 protein / Data 13 drug headers; COG-LL-319/332/394/402 → corresponding h labels, TMD-5 → TMD5. These are dataset-specific labels, not universal aliases of different culture derivatives.')
        st.subheader('Within-group correlations')
        st.dataframe(grouped,hide_index=True)
        st.caption('BH correction across eligible displayed groups for the selected method. Groups with fewer than 3 pairs omitted; pooled association may reflect group differences.')
        st.subheader('Group-adjusted association')
        st.caption('Controls categorical cell-line group. Primary analysis: global average ranks, then remove group means; permute within groups. Pearson/linear regression is a sensitivity analysis. Primary columns only.')
        try:
            adjusted,residuals=group_adjusted(clean.loc[~clean.replicate])
            computed_statistics['group_adjusted']=adjusted
            rank=adjusted['partial_Spearman']; linear=adjusted['linear_model']
            st.metric('Group-adjusted partial Spearman',f"{rank['coefficient']:.3f}")
            st.caption(f"n={adjusted['n']}; stratified permutation p={rank['permutation_p']:.4g} (9,999 permutations); stratified bootstrap 95% CI [{rank['bootstrap_95_ci'][0]:.3f}, {rank['bootstrap_95_ci'][1]:.3f}] (2,000 resamples); seed 20261022.")
            st.write(f"Adjusted Pearson r={adjusted['partial_Pearson']['coefficient']:.3f}; nominal t p={adjusted['partial_Pearson']['nominal_t_p']:.4g}. Common linear slope={linear['slope']:.3f}, HC3 95% CI [{linear['slope_HC3_95_ci'][0]:.3f}, {linear['slope_HC3_95_ci'][1]:.3f}], HC3 approximate p={linear['slope_HC3_t_approx_p']:.4g}.")
            st.caption(f"Outcome R²: group only={linear['group_only_R2']:.3f}; group + first variable={linear['full_R2']:.3f}; incremental R²={linear['incremental_R2']:.3f}. A common within-group slope is assumed; this is not proof of causal independence or adjustment for all confounders.")
            st.info('Interpretation: compare this adjusted coefficient with the pooled coefficient above. Residual association after removing group means is evidence beyond group mean differences under this model, not proof of causal independence. Separate group tests have smaller sample sizes and may be less precise. Select only BCP-ALL/T-ALL or turn off label harmonization to check sensitivity.')
            view=st.selectbox('Adjusted scatter',['Group-adjusted global ranks','Group-demeaned original values'])
            ax,ay=('x_rank_group_residual','y_rank_group_residual') if view=='Group-adjusted global ranks' else ('x_group_residual','y_group_residual')
            residual_fig=px.scatter(residuals,x=ax,y=ay,color='group',hover_name='sample',hover_data=['second_sample'],
                labels={ax:f'{term} · '+('rank residual' if 'rank' in ax else 'group-demeaned '+UNITS[kind]),
                        ay:f'{second} · '+('rank residual' if 'rank' in ay else 'group-demeaned '+UNITS[second_kind])})
            residual_fig.add_hline(y=0,line_dash='dot',line_color='gray');residual_fig.add_vline(x=0,line_dash='dot',line_color='gray')
            residual_fig.update_layout(template='plotly_white',font_family='Arial',height=500)
            st.plotly_chart(residual_fig,use_container_width=True)
            st.download_button('Download adjusted interactive HTML',residual_fig.to_html(include_plotlyjs=True),'forall-adjusted.html','text/html')
            st.download_button('Download adjusted data CSV',residuals.to_csv(index=False),'forall-adjusted.csv','text/csv')
            st.download_button('Download adjustment statistics',json.dumps(adjusted,indent=2),'forall-adjustment.json')
        except ValueError as error:
            st.info(str(error))
        params.update(second=second,second_layer=second_kind,method=method,label_harmonization=harmonize,statistics=computed_statistics)
    elif module=='Correlation scan':
        target_kind=st.selectbox('Scan against data layer',['protein','drug'],index=1)
        method=st.selectbox('Scan correlation',['Spearman','Pearson'])
        signature=(kind,term,tuple(selected),target_kind,method)
        if st.button('Compute correlation scan') or run_requested:
            with st.spinner('Computing pairwise-complete correlations across the full target layer…'):
                st.session_state['scan_result']=(signature,correlation_scan(frames,metadata,kind,term,selected,target_kind,method))
        stored=st.session_state.get('scan_result');result=stored[1] if stored and stored[0]==signature else pd.DataFrame()
        fig=None
        if result.empty:st.info('Choose the feature and target layer, then compute the scan.')
        else:
            plot=result.dropna(subset=['coefficient']).copy();plot['rank']=np.arange(1,len(plot)+1)
            fig=px.scatter(plot,x='rank',y='coefficient',hover_name='feature',hover_data=['n','nominal_p','q_BH'],render_mode='webgl')
            fig.update_layout(xaxis_title='Correlation rank (descending)',yaxis_title=f'{method} correlation')
            st.caption(f'{len(result)} target features; {len(plot)} testable. Pairwise-complete data; BH across all testable target features. Spearman p uses the asymptotic t approximation. Audited FORALL labels enabled; replicates follow the sidebar selection.')
        params.update(second_layer=target_kind,method=method)
    elif module=='Cell-line profile':
        if st.session_state.get('cell_line') not in selected:st.session_state['cell_line']='REH' if 'REH' in selected else selected[0]
        cell=st.selectbox('Cell line',selected,key='cell_line')
        result=pd.DataFrame({'feature':frame.index,'value':frame[cell].to_numpy()}).sort_values('value',ascending=False,na_position='last')
        fig=px.scatter(result.dropna(),x='feature',y='value',hover_name='feature')
        fig.update_layout(xaxis_title='Protein / drug (ranked by value)',yaxis_title=UNITS[kind],xaxis_showticklabels=False)
        st.caption(f'{cell}: all {len(result)} features retained in the table/export. No averaging or imputation.')
        params.update(cell_line=cell)
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
        if st.button('Compute embedding') or run_requested:
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
    session_key=None if identity else st.text_input('OpenLLM API key (optional if configured by host)',type='password',help='Used only in this browser session; not saved to disk.')
    st.caption('Your question and the displayed public-data summary are sent to SciLifeLab OpenLLM only when you click a button. Keys stay server-side.')
    prompt=st.text_area('Ask AI',placeholder='Plot BCL2 protein abundance for T-ALL, or explain the current analysis.',key='ask_ai_prompt')
    st.caption('Press Tab to complete the example prompt. Shift+Tab moves to the previous control.')
    components.html('<script>'+Path(__file__).with_name('prompt_completion.js').read_text(encoding='utf-8')+'</script>',height=0)
    if st.button('Run AI',type='primary',use_container_width=True,help='Read your prompt, apply validated analysis settings, and run the analysis.'):
        try:
            if not prompt.strip():
                st.warning('Enter a prompt, or press Tab in the input to use the example.');st.stop()
            group_options={k:metadata[k].group.unique().tolist() for k in metadata}
            instruction=GLOSSARY+' Translate into web settings only. Return JSON kind (protein/drug), module from '+json.dumps(MODULES)+', feature, groups(array). For Pairwise also second_kind, second_feature, method(Spearman/Pearson). For Cell-line profile optionally cell_line. No explanation or invented results. For all groups use all available groups in requested layer. Available groups: '+json.dumps(group_options)
            with st.spinner('Reading your prompt and preparing the analysis…'):
                with hosted_request(identity,prompt):
                    suggestion=chat([{'role':'system','content':instruction}, {'role':'user','content':prompt}],True,session_key)
                plan=validate_controls(suggestion,frames,metadata)
            st.session_state['suggestion']=plan
            st.session_state['pending_controls']=plan
            st.session_state['run_requested']=True
            st.rerun()
        except Exception:st.error('OpenLLM did not return valid controls. No analysis changed.')
    if 'suggestion' in st.session_state:
        plan=st.session_state['suggestion'];st.json(plan)
        if st.button('Apply suggested controls'):
            newkind=plan.get('kind',kind);newmodule=plan.get('module',module);newfeature=plan.get('feature')
            if newkind not in frames or newkind not in metadata:st.error('Unknown layer.');st.stop()
            if newmodule not in MODULES:st.error('Unknown module.');st.stop()
            validgroups=plan.get('groups') or metadata[newkind].group.unique().tolist()
            if any(g not in metadata[newkind].group.values for g in validgroups):st.error('Unknown groups.');st.stop()
            if newfeature and newfeature not in frames[newkind].index:st.error('Feature not found in original table.');st.stop()
            st.session_state['pending_controls']={**plan,'kind':newkind,'module':newmodule,'groups':validgroups}
            st.rerun()
    if st.button('Explain result',use_container_width=True,help='Explain the currently displayed computed result using your question.'):
        try:
            summary=result_context(result,params,sources,computed_statistics)
            with hosted_request(identity,prompt):
                answer=explain_result(summary,prompt,api_key=session_key)
            st.write(answer)
        except Exception:st.error('OpenLLM unavailable; your computed result remains available.')
    with st.expander('Data provenance'):st.json(sources)
    with st.expander('Trial scope and original FORALL coverage'):
        st.write('Independent AI-assisted protein/drug explorer; not a full replacement of original FORALL. Includes profiles, pairwise/group-adjusted correlations, full-layer correlation scans/ranks, cell-line profiles, exploratory group tests, sample-wise PCA/UMAP, metadata and exports.')
        st.write('Not yet included: RNA profiles/correlations, edgeR, DEqMS, fusion genes, GSEA/functional enrichment, correlation networks, subtype/consensus-cluster/custom sample grouping, drug-wise UMAP, or full original-method equivalence. AI suggests validated controls and explains computed results; some steps still require manual selection or Compute buttons.')
