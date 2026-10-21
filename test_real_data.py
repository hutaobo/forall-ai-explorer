"""Read-only integration validation against independently opened publisher workbooks."""
import json
import os
import numpy as np
from openpyxl import load_workbook
from core import load, samples, compare, embedding, bh, FILES, validate_controls

frames,meta,sources=load(os.environ['FORALL_DATA_ROOT'])
assert frames['protein'].shape==(12446,82)
assert frames['drug'].shape==(528,43)
assert len(samples(meta['protein'],meta['protein'].group.unique(),False))==51
assert len(samples(meta['protein'],meta['protein'].group.unique(),True))==82
assert frames['drug'].loc['Venetoclax','REH']==0.9
assert frames['drug'].loc['Venetoclax','SEM']==21.1
plan=validate_controls({'kind':'drug','module':'Profile','feature':'Venetoclax','groups':['BCP-ALL','T-ALL','B-ALL'],'explanation':'untrusted extra'},frames,meta)
assert set(plan)=={'kind','module','feature','groups'}
try:
    validate_controls({**plan,'groups':['EBV']},frames,meta)
except ValueError:
    pass
else:
    raise AssertionError('Unavailable drug group accepted')
assert np.isclose(frames['protein'].loc['BCL2','REH'],0.541117504490897)
for kind,feature,col,coordinate in [('drug','Venetoclax','REH','V8'),('protein','BCL2','REH','C1080')]:
    book=load_workbook(os.path.join(os.environ['FORALL_DATA_ROOT'],FILES[kind]),read_only=True,data_only=True)
    assert book.worksheets[0][coordinate].value==frames[kind].loc[feature,col]
    book.close()
assert np.allclose(bh([0.01,0.04,0.03]),[0.03,0.04,0.04])
primary=samples(meta['drug'],['BCP-ALL','T-ALL'],False)
out=compare(frames['drug'][primary],meta['drug'].loc[primary],'BCP-ALL','T-ALL','Welch t-test')
assert len(out)==528 and out.q_BH.dropna().between(0,1).all()
assert out.q_BH.notna().sum()>500
xy,details=embedding(frames['drug'][primary],'PCA',100,15,0.1,42)
assert xy.shape==(len(primary),2) and np.isfinite(xy.to_numpy()).all()
print(json.dumps({'real_workbook_validation':'passed','protein_shape':list(frames['protein'].shape),
    'drug_shape':list(frames['drug'].shape),'replicate_policy':'51 primary / 82 total columns',
    'independent_cell_checks':2,'BH_check':True,'comparison_features':len(out),'PCA_samples':len(xy)}))
