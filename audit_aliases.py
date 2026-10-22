"""Read publisher metadata and exact channel evidence for five discrepant labels."""
import os,json
from pathlib import Path
import pandas as pd
from core import load,FILES
root=Path(os.environ['FORALL_DATA_ROOT'])
f,m,s=load(root)
cells=f['cells']
cells.columns=[' '.join(str(c).split()) for c in cells.columns]
mask=cells['Cell Line Name'].astype(str).str.contains('COG-LL-319|COG-LL-332|COG-LL-394|COG-LL-402|TMD',regex=True)
cols=['Cell Line Name','TMT10 channel','HiRIEF Set','RNA-seq sample ID','DSRT','Leukemia Lineage','Collection','Reference']
print('columns',list(cells.columns))
print(cells.loc[mask,[c for c in cols if c in cells.columns]].to_json(orient='records'))
for kind in ('protein','drug'):
    raw=pd.read_excel(root/FILES[kind],header=None,keep_default_na=False)
    hits=[i for i in range(1,raw.shape[1]) if any(t in str(raw.iloc[2,i]) for t in ['COG-LL-319','COG-LL-332','COG-LL-394','COG-LL-402','TMD'])]
    print(kind, raw.iloc[:8,[0]+hits].to_json(orient='split'))
