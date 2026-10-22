"""UI-only regression with explicitly synthetic data and mocked OpenLLM.

No network request, credentials, or scientific conclusions are involved.
"""
import os
import pandas as pd
from unittest.mock import patch
from streamlit.testing.v1 import AppTest
import core

os.environ['FORALL_SEEKYOO_AUTH']='0'
names=['synthetic-B','synthetic-T1','synthetic-T2']
meta=pd.DataFrame({'sample':names,'group':['BCP-ALL','T-ALL','T-ALL'],
                   'replicate':[False]*3,'excel_column':['C','D','E']},index=names)
frames={'protein':pd.DataFrame([[1.,2.,3.]],index=['BCL2'],columns=names),
        'drug':pd.DataFrame([[4.,5.,6.]],index=['Venetoclax'],columns=names)}
plan={'kind':'protein','module':'Profile','feature':'BCL2','groups':['T-ALL']}
with patch.object(core,'load',return_value=(frames,{'protein':meta,'drug':meta},{})), patch.object(core,'chat',return_value=plan):
    app=AppTest.from_file('app.py').run()
    assert not app.exception
    app.text_area(key='ask_ai_prompt').set_value('Plot BCL2 for T-ALL')
    next(b for b in app.button if b.label=='Run AI').click().run()
    assert not app.exception
    assert app.multiselect(key='groups').value==['T-ALL']
    assert len(app.dataframe[0].value)==2
    assert app.text_area(key='ask_ai_prompt').value=='Plot BCL2 for T-ALL'
    assert any(b.label=='Explain result' for b in app.button)
print('Synthetic UI fixture passed: Run AI applies validated settings and redraws without an Apply click.')
