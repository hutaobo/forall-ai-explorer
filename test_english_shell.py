"""Synthetic session UI checks; no real accounts, credentials or network."""
import os,tempfile,time
from unittest.mock import patch
from streamlit.testing.v1 import AppTest

os.environ['FORALL_SEEKYOO_AUTH']='1'
script="from ui import theme\nfrom seekyoo_access import gate\nimport streamlit as st\nstyle=theme(st)\ngate(st)"
with tempfile.TemporaryDirectory() as folder:
    os.environ['FORALL_PRIVATE_ROOT']=folder
    app=AppTest.from_string(script).run()
    assert not app.exception
    assert any('You are not signed in' in x.value for x in app.info)
    assert len(app.get('link_button'))==2
    with patch('seekyoo_access.request_json',return_value={'valid':True,'email':'fictional@example.test','display_name':'<script>Fictional</script>','trial_unlimited':True}):
        app.session_state['seekyoo_session']={'token':'synthetic','claims':{'exp':int(time.time())+3600,'admin':False}}
        app.run()
        assert not app.exception
        text='\n'.join(x.value for x in app.markdown)
        assert 'SIGNED IN' in text and 'fictional@example.test' in text
        assert '&lt;script&gt;Fictional&lt;/script&gt;' in text
        assert '<script>Fictional</script>' not in text
        assert 'No daily app usage cap' in '\n'.join(x.value for x in app.caption)
print('Passed: English signed-out page, signed-in identity, unlimited-trial label and escaped display name.')
