"""Seekyoo member sessions and server-side credentials for the hosted explorer."""
import base64
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import time
import urllib.request
import uuid

PLATFORM = 'https://app.seekyoo.com'
AUDIENCE = 'https://forall-ai-explorer.serve.scilifelab.se'
PROTOCOL = 'seekyoo-forall-session/1'
AGENT = 'a92604df-ec7c-462c-a246-32e179d11b8f'

def private_root():
    root = Path(os.environ.get('FORALL_PRIVATE_ROOT', '/home/data/forall-private'))
    root.mkdir(mode=0o700, parents=True, exist_ok=True)
    return root

def request_json(path, body=None):
    url = PLATFORM + '/api/play/tools/forall/' + path
    req = urllib.request.Request(url, data=None if body is None else json.dumps(body).encode(),
                                 headers={'Content-Type':'application/json'})
    try:
        with urllib.request.urlopen(req, timeout=20) as response:
            if response.geturl() != url: raise ValueError('Unexpected platform redirect')
            return json.loads(response.read(16000))
    except Exception:
        raise ValueError('Seekyoo authorization unavailable or usage limit reached. Return to Seekyoo to retry.') from None

def decode_session(token, public_key, now=None):
    from cryptography.hazmat.primitives.serialization import load_pem_public_key
    from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PublicKey
    now = int(time.time() if now is None else now)
    if not isinstance(token,str) or len(token)>4096: raise ValueError('Invalid session')
    parts=token.split('.')
    if len(parts)!=3: raise ValueError('Invalid session')
    decode=lambda value:base64.urlsafe_b64decode(value+'='*(-len(value)%4))
    if json.loads(decode(parts[0])) != {'alg':'EdDSA','typ':'JWT'}: raise ValueError('Invalid session')
    key=load_pem_public_key(public_key.encode())
    if not isinstance(key,Ed25519PublicKey): raise ValueError('Invalid key type')
    key.verify(decode(parts[2]),(parts[0]+'.'+parts[1]).encode())
    c=json.loads(decode(parts[1]))
    if c.get('protocol')!=PROTOCOL or c.get('iss')!=PLATFORM or c.get('aud')!=AUDIENCE or c.get('agent_id')!=AGENT:
        raise ValueError('Invalid audience')
    if type(c.get('iat')) is not int or type(c.get('exp')) is not int or c['iat']>now+30 or c['exp']<=now or c['exp']-c['iat']!=3600:
        raise ValueError('Expired session')
    uuid.UUID(c['sub']);uuid.UUID(c['jti'])
    if type(c.get('admin')) is not bool: raise ValueError('Invalid role')
    return c

def redeem_session(token, public_key):
    c=decode_session(token,public_key)
    if int(time.time())-c['iat']>300: raise ValueError('Launch link expired')
    db=sqlite3.connect(private_root()/'sessions.sqlite', timeout=10)
    try:
        with db:
            db.execute('CREATE TABLE IF NOT EXISTS redeemed (id TEXT PRIMARY KEY, expires INTEGER NOT NULL)')
            db.execute('DELETE FROM redeemed WHERE expires < ?', (int(time.time()),))
            db.execute('INSERT INTO redeemed VALUES (?,?)',(c['jti'],c['exp']))
    finally: db.close()
    return c

def save_key(value):
    value=value.strip()
    if not 20<=len(value)<=1000 or any(ch.isspace() for ch in value): raise ValueError('Invalid key format')
    root=private_root(); target=root/'openllm-key'; temp=root/('key-'+uuid.uuid4().hex)
    fd=os.open(temp,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
    try:
        with os.fdopen(fd,'w') as stream: stream.write(value)
        os.replace(temp,target)
    finally:
        if temp.exists(): temp.unlink()

def gate(st):
    if os.environ.get('FORALL_SEEKYOO_AUTH')!='1': return None
    if 'seekyoo_ticket' in st.query_params:
        token=st.query_params['seekyoo_ticket']
        st.query_params.clear()
        try:
            public=request_json('key')
            if public.get('protocol')!=PROTOCOL: raise ValueError('Unknown issuer')
            claims=redeem_session(token,public['public_key'])
            st.session_state['seekyoo_session']={'token':token,'claims':claims}
        except Exception:
            st.error('登录链接已过期或不可用，请从 Seekyoo 重新打开。 / Reopen through Seekyoo.')
            st.stop()
    identity=st.session_state.get('seekyoo_session')
    if not identity or identity['claims']['exp']<=time.time():
        st.info('请登录 Seekyoo 后打开交互分析。 / Sign in through Seekyoo to use this tool.')
        st.link_button('Open through Seekyoo',PLATFORM+'/api/play/tools/forall/launch')
        st.stop()
    os.environ['OPENLLM_KEY_FILE']=str(private_root()/'openllm-key')
    st.caption('Seekyoo member session · AI provided by the administrator · session valid for up to one hour')
    if identity['claims']['admin']:
        with st.expander('Administrator: OpenLLM connection'):
            configured=(private_root()/'openllm-key').is_file()
            st.caption('Configured' if configured else 'An administrator must configure the shared OpenLLM connection.')
            with st.form('shared_openllm_setup',clear_on_submit=True):
                value=st.text_input('Administrator OpenLLM key',type='password')
                if st.form_submit_button('Save server-side key'):
                    try: save_key(value);st.success('Saved in private Serve storage. Members do not receive the key.')
                    except ValueError: st.error('Check the key format; it was not saved.')
    return identity

@contextmanager
def hosted_request(identity, question):
    if identity is None:
        yield
        return
    token=identity['token']
    reservation=request_json('calls',{'action':'reserve','token':token,'intent_id':str(uuid.uuid4()),
                                    'input_sha256':hashlib.sha256(question.encode()).hexdigest()})
    call_id=reservation['call_id']; success=False
    try:
        if not request_json('calls',{'action':'validate','token':token,'call_id':call_id}).get('valid'):
            raise ValueError('Permission changed')
        yield
        success=True
    finally:
        result=request_json('calls',{'action':'finish','token':token,'call_id':call_id,'success':success})
        if success and not result.get('finished'): raise ValueError('Permission changed; answer discarded')
