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
import urllib.error
import uuid
from ui import sign_in, account

PLATFORM = 'https://app.seekyoo.com'
AUDIENCE = 'https://forall-ai-explorer.serve.scilifelab.se'
PROTOCOL = 'seekyoo-forall-session/1'
AGENT = 'a92604df-ec7c-462c-a246-32e179d11b8f'

class AccessError(ValueError):
    """A safe, user-facing platform status without response or credential data."""

def private_root():
    root = Path(os.environ.get('FORALL_PRIVATE_ROOT', '/home/data/forall-private'))
    root.mkdir(mode=0o700, parents=True, exist_ok=True)
    os.chmod(root,0o700)
    return root

def request_json(path, body=None):
    url = PLATFORM + '/api/play/tools/forall/' + path
    req = urllib.request.Request(url, data=None if body is None else json.dumps(body).encode(),
                                 headers={'Content-Type':'application/json'})
    try:
        with urllib.request.urlopen(req, timeout=20) as response:
            if response.geturl() != url: raise ValueError('Unexpected platform redirect')
            return json.loads(response.read(16000))
    except urllib.error.HTTPError as error:
        messages={401:'Your Seekyoo session expired. Reopen FORALL through Seekyoo.',
                  403:'This Seekyoo account does not currently have access.',
                  409:'FORALL changed. Reopen it through Seekyoo.',
                  429:'Seekyoo is temporarily busy. Please retry shortly.'}
        raise AccessError(messages.get(error.code,f'Seekyoo account verification is unavailable (HTTP {error.code}). Please retry.')) from None
    except Exception:
        raise AccessError('Cannot reach Seekyoo to verify your account. Please retry shortly.') from None

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
            sign_in(st,expired=True)
            st.stop()
    identity=st.session_state.get('seekyoo_session')
    if not identity or identity['claims']['exp']<=time.time():
        sign_in(st,expired=bool(identity))
        st.stop()
    try:
        identity['account']=request_json('calls',{'action':'authorize','token':identity['token']})
    except AccessError as error:
        sign_in(st)
        st.error(str(error));st.stop()
    os.environ['OPENLLM_KEY_FILE']=str(private_root()/'openllm-key')
    if (private_root()/'openllm-key').is_file(): os.chmod(private_root()/'openllm-key',0o600)
    account(st,identity)
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
    # Independent explorer trial: live eligibility checks, no marketplace budget.
    # The original agent's reserve/validate/finish policy is unchanged.
    if identity is None:
        yield
        return
    request={'action':'authorize','token':identity['token']}
    request_json('calls',request)
    yield
    request_json('calls',request)
