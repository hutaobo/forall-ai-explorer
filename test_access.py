"""Security fixtures use ephemeral synthetic identities, never real member data."""
import base64,json,os,tempfile,time,uuid
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
from cryptography.hazmat.primitives.serialization import Encoding,PublicFormat
import seekyoo_access as access

key=Ed25519PrivateKey.generate()
public=key.public_key().public_bytes(Encoding.PEM,PublicFormat.SubjectPublicKeyInfo).decode()
encode=lambda x:base64.urlsafe_b64encode(json.dumps(x,separators=(',',':')).encode()).decode().rstrip('=')
def token(claims):
    body=encode({'alg':'EdDSA','typ':'JWT'})+'.'+encode(claims)
    return body+'.'+base64.urlsafe_b64encode(key.sign(body.encode())).decode().rstrip('=')
now=int(time.time())
claims={'protocol':access.PROTOCOL,'iss':access.PLATFORM,'aud':access.AUDIENCE,'sub':str(uuid.uuid4()),'jti':str(uuid.uuid4()),
        'iat':now,'exp':now+3600,'agent_id':access.AGENT,'manifest_sha256':'a'*64,'admin':False}
assert not access.decode_session(token(claims),public)['admin']
for invalid in [{**claims,'aud':'https://attacker.invalid'}, {**claims,'exp':now}, {**claims,'admin':'true'}]:
    try: access.decode_session(token(invalid),public)
    except Exception: pass
    else: raise AssertionError('Invalid session accepted')
with tempfile.TemporaryDirectory() as folder:
    os.environ['FORALL_PRIVATE_ROOT']=folder
    access.redeem_session(token(claims),public)
    try: access.redeem_session(token(claims),public)
    except Exception: pass
    else: raise AssertionError('Launch replay accepted')
    access.save_key('synthetic-test-key-never-a-real-credential')
    assert (access.private_root()/'openllm-key').stat().st_mode & 0o777 == 0o600
print('Member gate: signature, audience, expiry, role type, replay and owner-only credential storage passed.')
