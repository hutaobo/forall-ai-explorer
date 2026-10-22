"""Synthetic, network-free checks for the explorer's unlimited trial gate."""
from unittest.mock import patch
from urllib.error import HTTPError
from seekyoo_access import hosted_request, request_json, AccessError

for admin in [False,True]:
    identity={'token':'synthetic-session','claims':{'admin':admin}}
    with patch('seekyoo_access.request_json',return_value={'valid':True}) as endpoint:
        with hosted_request(identity,'synthetic prompt'):
            pass
        assert endpoint.call_count==2
        assert all(call.args==('calls',{'action':'authorize','token':'synthetic-session'}) for call in endpoint.call_args_list)
with patch('urllib.request.urlopen',side_effect=HTTPError('https://example.test',401,'test',{},None)):
    try:request_json('calls',{})
    except AccessError as error:assert 'session expired' in str(error) and 'limit' not in str(error)
    else:raise AssertionError('Expired session accepted')
with patch('seekyoo_access.request_json',side_effect=AccessError('Account unavailable')):
    executed=False
    try:
        with hosted_request(identity,'synthetic prompt'):executed=True
    except AccessError:pass
    assert not executed
print('Passed: owner/member checks before and after AI, no budget reservations, expiry is distinct, revoked access stops execution.')
