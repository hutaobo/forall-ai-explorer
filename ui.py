"""English research workspace styling, using Seekyoo's ivory/teal palette."""
import html
import time

def theme(st):
    st.markdown('''<style>
    :root { --ink:#23463b; --teal:#0e6f6a; --paper:#fafaf4; --line:#dfe5d7; }
    .stApp { background:var(--paper); color:var(--ink); }
    [data-testid="stHeader"] { background:var(--paper); }
    [data-testid="stSidebar"] { background:#f0f3e9; border-right:1px solid var(--line); }
    h1,h2,h3 { color:var(--ink); letter-spacing:-.035em; }
    a { color:var(--teal)!important; }
    .block-container { padding-top:2rem; max-width:1600px; }
    .sy-brand { display:flex; align-items:center; gap:.7rem; color:var(--teal); font-weight:750; font-size:1.2rem; }
    .sy-mark { width:32px; height:32px; display:inline-flex; align-items:center; justify-content:center; background:var(--teal); color:white; border-radius:10px; }
    .sy-tag { color:#627361; font-size:.75rem; text-transform:uppercase; letter-spacing:.14em; margin:1.6rem 0 .7rem; }
    .sy-hero { padding:2.5rem; border:1px solid var(--line); border-radius:22px; background:#fffef9; margin:1.25rem 0; }
    .sy-hero h1 { font-size:clamp(2.1rem,4vw,3.5rem); line-height:1.12; margin:.6rem 0 1rem; }
    .sy-hero p { color:#627361; max-width:700px; line-height:1.7; }
    .sy-account { padding:1rem; border:1px solid #cddcce; border-radius:14px; background:#fffef9; margin:1rem 0; overflow-wrap:anywhere; }
    .sy-account .status { font-size:.78rem; color:var(--teal); font-weight:650; margin-bottom:.65rem; }
    .sy-account .name { font-size:1rem; font-weight:700; }
    .sy-account .email,.sy-account .expiry { font-size:.8rem; color:#627361; margin-top:.3rem; }
    .sy-feature { border-top:1px solid var(--line); padding:1.2rem 0; }
    .sy-feature b { color:var(--ink); }
    .sy-feature p { color:#627361; margin:.4rem 0; }
    [data-testid="stBaseButton-primary"],button[kind="primary"] { background:var(--teal)!important; border-color:var(--teal)!important; border-radius:10px; min-height:44px; }
    [data-testid="stBaseButton-secondary"], [data-testid="stLinkButton"] a { border-color:#ccd8c9; border-radius:10px; }
    [data-testid="stTextArea"] textarea { background:#fffef9!important; border-radius:12px; }
    [data-testid="stExpander"] { background:#fffef9; border-radius:12px; }
    [data-testid="stPlotlyChart"] { border:1px solid var(--line); border-radius:16px; background:white; padding:.5rem; }
    [data-testid="stVerticalBlockBorderWrapper"] { border-color:var(--line); }
    @media(max-width:700px) { .sy-hero { padding:1.5rem; } .block-container { padding:1rem; } }
    </style>''',unsafe_allow_html=True)

def brand(st):
    st.markdown('<div class="sy-brand"><span class="sy-mark">s</span> seekyoo <span style="color:#98a78e;font-weight:400">/</span> research tools</div>',unsafe_allow_html=True)

def sign_in(st, expired=False, unavailable=False):
    brand(st)
    st.markdown('''<section class="sy-hero"><div class="sy-tag">FORALL · AI Explorer</div>
    <h1>From a question<br>to an interactive discovery.</h1>
    <p>Explore published ALL cell-line protein expression and drug sensitivity.
    Choose your own parameters, or ask AI to set up an analysis and explain the results.</p></section>''',unsafe_allow_html=True)
    left,right=st.columns([1.5,1],gap='large')
    with left:
        st.subheader('Your research workspace starts here')
        st.write('Sign in with your Seekyoo account, then open FORALL. The shared SciLifeLab OpenLLM connection is provided for you.')
        status='Account verification is temporarily unavailable. Refresh your session or retry shortly.' if unavailable else ('Session expired. Sign in again to continue.' if expired else 'You are not signed in to this FORALL workspace.')
        st.info(status)
        st.link_button('1 · Sign in to Seekyoo','https://app.seekyoo.com/?lang=en',type='primary',use_container_width=True)
        st.link_button('2 · Open FORALL with my account','https://app.seekyoo.com/api/play/tools/forall/launch',use_container_width=True)
        st.caption('Already signed in to Seekyoo? Use step 2. New members can register and verify their email in step 1. No OpenLLM key needed.')
    with right:
        st.markdown('''<div class="sy-account"><div class="status">MEMBER TRIAL</div>
        <div class="name">No daily app usage cap</div><p>Sign-in is required. Your account and access are verified before each AI request.</p></div>
        <div class="sy-feature"><b>Explore the published data</b><p>Interactive charts, correlations and group comparisons.</p></div>
        <div class="sy-feature"><b>Ask, run, understand</b><p>Use natural language to configure analysis and interpret computed results.</p></div>
        <div class="sy-feature"><b>Keep the evidence</b><p>Export data, interactive figures and analysis settings.</p></div>''',unsafe_allow_html=True)
    st.caption('Hosted on SciLifeLab Serve · Powered by SciLifeLab OpenLLM · Independent implementation; not endorsed by the original FORALL authors.')

def account(st, identity):
    info=identity['account']; claims=identity['claims']
    name=html.escape(info['display_name']); email=html.escape(info['email'])
    minutes=max(0,int((claims['exp']-time.time())/60))
    with st.sidebar:
        brand(st)
        st.markdown(f'''<div class="sy-account"><div class="status">● SIGNED IN · SEEKYOO</div>
        <div class="name">{name}</div><div class="email">{email}</div>
        <div class="expiry">Session expires in {minutes} min</div></div>''',unsafe_allow_html=True)
        st.caption('Member trial · No daily app usage cap')
        st.link_button('My Seekyoo account','https://app.seekyoo.com/app?lang=en',use_container_width=True)
        st.link_button('Refresh session','https://app.seekyoo.com/api/play/tools/forall/launch',use_container_width=True)
        st.divider()
