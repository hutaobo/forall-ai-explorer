"""Run this independent SSCB trial on one allowed logical CPU."""
import os
import sys
if hasattr(os,'sched_setaffinity'):
    os.sched_setaffinity(0,{max(os.sched_getaffinity(0))})
os.execv(sys.executable,[sys.executable,'-m','streamlit','run','app.py','--server.address=127.0.0.1','--server.port=8791','--browser.gatherUsageStats=false'])
