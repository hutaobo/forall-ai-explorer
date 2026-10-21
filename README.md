# FORALL AI Explorer

Independent interactive exploration of the CC BY 4.0 supplements to
[Leo et al. 2022](https://doi.org/10.1038/s41467-022-29224-5).
This is not the original FORALL application, is not endorsed by its authors,
and does not reproduce DEqMS, edgeR or all upstream Shiny analyses.

Features: protein/drug profiles, group filters, replicate controls, pairwise
correlations, exploratory group comparisons with BH adjustment, PCA and UMAP,
metadata tables, CSV/HTML/parameter exports, and SciLifeLab OpenLLM assistance.
No upstream R source is copied. Values are not imputed or averaged. UMAP/PCA
drop incomplete features, select variable features, and record parameters.
Protein units are log2 relative protein abundance; drug response is sDSS.
Nominal correlation tests with replicate columns are exploratory and do not
establish independent biological replication or causal effects.

## Seekyoo member deployment

The container requires a signed Seekyoo member launch session before showing the
application. Open it from the FORALL card in Seekyoo. Launch links expire for
redemption after five minutes, are single-use, and establish a session of at most
one hour. AI calls revalidate membership, the live FORALL revision and existing
FORALL quotas through the platform; questions and API keys are not sent to that
quota endpoint. The platform public verification key is fetched from a fixed
HTTPS endpoint. Anonymous visitors cannot reach charts or shared-key AI calls.

Attach private persistent project storage at `/home/data`. The verified owner
may configure the shared OpenLLM key through the administrator expander, which
saves an owner-only file in `/home/data/forall-private`. Ordinary members have
no key input and never receive the key. Do not put the key into this repository,
container build arguments, Docker image or a public app download. Expired OpenLLM
keys must be replaced by the administrator. Local development retains optional
session-key entry when `FORALL_SEEKYOO_AUTH` is unset; the published container
sets it to `1` and fails closed when platform authorization is unavailable.

## Run

Install requirements in an isolated environment. Set FORALL_DATA_ROOT to the
four publisher workbooks identified in sources.json, then run:

```
streamlit run app.py --server.address=127.0.0.1 --server.port=8791
```

OpenLLM uses the official https://openllm.scilifelab.se/api endpoint.
Set OPENLLM_KEY_FILE to a server-private credential file, or supply your own
key in the session password input. Credentials must not be included in source,
images, publications, environment dumps, or logs. AI suggestions are validated
and applied by a separate user click; arbitrary model-generated code is never run.

Container port: 8501. The build downloads only the four public workbooks and
verifies publisher byte counts and MD5. No key is part of the image.
For initial Serve testing choose Private visibility. Published source and image
contain only this independent explorer, not the Seekyoo platform or node runtime.

## Data sources and differences

The supplement has 12,446 protein rows and 82 sample columns (51 primary columns
including EBV controls, plus 31 replicates), and 528 drugs across 43 cell lines.
The original 2023 FORALL app paper reports 12,326 proteins; these versions are
not assumed equivalent. RNA data are not bundled in this release.

Data license: CC BY 4.0, see sources.json. Software: MIT (LICENSE).
