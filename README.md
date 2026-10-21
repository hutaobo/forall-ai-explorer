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
