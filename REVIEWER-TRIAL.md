# FORALL AI Explorer — reviewer trial

Entry: https://app.seekyoo.com/api/play/tools/forall/launch (register, verify email and sign in).
Hosted app: https://forall-ai-explorer.serve.scilifelab.se/ (requires a fresh Seekyoo launch).

Independent prototype using published CC BY 4.0 tables and SciLifeLab OpenLLM. It is not the original FORALL, an author-endorsed service, or a validated replacement. API credentials remain server-side. AI suggests validated controls and explains computed data; it does not invent statistics. Some steps still require Apply or Compute buttons. Users' requests and public-data context go to SciLifeLab OpenLLM when an AI button is clicked.

## Suggested trial

1. In Ask AI, request: “Show the Spearman correlation between BCL2 protein abundance and Venetoclax drug sensitivity across all available groups.” Click AI: suggest controls, then Apply suggested controls.
2. Inspect 43 audited pairs, within-group results and Group-adjusted association. Hover to see original cross-table labels. Compare global-rank residuals with group-demeaned values; download CSV, HTML and statistics.
3. Select Correlation scan with BCL2 against drug, then Compute correlation scan. Inspect all 528 drug correlations, pair counts, nominal p and BH q.
4. Select drug / Cell-line profile / REH to inspect all 528 drug scores for one cell line.
5. Ask AI to explain the current result, including uncertainty and the difference between pooled and group-adjusted associations.

## Comparison with original FORALL

Audit basis: live original module menu plus the 17 functions described in [Aswad & Jafari 2023](https://doi.org/10.1093/bioadv/vbad143). New-side evidence is source inspection, real-data tests and browser acceptance; this is not exhaustive numerical equivalence testing against every original module.

| Original function (paper item) | Trial status / difference |
|---|---|
| 1. Gene profiles | Partial: protein; RNA absent |
| 2. RNA–protein association | Absent |
| 3. Grouped profiles | Partial: lineage only |
| 4. Pair correlation by lineage | Protein/drug supported; RNA absent; group table and filter rather than identical layout |
| 5. Custom sample subsets | Partial: lineage filter; no arbitrary cell-name subset |
| 6. Feature-to-all correlations | Protein/drug scan added; RNA absent |
| 7. Correlation network | Absent |
| 8. Correlation ranks | Partial: interactive ranked scatter/table; original waterfall/highlighting not replicated |
| 9. Drug/cell sensitivity views | Drug across cells and cell across drugs supported; layout differs |
| 10. RNA differential testing | edgeR absent |
| 11. Protein differential testing | DEqMS absent; Welch/Mann–Whitney exploratory alternatives are not equivalent |
| 12. Drug group testing | Partial: Welch/Mann–Whitney + BH; original arbitrary groups not available |
| 13. Two-group plots | Partial: existing lineage groups in profile boxplots |
| 14. UMAP | Partial: sample-wise protein/drug; no RNA or drug-wise mode; PCA added |
| 15. Characteristics | Published cell/drug metadata available; full original annotation parity unverified |
| 16. Fusion genes | Absent |
| 17. Enrichment | GSEA/functional enrichment absent |

Exports: interactive HTML, CSV, provenance and Plotly SVG; original PDF workflow and per-module guides not reproduced. Original paper reports 12,326 protein features; the verified publisher workbook here has 12,446. Dataset/method parity must be resolved before claiming equivalence. Original RNA layer is not bundled.

Added capabilities: Seekyoo member access, administrator-managed OpenLLM connection, validated AI parameter suggestions, full selected profile/pair context, explicit label audit, group-adjusted rank/linear analysis and reproducibility exports. AI currently sets layer/module/first feature/groups, both pair variables/method, and optional cell line; arbitrary multi-step workflows and all numerical controls are not fully automated.

Feedback requested: scientific correctness, missing high-value original workflows, clarity of AI suggestions/uncertainty, labels and usability. The current prototype is suitable for feedback under these limits, not presentation as full original-function coverage.
