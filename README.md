# SIDISH-derived glioma prognostic framework

Study-specific code accompanying *Evidence-gated translation of SIDISH-derived single-cell risk states into a prognostic framework for diffuse glioma*.

This repository links single-cell risk-state discovery to a fixed nine-gene expression score,
external survival assessment and matched clinical evaluation. It supplies the original
downstream R workflow and the additional analyses prepared for manuscript revision.

## Start here

- Browse `original_sources/` for the seven original R scripts and three independent audit scripts.
- Browse `revision_sources/` for the revision analysis, verification and working-figure sources.
- `SOURCE_MAP.csv` maps every readable source to its execution path in the complete packages.
- The current execution layouts, retained inputs and final figure-rendering sources are attached to [release v1.0.1](https://github.com/Liren100290/SIDISH-glioma-prognostic-framework/releases/tag/v1.0.1).
  Run analyses from the extracted release packages. The readable source folders provide a browsing view of the original and revision workflows.
- [Release v1.0.0](https://github.com/Liren100290/SIDISH-glioma-prognostic-framework/releases/tag/v1.0.0) remains available as the historical initial release.
- `locked_nine_gene_coefficients.csv` and `locked_model.json` record the fixed formula and cutoff.

## Download and verify

Download all five `Supplementary_Code_1_v07_*.zip` assets for the revision replay,
`Supplementary_Code_2_original_workflow_v03.zip` for the original sources and retained-input audits,
and `Supplementary_Code_3.zip` for final Figure 5 and Supplementary Figures 2-4.
The release also supplies `README.txt` and `SHA256SUMS.txt`. Each of the seven ZIPs is smaller than 20 MB.

Run the following from the repository directory:

```console
python download_release.py --output ./download --extract
cd download/SIDISH_revision_replay_v07
python RUN_REVISION.py --check-only
python RUN_REVISION.py --rscript PATH_TO_STATISTICS_RSCRIPT --graphics-rscript PATH_TO_GRAPHICS_RSCRIPT
```

The downloader uses Python's standard library and is pinned to release v1.0.1 and its
checksum-manifest digest. It verifies the seven ZIPs and release README before extracting
ZIP files only; it refuses a mismatched manifest or differing existing file. Keep
`download_release.py` and `SHA256SUMS.txt` together from the same repository revision.
Install analytical dependencies separately
using the supplied Python requirements and R dependency/version notes. The runner installs
no packages and does not delete files. On Windows, use a short ASCII destination path.

For the original retained-input audits, start again from the repository directory:

```console
cd download/SIDISH_original_workflow_v03
python audit_tools/CHECK_SOURCE_AND_SCORES.py
python audit_tools/AUDIT_EXACT_SOURCE_GATE.py
Rscript audit_tools/AUDIT_MODEL_FITS.R .
```

For final figure rendering, start from the repository directory and run in this order:

```console
cd download/SIDISH_final_plot_rendering
Rscript Figure5/Figure5.R
Rscript Supplementary/Draw_S2_S3.R
Rscript Supplementary/FigureS4.R
```

Code 3 supplies frozen inputs and panel objects. The last command produces the final
Supplementary Figure 4; the preceding command also produces an earlier view of that figure.
Use the graphics environment and dependencies documented in the package README.

## Analyses and scope

The revision package covers donor/cell-type/assay expression sensitivities, full-transcriptome
Hallmark testing, fixed-score evaluation in GSE16011, REMBRANDT and GSE43378, matched
clinical models with paired bootstrap and repeated cross-validation, and measured PRISM
glioma drug-response analyses. The gene coefficients and development cutoff are unchanged.

The recorded version 05 relocated replay reproduced 77 numerical tables exactly and seven
working PNGs byte-for-byte. Public distribution version 07 keeps those numerical scripts and
payloads unchanged. Version 03 of the original-workflow archive retains the seven original
R sources and independent coefficient/selected-model checks. See `REPRODUCIBILITY.md`.

The complete historical raw preprocessing, all 245 development configurations, all 100
bootstrap selections and historical GDSC training have not been independently replayed.
The development direction gate preceded LOCO and used all three development cohorts;
internal estimates are not fully nested external validation. Fixed-label sensitivities do
not retrain SIDISH. No new laboratory experiment is represented in this release.
Code 1 working-figure references support numerical reproduction. Code 3 separately
supplies the final Figure 5 and Supplementary Figures 2-4 rendering sources; all four
600-dpi TIFFs were reproduced byte-for-byte from the release inputs. The score is a
research prognostic framework.

## Upstream software, data and attribution

The published [SIDISH implementation](https://github.com/mcgilldinglab/SIDISH) is the upstream
algorithm, separately archived as [v1.0.0](https://doi.org/10.5281/zenodo.17420265).
This repository supplies the study-specific downstream workflow; it does not claim ownership
of SIDISH or replace its original authors' software. Ten archived local core files match
commit `9e9e86e16a08f88f0d4942af0854335b5acbe71b` after line-ending normalization.

See `DATA_SOURCES.md` and `THIRD_PARTY_NOTICES.md` for source links and rights.
Public-source derived inputs remain subject to their source terms; they are not relicensed
as author-owned software. Manuscript drafts, reviewer reports and private correspondence
are not part of this repository.

## Citation and version

Use the author/title information in `CITATION.cff` and cite the exact release tag/commit
used. Release v1.0.1 contains the final revision copies of Code 1 distribution 07 and
Code 2 public archive 03, plus Code 3 final figure rendering. Release v1.0.0 is preserved
as the initial-release record. This study has no publication DOI assigned in this
repository. No new Zenodo DOI is asserted.

## License

Original study code and newly authored reproducibility helpers are available under the MIT
license in `LICENSE`, as authorized by the author on 16 October 2026. Third-party software,
data, publications and gene-set material retain their source licenses; see `THIRD_PARTY_NOTICES.md`.
