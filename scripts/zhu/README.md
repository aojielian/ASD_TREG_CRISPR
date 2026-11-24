# Zhu et al. candidate evidence

The portable extractor transcribes published gene-level summary statistics and measures candidate coverage in conventional CD4+ T-cell CRISPRi assays. It uses Python 3.11 or later and the standard library. Run the commands from the repository root. Set `DATA_DIR` to the extracted processed-data archive root.

## Verify the distributed tables

```sh
python scripts/zhu/extract_zhu_evidence.py verify \
  --derived-dir "$DATA_DIR/data/r8/zhu_derived" \
  --candidate-master "$DATA_DIR/data/r7/mechanism/CANDIDATE_DIRECTION_MECHANISM.tsv"
```

This offline check validates the frozen candidate order, condition matrix, missing entries, guide-quality totals, split assay denominators, and original ZMYND8 anchors. The project plots consume the processed tables. See `data/r8/zhu_derived/DATA_DICTIONARY.md` inside the data package for full schema and statistical definitions.

## Reconstruct from the study sources

The optional source reconstruction requires the eight exact files linked in the processed archive's `data/r8/zhu_derived/ZHU_SOURCE_FILES.tsv`. Use the author file links there and retain their basenames in one local directory. Retrieve the original file bytes (GitHub's raw download action). Reference the listed SHA256 digests to identify the snapshot if the source branch changes. Consult the sources' access and reuse terms before retrieving third-party data. Set `ZHU_SOURCE_DIR` to the directory containing those files.

```sh
python scripts/zhu/extract_zhu_evidence.py extract \
  --source-dir "$ZHU_SOURCE_DIR" \
  --candidate-master "$DATA_DIR/data/r7/mechanism/CANDIDATE_DIRECTION_MECHANISM.tsv" \
  --out-dir reconstructed-zhu
```

Choose a new output directory outside the source directory. The command hashes the source bytes and reconstructs the five TSVs. It keeps donor labels within memory long enough to determine matched coverage and writes only counts. Summary effect sizes and adjusted P values are sourced from the authors' table. The source assay comparison figures and publication are attributed in the bibliography.
