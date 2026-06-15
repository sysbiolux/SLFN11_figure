# SLFN11 Figure

This repository generates SLFN11 expression figures from TCGA and GTEx Toil/Xena RNA-seq data.

## Outputs

The publication figures are written to `results/`:

- `slfn11_pan_cancer_violin_log2TPM.png`
- `slfn11_pan_cancer_violin_tpm.png`
- `slfn11_tcga_tumor_vs_normal.png`
- `slfn11_gtex_tissues.png`
- `slfn11_matched_tcga_gtex.png`

Raw downloaded data are cached in `results/raw/` and intentionally excluded from git.

## Requirements

Install the Python dependencies:

```bash
pip install -r requirements.txt
```

On Windows, use `py -m pip install -r requirements.txt` if `python` is not on `PATH`.

## Usage

Run the analysis from the repository root:

```bash
python script.py
```

On Windows, use `py script.py` if needed.

On the first run, the script downloads TCGA and GTEx expression and phenotype tables into `results/raw/`. Later runs reuse those files if present.

## Data Sources

Expression and phenotype data are downloaded from the UCSC Xena Toil recompute datasets:

- TCGA TPM matrix: `tcga_RSEM_gene_tpm.gz`
- GTEx TPM matrix: `gtex_RSEM_gene_tpm.gz`
- TCGA/TARGET/GTEx phenotype table: `TcgaTargetGTEX_phenotype.txt.gz`
- GTEx phenotype table: `GTEX_phenotype.gz`
