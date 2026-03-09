"""
SLFN11 pan-cancer analysis

This script:
1. Downloads TCGA and GTEx expression data (Toil recompute)
2. Extracts SLFN11 expression
3. Generates a pan-cancer expression figure

Run simply with:

python slfn11_analysis.py

Output:
results/slfn11_pan_cancer.png
"""

import os
import gzip
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# CONFIG
# -----------------------------

GENE = "ENSG00000172716"
OUTDIR = "results"

TCGA_URL = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz"
GTEX_URL = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/gtex_RSEM_gene_tpm.gz"

TCGA_PHENO_URL = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz"
GTEX_PHENO_URL = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/GTEX_phenotype.gz"

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(f"{OUTDIR}/raw", exist_ok=True)


# -----------------------------
# DOWNLOAD
# -----------------------------

def download(url, outpath):

    if os.path.exists(outpath):
        return

    print("Downloading", url)

    r = requests.get(url, stream=True)
    r.raise_for_status()

    with open(outpath, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)


download(TCGA_URL, f"{OUTDIR}/raw/tcga_tpm.gz")
download(GTEX_URL, f"{OUTDIR}/raw/gtex_tpm.gz")
download(TCGA_PHENO_URL, f"{OUTDIR}/raw/tcga_pheno.gz")
download(GTEX_PHENO_URL, f"{OUTDIR}/raw/gtex_pheno.gz")


# -----------------------------
# EXTRACT GENE
# -----------------------------

def extract_gene(matrix_path, gene):

    print("Extracting", gene, "from", matrix_path)

    with gzip.open(matrix_path, "rt") as f:

        header = f.readline().strip().split("\t")
        samples = header[1:]

        for line in f:

            fields = line.split("\t")
            symbol_pre = fields[0].split("|")[0]
            symbol = symbol_pre.split(".")[0]
            print(symbol)

            if symbol == gene:

                values = np.array(fields[1:], dtype=float)

                return pd.Series(values, index=samples)

    raise Exception("Gene not found")


tcga_expr = extract_gene(f"{OUTDIR}/raw/tcga_tpm.gz", GENE)
gtex_expr = extract_gene(f"{OUTDIR}/raw/gtex_tpm.gz", GENE)


# -----------------------------
# LOAD PHENOTYPE
# -----------------------------

print("Loading phenotype tables")

tcga_pheno = pd.read_csv(
    f"{OUTDIR}/raw/tcga_pheno.gz",
    sep="\t",
    compression="gzip",
    encoding="latin1",
    low_memory=False
)

gtex_pheno = pd.read_csv(
    f"{OUTDIR}/raw/gtex_pheno.gz",
    sep="\t",
    compression="gzip",
    encoding="latin1",
    low_memory=False
)


# -----------------------------
# ANNOTATE TCGA
# -----------------------------

print("Preparing TCGA table")

df = pd.DataFrame({
    "sample": tcga_expr.index,
    "expr": tcga_expr.values
})

df["logexpr"] = np.log2(df["expr"] + 1)

df["sample_type"] = df["sample"].str[13:15]

df["type"] = np.where(df["sample_type"] == "01", "tumor",
             np.where(df["sample_type"] == "11", "normal", "other"))

df = df[df["type"] != "other"]


# -----------------------------
# TCGA CANCER TYPE
# -----------------------------

if "sample" in tcga_pheno.columns:
    tcga_pheno = tcga_pheno.set_index("sample")

df["cancer"] = tcga_pheno.loc[df["sample"], "primary disease or tissue"].values


# -----------------------------
# GTEX
# -----------------------------

print("Preparing GTEx table")

gtex_df = pd.DataFrame({
    "sample": gtex_expr.index,
    "expr": gtex_expr.values
})

gtex_df["logexpr"] = np.log2(gtex_df["expr"] + 1)

if "Sample" in gtex_pheno.columns:
    gtex_pheno = gtex_pheno.set_index("Sample")

gtex_df["tissue"] = gtex_pheno.loc[gtex_df["sample"], "body_site_detail (SMTSD)"].values


# -----------------------------
# FIGURE
# -----------------------------
print("Generating figure (violin plot, TPM units, TCGA acronyms)")

import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# CONFIG
# -----------------------------
PSEUDOCOUNT = 1e-1          # because your floor is ~log2(0.001) = -9.965784...
MIN_N = 30                  # minimum samples per cancer type to show
X_FONTSIZE = 6              # smaller x-axis labels
OUTPNG = f"{OUTDIR}/slfn11_pan_cancer_violin_tpm.png"

TCGA_ACRONYM = {
    "Adrenocortical Cancer": "ACC",
    "Bladder Urothelial Carcinoma": "BLCA",
    "Breast Invasive Carcinoma": "BRCA",
    "Cervical & Endocervical Cancer": "CESC",
    "Cholangiocarcinoma": "CHOL",
    "Colon Adenocarcinoma": "COAD",
    "Diffuse Large B-Cell Lymphoma": "DLBC",
    "Esophageal Carcinoma": "ESCA",
    "Glioblastoma Multiforme": "GBM",
    "Head & Neck Squamous Cell Carcinoma": "HNSC",
    "Kidney Chromophobe": "KICH",
    "Kidney Clear Cell Carcinoma": "KIRC",
    "Kidney Papillary Cell Carcinoma": "KIRP",
    "Liver Hepatocellular Carcinoma": "LIHC",
    "Lung Adenocarcinoma": "LUAD",
    "Lung Squamous Cell Carcinoma": "LUSC",
    "Mesothelioma": "MESO",
    "Ovarian Serous Cystadenocarcinoma": "OV",
    "Pancreatic Adenocarcinoma": "PAAD",
    "Pheochromocytoma & Paraganglioma": "PCPG",
    "Prostate Adenocarcinoma": "PRAD",
    "Rectum Adenocarcinoma": "READ",
    "Sarcoma": "SARC",
    "Skin Cutaneous Melanoma": "SKCM",
    "Stomach Adenocarcinoma": "STAD",
    "Testicular Germ Cell Tumor": "TGCT",
    "Thymoma": "THYM",
    "Thyroid Carcinoma": "THCA",
    "Uterine Carcinosarcoma": "UCS",
    "Uterine Corpus Endometrioid Carcinoma": "UCEC",
    "Uveal Melanoma": "UVM",
    "Brain Lower Grade Glioma": "LGG",
}

# -----------------------------
# PREP DATA
# -----------------------------
# Expect df to already contain:
#   - df["type"]   : "tumor"/"normal"/...
#   - df["cancer"] : cancer label (full name)
#   - df["logexpr"]: expression as provided by Toil/Xena (log2(TPM+0.001) style)
tumor = df[df["type"] == "tumor"].copy()
tumor = tumor.dropna(subset=["cancer", "logexpr"]).copy()

# Convert back to TPM: TPM = 2^x - 0.001
tumor["tpm"] = (2.0 ** tumor["logexpr"].astype(float)) - PSEUDOCOUNT
tumor.loc[tumor["tpm"] < 0, "tpm"] = 0.0
tumor = tumor[np.isfinite(tumor["tpm"].to_numpy(dtype=float))].copy()

# Filter cohorts by sample count
counts = tumor["cancer"].value_counts()
keep = counts[counts >= MIN_N].index
tumor = tumor[tumor["cancer"].isin(keep)].copy()

# Order by median TPM (descending)
order_full = (
    tumor.groupby("cancer")["tpm"]
    .median()
    .sort_values(ascending=False)
    .index
    .tolist()
)

# Build violin groups
groups = []
labels_full = []
for c in order_full:
    v = tumor.loc[tumor["cancer"] == c, "tpm"].to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    if v.size >= MIN_N:
        groups.append(v)
        labels_full.append(c)

# Convert labels to TCGA acronyms (fallback to full name if not in dict)
labels = [TCGA_ACRONYM.get(c, c) for c in labels_full]

print(f"Plotting {len(labels)} cancer types with n >= {MIN_N}")
missing = [c for c in labels_full if c not in TCGA_ACRONYM]
if missing:
    print("[WARN] No acronym mapping for:", ", ".join(missing))

# -----------------------------
# PLOT
# -----------------------------
plt.figure(figsize=(14, 8))

parts = plt.violinplot(
    groups,
    showmeans=False,
    showmedians=True,
    showextrema=False
)

plt.xticks(
    ticks=np.arange(1, len(labels) + 1),
    labels=labels,
    rotation=90,
    fontsize=X_FONTSIZE
)

plt.yscale("log")
plt.ylabel("SLFN11 expression (TPM)")
plt.title("SLFN11 expression across TCGA cancers (tumor samples)")

plt.tight_layout()
plt.savefig(OUTPNG, dpi=300)

print("Saved figure to", OUTPNG)

###################

print("Generating figure (seaborn violin, log2(TPM+1e-3) units, TCGA acronyms)")

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------------
# CONFIG
# -----------------------------
MIN_N = 30
X_FONTSIZE = 7
OUTPNG = f"{OUTDIR}/slfn11_pan_cancer_violin_log2TPM.png"

TCGA_ACRONYM = {
    "Adrenocortical Cancer": "ACC",
    "Bladder Urothelial Carcinoma": "BLCA",
    "Breast Invasive Carcinoma": "BRCA",
    "Cervical & Endocervical Cancer": "CESC",
    "Cholangiocarcinoma": "CHOL",
    "Colon Adenocarcinoma": "COAD",
    "Diffuse Large B-Cell Lymphoma": "DLBC",
    "Esophageal Carcinoma": "ESCA",
    "Glioblastoma Multiforme": "GBM",
    "Head & Neck Squamous Cell Carcinoma": "HNSC",
    "Kidney Chromophobe": "KICH",
    "Kidney Clear Cell Carcinoma": "KIRC",
    "Kidney Papillary Cell Carcinoma": "KIRP",
    "Liver Hepatocellular Carcinoma": "LIHC",
    "Lung Adenocarcinoma": "LUAD",
    "Lung Squamous Cell Carcinoma": "LUSC",
    "Mesothelioma": "MESO",
    "Ovarian Serous Cystadenocarcinoma": "OV",
    "Pancreatic Adenocarcinoma": "PAAD",
    "Pheochromocytoma & Paraganglioma": "PCPG",
    "Prostate Adenocarcinoma": "PRAD",
    "Rectum Adenocarcinoma": "READ",
    "Sarcoma": "SARC",
    "Skin Cutaneous Melanoma": "SKCM",
    "Stomach Adenocarcinoma": "STAD",
    "Testicular Germ Cell Tumor": "TGCT",
    "Thymoma": "THYM",
    "Thyroid Carcinoma": "THCA",
    "Uterine Carcinosarcoma": "UCS",
    "Uterine Corpus Endometrioid Carcinoma": "UCEC",
    "Uveal Melanoma": "UVM",
    "Brain Lower Grade Glioma": "LGG",
}

TCGA_TO_GTEX_TISSUE = {
    "Breast Invasive Carcinoma": "Breast",
    "Colon Adenocarcinoma": "Colon - Transverse",
    "Rectum Adenocarcinoma": "Colon - Sigmoid",
    "Esophageal Carcinoma": "Esophagus - Mucosa",
    "Liver Hepatocellular Carcinoma": "Liver",
    "Lung Adenocarcinoma": "Lung",
    "Lung Squamous Cell Carcinoma": "Lung",
    "Pancreatic Adenocarcinoma": "Pancreas",
    "Prostate Adenocarcinoma": "Prostate",
    "Stomach Adenocarcinoma": "Stomach",
    "Thyroid Carcinoma": "Thyroid",
    "Bladder Urothelial Carcinoma": "Bladder",
    "Kidney Clear Cell Carcinoma": "Kidney - Cortex",
    "Kidney Papillary Cell Carcinoma": "Kidney - Cortex",
    "Kidney Chromophobe": "Kidney - Cortex",
    "Head & Neck Squamous Cell Carcinoma": "Minor Salivary Gland",
    "Cervical & Endocervical Cancer": "Cervix - Ectocervix",
    "Uterine Corpus Endometrioid Carcinoma": "Uterus",
    "Ovarian Serous Cystadenocarcinoma": "Ovary",
    "Skin Cutaneous Melanoma": "Skin - Sun Exposed (Lower leg)",
    "Glioblastoma Multiforme": "Brain - Cortex",
    "Brain Lower Grade Glioma": "Brain - Cortex",
}

# -----------------------------
# PREP DATA
# -----------------------------
tumor = df[df["type"] == "tumor"].copy()
tumor = tumor.dropna(subset=["cancer", "logexpr"]).copy()

# Ensure numeric and finite
tumor["logexpr"] = tumor["logexpr"].astype(float)
tumor = tumor[np.isfinite(tumor["logexpr"].to_numpy())].copy()

# Filter cancers by count
counts = tumor["cancer"].value_counts()
keep = counts[counts >= MIN_N].index
tumor = tumor[tumor["cancer"].isin(keep)].copy()

# Order by median on the plotted scale (logexpr scale)
order_full = (
    tumor.groupby("cancer")["logexpr"]
    .median()
    .sort_values(ascending=False)
    .index
    .tolist()
)

# Add acronym labels for plotting
tumor["cancer_acr"] = tumor["cancer"].map(lambda x: TCGA_ACRONYM.get(x, x))

# Build ordered acronym list in the same order
order_acr = [TCGA_ACRONYM.get(c, c) for c in order_full]

# -----------------------------
# PLOT
# -----------------------------
plt.figure(figsize=(14, 8))

ax = sns.violinplot(
    data=tumor,
    x="cancer_acr",
    y="logexpr",
    order=order_acr,
    cut=0,            # key: don't extend KDE beyond data range (fixes many "weird" tails)
    inner="quartile", # quartile lines are more informative than just a median
    linewidth=0.8,
    scale="width"     # comparable violin widths
)

ax.set_xlabel("")
ax.set_ylabel("SLFN11 expression (log2(TPM + 1e-3))")
ax.set_title("SLFN11 expression across TCGA cancers (tumor samples)")

plt.xticks(rotation=90, fontsize=X_FONTSIZE)
plt.tight_layout()
plt.savefig(OUTPNG, dpi=300)

print("Saved figure to", OUTPNG)

####################################
print("Generating figure: TCGA tumor vs TCGA normal")

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUTPNG = f"{OUTDIR}/slfn11_tcga_tumor_vs_normal.png"
MIN_PER_GROUP = 10
X_FONTSIZE = 7

# Keep only tumor and normal
df_tn = df[df["type"].isin(["tumor", "normal"])].copy()
df_tn = df_tn.dropna(subset=["cancer", "logexpr"]).copy()

# Count per cancer and per type
counts = (
    df_tn.groupby(["cancer", "type"])
    .size()
    .unstack(fill_value=0)
)

# Keep only cancers that have enough tumor and enough normal samples
valid_cancers = counts[
    (counts.get("tumor", 0) >= MIN_PER_GROUP) &
    (counts.get("normal", 0) >= MIN_PER_GROUP)
].index

df_tn = df_tn[df_tn["cancer"].isin(valid_cancers)].copy()

# Order by tumor-normal median difference
medians = (
    df_tn.groupby(["cancer", "type"])["logexpr"]
    .median()
    .unstack()
)

medians["delta"] = medians["tumor"] - medians["normal"]
order_full = medians.sort_values("delta", ascending=False).index.tolist()

df_tn["cancer_acr"] = df_tn["cancer"].map(lambda x: TCGA_ACRONYM.get(x, x))
order_acr = [TCGA_ACRONYM.get(c, c) for c in order_full]

plt.figure(figsize=(14, 8))

ax = sns.violinplot(
    data=df_tn,
    x="cancer_acr",
    y="logexpr",
    hue="type",
    order=order_acr,
    cut=0,
    inner="quartile",
    linewidth=0.8,
    scale="width"
)

ax.set_xlabel("")
ax.set_ylabel("SLFN11 expression (log2 scale)")
ax.set_title("SLFN11 expression in TCGA tumor vs normal samples")

plt.xticks(rotation=90, fontsize=X_FONTSIZE)
plt.tight_layout()
plt.savefig(OUTPNG, dpi=300)

print("Saved figure to", OUTPNG)

##############################

print("Generating figure: GTEx tissue expression")

OUTPNG = f"{OUTDIR}/slfn11_gtex_tissues.png"
MIN_N = 20
X_FONTSIZE = 7

gtex_plot = gtex_df.dropna(subset=["tissue", "logexpr"]).copy()
gtex_plot["logexpr"] = gtex_plot["logexpr"].astype(float)
gtex_plot = gtex_plot[np.isfinite(gtex_plot["logexpr"].to_numpy())].copy()

counts = gtex_plot["tissue"].value_counts()
keep = counts[counts >= MIN_N].index
gtex_plot = gtex_plot[gtex_plot["tissue"].isin(keep)].copy()

order_tissue = (
    gtex_plot.groupby("tissue")["logexpr"]
    .median()
    .sort_values(ascending=False)
    .index
    .tolist()
)

plt.figure(figsize=(16, 8))

ax = sns.violinplot(
    data=gtex_plot,
    x="tissue",
    y="logexpr",
    order=order_tissue,
    cut=0,
    inner="quartile",
    linewidth=0.8,
    scale="width"
)

ax.set_xlabel("")
ax.set_ylabel("SLFN11 expression (log2 scale)")
ax.set_title("SLFN11 expression across GTEx normal tissues")

plt.xticks(rotation=90, fontsize=X_FONTSIZE)
plt.tight_layout()
plt.savefig(OUTPNG, dpi=300)

print("Saved figure to", OUTPNG)

######################################

print("Generating figure: matched TCGA tumor / TCGA normal / GTEx normal")

OUTPNG = f"{OUTDIR}/slfn11_matched_tcga_gtex.png"
MIN_TCGA_TUMOR = 20
MIN_TCGA_NORMAL = 5
MIN_GTEX = 20
X_FONTSIZE = 8

rows = []

for cancer_name, gtex_tissue in TCGA_TO_GTEX_TISSUE.items():

    # TCGA tumor
    tmp_tumor = df[
        (df["cancer"] == cancer_name) &
        (df["type"] == "tumor")
    ].dropna(subset=["logexpr"]).copy()

    # TCGA normal
    tmp_normal = df[
        (df["cancer"] == cancer_name) &
        (df["type"] == "normal")
    ].dropna(subset=["logexpr"]).copy()

    # GTEx normal
    tmp_gtex = gtex_df[
        (gtex_df["tissue"] == gtex_tissue)
    ].dropna(subset=["logexpr"]).copy()

    if len(tmp_tumor) >= MIN_TCGA_TUMOR and len(tmp_normal) >= MIN_TCGA_NORMAL and len(tmp_gtex) >= MIN_GTEX:

        cancer_acr = TCGA_ACRONYM.get(cancer_name, cancer_name)

        rows.append(pd.DataFrame({
            "group": cancer_acr,
            "source": "TCGA tumor",
            "logexpr": tmp_tumor["logexpr"].astype(float).values
        }))

        rows.append(pd.DataFrame({
            "group": cancer_acr,
            "source": "TCGA normal",
            "logexpr": tmp_normal["logexpr"].astype(float).values
        }))

        rows.append(pd.DataFrame({
            "group": cancer_acr,
            "source": "GTEx normal",
            "logexpr": tmp_gtex["logexpr"].astype(float).values
        }))

matched_df = pd.concat(rows, ignore_index=True)
matched_df = matched_df[np.isfinite(matched_df["logexpr"].to_numpy())].copy()

# Order by tumor median
order_groups = (
    matched_df[matched_df["source"] == "TCGA tumor"]
    .groupby("group")["logexpr"]
    .median()
    .sort_values(ascending=False)
    .index
    .tolist()
)

plt.figure(figsize=(15, 8))

ax = sns.violinplot(
    data=matched_df,
    x="group",
    y="logexpr",
    hue="source",
    order=order_groups,
    cut=0,
    inner="quartile",
    linewidth=0.8,
    scale="width"
)

ax.set_xlabel("")
ax.set_ylabel("SLFN11 expression (log2 scale)")
ax.set_title("SLFN11 expression in matched TCGA tumor, TCGA normal, and GTEx normal tissues")

plt.xticks(rotation=90, fontsize=X_FONTSIZE)
plt.tight_layout()
plt.savefig(OUTPNG, dpi=300)

print("Saved figure to", OUTPNG)
