"""
parse_gse21501_survival.py
Stream-parse the real GSE21501 family SOFT file to extract:
  1. Per-sample clinical fields (os_time, os_event, risk_group) from
     !Sample_characteristics_ch2.
  2. Per-sample expression values (VALUE = log2 ratio, tumor/reference) for the
     hypoxia and acinar_identity marker genes only, via the GPL4133 platform's
     probe -> GENE_SYMBOL mapping.

The family SOFT file is ~188MB gzipped / ~6M lines; this streams it once with
a plain line-by-line parser (no GEOparse) and keeps only the handful of probe
rows and characteristic lines needed, so it runs in well under a minute.

Outputs: data/processed/GSE21501_clinical_expression.tsv
  (one row per sample with os_time, os_event, risk_group, and one column per
   hypoxia/acinar marker gene with its log2 ratio value, where available)
"""

import gzip
import os
import re
import sys
import yaml

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(BASE_DIR, "config")
RAW_FILE = os.path.join(BASE_DIR, "data", "raw", "validation_extra", "GSE21501_family.soft.gz")
OUT_FILE = os.path.join(BASE_DIR, "data", "processed", "GSE21501_clinical_expression.tsv")


def load_genes_of_interest():
    with open(os.path.join(CONFIG_DIR, "gene_sets.yml")) as f:
        gene_sets = yaml.safe_load(f)
    genes = set(gene_sets["hypoxia"]) | set(gene_sets["acinar_identity"])
    return {g.upper() for g in genes}


def parse():
    genes_of_interest = load_genes_of_interest()
    print(f"Looking for {len(genes_of_interest)} genes: {sorted(genes_of_interest)}")

    probe_to_gene = {}   # probe_id (str) -> gene_symbol (only for genes_of_interest)
    samples = {}         # geo_accession -> {"os_time":..., "os_event":..., "risk_group":..., genes: {gene: [values]}}

    in_platform_table = False
    platform_id_col = platform_sym_col = None

    cur_sample = None
    in_sample_table = False
    sample_id_col = sample_val_col = None

    n_lines = 0
    with gzip.open(RAW_FILE, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            n_lines += 1
            line = line.rstrip("\n")

            if line.startswith("^SAMPLE"):
                cur_sample = line.split("=", 1)[1].strip()
                samples[cur_sample] = {"os_time": None, "os_event": None,
                                        "risk_group": None, "genes": {}}
                continue

            if cur_sample is not None and line.startswith("!Sample_characteristics_ch2"):
                val = line.split("=", 1)[1].strip()
                low = val.lower()
                if low.startswith("os time:"):
                    v = val.split(":", 1)[1].strip()
                    samples[cur_sample]["os_time"] = v if v else None
                elif low.startswith("os event:"):
                    v = val.split(":", 1)[1].strip()
                    samples[cur_sample]["os_event"] = v if v else None
                elif low.startswith("risk group:"):
                    v = val.split(":", 1)[1].strip()
                    samples[cur_sample]["risk_group"] = v if v else None
                continue

            if line.startswith("!platform_table_begin"):
                in_platform_table = True
                continue
            if line.startswith("!platform_table_end"):
                in_platform_table = False
                continue
            if in_platform_table:
                if platform_id_col is None:
                    header = line.split("\t")
                    platform_id_col = header.index("ID")
                    platform_sym_col = header.index("GENE_SYMBOL")
                    continue
                cols = line.split("\t")
                if len(cols) <= max(platform_id_col, platform_sym_col):
                    continue
                sym = cols[platform_sym_col].strip().upper()
                if sym in genes_of_interest:
                    probe_to_gene[cols[platform_id_col]] = sym
                continue

            if line.startswith("!sample_table_begin"):
                in_sample_table = True
                sample_id_col = sample_val_col = None
                continue
            if line.startswith("!sample_table_end"):
                in_sample_table = False
                continue
            if in_sample_table and cur_sample is not None:
                if sample_id_col is None:
                    header = line.split("\t")
                    sample_id_col = header.index("ID_REF")
                    sample_val_col = header.index("VALUE")
                    continue
                cols = line.split("\t")
                if len(cols) <= max(sample_id_col, sample_val_col):
                    continue
                probe_id = cols[sample_id_col]
                gene = probe_to_gene.get(probe_id)
                if gene is None:
                    continue
                val_str = cols[sample_val_col].strip()
                if not val_str:
                    continue
                try:
                    val = float(val_str)
                except ValueError:
                    continue
                samples[cur_sample]["genes"].setdefault(gene, []).append(val)

            if n_lines % 1_000_000 == 0:
                print(f"  ...{n_lines:,} lines processed, {len(probe_to_gene)} probes of interest mapped")

    print(f"Total lines: {n_lines:,}; probes mapped to genes of interest: {len(probe_to_gene)}")
    print(f"Total samples parsed: {len(samples)}")
    return samples, sorted(genes_of_interest)


def main():
    samples, genes = parse()

    rows = []
    for acc, d in samples.items():
        row = {"geo_accession": acc, "os_time": d["os_time"], "os_event": d["os_event"],
               "risk_group": d["risk_group"]}
        for g in genes:
            vals = d["genes"].get(g)
            row[g] = sum(vals) / len(vals) if vals else None
        rows.append(row)

    import pandas as pd
    df = pd.DataFrame(rows)
    n_with_clinical = df["os_time"].notna().sum()
    n_with_event = df["os_event"].notna().sum()
    print(f"Samples with os_time: {n_with_clinical}; with os_event: {n_with_event}")
    for g in genes:
        print(f"  {g}: {df[g].notna().sum()}/{len(df)} samples with a value")

    os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
    df.to_csv(OUT_FILE, sep="\t", index=False)
    print(f"Saved: {OUT_FILE}")


if __name__ == "__main__":
    main()
