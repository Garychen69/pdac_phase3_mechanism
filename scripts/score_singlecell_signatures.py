"""
score_singlecell_signatures.py
Score each cell for all 7 signatures using sc.tl.score_genes.
Saves per-cell scores TSV with required columns.
"""

import os
import sys
import numpy as np
import pandas as pd
import random
import yaml
import warnings
warnings.filterwarnings("ignore")

np.random.seed(1234)
random.seed(1234)

import scanpy as sc

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(BASE_DIR, "config")
sc.settings.verbosity = 0

SIGNATURES = [
    "hypoxia",
    "acinar_identity",
    "lipid_synthesis_srebp",
    "desaturation_elongation",
    "fatty_acid_uptake_oxidation",
    "caf",
    "emt",
]


def load_configs():
    with open(os.path.join(CONFIG_DIR, "singlecell_cohorts.yml")) as f:
        cohort_cfg = yaml.safe_load(f)
    with open(os.path.join(CONFIG_DIR, "gene_sets.yml")) as f:
        gene_sets = yaml.safe_load(f)
    return cohort_cfg, gene_sets


def score_cohort(cohort, gene_sets):
    name = cohort["name"]
    processed_file = os.path.join(BASE_DIR, cohort["processed_file"])
    annotations_file = os.path.join(BASE_DIR, cohort["annotations_file"])
    scores_file = os.path.join(BASE_DIR, cohort["scores_file"])

    if not os.path.exists(processed_file):
        print(f"  Preprocessed file not found: {processed_file}. Skipping.")
        return

    print(f"  Loading {name}...")
    adata = sc.read_h5ad(processed_file)

    # Load cell type annotations if available
    if os.path.exists(annotations_file):
        ann_df = pd.read_csv(annotations_file, sep="\t")
        ann_df = ann_df.set_index("cell_id")
        # Add cell_type to adata.obs
        common = adata.obs.index.intersection(ann_df.index)
        if "cell_type" not in adata.obs.columns:
            adata.obs["cell_type"] = "unknown"
        adata.obs.loc[common, "cell_type"] = ann_df.loc[common, "cell_type"]

    # Score each signature
    score_cols = {}
    for sig_name in SIGNATURES:
        gene_list = gene_sets.get(sig_name, [])
        available = [g for g in gene_list if g in adata.var_names]
        score_key = f"{sig_name}_score"

        if len(available) >= 1:
            np.random.seed(1234)
            sc.tl.score_genes(adata, available, score_name=score_key, random_state=1234)
            score_cols[score_key] = adata.obs[score_key].values
            print(f"    {sig_name}: scored {len(available)}/{len(gene_list)} genes")
        else:
            # All genes missing: assign zero score with warning
            print(f"    WARNING: no genes from {sig_name} found in {name}. Assigning zero score.")
            score_cols[score_key] = np.zeros(adata.n_obs)

    # Build output dataframe
    cell_id = adata.obs.index.values
    patient_id = adata.obs.get("patient_id", pd.Series(["unknown"] * adata.n_obs, index=adata.obs.index)).values
    cohort_col = np.array([name] * adata.n_obs)
    cell_type = adata.obs.get("cell_type", pd.Series(["unknown"] * adata.n_obs, index=adata.obs.index)).values
    is_simulated = np.array([adata.uns.get("simulated", False)] * adata.n_obs)

    out_df = pd.DataFrame({
        "cell_id": cell_id,
        "patient_id": patient_id,
        "cohort": cohort_col,
        "cell_type": cell_type,
        "hypoxia_score": score_cols["hypoxia_score"],
        "acinar_identity_score": score_cols["acinar_identity_score"],
        "lipid_synthesis_srebp_score": score_cols["lipid_synthesis_srebp_score"],
        "desaturation_elongation_score": score_cols["desaturation_elongation_score"],
        "fatty_acid_uptake_oxidation_score": score_cols["fatty_acid_uptake_oxidation_score"],
        "caf_score": score_cols["caf_score"],
        "emt_score": score_cols["emt_score"],
        "is_simulated": is_simulated,
    })

    os.makedirs(os.path.dirname(scores_file), exist_ok=True)
    out_df.to_csv(scores_file, sep="\t", index=False)
    print(f"  Scores saved: {scores_file}")

    # Also save scores back to adata
    for sk, vals in score_cols.items():
        adata.obs[sk] = vals
    adata.write_h5ad(processed_file)

    return out_df


def main():
    print("=== Scoring Single-Cell Signatures ===\n")
    cohort_cfg, gene_sets = load_configs()

    for cohort in cohort_cfg["cohorts"]:
        print(f"--- Scoring: {cohort['name']} ---")
        try:
            score_cohort(cohort, gene_sets)
        except Exception as e:
            import traceback
            print(f"  ERROR in {cohort['name']}: {e}")
            traceback.print_exc()
        print()

    print("=== Scoring Complete ===")


if __name__ == "__main__":
    main()
