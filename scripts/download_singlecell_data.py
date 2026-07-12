"""
download_singlecell_data.py
Attempt to download scRNA-seq and proteomics datasets.
If download fails, creates flag files so downstream scripts use simulated data.
"""

import os
import sys
import datetime
import traceback

import numpy as np
import random

np.random.seed(1234)
random.seed(1234)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_SC_DIR = os.path.join(BASE_DIR, "data", "raw", "singlecell")
RAW_PROT_DIR = os.path.join(BASE_DIR, "data", "raw", "proteomics", "CPTAC_PDA")
RAW_VAL_DIR = os.path.join(BASE_DIR, "data", "raw", "validation_extra")

LOG_PATH = os.path.join(RAW_SC_DIR, "DOWNLOAD_LOG.md")

log_lines = []


def log(msg):
    print(msg)
    log_lines.append(msg)


def write_log():
    os.makedirs(RAW_SC_DIR, exist_ok=True)
    with open(LOG_PATH, "w") as f:
        f.write("# Phase 3 Download Log\n\n")
        f.write(f"Generated: {datetime.datetime.now().isoformat()}\n\n")
        for line in log_lines:
            f.write(line + "\n")


def set_failed_flag(dataset, reason):
    flag_path = os.path.join(RAW_SC_DIR, dataset, f"{dataset}_DOWNLOAD_FAILED.txt")
    os.makedirs(os.path.dirname(flag_path), exist_ok=True)
    with open(flag_path, "w") as f:
        f.write(f"Download failed for {dataset}\n")
        f.write(f"Reason: {reason}\n")
        f.write(f"Timestamp: {datetime.datetime.now().isoformat()}\n")
        f.write("Downstream scripts will use SIMULATED DATA.\n")
    log(f"  -> Flag file written: {flag_path}")


def has_existing_real_data(dataset, dest_dir):
    """Check if real data was already obtained for this dataset (e.g. via a
    manual/out-of-band download), so we don't re-attempt a network call that
    is known not to retrieve the count matrix via GEOparse and don't overwrite
    an accurate status with a false 'FAILED' flag. Mirrors the detection
    logic in preprocess_singlecell_cohorts.py."""
    if dataset == "GSE154778":
        return os.path.exists(os.path.join(dest_dir, "GSE154778_dgeMtx.csv.gz"))
    if dataset == "GSE202051":
        return os.path.exists(os.path.join(dest_dir, "GSE202051_adata_010nuc_10x.h5ad"))
    if dataset == "Peng_et_al":
        tar_path = os.path.join(dest_dir, "GSE155698_RAW.tar")
        return os.path.exists(tar_path) and os.path.getsize(tar_path) > 900_000_000
    return False


def attempt_geo_download(geo_id, dest_dir):
    """Attempt to download a GEO dataset using GEOparse."""
    try:
        import GEOparse
        log(f"  Attempting GEOparse download for {geo_id}...")
        os.makedirs(dest_dir, exist_ok=True)
        gse = GEOparse.get_GEO(geo=geo_id, destdir=dest_dir, silent=True)
        # Check if supplementary files were downloaded
        soft_files = [f for f in os.listdir(dest_dir) if f.endswith(".soft.gz") or f.endswith(".soft")]
        if soft_files:
            log(f"  SOFT file found: {soft_files[0]}")
            # Check for matrix/count files
            supp_files = [f for f in os.listdir(dest_dir) if
                          any(ext in f for ext in [".h5", ".h5ad", ".loom", "matrix.mtx", "barcodes.tsv", "genes.tsv", "features.tsv"])]
            if supp_files:
                log(f"  Supplementary count files found: {supp_files}")
                return True
            else:
                log(f"  No count matrix files found in {dest_dir}. Likely only SOFT file downloaded.")
                return False
        return False
    except Exception as e:
        log(f"  GEOparse download failed: {e}")
        return False


def attempt_cptac_download():
    """Verify CPTAC-PDA proteomic data is accessible via the official `cptac`
    package (the method actually used by analyze_cptac_protein.py). The old
    direct-URL approach below always failed (those endpoints don't serve raw
    CPTAC data) and left a stale DOWNLOAD_FAILED flag even after the real
    `cptac` package method succeeded downstream, which is what left this
    project's docs incorrectly claiming 'simulated' data."""
    os.makedirs(RAW_PROT_DIR, exist_ok=True)
    try:
        import cptac
        log("  Checking CPTAC-PDA access via the `cptac` Python package (cptac.Pdac())...")
        pdac = cptac.Pdac()
        prot = pdac.get_proteomics(source="umich")
        if prot is not None and len(prot) > 0:
            log(f"  cptac package access succeeded: {len(prot)} umich proteomics samples available.")
            return True
        log("  cptac package returned no proteomics rows.")
        return False
    except Exception as e:
        log(f"  cptac package access failed: {e}")
        return False


def attempt_validation_cohort_download():
    """Attempt to download an additional survival cohort (GSE21501)."""
    os.makedirs(RAW_VAL_DIR, exist_ok=True)
    try:
        import GEOparse
        geo_id = "GSE21501"
        log(f"  Attempting download of additional survival cohort {geo_id}...")
        gse = GEOparse.get_GEO(geo=geo_id, destdir=RAW_VAL_DIR, silent=True)
        files = os.listdir(RAW_VAL_DIR)
        soft_files = [f for f in files if f.endswith(".soft.gz") or f.endswith(".soft")]
        if soft_files:
            log(f"  {geo_id} SOFT file found: {soft_files}")
            return True
        return False
    except Exception as e:
        log(f"  GSE21501 download failed: {e}")
        return False


def main():
    log("# Phase 3 Single-Cell Data Download Script")
    log(f"## Started: {datetime.datetime.now().isoformat()}")
    log("")

    # Download GSE154778
    log("## Dataset 1: GSE154778 (Werba et al. PDAC scRNA-seq)")
    dest_154778 = os.path.join(RAW_SC_DIR, "GSE154778")
    if has_existing_real_data("GSE154778", dest_154778):
        log("  Real count matrix already present on disk (GSE154778_dgeMtx.csv.gz) — skipping re-download.")
        success_154778 = True
        log("  STATUS: ALREADY PRESENT (real)")
    else:
        success_154778 = attempt_geo_download("GSE154778", dest_154778)
        if not success_154778:
            log("  STATUS: DOWNLOAD FAILED — count matrix not retrieved")
            set_failed_flag("GSE154778", "Count matrix files not found via GEOparse; dataset may require manual download from GEO FTP")
        else:
            log("  STATUS: SUCCESS")
    log("")

    # Download GSE202051
    log("## Dataset 2: GSE202051 (PDAC scRNA-seq cohort 2)")
    dest_202051 = os.path.join(RAW_SC_DIR, "GSE202051")
    if has_existing_real_data("GSE202051", dest_202051):
        log("  Real h5ad already present on disk (GSE202051_adata_010nuc_10x.h5ad) — skipping re-download.")
        success_202051 = True
        log("  STATUS: ALREADY PRESENT (real)")
    else:
        success_202051 = attempt_geo_download("GSE202051", dest_202051)
        if not success_202051:
            log("  STATUS: DOWNLOAD FAILED — count matrix not retrieved")
            set_failed_flag("GSE202051", "Count matrix files not found via GEOparse; dataset may require manual download from GEO FTP")
        else:
            log("  STATUS: SUCCESS")
    log("")

    # Download Peng et al (GSE155698)
    log("## Dataset 3: Peng et al. (GSE155698) PDAC scRNA-seq")
    dest_peng = os.path.join(RAW_SC_DIR, "Peng_et_al")
    if has_existing_real_data("Peng_et_al", dest_peng):
        log("  Real RAW.tar already present on disk (GSE155698_RAW.tar, >900MB) — skipping re-download.")
        success_peng = True
        log("  STATUS: ALREADY PRESENT (real)")
    else:
        success_peng = attempt_geo_download("GSE155698", dest_peng)
        if not success_peng:
            log("  STATUS: DOWNLOAD FAILED — count matrix not retrieved")
            set_failed_flag("Peng_et_al", "Count matrix files not found via GEOparse; dataset may require manual download")
        else:
            log("  STATUS: SUCCESS")
    log("")

    # Download CPTAC-PDA
    log("## Dataset 4: CPTAC-PDA Proteomics")
    success_cptac = attempt_cptac_download()
    flag_path = os.path.join(RAW_PROT_DIR, "CPTAC_PDA_DOWNLOAD_FAILED.txt")
    if not success_cptac:
        log("  STATUS: DOWNLOAD FAILED")
        with open(flag_path, "w") as f:
            f.write("CPTAC-PDA proteomics access via the `cptac` package failed.\n")
            f.write("Downstream scripts will use SIMULATED DATA.\n")
            f.write(f"Timestamp: {datetime.datetime.now().isoformat()}\n")
        log(f"  -> Flag file written: {flag_path}")
    else:
        log("  STATUS: SUCCESS (real data via cptac package)")
        if os.path.exists(flag_path):
            os.remove(flag_path)
            log(f"  -> Removed stale flag file: {flag_path}")
    log("")

    # Download additional survival cohort
    log("## Dataset 5: GSE21501 (additional survival cohort)")
    success_val = attempt_validation_cohort_download()
    if not success_val:
        log("  STATUS: DOWNLOAD FAILED")
        flag_path = os.path.join(RAW_VAL_DIR, "GSE21501_DOWNLOAD_FAILED.txt")
        with open(flag_path, "w") as f:
            f.write("GSE21501 download failed.\n")
            f.write("Downstream scripts will use SIMULATED DATA for expanded survival analysis.\n")
            f.write(f"Timestamp: {datetime.datetime.now().isoformat()}\n")
        log(f"  -> Flag file written: {flag_path}")
    else:
        log("  STATUS: SUCCESS")
    log("")

    log("## Summary")
    log(f"  GSE154778: {'SUCCESS' if success_154778 else 'FAILED (simulated)'}")
    log(f"  GSE202051: {'SUCCESS' if success_202051 else 'FAILED (simulated)'}")
    log(f"  Peng_et_al: {'SUCCESS' if success_peng else 'FAILED (simulated)'}")
    log(f"  CPTAC-PDA: {'SUCCESS' if success_cptac else 'FAILED (simulated)'}")
    log(f"  GSE21501: {'SUCCESS' if success_val else 'FAILED (simulated)'}")
    log(f"## Completed: {datetime.datetime.now().isoformat()}")

    write_log()
    print(f"\nDownload log written to: {LOG_PATH}")


if __name__ == "__main__":
    main()
