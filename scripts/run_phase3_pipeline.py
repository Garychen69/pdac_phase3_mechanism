"""
run_phase3_pipeline.py
Master runner for the Phase 3 PDAC analysis pipeline.
Runs all scripts in order, catching errors per step.
"""

import os
import sys
import subprocess
import time
import datetime
import numpy as np
import random

np.random.seed(1234)
random.seed(1234)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")
PYTHON = sys.executable

PIPELINE_STEPS = [
    ("download_singlecell_data.py", "Download scRNA-seq and CPTAC data"),
    ("preprocess_singlecell_cohorts.py", "Preprocess single-cell cohorts (or simulate)"),
    ("annotate_cell_types.py", "Annotate cell types via Leiden clustering"),
    ("score_singlecell_signatures.py", "Score gene signatures per cell"),
    ("analyze_hypoxia_acinar_cooccurrence.py", "Analyze hypoxia/acinar co-occurrence"),
    ("analyze_lipid_cell_of_origin.py", "Analyze lipid rewiring cell-of-origin"),
    ("analyze_caf_emt_cell_of_origin.py", "Analyze CAF/EMT cell-of-origin"),
    ("deconvolve_bulk_purity.py", "Purity deconvolution and adjustment"),
    ("analyze_cptac_protein.py", "CPTAC protein-level validation"),
    ("parse_gse21501_survival.py", "Parse real GSE21501 SOFT file (clinical + marker genes)"),
    ("run_expanded_survival_meta.py", "Expanded survival meta-analysis"),
    ("generate_phase3_report.py", "Generate final report"),
]


def run_step(script_name, description, step_num, total_steps):
    script_path = os.path.join(SCRIPTS_DIR, script_name)
    print(f"\n{'='*60}")
    print(f"Step {step_num}/{total_steps}: {description}")
    print(f"Script: {script_name}")
    print(f"{'='*60}")
    t_start = time.time()

    if not os.path.exists(script_path):
        print(f"ERROR: Script not found: {script_path}")
        return False, "Script not found"

    try:
        result = subprocess.run(
            [PYTHON, script_path],
            capture_output=False,
            text=True,
            cwd=BASE_DIR,
            timeout=600,
        )
        elapsed = time.time() - t_start
        if result.returncode == 0:
            print(f"\nStep {step_num} SUCCEEDED ({elapsed:.1f}s)")
            return True, None
        else:
            print(f"\nStep {step_num} FAILED (exit code {result.returncode}, {elapsed:.1f}s)")
            return False, f"Exit code {result.returncode}"
    except subprocess.TimeoutExpired:
        print(f"\nStep {step_num} TIMED OUT (>600s)")
        return False, "Timeout"
    except Exception as e:
        print(f"\nStep {step_num} ERROR: {e}")
        return False, str(e)


def main():
    print(f"{'='*60}")
    print(f"PDAC Phase 3 Pipeline — Master Runner")
    print(f"Started: {datetime.datetime.now().isoformat()}")
    print(f"Python: {PYTHON}")
    print(f"Base dir: {BASE_DIR}")
    print(f"{'='*60}")

    total = len(PIPELINE_STEPS)
    results = []

    for i, (script, desc) in enumerate(PIPELINE_STEPS, 1):
        success, error = run_step(script, desc, i, total)
        results.append((script, desc, success, error))

    # Summary
    print(f"\n{'='*60}")
    print(f"PIPELINE SUMMARY")
    print(f"{'='*60}")
    n_success = sum(1 for _, _, s, _ in results if s)
    n_failed = total - n_success
    for script, desc, success, error in results:
        status = "OK" if success else f"FAILED ({error})"
        print(f"  {'[OK]  ' if success else '[FAIL]'} {script}: {status}")

    print(f"\n{n_success}/{total} steps succeeded, {n_failed} failed.")
    print(f"Completed: {datetime.datetime.now().isoformat()}")

    if n_failed > 0:
        print("\nSome steps failed. Check output above for details.")
        sys.exit(1)
    else:
        print("\nAll steps completed successfully.")


if __name__ == "__main__":
    main()
