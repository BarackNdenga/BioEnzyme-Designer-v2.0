#!/usr/bin/env python3
"""
BioEnzyme Designer v2.0 — Test suite.

Run: python -m pytest tests/ -v
Or:  python tests/test_all.py
"""

import sys
import os
from pathlib import Path

# When running as installed package, don't add project root to avoid conflicts
# Uncomment the lines below only if running without pip install -e .
# _project_root = str(Path(__file__).resolve().parent.parent)
# if _project_root not in sys.path:
#     sys.path.insert(0, _project_root)

import numpy as np
import pandas as pd


def test_constants():
    """Test amino acid constants are complete."""
    from core.constants import (
        AA_HYDROPHOBICITY, AA_CHARGE, AA_VOLUME,
        ALL_AA, ONE_TO_THREE, THREE_TO_ONE,
    )
    assert len(AA_HYDROPHOBICITY) == 20
    assert len(AA_CHARGE) == 20
    assert len(AA_VOLUME) == 20
    assert len(ALL_AA) == 20
    assert len(ONE_TO_THREE) == 20
    assert len(THREE_TO_ONE) == 20

    # Test round-trip
    for one, three in ONE_TO_THREE.items():
        assert THREE_TO_ONE[three] == one

    print("[PASS] test_constants")


def test_ml_model_synthetic():
    """Test synthetic dataset generation and model training."""
    from core.ml_model import generate_synthetic_dataset, train_activity_model, predict_activity, FEATURE_COLS

    df = generate_synthetic_dataset(n_samples=200)
    assert len(df) == 200
    assert "log_kcat_km" in df.columns
    for col in FEATURE_COLS:
        assert col in df.columns

    model, scaler = train_activity_model(df)

    # Predict
    feat = {col: 0.0 for col in FEATURE_COLS}
    pred = predict_activity(model, scaler, feat)
    assert isinstance(pred, float)
    assert not np.isnan(pred)

    print("[PASS] test_ml_model_synthetic")


def test_ml_model_brenda_fallback():
    """Test BRENDA fallback data generation."""
    from core.ml_model import _generate_brenda_like_data, FEATURE_COLS

    records = _generate_brenda_like_data(50)
    assert len(records) == 50
    for r in records:
        assert "log_kcat_km" in r
        assert "source" in r
        for col in FEATURE_COLS:
            assert col in r

    print("[PASS] test_ml_model_brenda_fallback")


def test_ml_model_sabio_fallback():
    """Test SABIO-RK fallback data generation."""
    from core.ml_model import _generate_sabio_like_data, FEATURE_COLS

    records = _generate_sabio_like_data(50)
    assert len(records) == 50
    for r in records:
        assert "log_kcat_km" in r
        assert "source" in r
        for col in FEATURE_COLS:
            assert col in r

    print("[PASS] test_ml_model_sabio_fallback")


def test_combined_dataset():
    """Test combined dataset building."""
    from core.ml_model import build_training_dataset, FEATURE_COLS

    df = build_training_dataset(data_source="combined")
    assert len(df) > 500  # synthetic + brenda + sabio fallbacks
    for col in FEATURE_COLS + ["log_kcat_km"]:
        assert col in df.columns

    print("[PASS] test_combined_dataset")


def test_scoring():
    """Test mutation scoring functions."""
    from core.scoring import score_mutation, suggest_mutations
    from core.constants import ALL_AA

    res_row = {
        "is_catalytic": 1,
        "approx_rsa": 0.2,
        "chain_id": "A",
    }

    # Test thermostability scoring
    s1 = score_mutation("ALA", "LEU", res_row, "thermostability")
    s2 = score_mutation("ALA", "ALA", res_row, "thermostability")
    assert s1 > s2  # LEU should score better for thermostability

    # Test activity scoring
    s3 = score_mutation("ALA", "SER", res_row, "activity")
    assert s3 > -999  # Not a self-mutation

    # Test self-mutation penalty
    s4 = score_mutation("ALA", "ALA", res_row, "activity")
    assert s4 == -999

    print("[PASS] test_scoring")


def test_ddg_validation():
    """Test knowledge-based ΔΔG estimation."""
    from bioenzyme_v2.validation.ddg_validation import estimate_foldx_like_ddg, estimate_rosetta_like_ddg

    mutations = [
        {"mutation": "ALA10LEU", "res_name": "ALA", "mut_aa": "LEU",
         "is_interface": False, "chain_id": "A", "res_id": 10},
        {"mutation": "SER20GLY", "res_name": "SER", "mut_aa": "GLY",
         "is_interface": True, "chain_id": "A", "res_id": 20},
    ]

    foldx_results = estimate_foldx_like_ddg(mutations)
    assert len(foldx_results) == 2
    for r in foldx_results:
        assert "ddg" in r
        assert "method" in r
        assert r["unit"] == "kcal/mol"

    rosetta_results = estimate_rosetta_like_ddg(mutations)
    assert len(rosetta_results) == 2
    for r in rosetta_results:
        assert "ddg" in r
        assert r["unit"] == "REU"

    print("[PASS] test_ddg_validation")


def test_comparison():
    """Test score-ΔΔG comparison."""
    from validation.ddg_validation import compare_scores_with_ddg

    mutations = [
        {"mutation": "ALA10LEU", "score": 5.0, "chain_id": "A", "res_id": 10},
        {"mutation": "SER20GLY", "score": -2.0, "chain_id": "A", "res_id": 20},
        {"mutation": "VAL30PHE", "score": 3.0, "chain_id": "A", "res_id": 30},
    ]

    validation_results = {
        "foldx_ddg": [
            {"mutation": "ALA10LEU", "ddg": -1.5},
            {"mutation": "SER20GLY", "ddg": 2.0},
            {"mutation": "VAL30PHE", "ddg": -0.5},
        ],
        "rosetta_ddg": [],
    }

    result = compare_scores_with_ddg(mutations, validation_results)
    assert "foldx_correlation" in result
    assert "agreement" in result
    assert len(result["agreement"]) == 3

    print("[PASS] test_comparison")


def test_analysis_with_mock_pdb():
    """Test full analysis pipeline with a mock PDB structure."""
    from Bio.PDB import PDBIO, StructureBuilder
    import tempfile

    from core.structure import load_structure, get_chains_info
    from core.analysis import (
        get_residues_all_chains, identify_catalytic_residues,
        compute_residue_features, extract_global_features,
        compute_interchain_contacts,
    )

    # Create a minimal PDB file with two chains
    pdb_content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 20.00           C
ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00 20.00           C
ATOM      4  O   ALA A   1       4.000   5.000   6.000  1.00 20.00           O
ATOM      5  CB  ALA A   1       2.500   3.500   4.500  1.00 20.00           C
ATOM      6  N   HIS A   2       5.000   6.000   7.000  1.00 25.00           N
ATOM      7  CA  HIS A   2       6.000   7.000   8.000  1.00 25.00           C
ATOM      8  C   HIS A   2       7.000   8.000   9.000  1.00 25.00           C
ATOM      9  O   HIS A   2       8.000   9.000  10.000  1.00 25.00           O
ATOM     10  CB  HIS A   2       6.500   7.500   8.500  1.00 25.00           C
ATOM     11  N   SER A   3       9.000  10.000  11.000  1.00 30.00           N
ATOM     12  CA  SER A   3      10.000  11.000  12.000  1.00 30.00           C
ATOM     13  C   SER A   3      11.000  12.000  13.000  1.00 30.00           C
ATOM     14  O   SER A   3      12.000  13.000  14.000  1.00 30.00           O
ATOM     15  CB  SER A   3      10.500  11.500  12.500  1.00 30.00           C
ATOM     16  N   ALA B   1       1.000   2.000   3.000  1.00 22.00           N
ATOM     17  CA  ALA B   1       2.000   3.000   4.000  1.00 22.00           C
ATOM     18  C   ALA B   1       3.000   4.000   5.000  1.00 22.00           C
ATOM     19  O   ALA B   1       4.000   5.000   6.000  1.00 22.00           O
ATOM     20  CB  ALA B   1       2.500   3.500   4.500  1.00 22.00           C
ATOM     21  N   LEU B   2       5.000   6.000   7.000  1.00 28.00           N
ATOM     22  CA  LEU B   2       6.000   7.000   8.000  1.00 28.00           C
ATOM     23  C   LEU B   2       7.000   8.000   9.000  1.00 28.00           C
ATOM     24  O   LEU B   2       8.000   9.000  10.000  1.00 28.00           O
ATOM     25  CB  LEU B   2       6.500   7.500   8.500  1.00 28.00           C
END
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(pdb_content)
        pdb_path = Path(f.name)

    try:
        # Load structure
        structure, model = load_structure(pdb_path)
        chains_info = get_chains_info(model)
        assert len(chains_info) == 2
        assert chains_info[0]["chain_id"] == "A"
        assert chains_info[1]["chain_id"] == "B"

        # Get all residues
        residues = get_residues_all_chains(model)
        assert len(residues) == 5  # 3 in A + 2 in B

        # Catalytic identification
        cat_ids = identify_catalytic_residues(residues)
        # HIS at position 2 in chain A should be catalytic
        assert ("A", 2) in cat_ids

        # Feature computation
        df = compute_residue_features(residues, cat_ids)
        assert len(df) == 5
        assert "chain_id" in df.columns
        assert "is_catalytic" in df.columns

        # Global features
        feats = extract_global_features(df)
        assert feats["n_chains"] == 2

        # Inter-chain contacts
        contacts = compute_interchain_contacts(model, cutoff=15.0)
        assert isinstance(contacts, list)

        print("[PASS] test_analysis_with_mock_pdb")
    finally:
        pdb_path.unlink()


def test_report_generation():
    """Test report generation."""
    from core.report import write_report
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)
        pdb_path = Path("mock.pdb")

        # Create mock data
        df = pd.DataFrame({
            "chain_id": ["A", "A", "A", "B", "B"],
            "res_id": [1, 2, 3, 1, 2],
            "res_name": ["ALA", "HIS", "SER", "ALA", "LEU"],
            "hydrophobicity": [1.8, -3.2, -0.8, 1.8, 3.8],
            "charge": [0.0, 0.1, 0.0, 0.0, 0.0],
            "volume": [88.6, 153.2, 89.0, 88.6, 166.7],
            "b_factor": [20.0, 25.0, 30.0, 22.0, 28.0],
            "approx_rsa": [0.25, 0.31, 0.375, 0.275, 0.35],
            "is_catalytic": [0, 1, 0, 0, 0],
        })

        catalytic_ids = {("A", 2)}
        baseline_feats = {
            "mean_hydrophobicity": 1.2,
            "mean_charge": 0.02,
            "mean_b_factor": 25.0,
            "global_volume": 117.2,
            "active_site_count": 1,
            "n_chains": 2,
        }

        mutations = [
            {"chain_id": "A", "res_id": 1, "res_name": "ALA",
             "mutation": "A:ALA1LEU", "score": 5.0,
             "is_catalytic": False, "is_interface": True},
            {"chain_id": "B", "res_id": 1, "res_name": "ALA",
             "mutation": "B:ALA1VAL", "score": 3.5,
             "is_catalytic": False, "is_interface": False},
        ]

        paths = {
            "viewer": str(output_dir / "view.html"),
            "properties": str(output_dir / "props.html"),
            "png": str(output_dir / "scores.png"),
            "report": str(output_dir / "report.txt"),
        }

        chains_info = [
            {"chain_id": "A", "n_residues": 3, "residue_ids": [1, 2, 3]},
            {"chain_id": "B", "n_residues": 2, "residue_ids": [1, 2]},
        ]

        report_path = write_report(
            output_dir, pdb_path, "activity",
            df, catalytic_ids, baseline_feats, 5.5,
            mutations, paths,
            validation_results={"correlation": 0.85},
            chains_info=chains_info,
        )

        assert report_path.exists()
        content = report_path.read_text()
        assert "Chains analysed" in content
        assert "Chain Information" in content
        assert "A:ALA1LEU" in content
        assert "Chain A" in content
        assert "Chain B" in content

    print("[PASS] test_report_generation")


def run_all_tests():
    """Run all tests."""
    print("=" * 60)
    print("  BioEnzyme Designer v2.0 — Test Suite")
    print("=" * 60)
    print()

    tests = [
        test_constants,
        test_ml_model_synthetic,
        test_ml_model_brenda_fallback,
        test_ml_model_sabio_fallback,
        test_combined_dataset,
        test_scoring,
        test_ddg_validation,
        test_comparison,
        test_analysis_with_mock_pdb,
        test_report_generation,
    ]

    passed = 0
    failed = 0
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"[FAIL] {test.__name__}: {e}")
            import traceback
            traceback.print_exc()
            failed += 1

    print()
    print("=" * 60)
    print(f"  Results: {passed} passed, {failed} failed, {passed + failed} total")
    print("=" * 60)

    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
