#!/usr/bin/env python3
"""
BioEnzyme Designer v2.0 — Main CLI entry point.

Usage:
    # Basic usage (same as v1):
    python -m bioenzyme_v2 --pdb enzyme.pdb --improve activity

    # Multi-chain analysis:
    python -m bioenzyme_v2 --pdb enzyme.pdb --improve thermostability \
        --chains all --mutations 15

    # Specific chain:
    python -m bioenzyme_v2 --pdb enzyme.pdb --improve pH_optimum \
        --chains A,B --mutations 10

    # With BRENDA/SABIO-RK data:
    python -m bioenzyme_v2 --pdb enzyme.pdb --improve activity \
        --data combined --ec 1.1.1.1 --organism "Homo sapiens"

    # With FoldX/Rosetta validation:
    python -m bioenzyme_v2 --pdb enzyme.pdb --improve activity \
        --validate foldx --validate rosetta --foldx-bin /path/to/foldx

    # Launch GUI:
    python -m bioenzyme_v2 --gui

    # Launch API server:
    python -m bioenzyme_v2 --api --port 8000
"""

import argparse
import os
import sys
from pathlib import Path

# Ensure the package root is in the path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def parse_args():
    parser = argparse.ArgumentParser(
        description="BioEnzyme Designer v2.0 — enzyme engineering assistant",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m bioenzyme_v2 --pdb enzyme.pdb --improve activity
  python -m bioenzyme_v2 --seq MKVL... --improve thermostability --mutations 10
  python -m bioenzyme_v2 --pdb enzyme.pdb --improve activity --chains all --validate foldx
  python -m bioenzyme_v2 --pdb enzyme.pdb --data brenda --ec 1.1.1.1
  python -m bioenzyme_v2 --gui
  python -m bioenzyme_v2 --api --port 8000
        """,
    )

    # ── Structure input ───────────────────────────────────────────────────
    src = parser.add_mutually_exclusive_group()
    src.add_argument("--pdb", metavar="FILE", help="Path to a PDB file.")
    src.add_argument("--seq", metavar="SEQUENCE",
                     help="Amino acid sequence (single-letter code). "
                          "Structure fetched from AlphaFold/RCSB.")

    # ── GUI / API mode ────────────────────────────────────────────────────
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--gui", action="store_true",
                      help="Launch the Tkinter graphical interface.")
    mode.add_argument("--api", action="store_true",
                      help="Launch the REST API server.")

    # ── Analysis parameters ───────────────────────────────────────────────
    parser.add_argument("--improve",
                        choices=["activity", "thermostability", "pH_optimum"],
                        default="activity",
                        help="Property to optimise (default: activity).")
    parser.add_argument("--output", default="./output", metavar="DIR",
                        help="Output folder (default: ./output).")
    parser.add_argument("--mutations", type=int, default=10, metavar="N",
                        help="Number of top mutations to suggest (default: 10).")

    # ── Multi-chain support ───────────────────────────────────────────────
    parser.add_argument("--chains", default="all", metavar="CHAIN_IDS",
                        help="Chains to analyse: 'all' for every chain, or "
                             "comma-separated list (e.g. 'A,B'). Default: all")

    # ── Data source ───────────────────────────────────────────────────────
    parser.add_argument("--data", default="synthetic",
                        choices=["synthetic", "brenda", "sabio_rk", "combined"],
                        help="ML training data source (default: synthetic).")
    parser.add_argument("--ec", default=None,
                        help="EC number filter for BRENDA/SABIO-RK data.")
    parser.add_argument("--organism", default=None,
                        help="Organism filter for BRENDA/SABIO-RK data.")

    # ── Validation ────────────────────────────────────────────────────────
    parser.add_argument("--validate", action="append", default=[],
                        choices=["foldx", "rosetta"],
                        help="Run ΔΔG validation (can be specified multiple times).")
    parser.add_argument("--foldx-bin", default=None,
                        help="Path to FoldX binary.")
    parser.add_argument("--rosetta-dir", default=None,
                        help="Path to Rosetta installation directory.")

    # ── API server options ────────────────────────────────────────────────
    parser.add_argument("--port", type=int, default=8000,
                        help="Port for API server (default: 8000).")
    parser.add_argument("--host", default="0.0.0.0",
                        help="Host for API server (default: 0.0.0.0).")

    return parser.parse_args()


def run_gui():
    """Launch the Tkinter GUI."""
    from gui.tkinter_gui import launch_gui
    print("[INFO] Launching BioEnzyme Designer GUI...")
    launch_gui()


def run_api(args):
    """Launch the REST API server."""
    from api.rest_api import run_server
    print(f"[INFO] Starting BioEnzyme Designer API on {args.host}:{args.port}")
    print(f"       API docs: http://{args.host}:{args.port}/docs")
    run_server(host=args.host, port=args.port)


def run_analysis(args):
    """Run the full analysis pipeline."""
    from core.structure import (
        load_structure, fetch_alphafold_pdb, get_chains_info,
    )
    from core.analysis import (
        get_residues_all_chains, get_residues_single_chain,
        identify_catalytic_residues, compute_residue_features,
        extract_global_features, compute_interchain_contacts,
    )
    from core.ml_model import (
        build_training_dataset, train_activity_model, predict_activity,
    )
    from core.scoring import suggest_mutations
    from core.visualization import (
        make_3d_html, make_residue_properties_html, make_mutation_score_png,
    )
    from core.report import write_report
    from validation.ddg_validation import run_full_validation

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Obtain structure ──────────────────────────────────────────
    try:
        if args.pdb:
            pdb_path = Path(args.pdb)
            if not pdb_path.exists():
                print(f"[ERROR] PDB file not found: {pdb_path}", file=sys.stderr)
                sys.exit(1)
        else:
            pdb_path = fetch_alphafold_pdb(args.seq, output_dir)
    except RuntimeError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

    # ── Step 2: Load structure ────────────────────────────────────────────
    try:
        structure, model = load_structure(pdb_path)
        chains_info = get_chains_info(model)
    except Exception as e:
        print(f"[ERROR] Failed to parse PDB: {e}", file=sys.stderr)
        sys.exit(1)

    # ── Step 3: Residue analysis (multi-chain aware) ──────────────────────
    if args.chains == "all":
        residues_with_chain = get_residues_all_chains(model)
    else:
        # Parse comma-separated chain IDs
        chain_ids = [c.strip() for c in args.chains.split(",")]
        residues_with_chain = []
        for cid in chain_ids:
            residues_with_chain.extend(
                get_residues_single_chain(model, chain_id=cid)
            )

    if not residues_with_chain:
        print("[ERROR] No standard amino acid residues found.", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Analysing {len(residues_with_chain)} residues "
          f"across {model[0] and len(list(model.get_chains()))} chain(s) ...")

    catalytic_ids = identify_catalytic_residues(residues_with_chain)
    print(f"[INFO] Identified {len(catalytic_ids)} putative catalytic residue(s).")

    df_residues = compute_residue_features(residues_with_chain, catalytic_ids)
    baseline_feats = extract_global_features(df_residues)

    # ── Step 4: ML model ──────────────────────────────────────────────────
    print(f"\n[INFO] Building training dataset (source: {args.data})...")
    train_df = build_training_dataset(
        data_source=args.data,
        ec_number=args.ec,
        organism=args.organism,
    )
    model_rf, scaler = train_activity_model(train_df)
    baseline_activity = predict_activity(model_rf, scaler, baseline_feats)
    print(f"[INFO] Baseline predicted log10(kcat/Km) = {baseline_activity:.3f} "
          f"≈ {10**baseline_activity:.2e} M⁻¹s⁻¹")

    # ── Step 5: Inter-chain contacts ──────────────────────────────────────
    interchain_contacts = compute_interchain_contacts(model)
    print(f"[INFO] Found {len(interchain_contacts)} inter-chain contacts "
          f"(cutoff: 8.0 Å)")

    # ── Step 6: Mutation scoring ──────────────────────────────────────────
    top_mutations = suggest_mutations(
        df_residues, args.improve, args.mutations,
        interchain_contacts,
    )
    print(f"\n[RESULT] Top {len(top_mutations)} mutations for '{args.improve}':")
    for i, m in enumerate(top_mutations, 1):
        cat = " ★ catalytic" if m["is_catalytic"] else ""
        iface = " ⬡ interface" if m["is_interface"] else ""
        print(f"   {i}. {m['mutation']:<22} score={m['score']:.3f}{cat}{iface}")
    print()

    # ── Step 7: Validation ────────────────────────────────────────────────
    validation_results = {}
    if args.validate:
        validation_results = run_full_validation(
            pdb_path, top_mutations, output_dir,
            foldx_binary=args.foldx_bin if "foldx" in args.validate else "no_foldx",
            rosetta_dir=args.rosetta_dir if "rosetta" in args.validate else "no_rosetta",
        )
        print(f"\n[VALIDATION] Score–ΔΔG correlation: "
              f"{validation_results.get('correlation', 'N/A')}")
        if "foldx_correlation" in validation_results:
            print(f"           FoldX correlation: "
                  f"{validation_results['foldx_correlation']}")
        if "rosetta_correlation" in validation_results:
            print(f"           Rosetta correlation: "
                  f"{validation_results['rosetta_correlation']}")

    # ── Step 8: Visualisations ────────────────────────────────────────────
    viewer_path = make_3d_html(pdb_path, top_mutations, output_dir)
    props_path = make_residue_properties_html(df_residues, top_mutations, output_dir)
    png_path = make_mutation_score_png(top_mutations, baseline_activity, output_dir)

    paths = {
        "viewer": viewer_path,
        "properties": props_path,
        "png": png_path,
        "report": output_dir / "report.txt",
    }

    # ── Step 9: Report ────────────────────────────────────────────────────
    report_path = write_report(
        output_dir, pdb_path, args.improve,
        df_residues, catalytic_ids, baseline_feats, baseline_activity,
        top_mutations, paths, validation_results, chains_info,
    )

    print(f"\n[DONE] All outputs written to: {output_dir.resolve()}")
    print(f"       3D viewer      → {viewer_path}")
    print(f"       Properties     → {props_path}")
    print(f"       Score plot     → {png_path}")
    print(f"       Report         → {report_path}")


def main():
    args = parse_args()

    if args.gui:
        run_gui()
    elif args.api:
        run_api(args)
    else:
        if not args.pdb and not args.seq:
            print("[ERROR] Provide --pdb or --seq (or use --gui/--api).",
                  file=sys.stderr)
            sys.exit(1)
        run_analysis(args)


if __name__ == "__main__":
    main()
