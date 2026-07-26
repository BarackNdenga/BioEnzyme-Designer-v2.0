#!/usr/bin/env python3
"""
Text report generator — summarises the full analysis.
Chain-aware: reports on all chains of a multimer.
"""

from pathlib import Path


def write_report(output_dir: Path, pdb_path: Path, improve: str,
                 df_residues, catalytic_ids: set,
                 baseline_features: dict, baseline_activity: float,
                 mutations: list, paths: dict,
                 validation_results: dict = None,
                 chains_info: list = None) -> Path:
    """Write a comprehensive plain-text summary report."""
    print("[INFO] Writing report ...")
    n_chains = df_residues["chain_id"].nunique()

    lines = [
        "=" * 70,
        "  BioEnzyme Designer v2.0 — Analysis Report",
        "=" * 70,
        "",
        f"  PDB source      : {pdb_path}",
        f"  Optimise for    : {improve}",
        f"  Total residues  : {len(df_residues)}",
        f"  Chains analysed : {n_chains}",
        "",
    ]

    # Chains info
    if chains_info:
        lines.append("── Chain Information ────────────────────────────────────────────────")
        for ci in chains_info:
            lines.append(
                f"  Chain {ci['chain_id']:<4}: {ci['n_residues']:>5} residues  "
                f"(IDs: {ci['residue_ids'][0]}–{ci['residue_ids'][-1]})"
            )
        lines.append("")

    lines.append("── Global Properties ──────────────────────────────────────────────")
    lines.append(f"  Mean hydrophobicity : {baseline_features['mean_hydrophobicity']:.3f}")
    lines.append(f"  Mean charge/residue : {baseline_features['mean_charge']:.3f}")
    lines.append(f"  Mean B-factor       : {baseline_features['mean_b_factor']:.2f} Å²")
    lines.append(f"  Mean residue volume : {baseline_features['global_volume']:.1f} Å³")
    lines.append(f"  Catalytic residues  : {baseline_features['active_site_count']}")
    lines.append("")

    lines.append("── Predicted Catalytic Efficiency ─────────────────────────────────")
    lines.append(f"  log10(kcat/Km) [baseline] : {baseline_activity:.3f}")
    lines.append(f"  kcat/Km estimate          : {10**baseline_activity:.2e} M⁻¹s⁻¹")
    lines.append("")

    # Catalytic residues per chain
    lines.append("── Catalytic Residue Positions ────────────────────────────────────")
    for cid in sorted(catalytic_ids):
        chain_id, res_id = cid
        row = df_residues[(df_residues["chain_id"] == chain_id) &
                          (df_residues["res_id"] == res_id)]
        if not row.empty:
            r = row.iloc[0]
            lines.append(
                f"  Chain {chain_id}: {r['res_name']}{res_id:>5}  "
                f"hydro={r['hydrophobicity']:+.1f}  "
                f"charge={r['charge']:+.1f}  B={r['b_factor']:.1f}"
            )
    lines.append("")

    # Mutations
    lines.append("── Suggested Mutations ─────────────────────────────────────────────")
    lines.append(
        f"  {'Rank':<6} {'Mutation':<22} {'Score':>8}  "
        f"{'Catalytic':>10} {'Interface':>10}"
    )
    lines.append("  " + "-" * 62)
    for rank, m in enumerate(mutations, 1):
        cat_flag = "YES ★" if m.get("is_catalytic") else "no"
        if_flag = "YES" if m.get("is_interface") else "—"
        lines.append(
            f"  {rank:<6} {m['mutation']:<22} {m['score']:>8.3f}  "
            f"{cat_flag:>10} {if_flag:>10}"
        )
    lines.append("")

    # Validation results
    if validation_results:
        lines.append("── FoldX / Rosetta Validation ─────────────────────────────────────")
        vr = validation_results
        if "foldx_ddg" in vr:
            lines.append(f"  FoldX ΔΔG (kcal/mol):")
            for m in vr.get("foldx_ddg", []):
                lines.append(
                    f"    {m['mutation']:<22} ΔΔG = {m['ddg']:>+.3f} kcal/mol"
                )
        if "rosetta_ddg" in vr:
            lines.append(f"  Rosetta ΔΔG (REU):")
            for m in vr.get("rosetta_ddg", []):
                lines.append(
                    f"    {m['mutation']:<22} ΔΔG = {m['ddg']:>+.3f} REU"
                )
        if vr.get("correlation"):
            lines.append(f"  Score–ΔΔG correlation: {vr['correlation']:.3f}")
        lines.append("")

    # Output files
    lines.append("── Output Files ────────────────────────────────────────────────────")
    for key, label in [
        ("viewer", "3D viewer"),
        ("properties", "Residue properties"),
        ("png", "Mutation scores"),
        ("report", "This report"),
    ]:
        if key in paths:
            lines.append(f"  {label:<20}: {paths[key]}")
    lines.append("")

    lines.append("── Interpretation Notes ────────────────────────────────────────────")
    lines.append("  • Catalytic residue identification is heuristic (sequence position")
    lines.append("    + B-factor). Always validate with literature / experimental data.")
    lines.append("  • The ML model was trained on real data (BRENDA/SABIO-RK)")
    lines.append("    supplemented with synthetic data for robustness.")
    lines.append("  • Mutation scores are energy-like proxies. Use FoldX/Rosetta")
    lines.append("    ΔΔG calculations for higher-fidelity energetics.")
    lines.append("  • For multimeric enzymes, inter-chain interface mutations are")
    lines.append("    scored with additional bonuses for stabilising contacts.")
    lines.append("  • Cross-validate top candidates with molecular dynamics simulations.")
    lines.append("")
    lines.append("=" * 70)
    lines.append("  Generated by BioEnzyme Designer v2.0")
    lines.append("=" * 70)

    out_path = output_dir / "report.txt"
    out_path.write_text("\n".join(lines))
    print(f"[INFO] Report saved to {out_path}")
    return out_path
