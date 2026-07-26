#!/usr/bin/env python3
"""
Mutation scoring engine.

Evaluates all possible single point mutations for all residues across
all chains. Supports inter-chain contact bonuses for multimeric enzymes.
"""

try:
    from .constants import (
        AA_HYDROPHOBICITY, AA_CHARGE, AA_VOLUME,
        ALL_AA, THERMOSTABILITY_FAVOR, ACTIVITY_FAVOR, PH_FAVOR,
    )
except ImportError:
    from bioenzyme_v2.core.constants import (
        AA_HYDROPHOBICITY, AA_CHARGE, AA_VOLUME,
        ALL_AA, THERMOSTABILITY_FAVOR, ACTIVITY_FAVOR, PH_FAVOR,
    )


def score_mutation(orig_aa: str, mut_aa: str, res_row: dict, improve: str,
                   is_interface: bool = False) -> float:
    """
    Score a single point mutation at one residue position.

    Energy-like scoring function — higher score is *better*.
    Components:
      Δhydrophobicity  — change in Kyte-Doolittle score
      Δcharge          — change in formal charge
      Δvolume          — change in side-chain volume (steric clash proxy)
      active_site_bonus — reward mutations at catalytic residues for activity
      target_bonus      — reward preferred amino acids for the chosen target
      interface_bonus   — reward mutations at inter-chain contacts for multimers
    """
    delta_hydro  = AA_HYDROPHOBICITY.get(mut_aa, 0) - AA_HYDROPHOBICITY.get(orig_aa, 0)
    delta_charge = AA_CHARGE.get(mut_aa, 0)         - AA_CHARGE.get(orig_aa, 0)
    delta_vol    = AA_VOLUME.get(mut_aa, 100)        - AA_VOLUME.get(orig_aa, 100)

    is_catalytic = bool(res_row["is_catalytic"])
    approx_rsa   = float(res_row["approx_rsa"])

    score = 0.0

    if improve == "thermostability":
        # Hydrophobic mutations in buried residues increase thermal stability
        buried = approx_rsa < 0.3
        if buried:
            score += delta_hydro * 1.5          # bury hydrophobics → stable
        score -= abs(delta_charge) * 2.0        # charge changes destabilise core
        score -= abs(delta_vol)    * 0.02       # steric clash penalty
        if mut_aa in THERMOSTABILITY_FAVOR:
            score += 2.0
        # Interface residues: strengthening hydrophobic contacts helps
        if is_interface and mut_aa in THERMOSTABILITY_FAVOR:
            score += 1.5

    elif improve == "activity":
        # Active-site flexibility and catalytic residue identity matter most
        if is_catalytic:
            score += 3.0 if mut_aa in ACTIVITY_FAVOR else -1.0
        score += approx_rsa * 1.0               # exposed loops tolerate mutations
        score -= abs(delta_vol) * 0.01          # small steric penalty
        if mut_aa in ACTIVITY_FAVOR:
            score += 1.5
        # Interface mutations that improve substrate channeling
        if is_interface and mut_aa in ACTIVITY_FAVOR:
            score += 1.0

    elif improve == "pH_optimum":
        # Ionisable residues shift pKa environment of the active site
        score += abs(delta_charge) * 1.5        # charge changes shift pH optimum
        if is_catalytic and mut_aa in PH_FAVOR:
            score += 3.5
        if mut_aa in PH_FAVOR:
            score += 1.0

    # Universal penalty: avoid mutations to identical amino acid
    if orig_aa == mut_aa:
        score = -999.0

    return score


def suggest_mutations(df_residues, improve: str, n: int,
                      interchain_contacts: list = None):
    """
    Evaluate all possible single point mutations for all residues
    across ALL chains. Return the top-N mutations sorted by score.

    Parameters
    ----------
    df_residues : pd.DataFrame
        Per-residue features (includes chain_id column).
    improve : str
        One of "activity", "thermostability", "pH_optimum".
    n : int
        Number of top mutations to return.
    interchain_contacts : list, optional
        List of dicts with keys chain_a, res_a, chain_b, res_b.
        Residues involved in inter-chain contacts get an interface bonus.

    Returns
    -------
    list of dict
        Each dict: {chain_id, res_id, res_name, mutation, score,
                    is_catalytic, is_interface}
    """
    # Build set of interface residues from inter-chain contacts
    interface_residues = set()
    if interchain_contacts:
        for c in interchain_contacts:
            interface_residues.add((c["chain_a"], c["res_a"]))
            interface_residues.add((c["chain_b"], c["res_b"]))

    print(f"[INFO] Scoring mutations across {df_residues['chain_id'].nunique()} "
          f"chain(s), {len(df_residues)} residues ...")

    candidates = []
    for _, row in df_residues.iterrows():
        orig = row["res_name"]
        chain_id = row["chain_id"]
        is_interface = (chain_id, int(row["res_id"])) in interface_residues

        for mut_aa in ALL_AA:
            s = score_mutation(
                orig, mut_aa, row, improve,
                is_interface=is_interface,
            )
            mutation_label = f"{chain_id}:{orig}{int(row['res_id'])}{mut_aa}"
            candidates.append({
                "chain_id":     chain_id,
                "res_id":       int(row["res_id"]),
                "res_name":     orig,
                "mutation":     mutation_label,
                "mut_aa":       mut_aa,
                "score":        s,
                "is_catalytic": bool(row["is_catalytic"]),
                "is_interface": is_interface,
            })

    # Sort by score descending, take top N (skip self-mutations via score=-999)
    candidates.sort(key=lambda x: x["score"], reverse=True)
    top = [c for c in candidates if c["score"] > -900][:n]
    return top
