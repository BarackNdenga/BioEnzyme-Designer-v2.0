#!/usr/bin/env python3
"""
Residue analysis — supports multi-chain (multimeric) enzyme structures.

When analysing a multimer, all chains are processed. Residue features
include chain identity so that inter-chain interactions can be scored.
"""

import numpy as np
import pandas as pd

from Bio.PDB import is_aa

from .constants import (
    AA_HYDROPHOBICITY, AA_CHARGE, AA_VOLUME, CATALYTIC_CANDIDATES,
)


def get_residues_all_chains(model):
    """
    Return all standard amino acid residues from ALL chains.
    Each residue record includes the chain_id for multimeric tracking.

    Returns a list of (chain_id, residue) tuples.
    """
    all_residues = []
    for chain in model.get_chains():
        for res in chain.get_residues():
            if is_aa(res, standard=True):
                all_residues.append((chain.id, res))
    return all_residues


def get_residues_single_chain(model, chain_id=None):
    """
    Return standard amino acid residues from a single chain.
    If chain_id is None, returns the first chain (legacy behaviour).
    """
    residues = []
    for chain in model.get_chains():
        if chain_id is None or chain.id == chain_id:
            for res in chain.get_residues():
                if is_aa(res, standard=True):
                    residues.append((chain.id, res))
            break  # Only first matching chain
    return residues


def identify_catalytic_residues(residues_with_chain):
    """
    Heuristic identification of catalytic residues across all chains.
    A residue is considered catalytic if:
      - Its 3-letter code is in CATALYTIC_CANDIDATES.
      - It lies within the first two-thirds of the chain (active site proximity proxy)
        OR its B-factor is above average (high mobility → active site loops).

    Returns a set of (chain_id, res_seq_id) tuples.
    """
    # Compute B-factors per chain separately
    b_factors_per_chain = {}
    for chain_id, res in residues_with_chain:
        atoms = list(res.get_atoms())
        b = np.mean([a.bfactor for a in atoms]) if atoms else 0.0
        b_factors_per_chain.setdefault(chain_id, []).append(b)

    mean_b_per_chain = {
        c: np.mean(bs) if bs else 0.0
        for c, bs in b_factors_per_chain.items()
    }

    catalytic_ids = set()
    # Count per chain to know n for the 2/3 rule
    counts_per_chain = {c: len(bs) for c, bs in b_factors_per_chain.items()}

    idx_per_chain = {}
    for chain_id, res in residues_with_chain:
        idx = idx_per_chain.setdefault(chain_id, 0)
        n = counts_per_chain[chain_id]
        mean_b = mean_b_per_chain[chain_id]

        atoms = list(res.get_atoms())
        b = np.mean([a.bfactor for a in atoms]) if atoms else 0.0

        name = res.get_resname().upper()
        in_first_two_thirds = idx < (2 * n // 3)
        high_mobility = b > mean_b * 1.2

        if name in CATALYTIC_CANDIDATES and (in_first_two_thirds or high_mobility):
            catalytic_ids.add((chain_id, res.get_id()[1]))

        idx_per_chain[chain_id] = idx + 1

    return catalytic_ids


def compute_residue_features(residues_with_chain, catalytic_ids):
    """
    Compute a feature vector for each residue across ALL chains.
    Each row includes chain_id for multimeric tracking.

    Features:
      - hydrophobicity (Kyte-Doolittle)
      - charge (formal, approximate)
      - van der Waals volume
      - mean B-factor
      - approximate RSA (relative surface area, 0-1 range via B-factor proxy)
      - is_catalytic flag

    Returns a DataFrame with one row per residue.
    """
    rows = []
    for chain_id, res in residues_with_chain:
        name = res.get_resname().upper()
        atoms = list(res.get_atoms())
        b = np.mean([a.bfactor for a in atoms]) if atoms else 0.0

        is_cat = 1 if (chain_id, res.get_id()[1]) in catalytic_ids else 0

        rows.append({
            "chain_id":       chain_id,
            "res_id":         res.get_id()[1],
            "res_name":       name,
            "hydrophobicity": AA_HYDROPHOBICITY.get(name, 0.0),
            "charge":         AA_CHARGE.get(name, 0.0),
            "volume":         AA_VOLUME.get(name, 100.0),
            "b_factor":       b,
            "approx_rsa":     min(b / 80.0, 1.0),
            "is_catalytic":   is_cat,
        })

    return pd.DataFrame(rows)


def extract_global_features(df_residues: pd.DataFrame) -> dict:
    """Compute global enzyme descriptors from the per-residue DataFrame."""
    return {
        "mean_hydrophobicity": df_residues["hydrophobicity"].mean(),
        "mean_charge":         df_residues["charge"].mean(),
        "active_site_count":   int(df_residues["is_catalytic"].sum()),
        "mean_b_factor":       df_residues["b_factor"].mean(),
        "global_volume":       df_residues["volume"].mean(),
        "n_chains":            df_residues["chain_id"].nunique(),
    }


def compute_interchain_contacts(model, cutoff=8.0):
    """
    Identify inter-chain residue contacts within a distance cutoff.
    Returns a list of (chain_a, res_a, chain_b, res_b, distance) tuples.
    Useful for scoring interface mutations in multimers.
    """
    contacts = []
    chains = list(model.get_chains())

    for i, ca in enumerate(chains):
        for cb in chains[i + 1:]:
            for res_a in ca.get_residues():
                if not is_aa(res_a, standard=True):
                    continue
                for res_b in cb.get_residues():
                    if not is_aa(res_b, standard=True):
                        continue
                    # Compute centroid distance
                    atoms_a = list(res_a.get_atoms())
                    atoms_b = list(res_b.get_atoms())
                    if not atoms_a or not atoms_b:
                        continue
                    centroid_a = np.mean([a.coord for a in atoms_a], axis=0)
                    centroid_b = np.mean([b.coord for b in atoms_b], axis=0)
                    dist = np.linalg.norm(centroid_a - centroid_b)
                    if dist < cutoff:
                        contacts.append({
                            "chain_a": ca.id,
                            "res_a":   res_a.get_id()[1],
                            "chain_b": cb.id,
                            "res_b":   res_b.get_id()[1],
                            "distance": round(dist, 2),
                        })
    return contacts
