#!/usr/bin/env python3
"""
Structure loading and fetching.

Supports:
  - Local PDB / CIF files
  - AlphaFold structure fetching (via UniProt lookup)
  - RCSB PDB REST API fallback
  - Multi-chain (multimeric) structure handling
"""

import os
import re
import time
from pathlib import Path

import requests

from Bio.PDB import PDBParser, is_aa
try:
    from Bio.PDB.Polypeptide import three_to_one
except ImportError:
    # Biopython >=1.85: three_to_one moved to SeqUtils
    from Bio.SeqUtils import seq1 as _seq1
    from Bio.SeqUtils import seq3 as _seq3
    # Build a simple mapping from standard_aa_names
    from Bio.PDB.Polypeptide import standard_aa_names
    _AA_MAP = {name: _seq1(name) for name in standard_aa_names}
    def three_to_one(name):
        return _AA_MAP.get(name, 'X')


# ─────────────────────────────────────────────────────────────────────────────
# Structure fetching
# ─────────────────────────────────────────────────────────────────────────────

def fetch_alphafold_pdb(sequence: str, output_dir: Path) -> Path:
    """
    Attempt to download an AlphaFold structure for the given sequence.
    AlphaFold structures are keyed by UniProt accession, so we first do a
    BLAST-lite lookup via the NCBI E-utilities REST API to find the best
    UniProt match, then pull the CIF/PDB from the AlphaFold CDN.
    Falls back to a RCSB text-search approach if AlphaFold fails.
    """
    print("[INFO] Searching AlphaFold for sequence ...")
    # Step 1: NCBI BLAST (quick) to get a UniProt accession
    blast_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    params_put = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": "swissprot",
        "QUERY": sequence,
        "FORMAT_TYPE": "JSON2",
        "HITLIST_SIZE": 1,
    }
    try:
        r = requests.post(blast_url, data=params_put, timeout=30)
        r.raise_for_status()
        rid = None
        for line in r.text.splitlines():
            if "RID" in line and "=" in line:
                rid = line.strip().split("=")[-1].strip()
                break
        if rid:
            print(f"[INFO] BLAST RID={rid}. Waiting for results ...")
            for _ in range(12):
                time.sleep(10)
                params_get = {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2"}
                r2 = requests.get(blast_url, params=params_get, timeout=30)
                if "Status=WAITING" not in r2.text:
                    accs = re.findall(r"\bsp\|([A-Z][0-9][A-Z0-9]{3}[0-9])\|", r2.text)
                    if accs:
                        uniprot_acc = accs[0]
                        af_url = (
                            f"https://alphafold.ebi.ac.uk/files/"
                            f"AF-{uniprot_acc}-F1-model_v4.pdb"
                        )
                        print(f"[INFO] Trying AlphaFold URL: {af_url}")
                        pdb_r = requests.get(af_url, timeout=60)
                        if pdb_r.status_code == 200:
                            out_path = output_dir / "fetched_structure.pdb"
                            out_path.write_text(pdb_r.text)
                            print(f"[INFO] AlphaFold structure saved to {out_path}")
                            return out_path
                    break
    except Exception as e:
        print(f"[WARN] AlphaFold fetch via BLAST failed: {e}")

    # Fallback: RCSB text search
    return fetch_rcsb_fallback(sequence, output_dir)


def fetch_rcsb_fallback(sequence: str, output_dir: Path) -> Path:
    """
    Fallback: query RCSB PDB REST API by sequence BLAST to get a PDB ID,
    then download that structure.
    """
    print("[INFO] Trying RCSB sequence search ...")
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": 0.5,
                "sequence_type": "protein",
                "value": sequence,
            },
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 1}},
    }
    try:
        r = requests.post(search_url, json=query, timeout=30)
        r.raise_for_status()
        data = r.json()
        pdb_id = data["result_set"][0]["identifier"]
        print(f"[INFO] Best RCSB hit: {pdb_id}")
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        r2 = requests.get(pdb_url, timeout=60)
        r2.raise_for_status()
        out_path = output_dir / f"{pdb_id}.pdb"
        out_path.write_text(r2.text)
        print(f"[INFO] RCSB structure saved to {out_path}")
        return out_path
    except Exception as e:
        raise RuntimeError(
            f"Could not fetch any structure for the given sequence. "
            f"Last error: {e}\n"
            "Tip: provide a PDB file directly with --pdb."
        )


def load_structure(pdb_path: Path):
    """
    Parse a PDB file with Biopython's PDBParser.
    Returns the structure and the first model.
    Supports multi-chain (multimeric) structures.
    """
    print(f"[INFO] Loading structure from {pdb_path} ...")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("enzyme", str(pdb_path))
    model = structure[0]
    chains = list(model.get_chains())
    chain_ids = [c.id for c in chains]
    print(f"[INFO] Structure loaded. Chains: {chain_ids}")
    print(f"[INFO] Number of chains: {len(chain_ids)}")
    for c in chains:
        n_res = sum(1 for r in c.get_residues() if is_aa(r, standard=True))
        print(f"       Chain {c.id}: {n_res} standard residues")
    return structure, model


def get_chains_info(model):
    """
    Get information about all chains in the structure.
    Returns a list of dicts with chain info.
    """
    chains_info = []
    for chain in model.get_chains():
        residues = [r for r in chain.get_residues() if is_aa(r, standard=True)]
        chains_info.append({
            "chain_id": chain.id,
            "n_residues": len(residues),
            "residue_ids": [r.get_id()[1] for r in residues],
        })
    return chains_info
