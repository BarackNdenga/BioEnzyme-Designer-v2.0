#!/usr/bin/env python3
"""
Physicochemical constants for the 20 standard amino acids.

Sources:
  - Kyte & Doolittle (1982) — hydrophobicity
  - Standard biochemistry tables — formal charge at pH 7
  - van der Waals volumes (Tsai et al., 1999, J. Mol. Biol.)
"""

# Kyte-Doolittle hydrophobicity index
AA_HYDROPHOBICITY = {
    "ALA": 1.8,  "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5,
    "GLN": -3.5, "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5,
    "LEU": 3.8,  "LYS": -3.9, "MET": 1.9,  "PHE": 2.8,  "PRO": -1.6,
    "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
}

# Formal charge at pH 7 (approximate; HIS partial)
AA_CHARGE = {
    "ARG": 1, "HIS": 0.1, "LYS": 1,
    "ASP": -1, "GLU": -1,
    **{aa: 0 for aa in ["ALA", "ASN", "CYS", "GLN", "GLY", "ILE", "LEU",
                         "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]},
}

# van der Waals side-chain volume (Å³)
AA_VOLUME = {
    "ALA": 88.6,  "ARG": 173.4, "ASN": 114.1, "ASP": 111.1, "CYS": 108.5,
    "GLN": 143.8, "GLU": 138.4, "GLY": 60.1,  "HIS": 153.2, "ILE": 166.7,
    "LEU": 166.7, "LYS": 168.6, "MET": 162.9, "PHE": 189.9, "PRO": 112.7,
    "SER": 89.0,  "THR": 116.1, "TRP": 227.8, "TYR": 193.6, "VAL": 140.0,
}

# Residues commonly found in catalytic triads / active sites
CATALYTIC_CANDIDATES = {"HIS", "SER", "ASP", "CYS", "GLU", "TYR", "LYS"}

# Preferred amino acids for each optimisation target
THERMOSTABILITY_FAVOR = {"ILE", "LEU", "VAL", "PHE", "TRP"}
ACTIVITY_FAVOR        = {"GLY", "SER", "THR", "ASN"}
PH_FAVOR              = {"HIS", "LYS", "ARG", "ASP", "GLU"}

# All 20 standard amino acids
ALL_AA = list(AA_HYDROPHOBICITY.keys())

# Single-letter → three-letter mapping
ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}

THREE_TO_ONE = {v: k for k, v in ONE_TO_THREE.items()}
