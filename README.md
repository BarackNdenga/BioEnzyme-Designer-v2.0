# BioEnzyme Designer v2.0

**Computational enzyme engineering toolkit** with machine learning, multi-chain analysis, ΔΔG validation, graphical interface, REST API, and cloud deployment.

---

## What's New in v2.0

| Feature | Description |
|---------|-------------|
| **BRENDA / SABIO-RK** | Train ML models on real experimental enzyme kinetics data from public databases |
| **FoldX / Rosetta** | Validate mutation scores with rigorous ΔΔG calculations |
| **Multimeric support** | Analyse all chains of a homomer, with inter-chain contact scoring |
| **Tkinter GUI** | Point-and-click interface for non-CLI users |
| **REST API** | FastAPI server for integration into automated bioinformatics pipelines |
| **Cloud deployment** | Docker, AWS ECS Fargate, Google Cloud Run, CloudFormation |

---

## Quick Start

### Installation

```bash
# From source
git clone https://github.com/your-repo/bioenzyme_v2.git
cd bioenzyme_v2
pip install -e .

# With API support
pip install -e ".[api]"

# With all features
pip install -e ".[all]"
```

### CLI Usage

```bash
# Basic: analyse a PDB file
python -m bioenzyme_v2 --pdb enzyme.pdb --improve activity

# Multi-chain analysis
python -m bioenzyme_v2 --pdb enzyme.pdb --improve thermostability --chains all

# With real data from BRENDA/SABIO-RK
python -m bioenzyme_v2 --pdb enzyme.pdb --data combined --ec 1.1.1.1

# With ΔΔG validation
python -m bioenzyme_v2 --pdb enzyme.pdb --validate foldx --validate rosetta

# Specific chains only
python -m bioenzyme_v2 --pdb enzyme.pdb --chains A,C --mutations 15

# Fetch structure by sequence
python -m bioenzyme_v2 --seq MKVL... --improve pH_optimum
```

### GUI

```bash
python -m bioenzyme_v2 --gui
```

### REST API

```bash
# Start server
python -m bioenzyme_v2 --api --port 8000

# API docs at http://localhost:8000/docs

# Upload a PDB and analyse
curl -X POST http://localhost:8000/analyze \
  -F "file=@enzyme.pdb" \
  -F "improve=activity" \
  -F "data_source=combined"

# Check job status
curl http://localhost:8000/status/{job_id}

# Get results
curl http://localhost:8000/mutations/{job_id}
curl http://localhost:8000/report/{job_id}
```

### Cloud Deployment

```bash
# Docker (local or any cloud)
docker compose -f cloud/docker-compose.yml up

# AWS ECS Fargate
bash cloud/deploy_aws.sh

# Google Cloud Run
bash cloud/deploy_gcp.sh
```

---

## Architecture

```
bioenzyme_v2/
├── __init__.py
├── __main__.py          # Package entry point
├── cli.py               # CLI argument parser + dispatcher
├── setup.py             # Package installation
├── requirements.txt     # Dependencies
│
├── core/                # Core analysis engine
│   ├── __init__.py
│   ├── constants.py     # AA physicochemical properties
│   ├── structure.py     # PDB loading, AlphaFold/RCSB fetching
│   ├── analysis.py      # Multi-chain residue analysis
│   ├── scoring.py       # Mutation scoring engine
│   ├── ml_model.py      # ML model (synthetic + BRENDA + SABIO-RK)
│   ├── visualization.py # 3D viewer, Plotly, Matplotlib
│   └── report.py        # Text report generator
│
├── data_sources/        # Public database connectors
│   └── __init__.py
│
├── validation/          # ΔΔG validation
│   ├── __init__.py
│   └── ddg_validation.py  # FoldX, Rosetta, knowledge-based
│
├── gui/                 # Graphical interface
│   ├── __init__.py
│   └── tkinter_gui.py   # Full Tkinter GUI
│
├── api/                 # REST API
│   ├── __init__.py
│   └── rest_api.py      # FastAPI server + endpoints
│
├── cloud/               # Cloud deployment
│   ├── __init__.py
│   ├── Dockerfile       # Multi-stage Docker build
│   ├── docker-compose.yml
│   ├── deploy_aws.sh    # AWS ECS Fargate deploy
│   ├── deploy_gcp.sh    # Google Cloud Run deploy
│   └── cloudformation.yml  # AWS CloudFormation template
│
├── templates/           # HTML/JS templates
├── docs/                # Documentation
└── tests/               # Test suite
```

---

## Data Sources

### BRENDA
The BRENDA enzyme database provides experimentally measured kcat and Km values. The REST API is queried automatically; when unavailable, a realistic fallback dataset is generated from published BRENDA statistics.

### SABIO-RK
SABIO-RK (System for the Analysis of Biochemical Pathways - Reaction Kinetics) provides curated biochemical reaction kinetics. Integrated via REST API with automatic fallback.

### Combined Mode
The `--data combined` flag merges synthetic, BRENDA, and SABIO-RK data for maximum robustness. The RandomForest model is trained on all sources simultaneously.

---

## Validation Pipeline

When `--validate foldx` and/or `--validate rosetta` are specified, the tool:

1. Runs FoldX `BuildModel` to compute ΔΔG for each proposed mutation
2. Runs Rosetta `ddg_monomer` for physics-based ΔΔG estimation
3. Computes Pearson correlation between mutation scores and ΔΔG values
4. Reports agreement rate (high score + negative ΔΔG = stabilising)

If FoldX or Rosetta binaries are not installed, knowledge-based energy estimations are used as a fallback (within ±2 kcal/mol).

---

## Multimeric Analysis

By default, `--chains all` processes every chain in the structure:

- Residues are tagged with chain identifiers
- Catalytic residues are identified per-chain
- Inter-chain contacts (within 8 Å) are computed
- Interface residues receive additional scoring bonuses
- The 3D viewer displays chain-aware highlighting

For specific chains: `--chains A,B` analyses only chains A and B.

---

## License

MIT License. See LICENSE file for details.
