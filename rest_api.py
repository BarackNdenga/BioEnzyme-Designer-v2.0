#!/usr/bin/env python3
"""
BioEnzyme Designer — REST API.

FastAPI-based REST API for integration into automated bioinformatics pipelines.

Endpoints:
  POST /analyze          — Run a full enzyme analysis
  POST /analyze/sequence — Fetch structure by sequence and analyze
  GET  /chains/{job_id}  — Get chain information for a structure
  GET  /mutations/{job_id} — Get top mutations for a completed job
  GET  /report/{job_id}  — Download the full report
  GET  /health           — Health check
  GET  /status/{job_id}  — Check job status

Usage:
    python -m bioenzyme_v2.api.rest_api

Or with uvicorn:
    uvicorn bioenzyme_v2.api.rest_api:app --host 0.0.0.0 --port 8000
"""

import sys
import os
import json
import uuid
import time
import tempfile
import threading
from pathlib import Path
from typing import Optional, List, Dict, Any

from fastapi import FastAPI, UploadFile, File, HTTPException, BackgroundTasks
from fastapi.responses import JSONResponse, FileResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# ─────────────────────────────────────────────────────────────────────────────
# Import core modules
# ─────────────────────────────────────────────────────────────────────────────

from core.structure import load_structure, fetch_alphafold_pdb, get_chains_info
from core.analysis import (
    get_residues_all_chains, get_residues_single_chain,
    identify_catalytic_residues, compute_residue_features,
    extract_global_features, compute_interchain_contacts,
)
from core.ml_model import build_training_dataset, train_activity_model, predict_activity
from core.scoring import suggest_mutations
from core.visualization import (
    make_3d_html, make_residue_properties_html, make_mutation_score_png,
)
from core.report import write_report
from validation.ddg_validation import run_full_validation


# ─────────────────────────────────────────────────────────────────────────────
# Application setup
# ─────────────────────────────────────────────────────────────────────────────

app = FastAPI(
    title="BioEnzyme Designer API",
    description="REST API for computational enzyme engineering and mutation design.",
    version="2.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# In-memory job storage (use Redis/DB in production)
jobs: Dict[str, Dict[str, Any]] = {}
OUTPUT_BASE = Path(tempfile.gettempdir()) / "bioenzyme_api"
OUTPUT_BASE.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────────────────────────────────────────────────────
# Pydantic models
# ─────────────────────────────────────────────────────────────────────────────

class AnalyzeRequest(BaseModel):
    improve: str = "activity"
    data_source: str = "combined"
    ec_number: Optional[str] = None
    organism: Optional[str] = None
    chain_mode: str = "all"       # "all" or "specific"
    chain_id: Optional[str] = "A"
    n_mutations: int = 10
    validate_foldx: bool = True
    validate_rosetta: bool = False


class SequenceAnalyzeRequest(AnalyzeRequest):
    sequence: str


class ErrorResponse(BaseModel):
    error: str
    detail: Optional[str] = None


# ─────────────────────────────────────────────────────────────────────────────
# Background job processing
# ─────────────────────────────────────────────────────────────────────────────

def process_job(job_id: str, pdb_path: Path, config: dict):
    """Process an analysis job in the background."""
    job = jobs[job_id]
    job["status"] = "running"

    try:
        output_dir = job["output_dir"]
        output_dir.mkdir(parents=True, exist_ok=True)

        # Load structure
        structure, model = load_structure(pdb_path)
        chains_info = get_chains_info(model)
        job["chains"] = chains_info

        # Residue analysis
        if config["chain_mode"] == "all":
            residues_with_chain = get_residues_all_chains(model)
        else:
            residues_with_chain = get_residues_single_chain(model, config["chain_id"])

        catalytic_ids = identify_catalytic_residues(residues_with_chain)
        df_residues = compute_residue_features(residues_with_chain, catalytic_ids)
        baseline_feats = extract_global_features(df_residues)

        # ML model
        train_df = build_training_dataset(
            data_source=config["data_source"],
            ec_number=config.get("ec_number"),
            organism=config.get("organism"),
        )
        model_rf, scaler = train_activity_model(train_df)
        baseline_activity = predict_activity(model_rf, scaler, baseline_feats)

        # Inter-chain contacts
        interchain_contacts = compute_interchain_contacts(model)

        # Mutation scoring
        top_mutations = suggest_mutations(
            df_residues, config["improve"],
            config["n_mutations"], interchain_contacts,
        )

        # Validation
        validation_results = {}
        if config.get("validate_foldx") or config.get("validate_rosetta"):
            validation_results = run_full_validation(
                pdb_path, top_mutations, output_dir,
                foldx_binary=None if config.get("validate_foldx") else "no_foldx",
                rosetta_dir=None if config.get("validate_rosetta") else "no_rosetta",
            )

        # Visualisations
        viewer_path = make_3d_html(pdb_path, top_mutations, output_dir)
        props_path = make_residue_properties_html(df_residues, top_mutations, output_dir)
        png_path = make_mutation_score_png(top_mutations, baseline_activity, output_dir)

        paths = {
            "viewer": str(viewer_path),
            "properties": str(props_path),
            "png": str(png_path),
            "report": str(output_dir / "report.txt"),
        }

        report_path = write_report(
            output_dir, pdb_path, config["improve"],
            df_residues, catalytic_ids, baseline_feats, baseline_activity,
            top_mutations, paths, validation_results, chains_info,
        )

        # Store results
        job["results"] = {
            "baseline_activity": baseline_activity,
            "kcat_km_estimate": 10**baseline_activity,
            "catalytic_residues": [
                {"chain": c, "res_id": r} for c, r in catalytic_ids
            ],
            "top_mutations": top_mutations,
            "validation": {
                k: v for k, v in validation_results.items()
                if k not in ("agreement",)
            },
            "files": paths,
        }
        job["status"] = "completed"
        job["completed_at"] = time.time()

    except Exception as e:
        job["status"] = "failed"
        job["error"] = str(e)
        job["completed_at"] = time.time()


# ─────────────────────────────────────────────────────────────────────────────
# API Endpoints
# ─────────────────────────────────────────────────────────────────────────────

@app.get("/health")
def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "version": "2.0.0"}


@app.post("/analyze", response_class=JSONResponse)
async def analyze(
    file: UploadFile = File(...),
    background_tasks: BackgroundTasks = None,
    improve: str = "activity",
    data_source: str = "combined",
    ec_number: Optional[str] = None,
    organism: Optional[str] = None,
    chain_mode: str = "all",
    chain_id: str = "A",
    n_mutations: int = 10,
    validate_foldx: bool = True,
    validate_rosetta: bool = False,
):
    """
    Upload a PDB file and run a full enzyme analysis.

    Returns a job ID for async processing.
    """
    job_id = str(uuid.uuid4())[:8]
    output_dir = OUTPUT_BASE / job_id
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save uploaded PDB
    pdb_path = output_dir / "uploaded.pdb"
    content = await file.read()
    pdb_path.write_bytes(content)

    # Store job
    config = {
        "improve": improve,
        "data_source": data_source,
        "ec_number": ec_number,
        "organism": organism,
        "chain_mode": chain_mode,
        "chain_id": chain_id,
        "n_mutations": n_mutations,
        "validate_foldx": validate_foldx,
        "validate_rosetta": validate_rosetta,
    }

    jobs[job_id] = {
        "status": "queued",
        "output_dir": output_dir,
        "pdb_path": pdb_path,
        "config": config,
        "created_at": time.time(),
    }

    # Start background processing
    if background_tasks:
        background_tasks.add_task(process_job, job_id, pdb_path, config)
    else:
        # Synchronous for testing
        thread = threading.Thread(target=process_job, args=(job_id, pdb_path, config))
        thread.start()

    return {
        "job_id": job_id,
        "status": "queued",
        "message": "Analysis started. Check /status/{job_id} for progress.",
    }


@app.post("/analyze/sequence", response_class=JSONResponse)
async def analyze_sequence(
    request: SequenceAnalyzeRequest,
    background_tasks: BackgroundTasks = None,
):
    """
    Fetch a structure by amino acid sequence and run analysis.

    The tool will attempt AlphaFold first, then RCSB fallback.
    """
    job_id = str(uuid.uuid4())[:8]
    output_dir = OUTPUT_BASE / job_id
    output_dir.mkdir(parents=True, exist_ok=True)

    config = {
        "improve": request.improve,
        "data_source": request.data_source,
        "ec_number": request.ec_number,
        "organism": request.organism,
        "chain_mode": request.chain_mode,
        "chain_id": request.chain_id,
        "n_mutations": request.n_mutations,
        "validate_foldx": request.validate_foldx,
        "validate_rosetta": request.validate_rosetta,
    }

    jobs[job_id] = {
        "status": "queued",
        "output_dir": output_dir,
        "sequence": request.sequence,
        "config": config,
        "created_at": time.time(),
    }

    def _fetch_and_process():
        try:
            pdb_path = fetch_alphafold_pdb(request.sequence, output_dir)
            jobs[job_id]["pdb_path"] = pdb_path
            process_job(job_id, pdb_path, config)
        except Exception as e:
            jobs[job_id]["status"] = "failed"
            jobs[job_id]["error"] = str(e)

    if background_tasks:
        background_tasks.add_task(_fetch_and_process)
    else:
        thread = threading.Thread(target=_fetch_and_process)
        thread.start()

    return {
        "job_id": job_id,
        "status": "queued",
        "message": f"Fetching structure for sequence (length={len(request.sequence)})...",
    }


@app.get("/status/{job_id}")
def get_job_status(job_id: str):
    """Check the status of an analysis job."""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = jobs[job_id]
    return {
        "job_id": job_id,
        "status": job["status"],
        "created_at": job.get("created_at"),
        "completed_at": job.get("completed_at"),
        "elapsed": time.time() - job.get("created_at", time.time()),
    }


@app.get("/chains/{job_id}")
def get_chains(job_id: str):
    """Get chain information for a structure."""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = jobs[job_id]
    chains = job.get("chains", [])
    return {"job_id": job_id, "chains": chains}


@app.get("/mutations/{job_id}")
def get_mutations(job_id: str):
    """Get top mutations for a completed job."""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = jobs[job_id]
    if job["status"] != "completed":
        raise HTTPException(
            status_code=400,
            detail=f"Job not completed yet. Status: {job['status']}",
        )

    results = job.get("results", {})
    return {
        "job_id": job_id,
        "baseline_activity": results.get("baseline_activity"),
        "kcat_km_estimate": results.get("kcat_km_estimate"),
        "catalytic_residues": results.get("catalytic_residues"),
        "top_mutations": results.get("top_mutations"),
    }


@app.get("/report/{job_id}")
def get_report(job_id: str):
    """Download the full analysis report."""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = jobs[job_id]
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail="Job not completed yet")

    report_path = job["output_dir"] / "report.txt"
    if report_path.exists():
        return FileResponse(
            str(report_path),
            media_type="text/plain",
            filename=f"bioenzyme_report_{job_id}.txt",
        )
    raise HTTPException(status_code=404, detail="Report file not found")


@app.get("/viewer/{job_id}")
def get_viewer(job_id: str):
    """Download the 3D HTML viewer file."""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = jobs[job_id]
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail="Job not completed yet")

    viewer_path = job["output_dir"] / "enzyme_view.html"
    if viewer_path.exists():
        return FileResponse(
            str(viewer_path),
            media_type="text/html",
            filename=f"enzyme_view_{job_id}.html",
        )
    raise HTTPException(status_code=404, detail="Viewer file not found")


# ─────────────────────────────────────────────────────────────────────────────
# CLI entry point for running the API server
# ─────────────────────────────────────────────────────────────────────────────

def run_server(host: str = "0.0.0.0", port: int = 8000, workers: int = 1):
    """Run the FastAPI server with uvicorn."""
    import uvicorn
    uvicorn.run(
        "bioenzyme_v2.api.rest_api:app",
        host=host,
        port=port,
        workers=workers,
        log_level="info",
    )


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run BioEnzyme Designer REST API")
    parser.add_argument("--host", default="0.0.0.0", help="Host to bind to")
    parser.add_argument("--port", type=int, default=8000, help="Port to bind to")
    parser.add_argument("--workers", type=int, default=1, help="Number of workers")
    args = parser.parse_args()
    run_server(host=args.host, port=args.port, workers=args.workers)
