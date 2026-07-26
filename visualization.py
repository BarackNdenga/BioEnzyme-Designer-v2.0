#!/usr/bin/env python3
"""
Visualisation generators for BioEnzyme Designer.

Generates:
  - 3D HTML viewer (py3Dmol) — interactive molecular viewer
  - Residue properties plot (Plotly) — 3D scatter
  - Mutation score bar chart (Matplotlib) — PNG
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.express as px


def make_3d_html(pdb_path: Path, mutations: list, output_dir: Path) -> Path:
    """
    Generate an HTML file with a py3Dmol interactive 3D viewer.
    - Backbone shown as a cartoon, coloured by secondary structure.
    - Original mutation sites highlighted in red spheres.
    - Catalytic residues highlighted in blue sticks.
    - Chain-aware: uses chain-specific residue selectors.
    """
    print("[INFO] Generating 3D HTML viewer ...")
    pdb_text = pdb_path.read_text()

    # Collect residue IDs per chain for mutation sites
    mut_res_by_chain = {}
    for m in mutations:
        chain = m.get("chain_id", "A")
        mut_res_by_chain.setdefault(chain, set()).add(int(m["res_id"]))

    # Build chain-aware selectors
    style_blocks = []
    for chain_id, res_ids in mut_res_by_chain.items():
        resid_str = "+".join(str(r) for r in sorted(res_ids))
        sel = json.dumps({"chain": chain_id, "resi": resid_str})
        style_blocks.append(
            f'viewer.addStyle({sel}, '
            f'{{sphere: {{color: "red", opacity: 0.6, radius: 1.4}}}});'
        )

    # Catalytic residues (blue sticks) — same chain selector
    cat_style_blocks = []
    for chain_id, res_ids in mut_res_by_chain.items():
        resid_str = "+".join(str(r) for r in sorted(res_ids))
        sel = json.dumps({"chain": chain_id, "resi": resid_str})
        cat_style_blocks.append(
            f'viewer.addStyle({sel}, '
            f'{{stick: {{color: "#4488ff", radius: 0.3}}}});'
        )

    spheres_js = "\n    ".join(style_blocks)
    sticks_js  = "\n    ".join(cat_style_blocks)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <title>BioEnzyme Designer — 3D Viewer</title>
  <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
  <style>
    body {{ margin: 0; background: #1a1a2e; color: #eee; font-family: sans-serif; }}
    #viewer {{ width: 100vw; height: 85vh; position: relative; }}
    #legend {{ padding: 10px 20px; font-size: 14px; }}
    .dot {{ display: inline-block; width: 12px; height: 12px; border-radius: 50%;
             margin-right: 6px; vertical-align: middle; }}
  </style>
</head>
<body>
  <div id="viewer"></div>
  <div id="legend">
    <span class="dot" style="background:#4488ff"></span> Catalytic / mutation sites &nbsp;&nbsp;
    <span class="dot" style="background:#ff4444"></span> Suggested mutations &nbsp;&nbsp;
    <span class="dot" style="background:#aaaaaa"></span> Other residues
  </div>
  <script>
    var config = {{ backgroundColor: "#1a1a2e" }};
    var viewer = $3Dmol.createViewer("viewer", config);

    var pdbData = `{pdb_text.replace("`", "'")}`;
    viewer.addModel(pdbData, "pdb");

    // Full cartoon rendering
    viewer.setStyle({{}}, {{ cartoon: {{ color: "spectrum" }} }});

    // Highlight mutation sites — red transparent spheres
    {spheres_js}

    // Highlight catalytic residues — blue sticks
    {sticks_js}

    viewer.zoomTo();
    viewer.render();
  </script>
</body>
</html>"""

    out_path = output_dir / "enzyme_view.html"
    out_path.write_text(html)
    print(f"[INFO] 3D viewer saved to {out_path}")
    return out_path


def make_residue_properties_html(df_residues, mutations: list,
                                  output_dir: Path) -> Path:
    """
    Interactive 3D scatter plot (Plotly):
    x = residue position, y = B-factor, z = hydrophobicity
    Mutation sites are coloured red; others are coloured by charge.
    Chain-aware labelling.
    """
    print("[INFO] Generating residue properties plot ...")
    mut_ids = {(m.get("chain_id", "A"), m["res_id"]) for m in mutations}
    df = df_residues.copy()
    df["colour"] = df.apply(
        lambda r: "Mutation site"
        if (r["chain_id"], r["res_id"]) in mut_ids else "Other",
        axis=1,
    )
    df["label"] = df["chain_id"] + ":" + df["res_name"] + df["res_id"].astype(str)

    fig = px.scatter_3d(
        df,
        x="res_id",
        y="b_factor",
        z="hydrophobicity",
        color="colour",
        facet_col="chain_id" if df["chain_id"].nunique() > 1 else None,
        color_discrete_map={"Mutation site": "#ff4444", "Other": "#4488ff"},
        hover_name="label",
        hover_data={"charge": True, "volume": True, "approx_rsa": True},
        labels={
            "res_id": "Residue position",
            "b_factor": "B-factor (Å²)",
            "hydrophobicity": "Hydrophobicity (KD)",
        },
        title="Residue Properties — B-factor × Hydrophobicity × Position",
        template="plotly_dark",
    )
    fig.update_traces(marker=dict(size=4))
    out_path = output_dir / "residue_properties.html"
    fig.write_html(str(out_path))
    print(f"[INFO] Residue properties plot saved to {out_path}")
    return out_path


def make_mutation_score_png(mutations: list, baseline_activity: float,
                             output_dir: Path) -> Path:
    """
    Matplotlib bar chart:
    x = mutation rank, y = predicted activity improvement (Δlog kcat/Km proxy = score).
    Saves as PNG.
    """
    print("[INFO] Generating mutation score PNG ...")
    names   = [m["mutation"] for m in mutations]
    scores  = [m["score"] for m in mutations]
    colours = []
    for m in mutations:
        if m.get("is_catalytic"):
            colours.append("#ff4444")
        elif m.get("is_interface"):
            colours.append("#44ff88")
        else:
            colours.append("#4488ff")

    fig, ax = plt.subplots(figsize=(max(8, len(mutations) * 1.4), 5))
    bars = ax.barh(names[::-1], scores[::-1], color=colours[::-1])
    ax.set_xlabel("Mutation score (higher = better)", fontsize=12)
    ax.set_title("Top Suggested Mutations by Score", fontsize=14)
    ax.axvline(x=0, color="white", linewidth=0.8, linestyle="--")
    ax.set_facecolor("#1a1a2e")
    fig.patch.set_facecolor("#12122a")
    ax.tick_params(colors="white")
    ax.xaxis.label.set_color("white")
    ax.title.set_color("white")
    for spine in ax.spines.values():
        spine.set_edgecolor("#444")

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#ff4444", label="Catalytic residue"),
        Patch(facecolor="#44ff88", label="Interface residue"),
        Patch(facecolor="#4488ff", label="Other residue"),
    ]
    ax.legend(handles=legend_elements, facecolor="#1a1a2e", labelcolor="white",
              loc="lower right")

    plt.tight_layout()
    out_path = output_dir / "mutation_scores.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Mutation score PNG saved to {out_path}")
    return out_path
