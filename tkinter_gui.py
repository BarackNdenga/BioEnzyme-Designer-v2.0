#!/usr/bin/env python3
"""
BioEnzyme Designer — Tkinter GUI.

Provides a graphical interface for non-CLI users to:
  - Load PDB files or fetch structures by sequence
  - Choose optimisation target (activity, thermostability, pH)
  - Select data source (synthetic, BRENDA, SABIO-RK, combined)
  - Configure mutation count and chain selection
  - Enable/disable FoldX/Rosetta validation
  - View results in an integrated window

Usage:
    python -m bioenzyme_v2.gui.tkinter_gui
"""

import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from pathlib import Path

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


class BioEnzymeGUI:
    """Main Tkinter GUI for BioEnzyme Designer."""

    def __init__(self, root):
        self.root = root
        self.root.title("BioEnzyme Designer v2.0")
        self.root.geometry("900x700")
        self.root.minsize(800, 600)

        # Configure style
        style = ttk.Style()
        style.theme_use("clam")
        style.configure("Title.TLabel", font=("Helvetica", 16, "bold"))
        style.configure("Section.TLabelframe.Label", font=("Helvetica", 11, "bold"))

        self._build_ui()

    def _build_ui(self):
        """Build the complete GUI."""
        # Main container with scrollbar
        main_frame = ttk.Frame(self.root, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Title
        title_label = ttk.Label(main_frame, text="BioEnzyme Designer v2.0",
                                style="Title.TLabel")
        title_label.pack(pady=(0, 15))

        # ── Structure Input ───────────────────────────────────────────────
        struct_frame = ttk.LabelFrame(main_frame, text="Structure Input", padding=10)
        struct_frame.pack(fill=tk.X, pady=5)

        self.source_var = tk.StringVar(value="pdb")
        ttk.Radiobutton(struct_frame, text="PDB File", variable=self.source_var,
                        value="pdb").pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(struct_frame, text="Sequence (AlphaFold/RCSB fetch)",
                        variable=self.source_var, value="seq").pack(side=tk.LEFT, padx=5)

        self.pdb_path_var = tk.StringVar()
        ttk.Entry(struct_frame, textvariable=self.pdb_path_var, width=60).pack(
            side=tk.LEFT, padx=5, pady=2)
        ttk.Button(struct_frame, text="Browse",
                   command=self._browse_pdb).pack(side=tk.LEFT, padx=5)

        self.seq_entry = ttk.Entry(struct_frame, width=60)
        self.seq_entry.pack(side=tk.LEFT, padx=5, pady=2)

        # ── Optimisation Target ───────────────────────────────────────────
        opt_frame = ttk.LabelFrame(main_frame, text="Optimisation Target", padding=10)
        opt_frame.pack(fill=tk.X, pady=5)

        self.improve_var = tk.StringVar(value="activity")
        for opt in ["activity", "thermostability", "pH_optimum"]:
            ttk.Radiobutton(opt_frame, text=opt, variable=self.improve_var,
                            value=opt).pack(side=tk.LEFT, padx=10)

        # ── Data Source ───────────────────────────────────────────────────
        data_frame = ttk.LabelFrame(main_frame, text="ML Data Source", padding=10)
        data_frame.pack(fill=tk.X, pady=5)

        self.data_source_var = tk.StringVar(value="combined")
        for ds in ["synthetic", "brenda", "sabio_rk", "combined"]:
            ttk.Radiobutton(data_frame, text=ds, variable=self.data_source_var,
                            value=ds).pack(side=tk.LEFT, padx=10)

        # EC number and organism filters
        filter_frame = ttk.Frame(data_frame)
        filter_frame.pack(fill=tk.X, pady=5)
        ttk.Label(filter_frame, text="EC number:").pack(side=tk.LEFT)
        self.ec_var = tk.StringVar()
        ttk.Entry(filter_frame, textvariable=self.ec_var, width=15).pack(
            side=tk.LEFT, padx=5)
        ttk.Label(filter_frame, text="Organism:").pack(side=tk.LEFT)
        self.organism_var = tk.StringVar()
        ttk.Entry(filter_frame, textvariable=self.organism_var, width=20).pack(
            side=tk.LEFT, padx=5)

        # ── Chain Selection ───────────────────────────────────────────────
        chain_frame = ttk.LabelFrame(main_frame, text="Chain Selection", padding=10)
        chain_frame.pack(fill=tk.X, pady=5)

        self.chain_mode_var = tk.StringVar(value="all")
        ttk.Radiobutton(chain_frame, text="All chains (multimeric analysis)",
                        variable=self.chain_mode_var, value="all").pack(
            side=tk.LEFT, padx=10)
        ttk.Radiobutton(chain_frame, text="Specific chain:",
                        variable=self.chain_mode_var, value="specific").pack(
            side=tk.LEFT, padx=5)
        self.chain_id_var = tk.StringVar(value="A")
        ttk.Entry(chain_frame, textvariable=self.chain_id_var, width=5).pack(
            side=tk.LEFT, padx=5)

        # ── Options ───────────────────────────────────────────────────────
        options_frame = ttk.LabelFrame(main_frame, text="Options", padding=10)
        options_frame.pack(fill=tk.X, pady=5)

        opt_grid = ttk.Frame(options_frame)
        opt_grid.pack(fill=tk.X)

        ttk.Label(opt_grid, text="Number of mutations:").grid(row=0, column=0, sticky=tk.W)
        self.mutations_var = tk.IntVar(value=10)
        ttk.Spinbox(opt_grid, textvariable=self.mutations_var, from_=1, to=50,
                     width=5).grid(row=0, column=1, padx=10)

        ttk.Label(opt_grid, text="Output directory:").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.output_var = tk.StringVar(value="./output")
        ttk.Entry(opt_grid, textvariable=self.output_var, width=50).grid(
            row=1, column=1, padx=5)
        ttk.Button(opt_grid, text="Browse", command=self._browse_output).grid(
            row=1, column=2, padx=5)

        # Validation checkboxes
        self.validate_foldx_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(opt_grid, text="Run FoldX ΔΔG validation",
                        variable=self.validate_foldx_var).grid(
            row=2, column=0, columnspan=2, sticky=tk.W, pady=2)

        self.validate_rosetta_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(opt_grid, text="Run Rosetta ΔΔG validation",
                        variable=self.validate_rosetta_var).grid(
            row=3, column=0, columnspan=2, sticky=tk.W, pady=2)

        # ── Action Button ─────────────────────────────────────────────────
        run_frame = ttk.Frame(main_frame)
        run_frame.pack(pady=15)
        self.run_btn = ttk.Button(run_frame, text="▶  Run Analysis",
                                  command=self._run_analysis, style="Accent.TButton")
        self.run_btn.pack()

        # ── Progress / Log ────────────────────────────────────────────────
        log_frame = ttk.LabelFrame(main_frame, text="Log", padding=5)
        log_frame.pack(fill=tk.BOTH, expand=True)

        self.log_text = scrolledtext.ScrolledText(log_frame, height=15,
                                                   font=("Courier", 9))
        self.log_text.pack(fill=tk.BOTH, expand=True)
        self.log_text.config(state=tk.DISABLED)

    def _browse_pdb(self):
        path = filedialog.askopenfilename(
            title="Select PDB file",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
        )
        if path:
            self.pdb_path_var.set(path)

    def _browse_output(self):
        path = filedialog.askdirectory(title="Select output directory")
        if path:
            self.output_var.set(path)

    def _log(self, msg):
        """Append a message to the log window."""
        self.log_text.config(state=tk.NORMAL)
        self.log_text.insert(tk.END, msg + "\n")
        self.log_text.see(tk.END)
        self.log_text.config(state=tk.DISABLED)
        self.root.update_idletasks()

    def _run_analysis(self):
        """Run the full analysis pipeline in a background thread."""
        self.run_btn.config(state=tk.DISABLED)
        self.log_text.config(state=tk.NORMAL)
        self.log_text.delete("1.0", tk.END)
        self.log_text.config(state=tk.DISABLED)

        thread = threading.Thread(target=self._analysis_worker)
        thread.daemon = True
        thread.start()

    def _analysis_worker(self):
        """Worker thread that executes the analysis pipeline."""
        try:
            from core.structure import load_structure, fetch_alphafold_pdb
            from core.analysis import (
                get_residues_all_chains, get_residues_single_chain,
                identify_catalytic_residues, compute_residue_features,
                extract_global_features, compute_interchain_contacts,
                get_chains_info,
            )
            from core.ml_model import build_training_dataset, train_activity_model, predict_activity
            from core.scoring import suggest_mutations
            from core.visualization import (
                make_3d_html, make_residue_properties_html, make_mutation_score_png,
            )
            from core.report import write_report
            from validation.ddg_validation import run_full_validation
        except ImportError as e:
            self._log(f"[ERROR] Failed to import modules: {e}")
            self.run_btn.config(state=tk.NORMAL)
            return

        self._log("[START] BioEnzyme Designer analysis starting...")

        # Step 1: Get structure
        output_dir = Path(self.output_var.get())
        output_dir.mkdir(parents=True, exist_ok=True)

        try:
            if self.source_var.get() == "pdb":
                pdb_path = Path(self.pdb_path_var.get())
                if not pdb_path.exists():
                    raise FileNotFoundError(f"PDB file not found: {pdb_path}")
            else:
                seq = self.seq_entry.get().strip()
                if not seq:
                    raise ValueError("Please enter a sequence.")
                pdb_path = fetch_alphafold_pdb(seq, output_dir)
            self._log(f"[OK] Structure loaded: {pdb_path}")
        except Exception as e:
            self._log(f"[ERROR] {e}")
            self.run_btn.config(state=tk.NORMAL)
            return

        # Step 2: Load and analyse
        try:
            structure, model = load_structure(pdb_path)
            chains_info = get_chains_info(model)

            if self.chain_mode_var.get() == "all":
                residues_with_chain = get_residues_all_chains(model)
            else:
                residues_with_chain = get_residues_single_chain(
                    model, chain_id=self.chain_id_var.get())

            self._log(f"[OK] {len(residues_with_chain)} residues analysed "
                      f"({model[0] and len(list(model.get_chains()))} chains)")
        except Exception as e:
            self._log(f"[ERROR] Structure analysis failed: {e}")
            self.run_btn.config(state=tk.NORMAL)
            return

        # Step 3: Catalytic identification
        catalytic_ids = identify_catalytic_residues(residues_with_chain)
        self._log(f"[OK] {len(catalytic_ids)} catalytic residues identified")

        df_residues = compute_residue_features(residues_with_chain, catalytic_ids)
        baseline_feats = extract_global_features(df_residues)

        # Step 4: ML model
        try:
            self._log("[INFO] Building training dataset...")
            train_df = build_training_dataset(
                data_source=self.data_source_var.get(),
                ec_number=self.ec_var.get() or None,
                organism=self.organism_var.get() or None,
            )
            model_rf, scaler = train_activity_model(train_df)
            baseline_activity = predict_activity(model_rf, scaler, baseline_feats)
            self._log(f"[OK] Baseline log10(kcat/Km) = {baseline_activity:.3f}")
        except Exception as e:
            self._log(f"[ERROR] ML model failed: {e}")
            self.run_btn.config(state=tk.NORMAL)
            return

        # Step 5: Inter-chain contacts
        interchain_contacts = compute_interchain_contacts(model)
        self._log(f"[OK] {len(interchain_contacts)} inter-chain contacts found")

        # Step 6: Mutation scoring
        top_mutations = suggest_mutations(
            df_residues, self.improve_var.get(),
            self.mutations_var.get(), interchain_contacts,
        )
        self._log(f"[OK] Top {len(top_mutations)} mutations suggested")

        # Step 7: Validation
        validation_results = {}
        if self.validate_foldx_var.get() or self.validate_rosetta_var.get():
            self._log("[INFO] Running ΔΔG validation...")
            validation_results = run_full_validation(
                pdb_path, top_mutations, output_dir,
                foldx_binary=None if self.validate_foldx_var.get() else "no_foldx",
                rosetta_dir=None if self.validate_rosetta_var.get() else "no_rosetta",
            )
            self._log(f"[OK] Validation complete. Correlation: "
                      f"{validation_results.get('correlation', 'N/A')}")

        # Step 8: Visualisations
        try:
            viewer_path = make_3d_html(pdb_path, top_mutations, output_dir)
            props_path = make_residue_properties_html(df_residues, top_mutations, output_dir)
            png_path = make_mutation_score_png(top_mutations, baseline_activity, output_dir)
            self._log("[OK] Visualisations generated")
        except Exception as e:
            self._log(f"[WARN] Some visualisations failed: {e}")
            viewer_path = props_path = png_path = output_dir

        # Step 9: Report
        paths = {
            "viewer": viewer_path,
            "properties": props_path,
            "png": png_path,
            "report": output_dir / "report.txt",
        }
        write_report(
            output_dir, pdb_path, self.improve_var.get(),
            df_residues, catalytic_ids, baseline_feats, baseline_activity,
            top_mutations, paths, validation_results, chains_info,
        )
        self._log("[DONE] Analysis complete! All outputs in: " + str(output_dir))
        self._log("")
        self._log("═══ Top Mutations ═══")
        for i, m in enumerate(top_mutations, 1):
            self._log(f"  {i}. {m['mutation']} score={m['score']:.3f}")

        self.run_btn.config(state=tk.NORMAL)

    def _show_results_dialog(self, output_dir):
        """Show a dialog with top results."""
        dialog = tk.Toplevel(self.root)
        dialog.title("Analysis Results")
        dialog.geometry("600x400")

        report_path = output_dir / "report.txt"
        if report_path.exists():
            text = report_path.read_text()
            text_widget = scrolledtext.ScrolledText(dialog, wrap=tk.WORD)
            text_widget.pack(fill=tk.BOTH, expand=True)
            text_widget.insert("1.0", text)
            text_widget.config(state=tk.DISABLED)
        else:
            ttk.Label(dialog, text="No results available.").pack()


def launch_gui():
    """Launch the Tkinter GUI."""
    root = tk.Tk()
    app = BioEnzymeGUI(root)
    root.mainloop()


if __name__ == "__main__":
    launch_gui()
