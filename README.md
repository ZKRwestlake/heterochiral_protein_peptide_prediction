# Structure prediction of heterochiral protein‑peptide interactions with Chai‑1

This repository contains the scripts and demo samples for the manuscript  
**"Structure prediction of heterochiral protein‑peptide interactions with Chai‑1"**.

It has two parts:

- **Prediction and analysis** — `chai-1/` and `af3/`, with a small demo dataset and the
  DockQ-based analysis (Section 4).
- **Binder design** — `design_protocol_scripts/`, the Rosetta/RifDock protocol for designing
  heterochiral protein–peptide complexes (Section 5).

---

## 1. System Requirements

### Operating System
- Red Hat Enterprise Linux 8.5.0 (tested)
- Other modern Linux distributions are expected to work, but have not been formally validated.

### Job Scheduler & Container
- [SLURM](https://slurm.schedmd.com/) workload manager
- [Singularity](https://sylabs.io/singularity/) (for running AlphaFold 3)

### Hardware
- NVIDIA GPU (e.g., L40 or A40) – required for both Chai‑1 and AlphaFold 3.

---

## 2. Software Dependencies

The following external software must be installed and accessible.  
Follow each tool’s official documentation for installation.

| Software            | Version | Link |
|---------------------|---------|------|
| Chai‑1              | 0.5.2   | [chaidiscovery/chai-lab](https://github.com/chaidiscovery/chai-lab) |
| AlphaFold 3         | 3.0.1   | [google-deepmind/alphafold3](https://github.com/google-deepmind/alphafold3) |
| Rosetta Modeling Suite | Any (for relax) | [RosettaCommons](https://www.rosettacommons.org/) |

All steps use a single conda environment, referred to below as `CHAI_ENV`. Chai‑1 requires
Python ≥ 3.10. Install Chai‑1 in it, plus the packages the analysis steps need:

```bash
conda create -n chai1 python=3.11
conda activate chai1
pip install chai-lab==0.5.2    # Chai-1
pip install scipy              # DockQ / RMSD analysis steps
```

On some systems `scipy` additionally needs the environment's `libstdc++` on
`LD_LIBRARY_PATH`; the analysis scripts export it automatically.

**Important**: After installing Rosetta, note the path to the `rosetta_scripts` binary and its suffix (e.g., `.default.linuxgccrelease`).  
You will need to update these in the relaxation SLURM scripts.

The binder design protocol has its own, larger set of dependencies (a pinned Rosetta release, RifDock, PatchDock, PyRosetta and others); they are listed in Section 5.

---

## 3. Installation Guide

No additional compilation is required. Simply clone the repository and ensure the dependencies above are available on your system.

```bash
git clone https://github.com/ZKRwestlake/heterochiral_protein_peptide_prediction.git
cd heterochiral_protein_peptide_prediction
```

**Typical install time**: The only setup is downloading/cloning the repository (< 1 minute). Installation of Chai‑1, AlphaFold 3, and Rosetta should follow their respective guides; these are not part of this repository.

---

## 4. Demo (Run on Provided Samples)

Each tool is run through four numbered steps, in this order:

| Step | Chai‑1 | AlphaFold 3 | Purpose |
|------|--------|-------------|---------|
| 1 | `1.predict.slurm` | `1.af3.sh` | Predict 5 seeds × 5 diffusion samples per target |
| 2 | `2.relax.slurm` | `2.relax.slurm` | Relax every prediction with Rosetta |
| 3 | `3.collect_dockq.sh` | `3.collect_dockq.sh` | Compute RMSD, DockQ and interface metrics per model |
| 4 | `4.best_scored.sh` | `4.best_scored.sh` | Collect the best-scored and most-confident model per target |

Steps 3 and 4 compare each model against a reference structure in the repository-root
directory `ground_truth/<TARGET>.pdb`, so that directory is required to run them
(see [Reference structures](#reference-structures)).

---

### 4.1 Chai‑1 Demo

**Step 1 – Predict**

1. **Navigate** to the Chai‑1 directory:
   ```bash
   cd chai-1/
   ```
   Here you will find subdirectories containing MSA and FASTA files.
2. **Edit the prediction SLURM script** `1.predict.slurm`:
   - Replace `CHAI_ENV` with the name of your Chai‑1 conda environment.
   - Change the array job parameter `-a` to select which sample indices (listed in `dirs.list`) you wish to run.
3. **Run prediction**:
   ```bash
   sbatch 1.predict.slurm
   ```
   - Each prediction runs with 5 seeds, each generating 5 diffusion samples.
   - **Expected run time** (per sample, on a single L40 GPU): from a few minutes up to ~1 hour, depending on sequence length.
   - **Expected output**: results stored under `prediction/*/outputs/*/`.

**Step 2 – Relax**

1. **Prepare the output list**:
   ```bash
   ls -d prediction/*/outputs/*/ > outputs.list
   ```
2. **Edit** `2.relax.slurm`:
   - Set `$rosettabin` to your Rosetta binary directory (e.g., `/home/software/rosetta.source.release/main/source/bin/`).
   - Set `$rosettasuffix` to the appropriate suffix (e.g., `.default.linuxgccrelease`).
   - If any input contains the non-canonical residue `D0C`, also export `$ROSETTAROOT`
     (its `main/database` is used as the fallback database).
   - Adjust the array job `-a` parameter to select predictions listed in `outputs.list`.
3. **Run relaxation**:
   ```bash
   sbatch 2.relax.slurm
   ```
   - **Expected run time**: from a few minutes up to ~1 hour per complex.
   - The project-local Rosetta chemistry for `D0C` (`rosetta_params/`) is detected and
     loaded automatically for every prediction that contains `D0C`.

**Step 3 – Per-model RMSD and DockQ**

1. **Edit** `3.collect_dockq.sh`: replace `CHAI_ENV` with the name of your Chai‑1 conda environment.
2. **Run**:
   ```bash
   sbatch 3.collect_dockq.sh
   ```
   - Output: `prediction/*/rmsd_confidence_dockq.txt` (one row per model, with `aligned_rmsd`,
     `deviation`, `dockq`, `f_nat`, `i_rmsd`, `l_rmsd`, `r_rmsd` and `iptm`) and
     `prediction/*/combined_score_dockq.sc` (the Rosetta `ddg` scores merged with the same metrics).

**Step 4 – Best scored / most confident model per target**

1. **Edit** `4.best_scored.sh` to point at the same `CHAI_ENV`.
2. **Run**:
   ```bash
   bash 4.best_scored.sh
   ```
   - Output: `best_scored_models/ddg/*.pdb` (lowest `ddg`) and
     `best_scored_models/iptm/*.pdb` (highest `iptm`), with the metric value encoded in the file name.

---

### 4.2 AlphaFold 3 Demo

**Step 1 – Predict**

1. **Navigate** to the AlphaFold 3 directory:
   ```bash
   cd af3/
   ```
   The folder `af_input` contains the input JSON files.
2. **(Optional)** To run a single sample instead of all, uncomment and edit the last line of
   `1.af3.sh`:
   ```
   --json_path=/root/af_input/[YOUR_JSON]
   ```
   instead of `--input_dir=/root/af_input`.
3. **Run prediction**:
   ```bash
   sbatch 1.af3.sh
   ```
   - Each prediction uses 5 seeds, each with 5 diffusion samples.
   - **Expected run time**: up to several hours for the full set on a single L40 GPU.
   - **Expected output**: results stored under `af_output/*/*/`.

**Step 2 – Relax**

1. **Prepare the output list**:
   ```bash
   ls -d af_output/*/*/ > outputs.list
   ```
2. **Edit** `2.relax.slurm` as in the Chai‑1 section (set `CHAI_ENV`, Rosetta paths and
   array job index). The AlphaFold 3 CIF is converted to PDB first, which needs `biopython`.
3. **Run relaxation**:
   ```bash
   sbatch 2.relax.slurm
   ```
   - **Expected run time**: from a few minutes to ~1 hour per complex.

**Steps 3 & 4 – DockQ and best models**

Same scripts as for Chai‑1, driven from the AlphaFold 3 output layout:

```bash
sbatch 3.collect_dockq.sh     # -> af_output/*/rmsd_confidence_dockq.txt
                              #    af_output/*/combined_score_dockq.sc
bash 4.best_scored.sh         # -> best_scored_models/{iptm,ddg}/*.pdb
```

### Reference structures

Steps 3 and 4 need one reference structure per target:

```
ground_truth/<TARGET>.pdb
```

`<TARGET>` is the upper-case target name (for example `2R3C`), matching the directory
name under `prediction/` and `af_output/`. The reference structures for the demo targets
are included in the repository.

---

## 5. Binder design protocol

`design_protocol_scripts/` holds the Rosetta/RifDock protocol used to design heterochiral
protein–peptide complexes, described in
[s41422-024-01014-2](https://www.nature.com/articles/s41422-024-01014-2). The pipeline for
other D-peptides is largely identical, differing mainly in the hotspot residues selected
and the filtering thresholds, which are target-dependent.

**The step-by-step instructions are in
`design_protocol_scripts/Design_ptotocol_README.txt`**, which is the reference for running
the protocol. The summary below is only a map of it.

### Pipeline

| Stage | What it does | Files |
|-------|--------------|-------|
| Target preparation | Build the D-chirality target and reduce it to a single chain (chain A) | `invert_chiral.py` |
| Scaffold preparation | Convert the scaffolds to poly-valine | `mutate_polyXXX.py` |
| Docking | rifgen, then PatchDock, then RifDock | `rifgen.flag`, `rifdock.flag` |
| Pre-handling | Rename residues to their D-forms and renumber the chains | `change_L_name_to_D_name.sh`, `renum_chain.pl` |
| Round 1 design | Predictor (discard poor docks), pilot runs, FastDesign | `predictor.xml`, `design.xml` (plus their `.flags`) |
| Motif grafting | Extract interface motifs, cluster them and pick hotspots, then graft the motifs onto scaffolds | `extract_PPI_motifs.py`, `motif_clustering/`, `notebook/look_at_motifs_by_length_w_hbond.ipynb`, `get_path_and_hotspot.py`, `graft.xml` |
| Round 2 design | Predictor, pilot runs and FastDesign on the grafted backbones | `predictor2.xml`, `design2.xml`, `constraints.comp` |
| Final filter | Score the designs and select the final set | `filter.xml`, `notebook/final_filtering.ipynb` |

Each Rosetta stage is driven by a matching `.flags` file; the input/output paths and the
thresholds in those files are what has to be edited for a new target. The notebooks under
`notebook/` are used alongside these stages: `look_at_motifs_by_length_w_hbond.ipynb` for the
motif analysis, and `predictor.ipynb` and `final_filtering.ipynb` for selecting which designs
advance from one round to the next.

### Additional requirements

| Software | Used for |
|----------|----------|
| Rosetta 2020.03 bundle | every `rosetta_scripts` stage (the release pinned by the protocol) |
| [RifDock](https://github.com/rifdock/rifdock) | rifgen and rifdock |
| PatchDock | initial scaffold docking |
| PyRosetta | `invert_chiral.py`, `mutate_polyXXX.py` |
| [silent_tools](https://github.com/bcov77/silent_tools) | reading and writing silent files |
| g++ | building `motif_clustering/` |
| Jupyter with numpy, scipy, pandas, scikit-learn, matplotlib and seaborn | the notebooks |

`motif_clustering/` is a small C++ helper that clusters the extracted motifs. A prebuilt
`cluster` binary is included; `motif_clustering/README` gives the build commands, and the
vendored `pybind11/` supplies the headers they need.

---

## 6. Instructions for Use (On Your Own Data)

1. **Prepare your input data** in the same format as the demo samples:
   - For **Chai‑1**: FASTA files and MSAs placed in a subdirectory under `chai-1/`, and the
     target directory listed in `dirs.list`.
   - For **AlphaFold 3**: JSON files placed in `af3/af_input/` (the D-amino-acid chemistry is
     embedded in each JSON as `userCCD`).
2. **Provide a reference structure** for each target as `ground_truth/<TARGET>.pdb` if you
   want to run steps 3 and 4.
3. **Update the SLURM scripts** with your environment name, Rosetta paths, and array job ranges.
4. **Run the four steps** exactly as shown in Section 4.

To design binders with the protocol of Section 5, follow
`design_protocol_scripts/Design_ptotocol_README.txt`.

For batch processing or using different compute clusters, adapt the `#SBATCH` directives in the scripts accordingly.

---

## 7. License & Third-party Software

### Code License
The scripts in this repository are distributed under the [MIT License](LICENSE).

### Third-party Software Notices
The workflows use external tools. Please adhere to their respective licenses.

---


*For questions or support, please open an issue on the repository or contact lupeilong@westlake.edu.cn.*
