# Accurate structure prediction of heterochiral protein‑peptide interactions by Chai‑1

This repository contains the main scripts and demo samples for the manuscript  
**"Accurate structure prediction of heterochiral protein‑peptide interactions by Chai‑1"**.

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

**Important**: After installing Rosetta, note the path to the `rosetta_scripts` binary and its suffix (e.g., `.default.linuxgccrelease`).  
You will need to update these in the relaxation SLURM scripts.

---

## 3. Installation Guide

No additional compilation is required. Simply clone the repository and ensure the dependencies above are available on your system.

```bash
git clone ZKRwestlake/heterochiral_protein_peptide_prediction
cd heterochiral_protein_peptide_prediction
```

**Typical install time**: The only setup is downloading/cloning the repository (< 1 minute). Installation of Chai‑1, AlphaFold 3, and Rosetta should follow their respective guides; these are not part of this repository.

---

## 4. Demo (Run on Provided Samples)

This section walks you through running the prediction and relaxation pipeline on the included demo dataset.

---

### 4.1 Chai‑1 Demo

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

4. **Prepare relaxation**:
   ```bash
   ls -d prediction/*/outputs/*/ > outputs.list
   ```
   Edit `2.relax.slurm`:
   - Set `$rosettabin` to your Rosetta binary directory (e.g., `/home/software/rosetta.source.release/main/source/bin/`).
   - Set `$rosettasuffix` to the appropriate suffix (e.g., `.default.linuxgccrelease`).
   - Adjust the array job `-a` parameter to select predictions listed in `outputs.list`.

5. **Run relaxation**:
   ```bash
   sbatch 2.relax.slurm
   ```
   - **Expected run time**: from a few minutes up to ~1 hour per complex.

---

### 4.2 AlphaFold 3 Demo

1. **Navigate** to the AlphaFold 3 directory:
   ```bash
   cd af3/
   ```
   The folder `af_input` contains the input JSON files.

2. **(Optional)** To run a single sample instead of all, edit `1.af3.sh` and replace:
   ```
   --input_dir=/root/af_input
   ```
   with:
   ```
   --json_path=/root/af_input/[YOUR_JSON]
   ```

3. **Run prediction**:
   ```bash
   sbatch 1.af3.sh
   ```
   - Each prediction uses 5 seeds, each with 5 diffusion samples.
   - **Expected run time**: up to several hours for the full set on a single L40 GPU.
   - **Expected output**: results stored under `af_output/*/*/`.

4. **Prepare relaxation**:
   ```bash
   ls -d af_output/*/*/ > outputs.list
   ```
   Edit `2.relax.slurm` as in the Chai‑1 section (set Rosetta paths and array job index).

5. **Run relaxation**:
   ```bash
   sbatch 2.relax.slurm
   ```
   - **Expected run time**: from a few minutes to ~1 hour per complex.

---

## 5. Instructions for Use (On Your Own Data)

1. **Prepare your input data** in the same format as the demo samples:
   - For **Chai‑1**: FASTA files and MSAs placed in a subdirectory under `chai-1/`.
   - For **AlphaFold 3**: JSON files placed in `af3/af_input/`.

2. **Update the SLURM scripts** with your environment names, Rosetta paths, and array job ranges.

3. **Run the prediction and relaxation steps** exactly as shown in Section 4.

For batch processing or using different compute clusters, adapt the `#SBATCH` directives in the scripts accordingly.

---

## 6. License & Third-party Software

### Code License
The scripts in this repository are distributed under the [MIT License](LICENSE).

### Third-party Software Notices
The workflows use external tools. Please adhere to their respective licenses.

---


*For questions or support, please open an issue on the repository or contact lupeilong@westlake.edu.cn.*