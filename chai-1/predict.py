from sys import argv
import logging
from pathlib import Path

from chai_lab.chai1 import run_inference

logging.basicConfig(level=logging.INFO)

n_cycle = 10
n_step = 200

dirname = argv[1]
output_dir = Path(dirname)
fasta_path = output_dir.joinpath("example.fasta")

for seed in range(42, 47):
    output_subdir = output_dir.joinpath(f'outputs/{seed}')
    candidates = run_inference(
        fasta_file=fasta_path,
        output_dir=output_subdir,
        constraint_path=None,
        msa_directory=Path(f"alignment"),
        use_msa_server=False,
        # 'default' setup
        num_trunk_recycles=n_cycle,
        num_diffn_timesteps=n_step,
        seed=seed,
        device="cuda",
        use_esm_embeddings=False,
    )
