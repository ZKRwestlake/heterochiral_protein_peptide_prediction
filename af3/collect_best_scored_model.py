from sys import argv
import os
import pandas as pd

in_dir = argv[1].rstrip('/')
out_dir = 'best_scored_models'
rank_by, higher_better = 'ddg', False
score_file = f'{in_dir}/combined_score_dockq.sc'
if not os.path.isfile(score_file):
    print(f'{score_file} not found.')
    raise SystemExit(0)
df = pd.read_csv(score_file)
if rank_by in df and 'aligned_rmsd' in df and 'deviation' in df and len(df) and df[rank_by].notna().any():
    if higher_better:
        best_scored = df.iloc[df[rank_by].idxmax()]
    else:
        best_scored = df.iloc[df[rank_by].idxmin()]
    best_model = os.path.join(best_scored["directory"], f'{best_scored["description"]}_aligned_dockq.pdb')
    base_name = os.path.basename(os.path.dirname(os.path.dirname(best_scored['directory'])))
    new_name = (f"{base_name}_"
                f"{best_scored['description']}_"
                f"{rank_by}{best_scored[rank_by]}_"
                f"aligned{0.005+best_scored['aligned_rmsd']:.2f}_deviate{0.005+best_scored['deviation']:.2f}.pdb")

    os.makedirs(f'{out_dir}/{rank_by}', exist_ok=True)
    new_model = f'{out_dir}/{rank_by}/{new_name}'
    os.popen(f'cp {best_model} {new_model}')
