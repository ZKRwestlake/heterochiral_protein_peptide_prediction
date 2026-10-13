import os
from sys import argv
from glob import glob
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from heterochiral_rmsd_dockq import main as compute_rmsd

n_cores = 5
# AlphaFold 3 stores one relaxed model per seed/sample directory
n_models = 1


def read_score_table(score_path):
    """Read a Rosetta score file.

    The header is written from the first structure and filter values outside the
    reporting range are omitted, so a row can be wider or narrower than the header.
    Recover the description (always the final field) for every model, and read ddg
    only when the row matches the header width.
    """
    with open(score_path) as fh:
        score_lines = [line for line in fh if line.startswith('SCORE:')]
    header = score_lines[0].split()[1:]
    ddg_index = header.index('ddg') if 'ddg' in header else None
    records = []
    for line in score_lines[1:]:
        tokens = line.split()[1:]
        ddg = np.nan
        if ddg_index is not None and len(tokens) == len(header):
            try:
                ddg = float(tokens[ddg_index])
            except ValueError:
                pass
        records.append({'ddg': ddg, 'description': tokens[-1]})
    return pd.DataFrame(records)

def process_directory(directory, idx, ref, align_canonical_chain):
    """Process a single directory, handling variable numbers of models"""
    result = {
        'scores': {},
        'rmsds': {},
        'model_count': 0,  # Track actual number of models processed
    }
    
    score_path = glob(os.path.join(directory, '*score.sc'))
    if not score_path:
        print(f'No score file found in {directory}')
        return idx, result
    
    score_path = score_path[0]
    
    # Process score.sc
    if os.path.isfile(score_path):
        try:
            score_df = read_score_table(score_path)
            
            # Handle variable number of models
            actual_n_models = len(score_df)
            result['model_count'] = actual_n_models
            
            if actual_n_models != n_models:
                print(f'Warning: {directory} has {actual_n_models} models (expected {n_models})')

            n_missing_ddg = int(score_df['ddg'].isna().sum())
            if n_missing_ddg:
                print(f'Warning: ddg unavailable for {n_missing_ddg} model(s) in {directory}')
            
            # Store all available data
            for term in score_df.columns:
                result['scores'][term] = score_df[term]
            
            # Compute RMSD and DockQ for available models
            pdbs = [f'{directory}/{i}.pdb' for i in score_df['description']]
            dockq_results = []
            for pdb in pdbs:
                if os.path.exists(pdb):
                    try:
                        rmsd_deviation_dockq = compute_rmsd(
                            ref, pdb,
                            output_pdb=f"{pdb.rstrip('.pdb')}_aligned_dockq.pdb",
                            align_canonical_chain=align_canonical_chain
                        )
                        dockq_results.append(rmsd_deviation_dockq)
                    except Exception as e:
                        print(f'Error computing RMSD for {pdb}: {e}')
                        dockq_results.append([np.nan] * 7)
                else:
                    print(f'PDB file not found: {pdb}')
                    dockq_results.append([np.nan] * 7)
            
            if dockq_results:
                rmsd_deviation_dockq = np.array(dockq_results, dtype=float)
                result['rmsds']['aligned_rmsd'] = rmsd_deviation_dockq[:, 0]
                result['rmsds']['deviation'] = rmsd_deviation_dockq[:, 1]
                result['rmsds']['dockq'] = rmsd_deviation_dockq[:, 2]
                result['rmsds']['f_nat'] = rmsd_deviation_dockq[:, 3]
                result['rmsds']['i_rmsd'] = rmsd_deviation_dockq[:, 4]
                result['rmsds']['l_rmsd'] = rmsd_deviation_dockq[:, 5]
                result['rmsds']['r_rmsd'] = rmsd_deviation_dockq[:, 6]

        except Exception as e:
            print(f'Error processing {score_path}: {e}')

    return idx, result

def process_binder_directory(directory, idx):
    """Process binder directory with optimized array handling"""
    result = {}
    
    try:
        npz_files = sorted([
            os.path.join(directory, f)
            for f in os.listdir(directory)
            if f.endswith('.npz')
        ])
        n_files = len(npz_files)
        
        if n_files > 0:
            # Preallocate array
            ptm_vals = np.empty(n_files, dtype=np.float64)
            
            # Direct memory assignment
            for i, f in enumerate(npz_files):
                with np.load(f) as data:
                    ptm_vals[i] = data['ptm']
            
            result['binder_ptm'] = ptm_vals
    except Exception as e:
        print(f'Error processing binder directory {directory}: {e}')
    
    return idx, result

def main():
    in_dir = argv[1].rstrip('/')
    # AlphaFold 3 writes one directory per seed/sample
    dirs = sorted(glob(f'{in_dir}/*/'))
    if not dirs:
        print(f'No sample directories found under {in_dir}.')
        return
    n_out = len(dirs)
    align_canonical_chain = True

    name = os.path.basename(in_dir).rstrip('/')
    ref = f'../ground_truth/{name.upper()}.pdb'
    if not os.path.isfile(ref):
        print(f'{ref} not found.')
        return

    # Use lists instead of preallocated arrays to handle variable model counts
    scores_data = {
        'ddg': [],
        'description': [],
        'aligned_rmsd': [],
        'deviation': [],
        'dockq': [],
        'f_nat': [],
        'i_rmsd': [],
        'l_rmsd': [],
        'r_rmsd': [],
        'directory': [],
        'model_index': []  # Track which model this is within its directory
    }

    # Track model counts per directory for debugging
    model_counts = {}

    # Parallel processing
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        futures = []
        for idx, dir in enumerate(dirs):
            futures.append(executor.submit(
                partial(process_directory), dir, idx, ref, align_canonical_chain
            ))
        
        for future in as_completed(futures):
            idx, res = future.result()
            directory = dirs[idx]
            model_counts[directory] = res['model_count']
            
            # Handle case where no models were processed
            if res['model_count'] == 0:
                print(f'Warning: No models processed for {directory}')
                continue
            
            # Append data for each model in this directory
            n_models_actual = res['model_count']
            for model_idx in range(n_models_actual):
                scores_data['directory'].append(directory)
                scores_data['model_index'].append(model_idx)
                
                # Add score data if available
                for term in ['ddg', 'description']:
                    if term in res['scores'] and model_idx < len(res['scores'][term]):
                        scores_data[term].append(res['scores'][term].iloc[model_idx])
                    else:
                        scores_data[term].append(np.nan)
                
                # Add RMSD/DockQ data if available
                for rmsd_term in ['aligned_rmsd', 'deviation', 'dockq', 'f_nat', 'i_rmsd', 'l_rmsd', 'r_rmsd']:
                    if rmsd_term in res['rmsds'] and model_idx < len(res['rmsds'][rmsd_term]):
                        scores_data[rmsd_term].append(res['rmsds'][rmsd_term][model_idx])
                    else:
                        scores_data[rmsd_term].append(np.nan)

    # Create DataFrame
    df = pd.DataFrame(scores_data)
    
    # Print summary of model counts
    print("\nModel count summary per directory:")
    for directory, count in model_counts.items():
        print(f"  {os.path.basename(directory)}: {count} models")
    
    total_models = len(df)
    expected_models = n_out * n_models
    print(f"\nTotal models processed: {total_models} (expected: {expected_models})")
    
    # Save results
    output_file = f'{in_dir}/combined_score_dockq.sc'
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

if __name__ == '__main__':
    main()
