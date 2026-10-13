import os
from sys import argv
from glob import glob
import numpy as np
from heterochiral_rmsd_dockq import main as compute
import multiprocessing as mp
from functools import partial

def process_cif(cif, ref, align_canonical_chain, keys):
    """Process a single CIF file to compute RMSD, deviation, DockQ, and confidence scores."""
    directory = os.path.dirname(cif)
    basename = os.path.basename(cif)
    npz_path = f"{directory}/scores.{basename.split('.')[-2]}.npz"
    npz = np.load(npz_path)
    confidence = [str(npz[k].item()) for k in keys]
    rmsd, deviation, dockq, f_nat, i_rmsd, l_rmsd, rms = compute(
        ref, cif, output_pdb=f"{cif.rstrip('.cif')}_aligned_dockq.pdb",
        align_canonical_chain=align_canonical_chain
    )
    return [directory, basename, str(rmsd), str(deviation), str(dockq),
            str(f_nat), str(i_rmsd), str(l_rmsd), str(rms)] + confidence

def main():
    in_dir = argv[1].rstrip('/')
    name = os.path.basename(in_dir)
    align_canonical_chain = True
    ref = f'../ground_truth/{name.upper()}.pdb'
    if not os.path.isfile(ref):
        print(f'{ref} not found.')
        return

    cifs = sorted(glob(f'{in_dir}/outputs/*/*.cif'))
    if not cifs:
        print(f'No cif files found under {in_dir}/outputs.')
        return
    keys = ['iptm']

    # Set up multiprocessing
    num_processes = 5
    process_func = partial(process_cif, ref=ref,
                          align_canonical_chain=align_canonical_chain, keys=keys)

    # Process CIFs in parallel
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_func, cifs)

    # Prepare and save output
    output = np.array(results, dtype='<U100')
    title = 'directory description aligned_rmsd deviation dockq f_nat i_rmsd l_rmsd r_rmsd iptm'
    np.savetxt(f'{in_dir}/rmsd_confidence_dockq.txt', output, fmt='%s', header=title, comments='')

if __name__ == '__main__':
    main()
