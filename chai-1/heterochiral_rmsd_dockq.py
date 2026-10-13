import sys
import numpy as np
from Bio import Align
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Superimposer
from Bio.SVDSuperimposer import SVDSuperimposer
from scipy.optimize import linear_sum_assignment
from itertools import permutations, product

CANONICAL_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'CYX',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'ACE', 'NH2', 'MSE'
}
conversion_l = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'G': 'GLY', 'P': 'PRO', 'K': 'LYS',
                'T': 'THR', 'S': 'SER', 'N': 'ASN', 'Q': 'GLN', 'H': 'HIS', 'W': 'TRP', 'Y': 'TYR',
                'I': 'ILE', 'L': 'LEU', 'F': 'PHE', 'M': 'MET', 'R': 'ARG', 'V': 'VAL',
                }
inversion_l = {v: k for k, v in conversion_l.items()}
conversion_d = {'A': 'DAL', 'C': 'DCY', 'D': 'DAS', 'E': 'DGL', 'P': 'DPR', 'K': 'DLY',
                'T': 'DTH', 'S': 'DSN', 'N': 'DSG', 'Q': 'DGN', 'H': 'DHI', 'W': 'DTR', 'Y': 'DTY',
                'I': 'DIL', 'L': 'DLE', 'F': 'DPN', 'M': 'MED', 'R': 'DAR', 'V': 'DVA',
                }
additional_inversion = {'DCS': 'c', 'DGU': 'e', 'DSE': 's', 'DAN': 'n', 'DPH': 'f', 'DME': 'm', 'MSE': 'M', 'D0C': 'f'}
inversion_d = {v: k.lower() for k, v in conversion_d.items()}
inversion = {**inversion_l, **inversion_d, **additional_inversion}
# AlphaFold 3 writes CCD/ligand components in lower case (e.g. "gly", "dal"),
# so accept lower-case residue names as well.
inversion.update({k.lower(): v for k, v in inversion.items()})
alignment_map = {}

# DockQ constants
FNAT_THRESHOLD = 5.0
FNAT_THRESHOLD_PEPTIDE = 4.0
INTERFACE_THRESHOLD = 10.0
INTERFACE_THRESHOLD_PEPTIDE = 8.0
BACKBONE_ATOMS = ("CA", "C", "N", "O")

aligner = Align.PairwiseAligner(mode='local')
aligner.target_end_gap_score = 0.0
aligner.query_end_gap_score = 0.0
aligner.match_score = 1.0
aligner.mismatch_score = -1.0
aligner.gap_score = -1.0
aligner.open_gap_score = -1.0
aligner.extend_gap_score = -0.5


def parse_structure(file_path):
    if file_path.endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure('structure', file_path)


def is_canonical(chain):
    for residue in chain:
        if residue.get_resname() not in CANONICAL_AA:
            return False
    return True


def get_chain_sequence(chain):
    residues = []
    for residue in chain:
        resname = residue.get_resname()
        if resname in inversion:
            try:
                one_letter = inversion[resname]
                residues.append((residue.id, one_letter))
            except KeyError:
                continue
    residues = sorted(residues, key=lambda x: x[0][1])
    sequence = ''.join([r[1] for r in residues]) if residues else ''
    return sequence, residues


def align_sequences(seq_a, seq_b):
    if (seq_a, seq_b) in alignment_map:
        return alignment_map[(seq_a, seq_b)]
    if (seq_b, seq_a) in alignment_map:
        b_idx, a_idx, score = alignment_map[(seq_b, seq_a)]
        return a_idx, b_idx, score

    alignments = aligner.align(seq_a, seq_b)
    try:
        best_alignment = next(alignments)
        score = best_alignment.score
        aligned_target = best_alignment.aligned[0]
        aligned_query = best_alignment.aligned[1]
        a_idx, b_idx = [], []

        for target_seg, query_seg in zip(aligned_target, aligned_query):
            target_start, target_end = target_seg
            query_start, query_end = query_seg
            for offset in range(target_end - target_start):
                a_idx.append(target_start + offset)
                b_idx.append(query_start + offset)

        alignment_map[(seq_a, seq_b)] = a_idx, b_idx, score
        return a_idx, b_idx, score
    except StopIteration:
        return [], [], 0


def find_chain_correspondence(actual_chains, predicted_chains):
    if not actual_chains or not predicted_chains:
        return {}
    score_matrix = np.zeros((len(actual_chains), len(predicted_chains)))
    for i, a_chain in enumerate(actual_chains):
        a_seq, _ = get_chain_sequence(a_chain)
        for j, p_chain in enumerate(predicted_chains):
            p_seq, _ = get_chain_sequence(p_chain)
            _, _, score = align_sequences(a_seq, p_seq)
            score_matrix[i][j] = score
    row_ind, col_ind = linear_sum_assignment(-score_matrix)
    return {j: i for i, j in zip(row_ind, col_ind)}


def group_identical_chains(chains):
    seqs = [get_chain_sequence(chain)[0] for chain in chains]
    unique, group_indices = np.unique(seqs, return_inverse=True)
    return unique, group_indices


def get_ca_atoms_for_permutation(actual_chain, predicted_chain):
    a_seq, a_res = get_chain_sequence(actual_chain)
    p_seq, p_res = get_chain_sequence(predicted_chain)
    a_idx, p_idx, _ = align_sequences(a_seq, p_seq)
    atoms = []
    for ai, pi in zip(a_idx, p_idx):
        try:
            atoms.append((actual_chain[a_res[ai][0]]['CA'], predicted_chain[p_res[pi][0]]['CA']))
        except KeyError:
            continue
    return atoms


def build_groups_with_correspondence(actual_chains, predicted_chains):
    if not actual_chains or not predicted_chains:
        return []
    correspondence = find_chain_correspondence(actual_chains, predicted_chains)
    _, group_indices = group_identical_chains(predicted_chains)
    groups = []
    for group_idx in np.unique(group_indices):
        pred_indices = np.where(group_indices == group_idx)[0]
        pred_group = [predicted_chains[i] for i in pred_indices]
        actual_group = [actual_chains[correspondence[i]] for i in pred_indices]
        groups.append((actual_group, pred_group))
    return groups


def find_optimal_joint_permutation(align_groups, deviate_groups, align_weight=0.2, deviate_weight=0.8):
    def precompute_group_atoms(groups):
        group_perms = []
        for actual_chains, predicted_chains in groups:
            perms = []
            for perm in permutations(predicted_chains):
                perm_atoms = []
                for a_chain, p_chain in zip(actual_chains, perm):
                    perm_atoms.extend(get_ca_atoms_for_permutation(a_chain, p_chain))
                perms.append((perm, perm_atoms))
            group_perms.append(perms)
        return group_perms
    
    align_group_perms = precompute_group_atoms(align_groups)
    deviate_group_perms = precompute_group_atoms(deviate_groups)
    
    best_score = float('inf')
    best_align_perms = None
    best_deviate_perms = None
    best_superimposer = None
    best_rmsd = float('inf')
    best_deviation = float('inf')
    
    align_ranges = [range(len(p)) for p in align_group_perms]
    deviate_ranges = [range(len(p)) for p in deviate_group_perms]
    
    for align_indices in product(*align_ranges):
        total_align_atoms = []
        current_align_perms = []
        for g_idx, perm_idx in enumerate(align_indices):
            perm, atoms = align_group_perms[g_idx][perm_idx]
            current_align_perms.append(perm)
            total_align_atoms.extend(atoms)
        if not total_align_atoms:
            continue
        sup = Superimposer()
        sup.set_atoms(*zip(*total_align_atoms))
        # sup.run()
        rmsd = sup.rms
        
        for deviate_indices in product(*deviate_ranges):
            total_deviate_atoms = []
            current_deviate_perms = []
            for g_idx, perm_idx in enumerate(deviate_indices):
                perm, atoms = deviate_group_perms[g_idx][perm_idx]
                current_deviate_perms.append(perm)
                total_deviate_atoms.extend(atoms)
            if not total_deviate_atoms:
                continue
            deviation_sum = 0
            count = 0
            rot = sup.rotran[0]
            tran = sup.rotran[1]
            for a_atom, p_atom in total_deviate_atoms:
                p_coord_transformed = np.dot(p_atom.coord, rot) + tran
                deviation_sum += np.sum(np.square(a_atom.coord - p_coord_transformed))
                count += 1
            deviation = np.sqrt(deviation_sum / count) if count > 0 else 0
            combined_score = align_weight * rmsd + deviate_weight * deviation
            if combined_score < best_score:
                best_score = combined_score
                best_align_perms = current_align_perms
                best_deviate_perms = current_deviate_perms
                best_superimposer = sup
                best_rmsd = rmsd
                best_deviation = deviation
    
    if best_superimposer is None:
        return None, None, None, float('inf'), float('inf'), float('inf')
    return best_align_perms, best_deviate_perms, best_superimposer, best_score, best_rmsd, best_deviation


def is_peptide_chain(chain, threshold=30):
    return len([r for r in chain if r.get_resname() in inversion]) < threshold


def get_cb_or_ca(residue):
    if "CB" in residue:
        return residue["CB"]
    elif "CA" in residue:
        return residue["CA"]
    return None


def get_backbone_atoms_list(residue):
    atoms = []
    for atom_type in BACKBONE_ATOMS:
        if atom_type in residue:
            atoms.append(residue[atom_type])
    return atoms


def calculate_dockq(actual_align_chains, actual_deviate_chains,
                   pred_align_chains, pred_deviate_chains,
                   capri_peptide=False):
    """
    Calculate DockQ and return the receptor-based transformation for re-superposition.
    Returns: dockq, fnat, irms, lrms, rotation_matrix, translation_vector, rms
    """
    if not actual_align_chains or not actual_deviate_chains:
        print('No actual_align_chains or actual_deviate_chains')
        return -1.0, -1.0, -1.0, -1.0, None, None, -1.
    
    any_peptide = any(is_peptide_chain(c) for c in actual_align_chains + actual_deviate_chains)
    use_peptide_mode = capri_peptide or any_peptide
    
    # Different thresholds for peptide vs protein
    fnat_threshold = FNAT_THRESHOLD_PEPTIDE if use_peptide_mode else FNAT_THRESHOLD
    interface_threshold = INTERFACE_THRESHOLD_PEPTIDE if use_peptide_mode else INTERFACE_THRESHOLD
    
    def get_aligned_residue_map(actual_chain, pred_chain):
        a_seq, a_res = get_chain_sequence(actual_chain)
        p_seq, p_res = get_chain_sequence(pred_chain)
        if not a_seq or not p_seq:
            return {}
        a_idx, p_idx, _ = align_sequences(a_seq, p_seq)
        mapping = {}
        for ai, pi in zip(a_idx, p_idx):
            actual_res_id = a_res[ai][0]
            pred_res = pred_chain[p_res[pi][0]]
            mapping[actual_res_id] = pred_res
        return mapping
    
    align_maps = [get_aligned_residue_map(a, p) for a, p in zip(actual_align_chains, pred_align_chains)]
    deviate_maps = [get_aligned_residue_map(a, p) for a, p in zip(actual_deviate_chains, pred_deviate_chains)]
    
    def get_mapped_residue(actual_res):
        parent = actual_res.get_parent()
        if parent in actual_align_chains:
            idx = actual_align_chains.index(parent)
            return align_maps[idx].get(actual_res.id)
        elif parent in actual_deviate_chains:
            idx = actual_deviate_chains.index(parent)
            return deviate_maps[idx].get(actual_res.id)
        return None
    
    def get_protein_residues(chains):
        return [res for chain in chains for res in chain if res.get_resname() in inversion]
    
    actual_align_res = get_protein_residues(actual_align_chains)
    actual_deviate_res = get_protein_residues(actual_deviate_chains)
    
    if not actual_align_res or not actual_deviate_res:
        print('No actual_align_res or actual_deviate_res')
        return -1.0, -1.0, -1.0, -1.0, None, None, -1.
    
    def calc_residue_distances_all_atom(res_list1, res_list2):
        """Calculate minimum distance between any heavy atoms of residue pairs."""
        if not res_list1 or not res_list2:
            return np.array([]).reshape(0, 0)
        dists = np.zeros((len(res_list1), len(res_list2)))
        for i, res1 in enumerate(res_list1):
            coords1 = np.array([a.coord for a in res1.get_atoms() if a.element != 'H'])
            if len(coords1) == 0:
                dists[i, :] = 1e6
                continue
            for j, res2 in enumerate(res_list2):
                coords2 = np.array([a.coord for a in res2.get_atoms() if a.element != 'H'])
                if len(coords2) == 0:
                    dists[i, j] = 1e6
                    continue
                all_dists = np.sqrt(np.sum((coords1[:, None] - coords2[None, :])**2, axis=2))
                dists[i, j] = np.min(all_dists)
        return dists
    
    def calc_residue_distances_cb(res_list1, res_list2):
        """Calculate CB-CB (or CA-CA) distances for peptide mode interface detection."""
        if not res_list1 or not res_list2:
            return np.array([]).reshape(0, 0)
        dists = np.zeros((len(res_list1), len(res_list2)))
        for i, res1 in enumerate(res_list1):
            atom1 = get_cb_or_ca(res1)
            if atom1 is None:
                dists[i, :] = 1e6
                continue
            for j, res2 in enumerate(res_list2):
                atom2 = get_cb_or_ca(res2)
                if atom2 is None:
                    dists[i, j] = 1e6
                    continue
                dists[i, j] = np.linalg.norm(atom1.coord - atom2.coord)
        return dists
    
    # Calculate ALL-ATOM distances for contacts (Fnat) - always use heavy atoms regardless of peptide mode
    native_dists_contacts = calc_residue_distances_all_atom(actual_align_res, actual_deviate_res)
    if native_dists_contacts.size == 0:
        print('native_dists_contacts.size == 0')
        return -1.0, -1.0, -1.0, -1.0, None, None, -1.
    
    # Calculate interface distances - use CB for peptides, all-atom for proteins
    if use_peptide_mode:
        native_dists_interface = calc_residue_distances_cb(actual_align_res, actual_deviate_res)
    else:
        native_dists_interface = native_dists_contacts  # Same as contacts for proteins
    
    actual_align_filtered = [r for r in actual_align_res if get_mapped_residue(r) is not None]
    actual_deviate_filtered = [r for r in actual_deviate_res if get_mapped_residue(r) is not None]
    
    pred_align_filtered = [get_mapped_residue(r) for r in actual_align_filtered]
    pred_deviate_filtered = [get_mapped_residue(r) for r in actual_deviate_filtered]
    
    if not pred_align_filtered or not pred_deviate_filtered:
        print('No pred_align_filtered or pred_deviate_filtered')
        return -1.0, -1.0, -1.0, -1.0, None, None, -1.
    
    # Model distances for contacts (all-atom)
    model_dists_contacts = calc_residue_distances_all_atom(pred_align_filtered, pred_deviate_filtered)
    native_dists_contacts_filtered = calc_residue_distances_all_atom(actual_align_filtered, actual_deviate_filtered)
    
    if model_dists_contacts.shape != native_dists_contacts_filtered.shape:
        print('model_dists_contacts.shape != native_dists_contacts_filtered.shape')
        return -1.0, -1.0, -1.0, -1.0, None, None, -1.
    
    # Calculate Fnat using all-atom distances with appropriate threshold
    native_contacts = native_dists_contacts_filtered < fnat_threshold
    model_contacts = model_dists_contacts < fnat_threshold
    nat_total = np.sum(native_contacts)
    if nat_total == 0:
        print('nat_total == 0')
        return -1.0, -1.0, -1.0, -1.0, None, None, -1.
    
    nat_correct = np.sum(native_contacts & model_contacts)
    fnat = nat_correct / nat_total
    
    # Recalculate interface distances for filtered lists (using CB for peptides)
    if use_peptide_mode:
        native_dists_interface_filtered = calc_residue_distances_cb(actual_align_filtered, actual_deviate_filtered)
    else:
        native_dists_interface_filtered = native_dists_contacts_filtered
    
    # Interface detection using appropriate distance metric
    interface_align_mask = np.any(native_dists_interface_filtered < interface_threshold, axis=1)
    interface_deviate_mask = np.any(native_dists_interface_filtered < interface_threshold, axis=0)
    interface_align_indices = np.where(interface_align_mask)[0]
    interface_deviate_indices = np.where(interface_deviate_mask)[0]
    
    def collect_interface_atoms(res_list, indices):
        atoms = []
        for idx in indices:
            res = res_list[idx]
            atoms.extend(get_backbone_atoms_list(res))
        return np.array([a.coord for a in atoms]) if atoms else np.array([]).reshape(0, 3)
    
    def collect_mapped_interface_atoms(res_list, indices):
        atoms = []
        for idx in indices:
            actual_res = res_list[idx]
            mapped_res = get_mapped_residue(actual_res)
            if mapped_res is not None:
                atoms.extend(get_backbone_atoms_list(mapped_res))
        return np.array([a.coord for a in atoms]) if atoms else np.array([]).reshape(0, 3)
    
    ref_interface_atoms = np.vstack([
        collect_interface_atoms(actual_align_filtered, interface_align_indices),
        collect_interface_atoms(actual_deviate_filtered, interface_deviate_indices)
    ]) if len(interface_align_indices) > 0 and len(interface_deviate_indices) > 0 else np.array([]).reshape(0, 3)
    
    sample_interface_atoms = np.vstack([
        collect_mapped_interface_atoms(actual_align_filtered, interface_align_indices),
        collect_mapped_interface_atoms(actual_deviate_filtered, interface_deviate_indices)
    ]) if len(interface_align_indices) > 0 and len(interface_deviate_indices) > 0 else np.array([]).reshape(0, 3)
    
    if len(ref_interface_atoms) == 0 or len(sample_interface_atoms) == 0:
        print('len(ref_interface_atoms) == 0 or len(sample_interface_atoms) == 0')
        return -1.0, -1.0, -1.0, -1.0, None, None, -1.
    
    sup = SVDSuperimposer()
    sup.set(ref_interface_atoms, sample_interface_atoms)
    sup.run()
    irms = sup.get_rms()
    
    # Determine receptor/ligand by size
    align_size = len(actual_align_filtered)
    deviate_size = len(actual_deviate_filtered)
    
    if align_size >= deviate_size:
        receptor_actual = actual_align_filtered
        receptor_pred = pred_align_filtered
        ligand_actual = actual_deviate_filtered
        ligand_pred = pred_deviate_filtered
    else:
        receptor_actual = actual_deviate_filtered
        receptor_pred = pred_deviate_filtered
        ligand_actual = actual_align_filtered
        ligand_pred = pred_align_filtered
    
    def collect_all_backbone_coords(res_list):
        coords = []
        for res in res_list:
            coords.extend([a.coord for a in get_backbone_atoms_list(res)])
        return np.array(coords) if coords else np.array([]).reshape(0, 3)
    
    receptor_ref_coords = collect_all_backbone_coords(receptor_actual)
    receptor_sample_coords = collect_all_backbone_coords(receptor_pred)
    ligand_ref_coords = collect_all_backbone_coords(ligand_actual)
    ligand_sample_coords = collect_all_backbone_coords(ligand_pred)
    
    if len(receptor_ref_coords) == 0 or len(receptor_sample_coords) == 0:
        print('len(receptor_ref_coords) == 0 or len(receptor_sample_coords) == 0')
        return -1.0, -1.0, -1.0, -1.0, None, None, -1.
    
    # Receptor-based superposition
    sup_receptor = SVDSuperimposer()
    sup_receptor.set(receptor_ref_coords, receptor_sample_coords)
    sup_receptor.run()
    rot, tran = sup_receptor.get_rotran()
    rms = sup_receptor.get_rms()
    
    # LRMSD calculation
    if len(ligand_sample_coords) > 0:
        ligand_sample_rotated = np.dot(ligand_sample_coords, rot) + tran
        if len(ligand_ref_coords) > 0:
            lrms = np.sqrt(np.mean(np.sum((ligand_ref_coords - ligand_sample_rotated)**2, axis=1)))
        else:
            lrms = -1.0
    else:
        lrms = -1.0
    
    if irms < 0 or lrms < 0:
        dockq = -1.0
    else:
        dockq = (fnat + 
                 1.0 / (1.0 + (irms / 1.5)**2) + 
                 1.0 / (1.0 + (lrms / 8.5)**2)) / 3.0
    
    return dockq, fnat, irms, lrms, rot, tran, rms


def main(actual_file, predicted_file, output_pdb=None, align_canonical_chain=True, capri_peptide=False):
    actual = parse_structure(actual_file)
    predicted = parse_structure(predicted_file)
    actual_model = next(iter(actual))
    predicted_model = next(iter(predicted))

    actual_canonical = [c for c in actual_model if is_canonical(c)]
    actual_non_canonical = [c for c in actual_model if not is_canonical(c)]
    predicted_canonical = [c for c in predicted_model if is_canonical(c)]
    predicted_non_canonical = [c for c in predicted_model if not is_canonical(c)]

    if align_canonical_chain:
        align_actual = actual_canonical
        align_predicted = predicted_canonical
        deviate_actual = actual_non_canonical
        deviate_predicted = predicted_non_canonical
    else:
        align_actual = actual_non_canonical
        align_predicted = predicted_non_canonical
        deviate_actual = actual_canonical
        deviate_predicted = predicted_canonical

    correspondence = []
    dockq = f_nat = i_rmsd = l_rmsd = rms = -1.0
    
    align_groups = build_groups_with_correspondence(align_actual, align_predicted)
    deviate_groups = build_groups_with_correspondence(deviate_actual, deviate_predicted)

    # Store initial coordinates to apply DockQ transformation
    dockq_rot = None
    dockq_tran = None

    if align_groups and deviate_groups:
        # Step 1: Find optimal permutation using weighted score
        result = find_optimal_joint_permutation(
            align_groups, deviate_groups,
            align_weight=0.2, deviate_weight=0.8
        )
        align_perms, deviate_perms, superimposer, combined_score, rmsd, deviation = result
        
        if superimposer is not None:
            # Apply initial weighted-score superposition
            superimposer.apply(predicted_model.get_atoms())
            
            # Build flat lists for DockQ
            actual_aligned_flat = []
            actual_deviated_flat = []
            pred_aligned_flat = []
            pred_deviated_flat = []
            
            for (actual_group, _), perm in zip(align_groups, align_perms):
                for a_chain, p_chain in zip(actual_group, perm):
                    correspondence.append((a_chain.id, p_chain.id))
                    actual_aligned_flat.append(a_chain)
                    pred_aligned_flat.append(p_chain)
            
            for (actual_group, _), perm in zip(deviate_groups, deviate_perms):
                for a_chain, p_chain in zip(actual_group, perm):
                    correspondence.append((a_chain.id, p_chain.id))
                    actual_deviated_flat.append(a_chain)
                    pred_deviated_flat.append(p_chain)
            
            # Step 2: Calculate DockQ and get the receptor-based transformation
            dockq, f_nat, i_rmsd, l_rmsd, dockq_rot, dockq_tran, rms = calculate_dockq(
                actual_aligned_flat, actual_deviated_flat,
                pred_aligned_flat, pred_deviated_flat,
                capri_peptide=capri_peptide
            )
            
            # Step 3: Re-apply DockQ transformation (re-superpose based on receptor alignment)
            if dockq_rot is not None and dockq_tran is not None:
                for atom in predicted_model.get_atoms():
                    atom.coord = np.dot(atom.coord, dockq_rot) + dockq_tran
                # rmsd = rms  # Since we re-superposed, original RMSD no longer applies
        else:
            print(f'Cannot find valid superposition for {actual_file}, {predicted_file}.')
            rmsd = deviation = -1.
    else:
        print(f'Insufficient chains for joint optimization in {actual_file}, {predicted_file}.')
        rmsd = deviation = -1.

    if output_pdb is not None:
        io = PDBIO()
        io.set_structure(predicted)
        io.save(output_pdb)

    return rmsd, deviation, dockq, f_nat, i_rmsd, l_rmsd, rms


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py actual.pdb predicted.pdb output.pdb")
        sys.exit(1)
    
    rmsd, deviation, dockq, f_nat, i_rmsd, l_rmsd, rms = main(sys.argv[1], sys.argv[2], sys.argv[3])
    
    output_text = f"""superposed-rmsd deviation dockq f_nat i_rmsd l_rmsd r_rmsd
{rmsd:.3f} {deviation:.3f} {dockq:.3f} {f_nat:.3f} {i_rmsd:.3f} {l_rmsd:.3f} {rms:.3f}"""
    
    with open(f"{sys.argv[3][:-4]}_results.txt", 'w') as f:
        f.write(output_text)
    
    print(f"Results saved. RMSD: {rmsd:.3f}, Deviation: {deviation:.3f}, DockQ: {dockq:.3f}")
