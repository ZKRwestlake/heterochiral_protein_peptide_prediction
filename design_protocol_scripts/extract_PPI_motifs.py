#!/usr/bin/env python

import os
import sys
import json
import argparse
import numpy as np
import itertools

import distutils.spawn
sys.path.append(os.path.dirname(distutils.spawn.find_executable("silent_tools.py")))
import silent_tools


from pyrosetta import *
from pyrosetta.rosetta import *

init("-corrections::beta_nov16 -holes:dalphaball /work/tlinsky/Rosetta/main/source/external/DAlpahBall/DAlphaBall.macgcc " + 
    "-use_truncated_termini -in:file:silent_struct_type binary")


script_dir = os.path.dirname(os.path.realpath(__file__))
xml = script_dir + "/extract_PPI_motifs_helper.xml"

# objs = protocols.rosetta_scripts.XmlObjects.create_from_file(xml)
# scorefxn = objs.get_score_function("sfxn")
# remove_polars = objs.get_mover("delete_polar")

scorefxn = get_fa_scorefxn()
hb_scorefxn = core.scoring.ScoreFunctionFactory.create_score_function("none")
hb_scorefxn.set_weight(core.scoring.hbond_sc, 1)
hb_scorefxn.set_weight(core.scoring.hbond_sr_bb, 1)
hb_scorefxn.set_weight(core.scoring.hbond_lr_bb, 1)
hb_scorefxn.set_weight(core.scoring.hbond_bb_sc, 1)


def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)
fix_scorefxn(hb_scorefxn, True)


def compute_filter(pose, filter, filtername):
    print("protocols.rosetta_scripts.ParsedProtocol.REPORT: ============Begin report for " + filtername + "=======================")
    value = filter.compute(pose)
    print("============End report for " + filtername + "=======================")
    extra_value = None

    # interface_sc returns a weird type
    if ( isinstance(value, pyrosetta.rosetta.core.scoring.sc._RESULTS)):
        extra_value = value.distance
        value = value.sc

    return value, extra_value


filters_to_apply = [
"interface_buried_sasa",
"ddg_norepack", 
# "ddg_hydrophobic",    // we have to do this manually because it doesn't have compute
"interface_sc",
"score_per_res",
#"buns_heavy_ball_1.1",
"vbuns5.5_heavy_ball_1.1",
"hydrophobic_residue_contacts",
"contact_molecular_surface"
]

def move_chainA_far_away(pose):
    pose = pose.clone()
    sel = core.select.residue_selector.ChainSelector("A")
    subset = sel.apply(pose)

    x_unit = numeric.xyzVector_double_t(1, 0, 0)
    far_away = numeric.xyzVector_double_t(10000, 0, 0)

    protocols.toolbox.pose_manipulation.rigid_body_move(x_unit, 0, far_away, pose, subset)

    return pose

abego_man = core.sequence.ABEGOManager()
def get_abego(pose, seqpos):
    return abego_man.index2symbol(abego_man.torsion2index_level1( pose.phi(seqpos), pose.psi(seqpos), pose.omega(seqpos)))



def better_dssp_hack(pose, length=-1):
    if ( length < 0 ):
        length = pose.size()

    dssp = core.scoring.dssp.Dssp(pose)
    dssp.dssp_reduced()
    the_dssp = dssp.get_dssp_secstruct()

    my_dssp = ""

    for seqpos in range(1, length+1):
        abego = get_abego(pose, seqpos)
        this_dssp = the_dssp[seqpos-1]
        if ( the_dssp[seqpos-1] == "H" and abego != "A" ):
            # print("!!!!!!!!!! Dssp - abego mismatch: %i %s %s !!!!!!!!!!!!!!!"%(seqpos, the_dssp[seqpos], abego))

            # This is the Helix-turn-helix HHHH case. See the test_scaffs folder
            if ( abego == "B" ):
                this_dssp = "L"

        my_dssp += this_dssp

    return my_dssp





parser = argparse.ArgumentParser()
parser.add_argument("-ref_pdb", type=str, default="", help="the reference pdb file")
parser.add_argument("-ddg_threshold", type=float, default=-15, help="the ddg cutoff value of the motif")
parser.add_argument("-pocket_threshold", type=float, default=0, help="the ddg cutoff value of specific residues")
parser.add_argument("-pocket_res", type=str, default="", help="the specific residues that make up the pocket")
parser.add_argument("-multi_segs", type=bool, default=False, help="dump multi segments or not if they are connected by a loop")
parser.add_argument("-out_prefix", type=str, default="", help="prefix on out files")
parser.add_argument("-motif_size_helix", type=int, default=-1, help="maximum size for helix motifs")
parser.add_argument("-motif_size_strand", type=int, default=-1, help="maximum size for strand motifs")
parser.add_argument("-hbond_atoms", type=str, default="",
# cd86    "103/OD1,85/NE2,85/OE1,101/NE2,34/NE2,34/OE1,89/NE2,89/ND1,38/OD1,38/ND2,36/O,38/O,37/O,40/N,53/O,43/OE2,45/OH", 
    help="maximum size for strand motifs")
parser.add_argument("pdbs", type=str, nargs="*", help="Inputs")
parser.add_argument("-in:file:silent", type=str, default="")




args = parser.parse_args(sys.argv[1:])

silent = args.__getattribute__("in:file:silent")


def get_ss_elements(dssp):
    ss_elements = []

    offset = 0
    ilabel = -1
    for label, group in itertools.groupby(dssp):
        ilabel += 1
        this_len = sum(1 for _ in group)
        next_offset = offset + this_len

        ss_elements.append( (label, offset, next_offset-1))

        offset = next_offset
    return ss_elements


# this will return all currently active TwoBodyEnergies
# this thing returns two AnalyticEtableEnergy which is correct.
# they get initialized with different flags for intra and inter
def find_all_2body_methods(scorefxn):
    nonzero = scorefxn.get_nonzero_weighted_scoretypes()

    to_find = {}
    for scoretype in nonzero:
        to_find[scoretype] = False

    sm = pyrosetta.rosetta.core.scoring.ScoringManager.get_instance()

    found = []

    for scoretype in to_find:
        if (to_find[scoretype]):
            continue

        method = sm.energy_method( scoretype, scorefxn.energy_method_options())
        if (isinstance(method, pyrosetta.rosetta.core.scoring.methods.TwoBodyEnergy)):
            score_types = method.score_types()
            if (scoretype in score_types):
                for this_score_type in method.score_types():
                    if (this_score_type in to_find):
                        to_find[this_score_type] = True  
                found.append(method)

    return found


methods = find_all_2body_methods(scorefxn)

fa_elec = None
etable = None

for method in methods:
    if ( isinstance(method, core.scoring.etable.AnalyticEtableEnergy) ):
        etable = method
    if ( isinstance(method, core.scoring.elec.FA_ElecEnergy) ):
        fa_elec = method

assert( not fa_elec is None and not etable is None )

def atid(resnum, atno):
    return core.id.AtomID( atno, resnum )

def score_atom_pair(res1, res2, bad_map ):
    score = 0
    rep = 0

    coulomb = fa_elec.coulomb()

    seqpos1 = res1.seqpos()
    seqpos2 = res2.seqpos()

    for at1 in range(1, res1.natoms()+1):
        if ( bad_map.get( atid(seqpos1, at1))):
            continue
        atom1 = res1.atom(at1)
        xyz1 = res1.xyz(at1)
        charge1 = res1.atomic_charge(at1)
        for at2 in range(1, res2.natoms()+1):
            if ( bad_map.get( atid(seqpos2, at2))):
                continue

            pe = core.scoring.etable.AtomPairEnergy()
            etable.atom_pair_energy(atom1, res2.atom(at2), 1, pe)

            score += pe.attractive
            score += pe.solvation
            rep += pe.repulsive

            elec = coulomb.eval_atom_atom_fa_elecE( xyz1, charge1, res2.xyz(at2), res2.atomic_charge(at2), 0)
            score += elec


    return score + rep * 0.55

def get_my_interface_graph(pose, bad_map):
    energy_graph = pose.energies().energy_graph()

    graph = np.zeros((pose.size()+1, pose.size()+1), np.float)

    real_weights = scorefxn.weights()

    other_weights = core.scoring.EMapVector()
    other_weights.assign(scorefxn.weights())
    other_weights.set(core.scoring.fa_rep, 0)
    other_weights.set(core.scoring.fa_atr, 0)
    other_weights.set(core.scoring.fa_sol, 0)
    other_weights.set(core.scoring.fa_elec, 0)

    ok_res = [True]
    for seqpos in range(1, pose.size()+1):
        is_good = True
        for atno in range(1, pose.residue(seqpos).natoms()+1):
            if ( bad_map.get(atid(seqpos, atno))):
                is_good = False

        ok_res.append(is_good)


    it = energy_graph.const_edge_list_begin()
    while it.valid():
        edge = it.__mul__()
        it.plus_plus()
        assert( not edge is None )

        low_node = edge.get_first_node_ind()
        high_node = edge.get_second_node_ind()
        if ( pose.chain(low_node) == pose.chain(high_node ) ):
            continue

        assert( low_node < high_node )

        score = edge.dot( real_weights )
        if ( score != 0 and ( not ok_res[low_node] or not ok_res[high_node] ) ):
            score = edge.dot( other_weights )
            score += score_atom_pair( pose.residue(low_node), pose.residue(high_node), bad_map )

        graph[low_node, high_node] = score

    return graph


def get_num_hbonds(pose, hbset, seqpos, atno, subset):
    num_hbonds = 0
    for hb in hbset.atom_hbonds(atid(seqpos, atno), False):
        if ( not subset[ hb.acc_res()] and not subset[ hb.don_res() ] ):
            continue
        if ( hb.energy() < -0.25 ):
            num_hbonds += 1

    res = pose.residue(seqpos)
    if ( res.heavyatom_has_polar_hydrogens(atno) ):
        for hatno in range(res.attached_H_begin(atno), res.attached_H_end(atno)+1):
            for hb in hbset.atom_hbonds(atid(seqpos, hatno), False):
                if ( not subset[ hb.acc_res()] and not subset[ hb.don_res() ] ):
                    continue
                if ( hb.energy() < -0.25 ):
                    num_hbonds += 1

    return num_hbonds

def get_hbond_atoms(pose, monomer_size):
    hbond_atoms = []
    if ( args.hbond_atoms != "" ):
        for item in args.hbond_atoms.split(","):
            sp = item.split("/")
            hbond_atoms.append((int(sp[0])+monomer_size, sp[1]))
    else:
        pose = pose.split_by_chain()[pose.num_chains()]
        for i in range(1, pose.size()+1):
            res = pose.residue(i)
            for j in range(1, res.nheavyatoms()+1):
                if ( res.heavyatom_is_an_acceptor(j) or res.heavyatom_has_polar_hydrogens(j) ):
                    hbond_atoms.append((i+monomer_size, res.atom_name(j)))
        return hbond_atoms



def get_seg_hbonds(pose, is_unsat, seg):

    monomer_size = pose.conformation().chain_end(pose.num_chains()-1)
    atomic_depth = core.scoring.atomic_depth.AtomicDepth( pose, 1.6, False, 0.49 )

    hb_scorefxn(pose)
    hbset = core.scoring.hbonds.HBondSet()
    core.scoring.hbonds.fill_hbond_set(pose, False, hbset)
    hbset.hbond_options().bb_donor_acceptor_check(False)
    core.scoring.hbonds.fill_hbond_set(pose, False, hbset)

    hbond_atoms = list(get_hbond_atoms(pose, monomer_size))

    is_sat = []
    for seqpos, atom in hbond_atoms:
        is_sat.append(not is_unsat.get(atid(seqpos, pose.residue(seqpos).atom_index(atom))))

    type_set = pose.residue(1).type().atom_type_set()

    is_buried = []
    is_buried_unsat = []
    is_buried_sat = []
    for ihbat, _ in enumerate(hbond_atoms):
        seqpos, atom = hbond_atoms[ihbat]
        res = pose.residue(seqpos)
        depth = atomic_depth.calcdepth(res.atom(atom), type_set)
        buried = depth > 3.0
        is_buried.append(buried)
        is_buried_unsat.append( buried and not is_sat[ihbat] )
        is_buried_sat.append( buried and is_sat[ihbat] )

    seg['buried_crit'] = np.sum(is_buried)
    seg['buried_crit_sat'] = np.sum(is_buried_sat)
    seg['buried_crit_unsat'] = np.sum(is_buried_unsat)


    subset = core.select.residue_selector.ChainSelector("A").apply(pose)

    hbond_to_crit = 0
    hbond_to_buried_crit = 0

    for ihbat, _ in enumerate(hbond_atoms):
        seqpos, atom = hbond_atoms[ihbat]

        buried = is_buried[ihbat]

        num_hbonds = get_num_hbonds(pose, hbset, seqpos, pose.residue(seqpos).atom_index(atom), subset)

        hbond_to_crit += num_hbonds
        if ( is_buried[ihbat] ):
            hbond_to_buried_crit += num_hbonds

    seg['hbond_to_crit'] = hbond_to_crit
    seg['hbond_to_buried_crit'] = hbond_to_buried_crit







ref_fname = args.ref_pdb
save_ref_pose = pose_from_file(ref_fname)

if ( silent != "" ):
    sfd_in = rosetta.core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
    sfd_in.read_file(silent)

    fnames = silent_tools.get_silent_index(silent)["tags"]
else:
    fnames = args.pdbs

for fname in fnames:
    # try:
    for k in [1]:
        basename = os.path.basename(fname)

        ftag = ""
        if (fname.endswith(".pdb")):
            ftag = basename[:-4]
        if (fname.endswith(".pdb.gz")):
            ftag = basename[:-7]
        if ( not ftag ):
            ftag = basename
        assert(ftag)
        ftag = args.out_prefix + ftag

        print("Processing pdb file %s"%fname )

        ddg_threshold = args.ddg_threshold
        multi_segs = args.multi_segs

        if ( silent == "" ):
            pose = pose_from_file(fname)
        else:
            pose = Pose()
            sfd_in.get_structure(fname).fill_pose(pose)
        ref_pose = save_ref_pose.clone()

        align = protocols.simple_moves.AlignChainMover()
        align.pose( ref_pose )
        align.source_chain( 2)
        align.target_chain( 1)
        align.apply(pose)




        binder_target = pose.split_by_chain()
        binder = binder_target[1]
        target = binder_target[2]
        chainA_len = binder.size()
        print("length of the binder %i"%chainA_len)


####################################################################################
####### satisfaction checking


        hb_scorefxn(pose)
        hbset = core.scoring.hbonds.HBondSet()
        core.scoring.hbonds.fill_hbond_set(pose, False, hbset)
        hbset.hbond_options().bb_donor_acceptor_check(False)
        core.scoring.hbonds.fill_hbond_set(pose, False, hbset)


        is_unsat = core.id.AtomID_Map_bool_t()
        is_unsat.resize(pose.size())
        is_next_to_unsat = core.id.AtomID_Map_bool_t()
        is_next_to_unsat.resize(pose.size())
        for i in range(1, pose.size()+1):
            is_unsat.resize( i, pose.residue(i).natoms(), False)
            is_next_to_unsat.resize( i, pose.residue(i).natoms(), False)



        subset = core.select.residue_selector.TrueResidueSelector().apply(pose)


        for seqpos in range(1, pose.size()+1):
            res = pose.residue(seqpos)
            for atno in range(1, res.nheavyatoms()+1):
                if ( res.heavyatom_is_an_acceptor(atno) or res.heavyatom_has_polar_hydrogens(atno) ):
                    if ( get_num_hbonds( pose, hbset, seqpos, atno, subset ) == 0 ):
                        is_unsat.set(atid(seqpos, atno), True)

                        is_next_to_unsat.set(atid(seqpos, atno), True)
                        is_next_to_unsat.set(atid(seqpos, res.atom_base(atno)), True)
                        if (res.heavyatom_has_polar_hydrogens(atno)):
                            for hatno in range(res.attached_H_begin(atno), res.attached_H_end(atno)+1):
                                is_next_to_unsat.set(atid(seqpos, hatno), True)




####################################################################################


        score = scorefxn(pose)
        print("total score of the complex: %.3f"%score)

        graph = get_my_interface_graph(pose, is_next_to_unsat)


        separate_pose = move_chainA_far_away(pose)
        scorefxn(separate_pose)


        per_res_ddg = utility.vector1_double()


        other_weights = core.scoring.EMapVector()
        other_weights.assign(scorefxn.weights())
        other_weights.set(core.scoring.fa_rep, 0)
        other_weights.set(core.scoring.fa_atr, 0)
        other_weights.set(core.scoring.fa_sol, 0)
        other_weights.set(core.scoring.fa_elec, 0)





        # per_res_score = utility.vector1_double()
        # per_res_bound = utility.vector1_double()
        for i in range(1, pose.size()+1):
            # ddgother = 2* ( pose.energies().residue_total_energies(i).dot(other_weights) \
            #                     - separate_pose.energies().residue_total_energies(i).dot(other_weights) )
            # print(ddgother)
            # ddg1 = 2* ( pose.energies().residue_total_energy(i) - separate_pose.energies().residue_total_energy(i) )
            ddg = np.sum( graph[:,i]) + np.sum( graph[i,:] )
            # print(ddg, ddg1)
            # assert(abs(ddg - ddg1) < 0.1)
            per_res_ddg.append(ddg)
            # per_res_score.append(separate_pose.energies().residue_total_energy(i))
            # per_res_bound.append(pose.energies().residue_total_energy(i))

        per_res_ddg0 = list(per_res_ddg)
        per_res_ddg = [0] + list(per_res_ddg)

############################################

        all_positions = utility.vector1_unsigned_long()
        for i in range(1, binder.size()+1):
            all_positions.append(i)

        poly_ala_binder_pose = pose.clone()
        protocols.toolbox.pose_manipulation.construct_poly_XXX_pose("GLY", poly_ala_binder_pose, all_positions, True, True, True)
        scorefxn(poly_ala_binder_pose)

        ala_graph = get_my_interface_graph(poly_ala_binder_pose, is_next_to_unsat)


        poly_ala_binder_separate_pose = move_chainA_far_away(poly_ala_binder_pose)
        scorefxn(poly_ala_binder_separate_pose)

        per_res_ddg_ala = utility.vector1_double()
        per_res_ala_scan = utility.vector1_double()
        # per_res_score = utility.vector1_double()
        # per_res_bound = utility.vector1_double()
        for i in range(1, pose.size()+1):
            # ddg = 2* ( poly_ala_binder_pose.energies().residue_total_energy(i) - poly_ala_binder_separate_pose.energies().residue_total_energy(i) )
            ddg = np.sum(ala_graph[:,i]) + np.sum(ala_graph[i,:])
            per_res_ddg_ala.append(ddg)
            per_res_ala_scan.append(per_res_ddg[i] - ddg)
            # per_res_score.append(separate_pose.energies().residue_total_energy(i))
            # per_res_bound.append(pose.energies().residue_total_energy(i))

        per_res_ddg0_ala = list(per_res_ddg_ala)


###############################################
        # pocket ddg

        if ( args.pocket_res != "" ):
            orig_pocket_res = [int(x) for x in args.pocket_res.split(",")]
            # 154,157,118,115,205
        

            real_pocket_res = []
            for old_res in orig_pocket_res:
                old_res += chainA_len
                real_pocket_res.append(old_res)
            pocket_sel = core.select.residue_selector.ResidueIndexSelector( ",".join([str(x) for x in real_pocket_res ] ) )
            target_sel = core.select.residue_selector.ChainSelector("B")

            not_pocket = core.select.residue_selector.NotResidueSelector( pocket_sel )
            to_delete = core.select.residue_selector.AndResidueSelector( not_pocket, target_sel )
            delete_mover = protocols.grafting.simple_movers.DeleteRegionMover()
            delete_mover.set_residue_selector( to_delete )

            pocket_binder = pose.clone()
            delete_mover.apply( pocket_binder )
            seperate_pocker_binder = move_chainA_far_away( pocket_binder )


            scorefxn(pocket_binder)
            scorefxn(seperate_pocker_binder)

            per_res_pocket_fa_atr = utility.vector1_double()
            # per_res_score = utility.vector1_double()
            # per_res_bound = utility.vector1_double()
            for i in range(1, pocket_binder.size()+1):
                ddg = 2* ( pocket_binder.energies().residue_total_energies(i)[core.scoring.fa_atr] 
                                - seperate_pocker_binder.energies().residue_total_energies(i)[core.scoring.fa_atr] )
                per_res_pocket_fa_atr.append(ddg)
                # per_res_score.append(separate_pose.energies().residue_total_energy(i))
                # per_res_bound.append(pose.energies().residue_total_energy(i))

            # pocket_binder.dump_pdb(fname + "_pocket_binder.pdb")
            per_res_pocket_fa_atr0 = list(per_res_pocket_fa_atr)
            # print(per_res_pocket_fa_atr0)


######################################
        dssp_str = better_dssp_hack(binder)

        if ( dssp_str[0] == "L" and dssp_str[1] != "L"):
            tmp = list(dssp_str)
            tmp[0] = tmp[1]
            dssp_str = ''.join(tmp) 
        if (dssp_str[-1] == "L" and dssp_str[-2] != "L"):
            tmp = list(dssp_str)
            tmp[-1] = tmp[-2]
            dssp_str = ''.join(tmp) 

        segs = []
        seg_temp = {}

        for ii in range(len(dssp_str)):
            if (ii == 0):
                seg_temp["start"] = 1
                seg_temp["sec_type"] = dssp_str[ii]

            if ( dssp_str[ii] != seg_temp["sec_type"]):
                seg_temp["end"] = ii
                segs.append(seg_temp)
                seg_temp = {}
                seg_temp["start"] = ii + 1
                seg_temp["sec_type"] = dssp_str[ii]

            if ( ii == len(dssp_str)-1 ):
                seg_temp["end"] = ii + 1
                segs.append(seg_temp)

        min_seg_length = 4
        new_segs = []
        for seg in segs:
            start = seg["start"]
            end = seg["end"]

            for i in range(start, end+1):
                for j in range(i + min_seg_length-1, end+1):
                    length = j - i + 1
                    if (length < min_seg_length):
                        continue

                    new_seg = {}
                    new_seg["start"] = i
                    new_seg["end"] = j
                    new_seg["sec_type"] = seg["sec_type"]
                    if ( per_res_ala_scan[i] > -0.8 or per_res_ala_scan[j] > -0.8 ): # must start and end with hotspot
                        continue
                    new_segs.append(new_seg)
        segs = new_segs

##########################################################

        dssp = core.scoring.dssp.Dssp(pose.split_by_chain()[1])

        pairing_set = dssp.strand_pairing_set()
        elements = get_ss_elements("x" + dssp_str)

        for ipairing in range(1, pairing_set.size()+1):
            pairing = pairing_set.strand_pairing(ipairing)

            ros_pairs = utility.vector1_core_scoring_dssp_Pairing()
            pairing.get_beta_pairs(ros_pairs)

            pairs = []
            for pair in ros_pairs:
                pairs.append((pair.Pos1(), pair.Pos2()))

            rep1 = pairs[0][0]
            rep2 = pairs[0][1]

            strand1 = None
            strand2 = None

            for element in elements:
                if ( element[0] != 'E' ):
                    continue
                if ( rep1 in range(element[1], element[2]+1) ):
                    assert(strand1 is None)
                    strand1 = element
                if ( rep2 in range(element[1], element[2]+1) ):
                    assert(strand2 is None)
                    strand2 = element
            assert(not strand1 is None and not strand2 is None)

            for start1 in range(strand1[1], strand1[2]+1):
                for end1 in range(start1, strand1[2]+1):
                    if ( end1 - start1 < 1 ):
                        continue
                    for start2 in range(strand2[1], strand2[2]+1):
                        for end2 in range(start2, strand2[2]+1):
                            if ( end2 - start2 < 1 ):
                                continue

                            any_pairs = False

                            two_range = range(start2, end2+1)
                            for one in range(start1, end1+1):
                                my_pair = None
                                for pair in pairs:
                                    if ( pair[0] == one ):
                                        my_pair = pair
                                        break
                                if ( my_pair is None ):
                                    continue
                                if ( pair[1] in two_range ):
                                    any_pairs = True

                            # print("Trying: %i-%i  %i-%i"%(start1, end1, start2, end2) + str(pairs))
                            if ( not any_pairs ):
                                # print("not paired")
                                continue

                            if ( per_res_ala_scan[start1] > -0.5 ): continue
                            if ( per_res_ala_scan[end1] > -0.5 ): continue
                            if ( per_res_ala_scan[start2] > -0.5 ): continue
                            if ( per_res_ala_scan[end2] > -0.5 ): continue
                            # print("Survived tails")

                            assert( start1 < start2 )

                            seg = {}
                            seg['start1'] = start1
                            seg['start2'] = start2
                            seg['end1'] = end1
                            seg['end2'] = end2
                            lengths = sorted([end1-start1+1, end2-start2+1])
                            seg['short_length'] = lengths[0]
                            seg['long_length'] = lengths[1]

                            seg['sec_type'] = 'P'

                            segs.append(seg)





        for seg in segs:
            seg["ddg"] = 0
            if ( "start" in seg ):
                for ires in range(seg["start"], seg["end"]+1):    
                    seg["ddg"] += per_res_ddg[ires]

            if ( "start1" in seg ):
                for ires in range(seg["start1"], seg["end1"]+1):    
                    seg["ddg"] += per_res_ddg[ires]
                for ires in range(seg["start2"], seg["end2"]+1):    
                    seg["ddg"] += per_res_ddg[ires]


        if (args.pocket_res != ""):
            if ( "start" in seg ):
                for seg in segs:
                    seg["fa_atr_pocket"] = 0
                    for ires in range(seg["start"], seg["end"]+1):    
                        seg["fa_atr_pocket"] += per_res_pocket_fa_atr[ires]

        good_motifs = []

        for seg in segs:
            if ( seg["ddg"] > ddg_threshold ):
                print("Reject ddg")
                continue
            if ( args.pocket_res != "" and seg["fa_atr_pocket"] > args.pocket_threshold ):
                print("Reject pocket")
                continue



            # good_motifs.append(seg)

            motif_res = utility.vector1_unsigned_long()

            if ( "start" in seg ):
                for i in range(seg["start"], seg["end"] + 1 ):
                    motif_res.append(i)
            if ( "start1" in seg ):
                for i in range(seg["start1"], seg["end1"] + 1 ):
                    motif_res.append(i)
                for i in range(seg["start2"], seg["end2"] + 1 ):
                    motif_res.append(i)


            motif_alone = core.pose.Pose()
            core.pose.pdbslice(motif_alone, binder, motif_res)

            motif_alone.pdb_info(core.pose.PDBInfo(motif_alone))
            for i in range(1, motif_alone.size()+1):
                motif_alone.pdb_info().chain(i, "A")

            print(core.select.residue_selector.ChainSelector("A").apply(motif_alone))


            motif_complex = motif_alone.clone()
            motif_complex.append_pose_by_jump(target, 1)

            motif_complex.pdb_info(core.pose.PDBInfo(motif_complex))
            for i in range(1, motif_complex.size()+1):
                if ( i <= motif_alone.size() ):
                    motif_complex.pdb_info().chain(i, "A")
                else:
                    motif_complex.pdb_info().chain(i, "B")

            # scorefxn(motif_complex)

            get_seg_hbonds(motif_complex, is_unsat, seg)

            # if ( seg['hbond_to_buried_crit'] < 1 ):
            #     continue
            # if ( seg['hbond_to_crit'] < 3 ):
            #     continue
            # if ( seg['buried_crit_unsat'] > 2 ):
            #     continue

            
            # for filt in filters_to_apply:
            #     the_filter = objs.get_filter(filt)

            #     # Get rid of stochastic filter
            #     if ( isinstance(the_filter, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
            #         the_filter = the_filter.subfilter()

            #     result, extra = compute_filter(motif_complex, the_filter, filt)
            #     seg[filt] = result
            #     if (not extra is None):
            #         seg["interface_sc_median_dist"] = extra

            # ddg_hydrophobic_pre = objs.get_filter("ddg_hydrophobic_pre").subfilter()
            # tmp = pose.clone()
            # remove_polars.apply(tmp)
            # seg["ddg_hydrophobic"] = ddg_hydrophobic_pre.compute(tmp)


            if ( "start" in seg ):
                seg["per_res_ddg"] = per_res_ddg0[seg["start"]-1:seg["end"]+1-1]
                seg["per_res_ddg_ala"] = per_res_ddg0_ala[seg["start"]-1:seg["end"]+1-1]

            if ( "start1" in seg ):
                seg["per_res_ddg1"] = per_res_ddg0[seg["start1"]-1:seg["end1"]+1-1]
                seg["per_res_ddg_ala1"] = per_res_ddg0_ala[seg["start1"]-1:seg["end1"]+1-1]

                seg["per_res_ddg2"] = per_res_ddg0[seg["start2"]-1:seg["end2"]+1-1]
                seg["per_res_ddg_ala2"] = per_res_ddg0_ala[seg["start2"]-1:seg["end2"]+1-1]

            if ( "start" in seg ):
                fout = ftag + "_%i_%i_%s.pdb.gz"%(seg["start"], seg["end"], seg["sec_type"])
            if ( "start1" in seg ):
                fout = ftag + "_%i_%i_%i_%i_%s.pdb.gz"%(seg["start1"], seg["end1"], seg["start2"], seg["end2"], seg["sec_type"])

            motif_alone.dump_pdb(fout)
            motif_complex.dump_pdb(fout.replace(".pdb.gz", "_og.pdb.gz"))

            def default(o):
                if isinstance(o, np.int64): return int(o)  
                raise TypeError

            fname = fout[:-7]
            seg["name"] = fname
            f = open(fname + ".json", "w")
            f.write(json.dumps(seg, indent=4, default=default))
            f.close()






        print("Job done!")

    # except Exception as e:
    #     print("Error!!!")
    #     print(e)



