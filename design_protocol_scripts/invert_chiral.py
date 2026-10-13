#!/usr/bin/env python
# coding: utf-8

import os
import sys
from pyrosetta import *
init()

infile = sys.argv[1] 
ofile = sys.argv[2]

pose = pose_from_pdb(F"{infile}")
for t in range(0, len(pose.residues)):
    for i in range(0, len(pose.residue(t+1).atoms())):
        pose.residue(t+1).set_xyz(i+1,pose.residue(t+1).xyz(pose.residue(t+1).atom_name(i+1)).negate())

pose.dump_pdb(F"{ofile}")
