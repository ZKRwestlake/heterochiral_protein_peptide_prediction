#!/usr/bin/env python

import os
import sys
import json
import numpy as np

files = sys.argv[1:]


hotspot_threshold = -2

pdbs = []
for file in files:
    with open(file) as f:
        for line in f:
            line = line.strip()
            if (len(line) == 0):
                continue
            pdbs.append(line)


for pdb in pdbs:
    js_file = pdb.replace(".pdb.gz", ".json")
    with open(js_file) as f:
        js = json.loads(f.read())

    stuffs = []
    names = ["", "1", "2"]

    for name in names:
        if ( 'per_res_ddg'+name not in js ):
            continue
        hotspot_pos = np.where(np.array(js['per_res_ddg'+name]) - np.array(js['per_res_ddg_ala'+name])  < hotspot_threshold)[0] + 1
        hotspots = ":".join(str(x) for x in hotspot_pos)

        stuffs.append(hotspots)

    print(pdb, ",".join(stuffs))













