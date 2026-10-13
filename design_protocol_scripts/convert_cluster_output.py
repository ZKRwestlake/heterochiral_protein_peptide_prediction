#!/usr/bin/env python

import os
import sys
import pandas as pd

file = sys.argv[1]


records = []


with open(file) as f:
    while (True):
        try:
            line1 = next(f)
            line2 = next(f)
        except:
            break

        the_map = {}
        sp = line1.split()
        for i in range(int(len(sp)/2)):
            key = sp[i*2+0].lower().replace(":", "")
            value = sp[i*2+1]
            try:
                value = float(value)
            except:
                pass
            the_map[key] = value

        pdb = line2.split()[0]
        og = pdb.replace(".pdb", "_og.pdb")


        the_map['path'] = pdb
        the_map['og'] = og
        the_map['description'] = os.path.basename(pdb)
        records.append( the_map )



df = pd.DataFrame(records)
df.to_csv("tmp", sep=" ", index=False)

with open("tmp") as f:
    for line in f:
        print(line.strip())


os.remove("tmp")

