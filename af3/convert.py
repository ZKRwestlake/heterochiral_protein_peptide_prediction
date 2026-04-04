# from sys import argv
from Bio.PDB import MMCIFParser, PDBIO 

parser = MMCIFParser(QUIET=True)
s = parser.get_structure('structure', 'model.cif')

for model in s:
    for chain in model:
        for residue in chain:
            residue.resname = residue.resname.upper()

io = PDBIO()
io.set_structure(s)
io.save('model.pdb')
