import pymol
from pymol import cmd

########### EDIT ONLY THESE LINES, USE FULL PATHS ############
#Your ESMfold protein
pdb_file1 = './pdb/mystery_protein.pdb'
text_file1 = './myst_prot/mystery_protein_100.txt'
#Your protein from RCSB
pdb_file2 = '/Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Assignment3/pdb/Calmodulin.pdb'
text_file2 = '/Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Assignment3/calmodulin/Calmodulin_100.txt'
##############################################################

# Load the PDB files
def load_and_align_pdb(pdb1, pdb2):
    cmd.load(pdb1, 'esmfold')
    cmd.load(pdb2, 'rcsb')
    
# Color top 10 ranked calcium binding sites
def color_atoms(protein_name, coordinates):
    for i, atom_coords in enumerate(coordinates[:10]):
        x, y, z = float(atom_coords[0]), float(atom_coords[1]), float(atom_coords[2])

        if protein_name == 'esmfold':
            cmd.pseudoatom(object=protein_name, pos=(x, y, z), color='red')
        if protein_name == 'rcsb':
            cmd.pseudoatom(object=protein_name, pos=(x, y, z), color='blue')

# Read coordinates from the input text file
def read_coordinates_from_text_file(text_file):
    coordinates = []
    with open(text_file, 'r') as file:
        for line in file:
            values = line.split()
            if len(values) == 5:
                x, y, z = float(values[1]), float(values[2]), float(values[3])
                coordinates.append((x, y, z))
    return coordinates

load_and_align_pdb(pdb_file1, pdb_file2)
coordinates1 = read_coordinates_from_text_file(text_file1)
coordinates2 = read_coordinates_from_text_file(text_file2)
color_atoms('esmfold', coordinates1)
color_atoms('rcsb', coordinates2)

cmd.color("green", "esmfold")
cmd.color("pink", "rcsb")
cmd.select("rcsb_calcium", "rcsb and elem Ca")
cmd.color("orange", "rcsb_calcium")
cmd.set("sphere_scale", 0.5, "rcsb_calcium")
cmd.super('esmfold', 'rcsb')
cmd.remove("(elem H)")
cmd.remove("(resn HOH)")
cmd.show('cartoon')
cmd.zoom('all')
