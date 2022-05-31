# Header section
__author__ = "Kirsty Ng and Pablo M. Scrosati"
__version__ = "dev_0.1"
__date__ = "May 30, 2022"
__title__ = "Program Name Pending"

# Global variables
coarse_cutoff = 10
short_cutoff = 2

# Change

# Import dependencies
import glob, os, multiprocessing

# Dictionary for residue heavy atoms
heavyAtoms = {
    "GUAN": "C",
    "LYS": "CE",
    "GLU": "CD",
    "PRO": "CG",
    "ASP": "CG",
    "PHE": "CZ",
    "GLN": "CD",
    "THR": "CG2",
    "ARG": "CZ",
    "ASN": "CG",
    "TYR": "CZ",
    "TRP": "CE2",
    "HIS": "CE1",
    "ALA": "CB",
    "GLY": "CA",
    "SER": "CB",
    "VAL": "CB",
    "LEU": "CG",
    "ILE": "CD",
    "MET": "CE"
}

# Functions
# Generate list of filenames matching pattern
def file_pattern(pattern, folder, files=[]):
    if not os.path.isdir(os.path.join(os.curdir, folder)):
        print("Folder path supplied does not exist!")
        print("Check that you have the correct relative path!")
        exit(1)
    for f in glob.glob(os.path.join(os.curdir, folder, (pattern + "*.gro"))):
        files.append(f)
    return files

# Open file and return as list
def file_open(file):
    with open(file) as f:
        config = [line.rstrip() for line in f]
    return config

# Parse list in GRO file format and return components
def gro_parse(gro_file, resid=[], resname=[], atomname=[], atomnum=[], x_c=[], y_c=[], z_c=[]):
    title, num_atoms, box_size = gro_file[0], gro_file[1], gro_file[-1]
    for i in gro_file:
        if i == title or i == num_atoms or i == box_size:
            pass
        else:
            resid.append(i[0:5].strip())
            resname.append(i[5:10].strip())
            atomname.append(i[10:15].strip())
            atomnum.append(i[15:20].strip())
            x_c.append(i[20:28].strip())
            y_c.append(i[28:36].strip())
            z_c.append(i[36:44].strip())
    return title.strip(), num_atoms.strip(), box_size.strip(), resid, resname, atomname, atomnum, x_c, y_c, z_c

if __name__ == "__main__":
    gro_file = file_open("4_1mguhcl_125.gro")
    resname = gro_parse(gro_file)[4]
    unique_list = list(set(resname))
    print(unique_list)