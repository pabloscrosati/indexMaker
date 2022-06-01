# Header section
__author__ = "Kirsty Ng and Pablo M. Scrosati"
__version__ = "dev_0.2"
__update__ = "May 31, 2022"
__title__ = "Program Name Pending"
__description__ = """
Automated extraction of Protein-Gdm contacts.
"""

# Global variables
coarse_cutoff = 2
short_cutoff = 0.5

# Import dependencies
import glob, os, multiprocessing, textwrap, numpy, time

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
# Print header information on launch of program
def program_launch(author=__author__, version=__version__, title=__title__, update=__update__):
    header_length = 45 # Full length of header is header_length + 2 characters
    print("|" + ("-" * header_length) + "|")
    print("|" + title.center(header_length) + "|")
    print("|" + ("Written by: " + author).center(header_length) + "|")
    print("|" + ("Version: " + version).center(header_length) + "|")
    print("|" + ("Last updated: " + update).center(header_length) + "|")
    print("|" + ("-" * header_length) + "|")
    return

def cmd_parse(argv=None, prog_description=__description__):
    parser = argparse.ArgumentParser(description=prog_description)
    parser.add_argument('-f', '--file-pattern', metavar='', help='filename pattern', required=True)
    parser.add_argument('-d', '--directory', metavar='', help='directory containing files', required=True)
    parser.add_argument('-o', '--output', metavar='', help='output file')
    return parser.parse_args(argv)

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
def gro_parse(gro_file):
    resid, resname, atomname, atomnum, x_c, y_c, z_c = [], [], [], [], [] ,[] ,[]
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

# Will find unique protein chains based on pattern matching
def find_protein_chain(resname, thinned_list=[]):
    for i in resname:
        if i in heavyAtoms and i != "GUAN":
            thinned_list.append(i)
    principal_string = ''.join(thinned_list)
    residue_pattern = (principal_string + principal_string).find(principal_string, 1, -1)
    if residue_pattern == -1:
        print("Multiple proteins not found!")
        return [len(thinned_list) - 1]
    else:
        protein_chain = textwrap.wrap(principal_string[:residue_pattern], 3)
        number_of_chains = principal_string.count(principal_string[:residue_pattern])
        print("%s protein chains found!" % number_of_chains)
        indices = []
        for i in range(number_of_chains):
            indices.append((len(protein_chain) - 1) + (i * len(protein_chain)))
        return indices, number_of_chains

# Average protein position
def average_position(x_coord, y_coord, z_coord):
    ave_x = sum(float(x) for x in x_coord) / len(x_coord)
    ave_y = sum(float(y) for y in y_coord) / len(y_coord)
    ave_z = sum(float(z) for z in z_coord) / len(z_coord)
    return ave_x, ave_y, ave_z

# Step 1 of the prefilter; finds max distance and sets cutoff
def max_distance_prefilter(x_list, y_list, z_list, x_ave, y_ave, z_ave, cutoff=coarse_cutoff):
    distance_log = []
    x_array = numpy.array(x_list).astype(float)
    y_array = numpy.array(y_list).astype(float)
    z_array = numpy.array(z_list).astype(float)
    x_ave_array = numpy.array(x_ave).astype(float)
    y_ave_array = numpy.array(y_ave).astype(float)
    z_ave_array = numpy.array(z_ave).astype(float)
    coordinate_array = numpy.stack((x_array, y_array, z_array), axis=1)
    average_array = numpy.stack((x_ave_array, y_ave_array, z_ave_array))
    twod_array = numpy.expand_dims(average_array, axis=0)
    for i in range(len(coordinate_array)):
        distance = numpy.linalg.norm(coordinate_array[i]-twod_array)
        distance_log.append(distance)
    max_distance = max(distance_log)
    buffered_max_distance = max_distance + float(cutoff)
    return buffered_max_distance

# Find the indices of all GUAN residues
def find_guan_indices(resname):
    indices = [i for i, x in enumerate(resname) if x == "GUAN"]
    return min(indices), max(indices)

# Find GUAN residues within cutoff distance
def guan_index_prefilter(x_ave, y_ave, z_ave, coord_x, coord_y, coord_z, atomname, guan_index, cutoff):
    success_guan_list = []
    guanonly = atomname[guan_index[0]:guan_index[1]+1]

    x_ave_array = numpy.array(x_ave).astype(float)
    y_ave_array = numpy.array(y_ave).astype(float)
    z_ave_array = numpy.array(z_ave).astype(float)
    average_array = numpy.stack((x_ave_array, y_ave_array, z_ave_array))
    twod_array = numpy.expand_dims(average_array, axis=0)
    guan_x_array = numpy.array(coord_x[guan_index[0]:guan_index[1]+1]).astype(float)
    guan_y_array = numpy.array(coord_y[guan_index[0]:guan_index[1]+1]).astype(float)
    guan_z_array = numpy.array(coord_z[guan_index[0]:guan_index[1]+1]).astype(float)
    guan_array = numpy.stack((guan_x_array, guan_y_array, guan_z_array), axis=1)
    for i in range(len(guan_array)):
        if guanonly[i] == heavyAtoms["GUAN"]:
            distance = numpy.linalg.norm(guan_array[i]-twod_array)
            if distance < cutoff:
                success_guan_list.append(i + guan_index[0])
    return success_guan_list

def residue_guan_pairs(coord_x, coord_y, coord_z, prot_index, guan_index, resid, resname, atomname, success_guan, cutoff=short_cutoff):
    prot_x_array = numpy.array(coord_x[prot_index[0]:prot_index[1]]).astype(float)
    prot_y_array = numpy.array(coord_y[prot_index[0]:prot_index[1]]).astype(float)
    prot_z_array = numpy.array(coord_z[prot_index[0]:prot_index[1]]).astype(float)
    prot_array = numpy.stack((prot_x_array, prot_y_array, prot_z_array), axis=1)
    guan_x_array = numpy.array(coord_x[guan_index[0]:guan_index[1]+1]).astype(float)
    guan_y_array = numpy.array(coord_y[guan_index[0]:guan_index[1]+1]).astype(float)
    guan_z_array = numpy.array(coord_z[guan_index[0]:guan_index[1]+1]).astype(float)
    guan_array = numpy.stack((guan_x_array, guan_y_array, guan_z_array), axis=1)
    distance_pairs = ["residue\tguan\tdistance"]
    guan_gro_string = ""
    for i in range(len(prot_array)):
        if resname[i] in heavyAtoms and atomname[i] == heavyAtoms[resname[i]]:
            for n in success_guan:
                distance = numpy.linalg.norm(prot_array[i] - guan_array[n - guan_index[0]])
                if distance < cutoff:
                    distance_pairs.append(resid[i] + resname[i] + '\t' + resid[n] + resname[n] + "\t" + "%.4f" % distance)
                    guan_gro_string = guan_gro_string + str(n+1).rjust(6) + str(n+2).rjust(6) + str(n+3).rjust(6) +\
                                     str(n+4).rjust(6) + str(n+5).rjust(6) + str(n+6).rjust(6) + str(n+7).rjust(6) +\
                                     str(n+8).rjust(6) + str(n+9).rjust(6) + str(n+10).rjust(6)
    gromacs_index = textwrap.wrap(guan_gro_string, 60)
    gromacs_index.insert(0, "[ Gdm_Index ]")
    return distance_pairs, gromacs_index


if __name__ == "__main__":
    execution_start = time.perf_counter()
    program_launch()

    # Open file, interpret and extract contents
    gro_file = file_open("4_1mguhcl_125.gro")
    _, _, _, resid, resname, atomname, _, x_c, y_c, z_c = gro_parse(gro_file)

    # Find multiple protein chains
    protein_indices, number_of_chains = find_protein_chain(resname)
    if number_of_chains > 1:
        for chain in range(number_of_chains):
            print("ell")

    # THIS IS WHERE I AM, ADDING MULTIPLE PROTEIN LOOP >> SINGLE LOOP WORKS

    average_pos = average_position(x_c[0:protein_indices[0]], y_c[0:protein_indices[0]], z_c[0:protein_indices[0]])

    cutoff_distance = max_distance_prefilter(x_c[0:protein_indices[0]], y_c[0:protein_indices[0]],
                                             z_c[0:protein_indices[0]], average_pos[0], average_pos[1], average_pos[2])

    guan_indices = find_guan_indices(resname)

    success_guan = close_guan_indices = guan_index_prefilter(average_pos[0], average_pos[1], average_pos[2], x_c, y_c,
                                                             z_c, atomname, guan_indices, cutoff_distance)

    pair_list = residue_guan_pairs(x_c, y_c, z_c, (0, protein_indices[0]), guan_indices, resid, resname, atomname, success_guan)

    execution_end = time.perf_counter()
    print(f"Program executed in {round(execution_end - execution_start, 5)} seconds.")
