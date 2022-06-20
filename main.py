# Header section
__author__ = "Pablo M. Scrosati"
__version__ = "0.9.3 beta"
__update__ = "June 20, 2022"
__title__ = "Auto Guanidinium Contact Detection"
__description__ = """
Automatic protein chain detection with residue-specific detection of guanidinium contacts.
Distances, residue ID, and guanidinium IDs are stored and exported.
Exported data is compatible with additional Excel resources.
Features coarse cut-off to speed initial contact discovery and multi-process support.
"""

# Debug string to print to monitor flow control
# REMOVE THIS IN FINAL VERSION!!
debug = "debug"

# Global variables
# Will be added as defaults, or overridden by command-line options
coarse_cutoff = 2
short_cutoff = 0.5

# Import dependencies
import glob, os, textwrap, numpy, time, argparse, sys, concurrent.futures

# Dictionary for residue heavy atoms used for identification
# Heavy atoms are corresponding to only residues contained within PDB 1WLA setup for pH 2
# This can be expanded manually, or added in as an import to include additional dictionary terms for different atoms
# The import functionality is not currently planned, but can be added with ease if needed
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
    """
    Prints header information on program launch
    :return: None
    """
    header_length = 45  # Full length of header is header_length + 2 characters
    print("|" + ("-" * header_length) + "|")
    print("|" + title.center(header_length) + "|")
    print("|" + ("Written by: " + author).center(header_length) + "|")
    print("|" + ("Version: " + version).center(header_length) + "|")
    print("|" + ("Last updated: " + update).center(header_length) + "|")
    print("|" + ("-" * header_length) + "|")
    return


# Parses command line arguments and returns list containing arguments and options
# This function is not finished as it required proper output file handling and other optional arguments
def cmd_parse(argv=None, prog_description=__description__):
    """
    Read command line arguments and parse
    :return: List of command line arguments
    """
    if sys.version_info.minor < 8:
        parser = argparse.ArgumentParser(description=prog_description)
        parser.add_argument('-f', '--file-pattern', help='filename pattern for input files', required=True)
        parser.add_argument('-d', '--directory', help='directory containing files to import', required=True)
        parser.add_argument('-o', '--output', help='output file', required=True)
        parser.add_argument('--csv', help='output CSV file alongside tab-delimited file', action='store_true')
        parser.add_argument('--st', help='run program in single-threaded mode', action='store_true')
        parser.add_argument('--coarse-cutoff', help='define coarse cutoff distance')
        parser.add_argument('--short-cutoff', help='define short cutoff distance')
        return parser.parse_args(argv)
    else:
        parser = argparse.ArgumentParser(description=prog_description)
        parser.add_argument('-f', '--file-pattern', metavar='', help='filename pattern for input files', required=True)
        parser.add_argument('-d', '--directory', metavar='', help='directory containing files to import', required=True)
        parser.add_argument('-o', '--output', metavar='', help='output file', required=True)
        parser.add_argument('--csv', help='output CSV file alongside tab-delimited file', action='store_true')
        parser.add_argument('--st', help='run program in single-threaded mode', action='store_true')
        parser.add_argument('--coarse-cutoff', metavar='', help='define coarse cutoff distance')
        parser.add_argument('--short-cutoff', metavar='', help='define short cutoff distance')
        return parser.parse_args(argv)


# Generate list of filenames in a folder matching pattern and return list of file names
# Folder parameter is required by argparse, but can be assigned default behavior in an update
def file_pattern(pattern, folder, files=[]):
    """
    Find files in a folder matching the pattern and return list of file names
    :param pattern: File name pattern to match, i.e. test_01.gro pattern is test_
    :param folder: Folder containing GRO files to find
    :param files: Empty list to be filled with file names
    :return: List of file names
    """
    if not os.path.isdir(os.path.join(os.curdir, folder)):
        print("Folder path supplied does not exist!")
        print("Check that you have the correct relative path!")
        exit(1)
    for f in glob.glob(os.path.join(os.curdir, folder, (pattern + "*.gro"))):
        files.append(f)
    return files


# Open file and return list containing file contents
def file_open(file):
    """
    Open a file and read contents into list
    :param file: String of the file name to be opened
    :return: List containing file contents
    """
    with open(file) as f:
        config = [line.rstrip() for line in f]
    return config


# Parse list in GRO file format and return components as lists
def gro_parse(gro_file):
    """
    Parse list containing GRO file contents
    :param gro_file: Input list in original GRO file format
    :return: Lists of components of the GRO file
    """
    resid, resname, atomname, atomnum, x_c, y_c, z_c = [], [], [], [], [], [], []
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
# This function relies on the dictionary as defined above to thin list and avoid problems with mixed charge states
def find_protein_chain(resname, file, silent=True):
    """
    Finds the longest repeating chain of atoms to detect protein residues; List is thinned using dictionary
    :param resname: List containing atom names (see gro_parse)
    :param file: File name corresponding to the list; Allows verbose and debugging functionality if needed
    :return: List containing indices of protein chain start and stop; number of chains detected
    """
    thinned_list = []
    for i in resname:
        if i in heavyAtoms and i != "GUAN":
            thinned_list.append(i)
    principal_string = ''.join(thinned_list)
    residue_pattern = (principal_string + principal_string).find(principal_string, 1, -1)
    if residue_pattern == -1:
        if silent is True:
            print("Multiple proteins not found!")
        return [len(thinned_list) - 1], 1
    else:
        protein_chain = textwrap.wrap(principal_string[:residue_pattern], 3)
        number_of_chains = principal_string.count(principal_string[:residue_pattern])
        if silent is True:
            print("%s protein chains found in %s" % (number_of_chains, os.path.split(file)[-1]))
        indices = [0]
        for i in range(number_of_chains):
            indices.append((len(protein_chain) - 1) + (i * len(protein_chain)))
        return indices, number_of_chains


# Find average protein coordinates
def average_position(x_coord, y_coord, z_coord):
    """
    Find average protein position from discrete atom positions
    :param x_coord: List containing protein x coordinates (indices are handled at the time the function is called)
    :param y_coord: List containing protein y coordinates (indices are handled at the time the function is called)
    :param z_coord: List containing protein z coordinates (indices are handled at the time the function is called)
    :return: Floating point values of average x,y, and z components
    """
    ave_x = sum(float(x) for x in x_coord) / len(x_coord)
    ave_y = sum(float(y) for y in y_coord) / len(y_coord)
    ave_z = sum(float(z) for z in z_coord) / len(z_coord)
    return ave_x, ave_y, ave_z


# Initial prefilter step
# Find maximum distance between protein residues and protein center, then applies coarse_cutoff for cutoff distance
# Cutoff is used to exclude residues during distance calculations to speed up residue-level calculations
def max_distance_prefilter(x_list, y_list, z_list, x_ave, y_ave, z_ave, cutoff=coarse_cutoff):
    """
    Find maximum distance within a protein and apply coarse_cutoff, then return final cutoff distance
    :param x_list: List containing protein x coordinates (indices are handled at the time the function is called)
    :param y_list: List containing protein y coordinates (indices are handled at the time the function is called)
    :param z_list: List containing protein z coordinates (indices are handled at the time the function is called)
    :param x_ave: Float value containing average x coordinate, typically protein
    :param y_ave: Float value containing average y coordinate, typically protein
    :param z_ave: Float value containing average z coordinate, typically protein
    :param cutoff: User-defined (or default value) float value to define coarse cutoff
    :return: Adjusted cutoff value for excluding guanidinium from distance calculations
    """
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
        distance = numpy.linalg.norm(coordinate_array[i] - twod_array)
        distance_log.append(distance)
    max_distance = max(distance_log)
    buffered_max_distance = max_distance + float(cutoff)
    return buffered_max_distance


# Find the first and last index of all guanidinium residues
def find_guan_indices(resname):
    """
    Find the indices of the start and end of guanidinium residues
    :param resname: Full list of residue names
    :return: Start and end indices denoting the region where guanidinium residues are found
    """
    indices = [i for i, x in enumerate(resname) if x == "GUAN"]
    return min(indices), max(indices)


# Find GUAN residues within cutoff distance
def guan_index_prefilter(x_ave, y_ave, z_ave, coord_x, coord_y, coord_z, atomname, guan_index, cutoff):
    """
    Sets maximum distance to solve residue-guanidinium distances
    :param x_ave: Protein average x value (float)
    :param y_ave: Protein average y value (float)
    :param z_ave: Protein average z value (float)
    :param coord_x: Coordinate list containing only x coordinates, indices are handled within this function
    :param coord_y: Coordinate list containing only y coordinates, indices are handled within this function
    :param coord_z: Coordinate list containing only z coordinates, indices are handled within this function
    :param atomname: List containing atom names
    :param guan_index: Indices containing start and stop positions to find guanidinium residues
    :param cutoff: Cutoff for solving, found in previous function
    :return: All possible guanidinium residues within cutoff is list
    """
    success_guan_list = []
    guanonly = atomname[guan_index[0]:guan_index[1] + 1]
    x_ave_array = numpy.array(x_ave).astype(float)
    y_ave_array = numpy.array(y_ave).astype(float)
    z_ave_array = numpy.array(z_ave).astype(float)
    average_array = numpy.stack((x_ave_array, y_ave_array, z_ave_array))
    twod_array = numpy.expand_dims(average_array, axis=0)
    guan_x_array = numpy.array(coord_x[guan_index[0]:guan_index[1] + 1]).astype(float)
    guan_y_array = numpy.array(coord_y[guan_index[0]:guan_index[1] + 1]).astype(float)
    guan_z_array = numpy.array(coord_z[guan_index[0]:guan_index[1] + 1]).astype(float)
    guan_array = numpy.stack((guan_x_array, guan_y_array, guan_z_array), axis=1)
    for i in range(len(guan_array)):
        if guanonly[i] == heavyAtoms["GUAN"]:
            distance = numpy.linalg.norm(guan_array[i] - twod_array)
            if distance < cutoff:
                success_guan_list.append(i + guan_index[0])
    return success_guan_list


# Find residue-guanidinium pairs within fine cutoff
def residue_guan_pairs(coord_x, coord_y, coord_z, prot_index, guan_index, resid, resname, atomname, success_guan, file,
                       chain, cutoff=short_cutoff):
    """
    Explicit solving of residue-guanidinium contact distances
    :param coord_x: Coordinate list containing only x coordinates, indices are handled within this function
    :param coord_y: Coordinate list containing only y coordinates, indices are handled within this function
    :param coord_z: Coordinate list containing only z coordinates, indices are handled within this function
    :param prot_index: Indices containing start and stop positions for protein chain
    :param guan_index: Indices containing start and stop positions for guanidinium residues
    :param resid: List containing residue ID numbers
    :param resname: List containing residues names
    :param atomname: List containing atom names
    :param success_guan: List containing indices for guanidinium residues within coarse cutoff from previous function
    :param file: File name for output formatting
    :param chain: Chain identifier for output formatting
    :param cutoff: Fine cutoff as a user-defined parameter
    :return: Formatted list containing residue-guanidinium contacts and GROMACS-compatible index list
    """
    prot_x_array = numpy.array(coord_x[prot_index[0]:prot_index[1]]).astype(float)
    prot_y_array = numpy.array(coord_y[prot_index[0]:prot_index[1]]).astype(float)
    prot_z_array = numpy.array(coord_z[prot_index[0]:prot_index[1]]).astype(float)
    prot_array = numpy.stack((prot_x_array, prot_y_array, prot_z_array), axis=1)
    guan_x_array = numpy.array(coord_x[guan_index[0]:guan_index[1] + 1]).astype(float)
    guan_y_array = numpy.array(coord_y[guan_index[0]:guan_index[1] + 1]).astype(float)
    guan_z_array = numpy.array(coord_z[guan_index[0]:guan_index[1] + 1]).astype(float)
    guan_array = numpy.stack((guan_x_array, guan_y_array, guan_z_array), axis=1)
    distance_pairs = ["%s\t%s\t" % (os.path.split(file)[-1], chain + 1), "residue\tguan\tdistance"]
    guan_gro_string = ''
    for i in range(len(prot_array)):
        if resname[i] in heavyAtoms and atomname[i] == heavyAtoms[resname[i]]:
            for n in success_guan:
                distance = numpy.linalg.norm(prot_array[i] - guan_array[n - guan_index[0]])
                if distance < cutoff:
                    distance_pairs.append(resid[i] + resname[i] + '\t' + resid[n] + resname[n] + "\t" + "%.4f" %
                                          distance)
                    guan_gro_string = guan_gro_string + str(n + 1).rjust(6) + str(n + 2).rjust(6) + str(n + 3).rjust(
                        6) + str(n + 4).rjust(6) + str(n + 5).rjust(6) + str(n + 6).rjust(6) + str(
                        n + 7).rjust(6) + str(n + 8).rjust(6) + str(n + 9).rjust(6) + str(n + 10).rjust(6)
    gromacs_index_temp = textwrap.wrap(guan_gro_string, 60)
    gromacs_index = [text.rjust(60) for text in gromacs_index_temp]
    gromacs_index.append('')
    gromacs_index.insert(0, "[ Gdm_Index %s %s ]" % (os.path.split(file)[-1], chain + 1))
    return distance_pairs, gromacs_index


# Write list to file
def write_file(output_list, filename):
    """
    Writes list to file
    :param output_list: List to be written to file
    :param filename: Name of file to be written
    :return: None
    """
    with open(filename, 'w') as f:
        for item in output_list:
            f.write("%s\n" % item)
    return


# Main program execution
def main_single_thread(file):
    """
    Single-threaded full program, optimized for sequential execution, but capable in multi-threaded mode
    :param file: File to work on, file name as present in directory
    :return: Guanidinium-residue pair list and GROMACS-compatible index list
    """
    # Open file, interpret and extract contents
    gro_file = file_open(file)
    _, _, _, resid, resname, atomname, _, x_c, y_c, z_c = gro_parse(gro_file)

    # Find multiple protein chains
    protein_indices, number_of_chains = find_protein_chain(resname, file)

    # Find where guanidinium residues start and end
    guan_indices = find_guan_indices(resname)

    # Execution handler for multiple protein chains
    export_pair_list, export_gro_index = [], []
    if number_of_chains > 1:
        for chain in range(number_of_chains):
            if chain == 0:
                modifier = 0
            else:
                modifier = 1
            average_pos = average_position(x_c[(protein_indices[chain] + modifier):protein_indices[chain + 1]],
                                           y_c[(protein_indices[chain] + modifier):protein_indices[chain + 1]],
                                           z_c[(protein_indices[chain] + modifier):protein_indices[chain + 1]])
            cutoff_distance = max_distance_prefilter(
                x_c[(protein_indices[chain] + modifier):protein_indices[chain + 1]],
                y_c[(protein_indices[chain] + modifier):protein_indices[chain + 1]],
                z_c[(protein_indices[chain] + modifier):protein_indices[chain + 1]],
                average_pos[0], average_pos[1], average_pos[2])
            success_guan = guan_index_prefilter(average_pos[0], average_pos[1], average_pos[2], x_c, y_c, z_c, atomname,
                                                guan_indices, cutoff_distance)
            pair_list, gro_index = residue_guan_pairs(x_c, y_c, z_c,
                                                      (protein_indices[chain] + modifier, protein_indices[chain + 1]),
                                                      guan_indices, resid, resname, atomname, success_guan, file, chain)
            export_pair_list.append(pair_list)
            export_gro_index.extend(gro_index)
            write_file(pair_list, "%s_output_protein_chain_%s.txt" % (
                os.path.join(os.path.split(file)[0], os.path.split(file)[-1].split(".")[0]), chain + 1))

        max_length_pair_list = max(len(x) for x in export_pair_list)
        for i in range(len(export_pair_list)):
            if len(export_pair_list[i]) < max_length_pair_list:
                export_pair_list[i].extend([''] * (max_length_pair_list - len(export_pair_list[i])))

        final_export_pair_list = []
        final_export_pair_list.extend([''] * max_length_pair_list)
        for n in range(len(export_pair_list)):
            for i in range(max_length_pair_list):
                if export_pair_list[n][i] == '':
                    final_export_pair_list[i] = final_export_pair_list[i] + "\t\t" + "\t"
                else:
                    final_export_pair_list[i] = final_export_pair_list[i] + export_pair_list[n][i] + "\t"

        return final_export_pair_list, export_gro_index

    # Single protein chain logic
    else:
        average_pos = average_position(x_c[0:protein_indices[0]], y_c[0:protein_indices[0]], z_c[0:protein_indices[0]])
        cutoff_distance = max_distance_prefilter(x_c[0:protein_indices[0]], y_c[0:protein_indices[0]],
                                                 z_c[0:protein_indices[0]], average_pos[0], average_pos[1],
                                                 average_pos[2])
        success_guan = guan_index_prefilter(average_pos[0], average_pos[1], average_pos[2], x_c, y_c, z_c, atomname,
                                            guan_indices, cutoff_distance)
        pair_list, gro_index = residue_guan_pairs(x_c, y_c, z_c, (0, protein_indices[0]), guan_indices, resid, resname,
                                                  atomname, success_guan, file, chain=0)
        write_file(pair_list,
                   "%s_output_protein_chain_1.txt" % (os.path.join(os.path.split(file)[0], os.path.split(file)[-1].split(".")[0])))

        export_pair_list.extend(pair_list)
        export_gro_index.extend(gro_index)

        return export_pair_list, export_gro_index


# Generate GROMACS-compatible protein indices to concatenate with guanidinium index list
def protein_gromacs_index(file):
    """
    Generate GROMACS-compatible index list for final index list export
    :param file: File name to be used, requires a file matched from original pattern
    :return: List in GROMACS-compatible form
    """
    gro_file = file_open(file)
    _, _, _, _, resname, _, _, _, _, _ = gro_parse(gro_file)
    prot_indices, number_of_prot = find_protein_chain(resname, file, silent=False)
    prot_index_list = []

    if len(prot_indices) == 1:
        prot_indices.insert(0, 0)

    for n in range(number_of_prot):
        index_string = ''
        prot_index_list.append("[ Protein_Chain_%s ]" % (n + 1))

        if n == 0:
            for i in range(prot_indices[n], prot_indices[n + 1] + 1):
                index_string = index_string + str(i + 1).rjust(6)

        else:
            for i in range(prot_indices[n] + 1, prot_indices[n + 1] + 1):
                index_string = index_string + str(i + 1).rjust(6)
        protein_list_temp = textwrap.wrap(index_string, 60)
        protein_list_temp2 = [text.ljust(60) for text in protein_list_temp]
        prot_index_list.extend(protein_list_temp2)
        prot_index_list.append('')

    return prot_index_list

# Program launch
if __name__ == "__main__":
    execution_start = time.perf_counter()
    program_launch()

    # Command-line argument parser
    args = cmd_parse(sys.argv[1:])
    print("Command line arguments:", ' '.join(sys.argv[0:]))

    # Override default cutoff information if supplied
    if args.coarse_cutoff:
        coarse_cutoff = float(args.coarse_cutoff)
    if args.short_cutoff:
        short_cutoff = float(args.short_cutoff)

    # Print cutoff information
    print("Using coarse cutoff value of: %s nm" % coarse_cutoff)
    print("Using residue cutoff value of: %s nm" % short_cutoff)

    # Find all GRO files contained within folder
    gro_files = file_pattern(args.file_pattern, args.directory)

    # Primary output data lists
    time_list = []
    gro_index_list = []

    # Run in single-threaded or multi-threaded mode
    if args.st:
        for file in gro_files:
            data_single_thread = main_single_thread(file)
            time_list.append(data_single_thread[0])
            gro_index_list.extend(data_single_thread[1])
    else:
        # Execute primary program with multiple processes handling individual files
        with concurrent.futures.ProcessPoolExecutor() as executor:
            executed_process = [executor.submit(main_single_thread, process_file) for process_file in gro_files]

        for f in concurrent.futures.as_completed(executed_process):
            time_list.append(f.result()[0])
            gro_index_list.extend(f.result()[1])

    # Formatter function for output list
    max_list_length = max(len(x) for x in time_list)
    for i in range(len(time_list)):
        if len(time_list[i]) < max_list_length:
            time_list[i].extend([''] * (max_list_length - len(time_list[i])))

    # Final list formatter
    final_list = []
    final_list.extend([''] * max_list_length)
    for i in range(len(time_list)):
        for n in range(max_list_length):
            if time_list[i][n] == '':
                final_list[n] = final_list[n] + "\t\t\t\t\t\t" + "\t"
            else:
                final_list[n] = final_list[n] + time_list[i][n] + "\t"

    # Check if output directory exists, and if not, create it
    if not os.path.isdir("output"):
        os.mkdir("output")

    # Build final index file, inlcuding protein chains
    final_index_list = protein_gromacs_index(gro_files[0])
    final_index_list.extend(gro_index_list)

    # Write GROMACS-compatible index list to file
    write_file(final_index_list, os.path.join("output", args.output.split('.')[0] + '.ndx'))

    # Write final list to file
    write_file(final_list, os.path.join("output", args.output.split('.')[0] + '.txt'))
    if args.csv:
        csv_list = [i.replace('\t', ',') for i in final_list]
        write_file(csv_list, os.path.join("output", args.output.split('.')[0] + '.csv'))

    # Print performance data
    execution_end = time.perf_counter()
    print("Program executed in %s seconds." % round(execution_end - execution_start, 5))
