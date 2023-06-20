"""
qmrebind.py

Qmrebind reparametrizes the partial charges on atoms using a QMMM ONIOM 
calculation on the ligand in the bound state of the receptor.
"""

# Standard library imports
from itertools import zip_longest
from sys import stdout
import itertools
import shutil
import re
import os

# Related third party imports
from matplotlib.backends.backend_pdf import PdfPages
from biopandas.pdb import PandasPdb
import matplotlib.pyplot as plt
from PyPDF2 import PdfMerger
import pandas as pd
import numpy as np
import parmed
import simtk

# Local application-specific imports
import qmrebind.defaults as defaults

# TODO: set overwrite to False for production code
#def make_work_dir(work_dir=None, overwrite=False):
def make_work_dir(files_to_copy, work_dir=None, overwrite=True):
    """
    Make a working directory to use for temporary files, log files, etc.
    
    Parameters
    ----------
    files_to_copy : list
        A list of files to copy to the working directory for the calculation.
        For instance, any forcefield or structure files.

    work_dir : str or None, Default None
        The working directory to create and run calculations in. If 'None',
        will default to 'work_dir_qmrebind'.
        
    overwrite : bool, Default True
        If False, an existing work directory of the same name as 'work_dir'
        will raise an error. If True, any existing directory will be removed
        and overwritten.
        
    """
    if work_dir is None:
        work_dir = "work_dir_qmrebind"
    work_dir_abs = os.path.abspath(work_dir)
    if overwrite and os.path.exists(work_dir_abs):
        shutil.rmtree(work_dir_abs)
    
    assert not os.path.exists(work_dir_abs), \
        f"Work directory {work_dir} already exists. Please remove directory, "\
        "choose a different work directory, or enable overwrite flag in this "\
        "function."
    os.mkdir(work_dir_abs)
    for myfile in files_to_copy:
        new_file_name = os.path.join(work_dir_abs, os.path.basename(myfile))
        print(f"Copying file '{myfile}' to '{work_dir}'.")
        shutil.copyfile(myfile, new_file_name)
        
    print(f"Moving to directory: {work_dir}.")
    os.chdir(work_dir_abs)
    return

def delete_files(to_delete):
    """
    Delete a list of files
    """
    for myfile in to_delete:
        if os.path.exists(myfile):
            os.remove(myfile)
    
    return

def get_anchor_list():

    """
    The function creates a list of all the directories
    , which starts with the word "anchor_". These are the
    directories where the SEEKR2 simulations are carried out.

    """
    current_pwd = os.getcwd()
    anchor_list = []
    for anchor in os.listdir(current_pwd):
        if anchor.startswith("anchor_"):
            anchor_list.append(anchor)
    return anchor_list


def group_str(iterable, n, fillvalue=None):

    """
    The function collects the data into fixed-length
    chunks or blocks.

    Parameters
    ----------
    lst : list
        Input list.

    n : int
        Length of the block.

    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def list_to_dict(lst):

    """
    The function converts an input list with mapped characters
    (every odd entry is the key of the dictionary, and every even
    entry adjacent to the odd entry is its corresponding
    value) to a dictionary.

    Parameters
    ----------
    lst : list
        Input list.

    Returns
    -------
    res_dct : dict
        A dictionary with every element mapped with
        its successive element starting from index 0.

    """
    res_dct = {lst[i]: lst[i + 1] for i in range(0, len(lst), 2)}
    return res_dct


def extract_charges_str(str_):

    """
    The function extracts the charges from the topology
    files and returns a list of these charges in multiples
    of five elements in each list.

    Parameters
    ----------
    str_ : str
        String in each line of the file.

    Returns
    -------
    extract : list
        A list of charges in multiples of five elements in
        each list.

    """
    extract = []
    for elem in str_.split():
        try:
            extract.append(float(elem))
        except ValueError:
            pass
    return extract

def prepare_pdb(input_pdb):

    """
    The function uses the pdb4amber module of AMBER to remove
    any "TER" and "CONECT" keyword in the PDB file. The function
    removes any solvent and ions present in the PDB file and
    replaces any "HETATM" keyword with the "ATOM" keyword to
    maintain homogeneity and backs up the original PDB file
    with an "_before_qmmm" extension.

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file.

    """
    input_pdb_base = os.path.splitext(input_pdb)[0]
    input_before_qmmm = input_pdb_base  + "_before_qmmm.pdb"
    shutil.copyfile(input_pdb, input_before_qmmm)
    intermediate_file_I = input_pdb_base + "_intermediate_I.pdb"
    intermediate_file_II = input_pdb_base + "_intermediate_II.pdb"
    intermediate_file_II_base = input_pdb_base + "_intermediate_II"
    with open(input_pdb) as f1, open(intermediate_file_I, "w") as f2:
        for line in f1:
            if not any(ion in line for ion in defaults.IONS):
                f2.write(line)
    command = (
        f"pdb4amber -i {intermediate_file_I} -o {intermediate_file_II}" \
        " --noter --no-conect --dry"
    )
    os.system(command)
    to_delete = (
        intermediate_file_II_base + "_nonprot.pdb",
        intermediate_file_II_base + "_renum.txt",
        intermediate_file_II_base + "_sslink",
        intermediate_file_II_base + "_water.pdb",
    )
    delete_files(to_delete)
    with open(intermediate_file_II, "r") as f1:
        filedata = f1.readlines()
    filedata2 = []
    for line in filedata:
        line2 = line.replace("HETATM", "ATOM  ")
        if line.startswith("CONECT"):
            continue
        filedata2.append(line2)
        
    with open(input_pdb, "w") as f2:
        f2.writelines(filedata2)
    delete_files([intermediate_file_I, intermediate_file_II])
    return

def strip_topology(forcefield_file):

    """
    The function removes the parameters for the solvent
    and the ions in the topology file. It does this using
    cpptraj

    Parameters
    ----------
    forcefield_file: str
        User-defined topology file (prmtop/parm7 file)

    """
    ff_base = os.path.splitext(forcefield_file)[0]
    ff_ext = os.path.splitext(forcefield_file)[1]
    ff_before_qmmm = ff_base  + f"_before_qmmm{ff_ext}"
    shutil.copyfile(forcefield_file, ff_before_qmmm)
    stripped_parm_file = f"{ff_base}_stripped{ff_ext}"
    cpptraj_input_filename = "strip_topology.cpptraj"
    with open(cpptraj_input_filename, "w") as f:
        f.write(f"parm {forcefield_file}\n")
        f.write("parmstrip :WAT\n")
        for i in defaults.IONS:
            f.write(f"parmstrip :{i}\n")
        f.write(f"parmwrite out {stripped_parm_file}\n")
        f.write("run")
    command = f"cpptraj -i {cpptraj_input_filename}"
    os.system(command)
    os.rename(stripped_parm_file, forcefield_file)
    delete_files([stripped_parm_file, cpptraj_input_filename])
    return

def get_ligand_pdb(input_pdb, ligand_pdb, ligand_resname):

    """
    The function reads the PDB file, extracts the coordinate
    information for the ligand molecule and saves it into a
    new PDB file.

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file.

    ligand_pdb: str
        User-defined ligand PDB file.

    ligand_resname: str
        Three-letter name for the ligand residue.

    """
    with open(input_pdb) as f1, open(ligand_pdb, "w") as f2:
        for line in f1:
            if ligand_resname in line:
                f2.write(line)

def get_receptor_pdb(input_pdb, receptor_pdb, ligand_resname):
    """
    The function reads the PDB file, extracts the coordinate
    information for the receptor and saves it into a PDB file.

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file.

    receptor_pdb: str
        User-defined receptor PDB coordinate file.

    ligand_resname: str
        Three-letter name for the ligand residue.

    """
    with open(input_pdb) as f1, open(receptor_pdb, "w") as f2:
        for line in f1:
            if not re.search(f"{ligand_resname}|CRYST|WAT|HOH", line) \
                    and not any(ion in line for ion in defaults.IONS):
                f2.write(line)
    
    return

def get_number_pdb_atoms(input_pdb):

    """
    The function counts the total number of atoms in the system,
    including the solvent and the ions, if any.

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file.

    """
    ppdb = PandasPdb().read_pdb(input_pdb)
    return len(ppdb.df["ATOM"]) + len(ppdb.df["HETATM"])

def get_indices_qm_region(input_pdb, ligand_resname):

    """
    The function returns the atom indices of the QM region, i.e.,
    the ligand molecule (beginning from 0).

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file.

    ligand_resname: str
        Three-letter name for the ligand residue.

    """
    ppdb = PandasPdb()
    ppdb.read_pdb(input_pdb)
    df = ppdb.df["ATOM"][["atom_number", "residue_number", "residue_name"]]
    indices = list(np.where(df["residue_name"] == ligand_resname))
    atom_indices = list(indices[0])
    return atom_indices

def get_indices_qm2_region(ligand_pdb, receptor_pdb, cut_off_distance):

    """
    The function returns the lists of residues and its indices (beginning
    from 0) that surrounds the QM region within the user-specified cut-off
    distance.

    Parameters
    ----------
    ligand_pdb: str
        User-defined ligand PDB file.

    receptor_pdb: str
        User-defined receptor PDB coordinate file.

    cut_off_distance: int
        Cut-off distance for the QM2 region within the
        vicinity of the QM region.

    """
    ppdb_lig = PandasPdb()
    ppdb_lig.read_pdb(ligand_pdb)
    coords_lig = ppdb_lig.df["ATOM"][["x_coord", "y_coord", "z_coord"]]
    ligand_coords = np.array(coords_lig.values.tolist())
    ppdb_rec = PandasPdb()
    ppdb_rec.read_pdb(receptor_pdb)
    
    receptor_atom_list = []
    for i in range(len(ligand_coords)):
        reference_point = ligand_coords[i]
        distances = ppdb_rec.distance(xyz=reference_point, records=("ATOM"))
        all_within_distance = ppdb_rec\
            .df["ATOM"][distances < float(cut_off_distance)]
        receptor_df = all_within_distance["atom_number"]
        receptor_list = receptor_df.values.tolist()
        receptor_atom_list.append(receptor_list)
    receptor_atom_list = list(set(list(itertools.chain(*receptor_atom_list))))
    receptor_atom_list.sort()
    df = ppdb_rec.df["ATOM"][["atom_number", "residue_number", "residue_name"]]
    index_list = []
    for i in receptor_atom_list:
        indices = np.where(df["atom_number"] == i)
        indices = list(indices)[0]
        indices = list(indices)
        index_list.append(indices)
    index_list = list(itertools.chain.from_iterable(index_list))
    df1 = df.iloc[index_list]
    resid_num = list(df1.residue_number.unique())
    atom_index_list = []
    for i in resid_num:
        atom_indices = list(np.where(df["residue_number"] == i))
        atom_indices = list(atom_indices[0])
        atom_index_list.append(atom_indices)
    receptor_atom_index_list = list(
        itertools.chain.from_iterable(atom_index_list))
    return (resid_num, receptor_atom_index_list)

# TODO: return to this function to determine what optimizations, if any, may
# be done to it.
def prepare_orca_pdb(
    input_pdb, ligand_pdb, orca_pdb, ligand_resname, receptor_pdb, 
    cut_off_distance):

    """
    ORCA QM/QM2/MM input file reads the PDB file. The occupancy
    column of the PDB file needs to be readjusted according to
    the definition of the QM/QM2/MM layer. An occupancy value of
    1.00 is defined for the QM region, 2.00 for the QM2 region,
    and 0.00 for the MM region. The function prepares the PDB file
    needed for the QM/QM2/MM calculations by saving a PDB file with
    the above considerations.

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file.

    ligand_pdb: str
        User-defined ligand PDB file.

    orca_pdb: str
        User-defined ORCA PDB file.

    ligand_resname: str
        Three-letter name for the ligand residue.

    recepor_pdb: str
        User-defined receptor PDB coordinate file.

    cut_off_distance: int
        Cut-off distance for the QM2 region within the
        vicinity of the QM region.

    """
    ppdb = PandasPdb().read_pdb(input_pdb)
    # Use occupancy value of 1.00 for the QM region and 2.00 for the QM2 region
    ligand_indices = get_indices_qm_region(
        input_pdb=input_pdb, ligand_resname=ligand_resname
    )
    ppdb.df["ATOM"].loc[ligand_indices[0] : ligand_indices[-1], "occupancy"] \
        = 1.00
    receptor_residues, receptor_indices = get_indices_qm2_region(
        ligand_pdb=ligand_pdb, receptor_pdb=receptor_pdb, 
        cut_off_distance=cut_off_distance
    )
    ranges = sum(
        (list(t) for t in zip(receptor_indices, receptor_indices[1:]) if t[0] + 1 != t[1]), []
    )
    iranges = iter(receptor_indices[0:1] + ranges + receptor_indices[-1:])
    receptor_input_indices = " ".join([str(n) + ":" + str(next(iranges)) for n in iranges])
    receptor_input_indices_list = [
        int(s) for s in re.findall(r"\b\d+\b", receptor_input_indices)
    ]
    split_list = [
        receptor_input_indices_list[i : i + 2]
        for i in range(0, len(receptor_input_indices_list), 2)
    ]
    for i in split_list:
        ppdb.df["ATOM"].loc[i[0] : i[1], "occupancy"] = 2.00
    ppdb.to_pdb(path=orca_pdb, records=None, gz=False, append_newline=True)
    return

def get_amber_to_orca_prms(forcefield_file):

    """
    The function uses the ORCA built-in command to convert
    the AMBER topology file, i.e., parm7/prmtop files to the
    ORCA forcefield file.

    Parameters
    ----------
    forcefield_file: str
        User-defined topology file (prmtop/parm7 file)

    """
    command = f"orca_mm -convff -AMBER {forcefield_file}"
    os.system(command)
    return


def get_orca_input(
    nprocs,
    maxiter,
    qm_method,
    qm_basis_set,
    qm2_method,
    optimization,
    frequency_calculation,
    qm_charge_scheme,
    qm2_charge_scheme,
    qm2_charge,
    qm2_mult,
    forcefield_file,
    input_pdb,
    ligand_resname,
    orca_pdb,
    orca_input_file,
    ligand_pdb,
    receptor_pdb,
    cut_off_distance,
    qm_charge,
    qm_mult,
):

    """
    The function creates an input file for the ORCA QM/QM2/MM
    calculation. It takes into account all the keywords required
    for performing the multiscale calculation.

    Parameters
    ----------
    nprocs: int
        Number of processors to be used for the MPI
        enabled calculation.

    maxiter: int
        Maximum number of iterations needed for the self
        consistent field to converge.

    qm_method: str
        QM theory to be implemented for the high layer,
        i.e., the QM layer. The options for the QM method
        include the perturbation method, reference method,
        local gradient corrected method, hybrid method,
        and meta GGA/hybrid meta GGA method.

    qm_basis_set: str
        Basis set for the QM level of theory. The options
        include pople basis set, pople polarized basis set,
        pople polarized diffused basis set, correlation
        consistent basis set, DEF2 basis set, DEF2 diffused
        basis set, auxiliary coulomb fitted basis_set, auxiliary
        coulomb fitted and exchange basis set, and auxiliary
        correlation consistent basis set.

    qm2_method: str
        QM theory to be implemented for the lower layer,
        i.e., the QM2 layer. The options for the QM2 method
        include the semiempirical method, the tight binding DFT
        method, and the composite method.

    optimization: str
        The string is set to "True" if the user chooses to
        optimize the geometry of the QM region. If not,
        the string is set to "False".

    frequency_calculation: str
        The string is set to "True" if the user chooses
        a frequency calculation for the QM region. If not,
        the string is set to "False".

    qm_charge_scheme: str
        Charge scheme for the calculation of QM charges for
        the QM region. The options include Hirshfeld,
        CHELPG, Mulliken, and the Loewdin charge calculation
        methods.

    qm2_charge_scheme: str
        Charge scheme for the calculation of QM charges for
        the QM2 region. The options include Hirshfeld,
        CHELPG, Mulliken, and the Loewdin charge calculation
        methods.

    qm2_charge: int
        Charge of the QM2 (low) region. The user does not have
        to provide the charge of the QM region since it is
        calculated by ORCA itself. The user also does not
        have to provide the charge of the MM region since it is
        determined by the ORCA force field file.

    qm2_mult: int
        The multiplicity of the QM2 (low) region. The user does
        not have to provide the multiplicity of the QM region
        since ORCA itself calculates it. The user also
        does not have to provide the multiplicity of the MM
        region since the ORCA forcefield file determines it.

    forcefield_file: str
        User-defined topology file (prmtop/parm7 file).

    input_pdb: str
        User-defined PDB file.

    ligand_resname: str
        Three-letter name for the ligand residue.

    orca_pdb: str
        User-defined ORCA PDB file.

    orca_input_file: str
        User-defined ORCA input file.

    ligand_pdb: str
        User-defined ligand PDB file.

    receptor_pdb: str
        User-defined receptor PDB coordinate file.

    cut_off_distance: int
        Cut-off distance for the QM2 region within the
        vicinity of the QM region.

    qm_charge: int
        Charge of the QM region.

    qm_mult: int
        Multiplicity of the QM region.

    implicit_solv: str
        The Universal Solvation Model (SMD) method, an implicit
        solvation model an improvement over the Conductor-like
        Polarizable Continuum Model (CPCM), since it uses the
        full solute electron density to compute the cavity-dispersion
        contribution instead of the area only.

    """
    line_0 = f"%PAL NPROCS {nprocs} END"
    line_1 = f"%scf MaxIter {maxiter} END"
    line_2 = f"!QM/{qm2_method}/MM {qm_method} {qm_basis_set}"
    if optimization:
        line_2 += " OPT"
        
    if frequency_calculation:
        line_2 += " NUMFREQ"
        
    line_3 = f"!{qm_charge_scheme}"
    line_4 = "%QMMM"
    line_5 = "PRINTLEVEL 5"
    line_6 = f"ORCAFFFilename \"{forcefield_file}.ORCAFF.prms\""
    ligand_indices = get_indices_qm_region(
        input_pdb=input_pdb, ligand_resname=ligand_resname
    )
    receptor_residues, receptor_indices = get_indices_qm2_region(
        ligand_pdb=ligand_pdb, receptor_pdb=receptor_pdb, 
        cut_off_distance=cut_off_distance
    )
    ranges = sum(
        (list(t) for t in zip(receptor_indices, receptor_indices[1:]) \
            if t[0] + 1 != t[1]), []
    )
    iranges = iter(receptor_indices[0:1] + ranges + receptor_indices[-1:])
    receptor_input_indices = " ".join([str(n) + ":" + str(next(iranges)) \
                                       for n in iranges])
    line_7 = f"QMATOMS {{{ligand_indices[0]}:{ligand_indices[-1]}}} END"
    line_8 = f"QM2ATOMS {{{receptor_input_indices}}} END"
    line_9 = f"CHARGE_MEDIUM {qm2_charge}"
    line_10 = f"MULT_MEDIUM {qm2_mult}"
    line_11 = f"CHARGE_METHOD {qm2_charge_scheme}"
    line_12 = "Use_QM_InfoFromPDB TRUE"
    line_13 = "Use_Active_InfoFromPDB TRUE END"
    line_14 = f"*PDBFILE {qm_charge} {qm_mult} {orca_pdb} "
    # TODO: do this more efficiently
    if nprocs == 0:
        cmd_list = []
    else:
        cmd_list = [line_0]
    
    cmd_list += [line_1, line_2, line_3, line_4, line_5, line_6, line_7,
                 line_8, line_9, line_10, line_11, line_12, line_13, line_14]
    commands = "\n".join(cmd_list)
    print("Writing ORCA input file:", orca_input_file)
    with open(orca_input_file, "w") as f:
        f.write(commands)
    return

def add_xtb_inputs(
    etemp,
    solvation,
    solvent,
    accuracy,
    xtb_memory,
    xtb_nprocs,
    orca_input_file,
    XTB_add_inputs,
):

    """
    The function adds additional XTB commands in the ORCA input file
    for QM2 calculations.

    Parameters
    ----------
    etemp: int
        Electronic temperature of the system.

    solvation: str
        Implicit solvation of the system.

    solvent: str
        Solvent for the implicit solvation of the
        system.

    accuracy: int
        Accuracy value for XTB optimization, default
        is ORCA’s accuracy x 1.e6.

    xtb_memory: int
        Memory in MB reserved for XTB calculation.

    xtb_nprocs: int
        Number of processors used for running XTB
        calculation.

    orca_input_file: str
        User-defined ORCA input file.

    XTB_add_inputs: bool
        if True, the function will add additional XTB 
        commands in the ORCA input file, and if False,
        the ORCA input file will remain unchanged.

    """
    line_1 = f"%XTB\n"
    line_2 = f"ETEMP        {etemp}\n"
    line_3 = f"DOALPB       {solvation}\n"
    line_4 = f"ALPBSOLVENT  \"{solvent}\"\n"
    line_5 = f"ACCURACY     {accuracy}\n"
    line_6 = f"MAXCORE      {xtb_memory}\n"
    line_7 = f"NPROCS       {xtb_nprocs}\n"
    line_8 = f"END\n"
    with open(orca_input_file, "r") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if "*PDBFILE" in line:
            to_begin = int(i)
    last_line_index = len(lines)-1
    lines_I = lines[: to_begin]
    lines_III = lines[to_begin : last_line_index + 1]
    lines_list = [line_1, line_2, line_3, line_4, line_5, line_6, line_7, 
                  line_8]
    all_lines = "".join(lines_list)
    if XTB_add_inputs:
        with open(orca_input_file, "w") as f:
            for line in lines_I:
                f.write(line)
            f.write(all_lines)
            for line in lines_III:
                f.write(line)
    else:
        with open(orca_input_file, "w") as f:
            for line in lines_I:
                f.write(line)
            for line in lines_III:
                f.write(line)
    
    return

def run_orca_qmmm(orca_dir_pwd, orca_input_file, orca_out_file):

    """
    The function runs the ORCA QM/QM2/MM calculations
    for the system.

    Parameters
    ----------
    orca_dir_pwd: str
        PWD of the directory where ORCA is installed.

    orca_input_file: str
        User-defined ORCA input file.

    orca_out_file: str
        User-defined ORCA output file.

    """
    orca_cmd = os.path.join(orca_dir_pwd, "orca")
    command = f"{orca_cmd} {orca_input_file} > {orca_out_file}"
    print("Running command:", command)
    os.system(command)
    return

def get_qm_charges(
    orca_out_file, qm_charge_file, input_pdb, ligand_resname, qm_charge_scheme
):

    """
    The function extracts the charges of the QM region of the system
    by parsing the output file of the ORCA QM/QM2/MMM calculations.


    Parameters
    ----------
    orca_out_file: str
        User-defined ORCA output file.

    qm_charge_file: str
        User-defined text file where the charges for the
        atoms of the QM region are saved.

    input_pdb: str
        User-defined PDB file.

    ligand_resname: str
        Three-letter name for the ligand residue.

    qm_charge_scheme: str
        Charge scheme for the calculation of QM charges for
        the QM region. The options include Hirshfeld,
        CHELPG, Mulliken, and the Loewdin charge calculation
        methods.

    """
    # TODO: make more efficient and apply DRY
    with open(orca_out_file, "r") as f:
        lines = f.readlines()
    if qm_charge_scheme == "HIRSHFELD":
        for i in range(len(lines)):
            if "HIRSHFELD ANALYSIS" in lines[i]:
                to_begin = int(i)
                to_end = len(
                    get_indices_qm_region(
                        input_pdb=input_pdb, ligand_resname=ligand_resname
                    )
                )
        charges = lines[to_begin + 7 : to_begin + 7 + to_end]
        with open("temp.txt", "w") as f:
            for i in charges:
                f.write(i)
        df_temp = pd.read_csv("temp.txt", delimiter=r"\s+", header=None)
        df_temp.columns = ["Index", "Atom", "Charge", "Spin"]
        df_temp["Charge"].to_csv(qm_charge_file, index=False, header=False, sep=" ")
        os.system("rm -rf temp.txt")
    if qm_charge_scheme == "CHELPG":
        for i in range(len(lines)):
            if "CHELPG Charges       " in lines[i]:
                to_begin = int(i)
            if "CHELPG charges calculated" in lines[i]:
                to_end = int(i)
        charges = lines[to_begin + 2 : to_end - 4]
        charge_list = []
        for i in range(len(charges)):
            charge_list.append(charges[i].strip().split())
        charge_list_value = []
        for i in range(len(charge_list)):
            charge_list_value.append(float(charge_list[i][3]))
        data = list(charge_list_value)
        df_charge = pd.DataFrame(data, columns=["Charge"])
        df_charge.to_csv(qm_charge_file, index=False, header=False, sep=" ")
    if qm_charge_scheme == "MULLIKEN":
        for i in range(len(lines)):
            if "MULLIKEN ATOMIC CHARGES" in lines[i]:
                to_begin = int(i)
                to_end = len(
                    get_indices_qm_region(
                        input_pdb=input_pdb, ligand_resname=ligand_resname
                    )
                )
        charges = lines[to_begin + 2 : to_begin + 2 + to_end]
        with open("temp.txt", "w") as f:
            for i in charges:
                f.write(i)
        df_temp = pd.read_csv("temp.txt", sep=":", header=None)
        df_temp.columns = ["Index_Atom", "Charge"]
        df_temp["Charge"].to_csv(qm_charge_file, index=False, header=False, sep=" ")
        os.system("rm -rf temp.txt")
    if qm_charge_scheme == "LOEWDIN":
        for i in range(len(lines)):
            if "LOEWDIN ATOMIC CHARGES" in lines[i]:
                to_begin = int(i)
                to_end = len(
                    get_indices_qm_region(
                        input_pdb=input_pdb, ligand_resname=ligand_resname
                    )
                )
        charges = lines[to_begin + 2 : to_begin + 2 + to_end]
        with open("temp.txt", "w") as f:
            for i in charges:
                f.write(i)
        df_temp = pd.read_csv("temp.txt", sep=":", header=None)
        df_temp.columns = ["Index_Atom", "Charge"]
        df_temp["Charge"].to_csv(qm_charge_file, index=False, header=False, sep=" ")
        os.system("rm -rf temp.txt")


def get_ff_charges(forcefield_file, ff_charges_file, input_pdb):

    """
    The function extracts the charges of the QM region of the system
    by parsing the output file of the ORCA QM/QM2/MM calculations.


    Parameters
    ----------
    forcefield_file: str
        User-defined topology file (prmtop/parm7 file).

    ff_charges_file: str
        User-defined file where the charges defined in the
        topology file (prmtop/parm7 file) is stored in
        a text file.

    input_pdb: str
        User-defined PDB file.

    """
    with open(forcefield_file, "r") as f:
        lines = f.readlines()
    no_atoms = get_number_pdb_atoms(input_pdb=input_pdb)
    if no_atoms % 5 == 0:
        lines_to_select = int(no_atoms / 5)
    else:
        lines_to_select = int(no_atoms // 5 + 1)
    for i in range(len(lines)):
        if "%FLAG CHARGE" in lines[i]:
            to_begin = int(i)
            to_end = int(i) + lines_to_select
    charges = lines[to_begin + 2 : to_end + 2]
    charge_list = []
    for i in range(len(charges)):
        charge_list.append(charges[i].strip().split())
    charge_list = [item for sublist in charge_list for item in sublist]
    df_charge = pd.DataFrame(charge_list, columns=["Charge"])
    df_charge.to_csv(ff_charges_file, index=False, header=False, sep=" ")


def get_ff_qm_charges(
    qm_charge_file,
    ff_charges_file,
    ff_charges_qm_fmt_file,
    input_pdb,
    ligand_resname,
):

    """
    The function uses the previous file that stored the charges of the
    the system obtained from the AMBER force field file, converts
    it into AMBER format, i.e. %FORMAT(5E16.8), and then replaces
    the charges of the QM region with the charges of the QM region
    obtained from the ORCA QM/QM2/MM calculations.

    Parameters
    ----------
    qm_charge_file: str
        User-defined text file where the charges for the
        atoms of the QM region are saved.

    ff_charges_file: str
        User-defined file where the charges defined in the
        topology file (prmtop/parm7 file) is stored in
        a text file.

    ff_charges_qm_fmt_file: str
        User-defined file where the replaced charges in the
        AMBER format is stored.

    input_pdb: str
        User-defined PDB file.

    ligand_resname: str
        Three-letter name for the ligand residue.

    """
    df_qm_charges = pd.read_csv(qm_charge_file, header=None, delimiter=r"\s+")
    df_qm_charges.columns = ["Charge"]
    df_qm_charges["Index"] = get_indices_qm_region(
        input_pdb=input_pdb, ligand_resname=ligand_resname
    )
    df_qm_charges = df_qm_charges[["Index", "Charge"]]
    qm_charge_list = df_qm_charges["Charge"].values.tolist()
    qm_charge_list = [
        i * 18.2223 for i in qm_charge_list
    ]  # charge conversion factor=18.2223
    # print(qm_charge_list)
    qm_charge_index_list = df_qm_charges["Index"].values.tolist()
    # print(qm_charge_index_list)
    dict_qm = {
        qm_charge_index_list[i]: qm_charge_list[i]
        for i in range(len(qm_charge_index_list))
    }
    df_ff_charges = pd.read_csv(ff_charges_file, header=None, delimiter=r"\s+")
    df_ff_charges.columns = ["Charge"]
    ff_charges = df_ff_charges["Charge"].values.tolist()
    for i in qm_charge_index_list:
        ff_charges[i] = dict_qm[i]
    ff_charges_fmt = []
    for group in group_str(ff_charges, 5):
        ff_charges_fmt.append(
            " ".join("{: .8E}".format(x) for x in group if x is not None) + "\n"
        )
    with open(ff_charges_qm_fmt_file, "w") as f:
        for i in ff_charges_fmt:
            f.write(" " + i)


def get_qmrebind_parm(forcefield_file, input_pdb, ff_charges_qm_fmt_file):

    """
    The function replaces the defined charges in the topology file
    (parm7/prmtop) file with the new set of charges obtained from
    the ORCA QM/QM2/MM calculations.

    Parameters
    ----------
    forcefield_file: str
        User-defined topology file (prmtop/parm7 file).

    input_pdb: str
        User-defined PDB file.

    ff_charges_qm_fmt_file: str
        User-defined file where the replaced charges in the
        AMBER format is stored.

    """

    if forcefield_file[-6:] == "prmtop":
        command = (
            "cp -r "
            + forcefield_file
            + " "
            + forcefield_file[:-7]
            + "_before_charge_replacement.prmtop"
        )
        os.system(command)
    if forcefield_file[-6:] == ".parm7":
        command = (
            "cp -r "
            + forcefield_file
            + " "
            + forcefield_file[:-6]
            + "_before_charge_replacement.parm7"
        )
        os.system(command)
    with open(forcefield_file, "r") as f1:
        ff_lines = f1.readlines()
    for i in range(len(ff_lines)):
        if "%FLAG CHARGE" in ff_lines[i]:
            to_begin = 0
            to_end = int(i)
    parm_lines_a = ff_lines[0 : to_end + 2]
    no_atoms = get_number_pdb_atoms(input_pdb=input_pdb)
    if no_atoms % 5 == 0:
        ff_lines_to_select = int(no_atoms / 5)
    else:
        ff_lines_to_select = int(no_atoms // 5 + 1)
    parm_lines_b = ff_lines[to_end + ff_lines_to_select + 2 :]
    with open(ff_charges_qm_fmt_file, "r") as f2:
        qm_ff_lines = f2.readlines()
    command = "rm -rf " + forcefield_file
    os.system(command)
    with open(forcefield_file, "w") as f:
        for i in parm_lines_a:
            f.write(i)
        for j in qm_ff_lines:
            f.write(j)
        for k in parm_lines_b:
            f.write(k)


def get_qmrebind_parm_solvent(input_pdb, forcefield_file, ff_charges_file):

    """
    The function replaces the defined charges in the original topology
    file (parm7/prmtop) file with the new set of charges obtained from
    the ORCA QM/QM2/MM calculations as defined in the topology file
    with no solvent.

    Parameters
    ----------

    input_pdb: str
        User-defined PDB file.

    forcefield_file: str
        User-defined topology file (prmtop/parm7 file).

    ff_charges_file: str
        User-defined file where the charges defined in the
        topology file (prmtop/parm7 file) is stored in
        a text file.

    """
    receptorligand_before_qmmm_pdb = input_pdb[:-4] + "_before_qmmm.pdb"
    ff_charges_file_before_qmmm = ff_charges_file[:-4] + "_before_qmmm.txt"
    ff_charges_file_after_qmmm = ff_charges_file[:-4] + "_after_qmmm.txt"
    if forcefield_file[-6:] == "prmtop":
        forcefield_file_before_qmmm = forcefield_file[:-7] + "_before_qmmm.prmtop"
    if forcefield_file[-6:] == ".parm7":
        forcefield_file_before_qmmm = forcefield_file[:-6] + "_before_qmmm.parm7"
    len_pdb_before_qmmm = get_number_pdb_atoms(input_pdb=receptorligand_before_qmmm_pdb)
    len_pdb_after_qmmm = get_number_pdb_atoms(input_pdb=input_pdb)
    get_ff_charges(
        forcefield_file=forcefield_file_before_qmmm,
        ff_charges_file=ff_charges_file_before_qmmm,
        input_pdb=receptorligand_before_qmmm_pdb,
    )
    get_ff_charges(
        forcefield_file=forcefield_file,
        ff_charges_file=ff_charges_file_after_qmmm,
        input_pdb=input_pdb,
    )
    list_ff_charges_file_before_qmmm = list(np.loadtxt(ff_charges_file_before_qmmm))
    list_ff_charges_file_after_qmmm = list(np.loadtxt(ff_charges_file_after_qmmm))
    for i in range(len(list_ff_charges_file_after_qmmm)):
        list_ff_charges_file_before_qmmm[i] = list_ff_charges_file_after_qmmm[i]
    list_ff_charges_file_before_qmmm_8e_format = []
    for i in list_ff_charges_file_before_qmmm:
        list_ff_charges_file_before_qmmm_8e_format.append("{:.8E}".format(i))
    df_charge = pd.DataFrame(
        list_ff_charges_file_before_qmmm_8e_format, columns=["Charge"]
    )
    ff_charges_file_after_qmmm_8e_format = (
        ff_charges_file[:-4] + "_after_qmmm_8e_format.txt"
    )
    df_charge.to_csv(
        ff_charges_file_after_qmmm_8e_format, index=False, header=False, sep=" "
    )
    df_ff_charges = pd.read_csv(
        ff_charges_file_after_qmmm_8e_format, header=None, delimiter=r"\s+"
    )
    df_ff_charges.columns = ["Charge"]
    ff_charges = df_ff_charges["Charge"].values.tolist()
    ff_charges_fmt = []
    for group in group_str(ff_charges, 5):
        ff_charges_fmt.append(
            " ".join("{: .8E}".format(x) for x in group if x is not None) + "\n"
        )
    ff_charges_file_after_qmmm_5_8e_format = (
        ff_charges_file[:-4] + "_after_qmmm_5_8e_format.txt"
    )
    with open(ff_charges_file_after_qmmm_5_8e_format, "w") as f:
        for i in ff_charges_fmt:
            f.write(" " + i)
    if forcefield_file[-6:] == "prmtop":
        command = (
            "cp -r "
            + forcefield_file
            + " "
            + forcefield_file[:-7]
            + "_no_solvent.prmtop"
        )
        os.system(command)
    if forcefield_file[-6:] == ".parm7":
        command = (
            "cp -r "
            + forcefield_file
            + " "
            + forcefield_file[:-6]
            + "_no_solvent.parm7"
        )
        os.system(command)
    with open(forcefield_file_before_qmmm, "r") as f1:
        ff_lines = f1.readlines()
    for i in range(len(ff_lines)):
        if "%FLAG CHARGE" in ff_lines[i]:
            to_begin = 0
            to_end = int(i)
    parm_lines_a = ff_lines[0 : to_end + 2]
    no_atoms = get_number_pdb_atoms(input_pdb=receptorligand_before_qmmm_pdb)
    if no_atoms % 5 == 0:
        ff_lines_to_select = int(no_atoms / 5)
    else:
        ff_lines_to_select = int(no_atoms // 5 + 1)
    parm_lines_b = ff_lines[to_end + ff_lines_to_select + 2 :]
    with open(ff_charges_file_after_qmmm_5_8e_format, "r") as f2:
        qm_ff_lines = f2.readlines()
    command = "rm -rf " + forcefield_file
    os.system(command)
    with open(forcefield_file, "w") as f:
        for i in parm_lines_a:
            f.write(i)
        for j in qm_ff_lines:
            f.write(j)
        for k in parm_lines_b:
            f.write(k)


def get_energy_diff_no_solvent(forcefield_file, input_pdb):

    """
    The function used the PARMED module to compare the total
    energy of the system before and after parameterization.

    Parameters
    ----------
    forcefield_file: str
        User-defined topology file (prmtop/parm7 file).

    input_pdb: str
        User-defined PDB file.

    """
    if forcefield_file[-6:] == "prmtop":
        non_qm_forcefield_file = (
            forcefield_file[:-7] + "_before_charge_replacement.prmtop"
        )
    if forcefield_file[-6:] == ".parm7":
        non_qm_forcefield_file = (
            forcefield_file[:-6] + "_before_charge_replacement.parm7"
        )
    parm_non_params = parmed.load_file(non_qm_forcefield_file, input_pdb)
    prmtop_energy_decomposition_non_params = parmed.openmm.energy_decomposition_system(
        parm_non_params, parm_non_params.createSystem()
    )
    prmtop_energy_decomposition_non_params_value = [
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_non_params
                ]
                for item in sublist
            ]
        ).get("HarmonicBondForce"),
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_non_params
                ]
                for item in sublist
            ]
        ).get("HarmonicAngleForce"),
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_non_params
                ]
                for item in sublist
            ]
        ).get("PeriodicTorsionForce"),
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_non_params
                ]
                for item in sublist
            ]
        ).get("NonbondedForce"),
    ]
    prmtop_energy_decomposition_non_params_list = [
        "HarmonicBondForce",
        "HarmonicAngleForce",
        "PeriodicTorsionForce",
        "NonbondedForce",
    ]
    df_energy_non_params = pd.DataFrame(
        list(
            zip(
                prmtop_energy_decomposition_non_params_list,
                prmtop_energy_decomposition_non_params_value,
            )
        ),
        columns=["Energy_term", "Energy_non_params"],
    )
    df_energy_non_params = df_energy_non_params.set_index("Energy_term")
    # print(df_energy_non_params)

    if forcefield_file[-6:] == "prmtop":
        qm_forcefield_file = forcefield_file[:-7] + "_no_solvent.prmtop"
    if forcefield_file[-6:] == ".parm7":
        qm_forcefield_file = forcefield_file[:-6] + "_no_solvent.parm7"

    parm_params = parmed.load_file(qm_forcefield_file, input_pdb)
    prmtop_energy_decomposition_params = parmed.openmm.energy_decomposition_system(
        parm_params, parm_params.createSystem()
    )
    prmtop_energy_decomposition_params_value = [
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_params
                ]
                for item in sublist
            ]
        ).get("HarmonicBondForce"),
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_params
                ]
                for item in sublist
            ]
        ).get("HarmonicAngleForce"),
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_params
                ]
                for item in sublist
            ]
        ).get("PeriodicTorsionForce"),
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_params
                ]
                for item in sublist
            ]
        ).get("NonbondedForce"),
    ]
    prmtop_energy_decomposition_params_list = [
        "HarmonicBondForce",
        "HarmonicAngleForce",
        "PeriodicTorsionForce",
        "NonbondedForce",
    ]
    df_energy_params = pd.DataFrame(
        list(
            zip(
                prmtop_energy_decomposition_params_list,
                prmtop_energy_decomposition_params_value,
            )
        ),
        columns=["Energy_term", "Energy_params"],
    )
    df_energy_params = df_energy_params.set_index("Energy_term")
    # print(df_energy_params)
    df_compare = pd.concat([df_energy_non_params, df_energy_params], axis=1)
    df_compare["Energy_difference"] = df_compare["Energy_params"].sub(
        df_compare["Energy_non_params"], axis=0
    )
    print(df_compare)


def get_energy_diff_solvent(forcefield_file, input_pdb):

    """
    The function used the PARMED module to compare the total
    energy of the system before and after parameterization.

    Parameters
    ----------
    forcefield_file: str
        User-defined topology file (prmtop/parm7 file).

    input_pdb: str
        User-defined PDB file.

    """
    if forcefield_file[-6:] == "prmtop":
        non_qm_forcefield_file = forcefield_file[:-7] + "_before_qmmm.prmtop"
    if forcefield_file[-6:] == ".parm7":
        non_qm_forcefield_file = forcefield_file[:-6] + "_before_qmmm.parm7"
    pdb_file = receptorligand_before_qmmm_pdb = input_pdb[:-4] + "_before_qmmm.pdb"
    parm_non_params = parmed.load_file(non_qm_forcefield_file, pdb_file)
    prmtop_energy_decomposition_non_params = parmed.openmm.energy_decomposition_system(
        parm_non_params, parm_non_params.createSystem()
    )
    prmtop_energy_decomposition_non_params_value = [
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_non_params
                ]
                for item in sublist
            ]
        ).get("HarmonicBondForce"),
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_non_params
                ]
                for item in sublist
            ]
        ).get("HarmonicAngleForce"),
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_non_params
                ]
                for item in sublist
            ]
        ).get("PeriodicTorsionForce"),
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_non_params
                ]
                for item in sublist
            ]
        ).get("NonbondedForce"),
    ]
    prmtop_energy_decomposition_non_params_list = [
        "HarmonicBondForce",
        "HarmonicAngleForce",
        "PeriodicTorsionForce",
        "NonbondedForce",
    ]
    df_energy_non_params = pd.DataFrame(
        list(
            zip(
                prmtop_energy_decomposition_non_params_list,
                prmtop_energy_decomposition_non_params_value,
            )
        ),
        columns=["Energy_term", "Energy_non_params"],
    )
    df_energy_non_params = df_energy_non_params.set_index("Energy_term")
    # print(df_energy_non_params)

    pdb_file = receptorligand_before_qmmm_pdb = input_pdb[:-4] + "_before_qmmm.pdb"
    parm_params = parmed.load_file(forcefield_file, pdb_file)
    prmtop_energy_decomposition_params = parmed.openmm.energy_decomposition_system(
        parm_params, parm_params.createSystem()
    )
    prmtop_energy_decomposition_params_value = [
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_params
                ]
                for item in sublist
            ]
        ).get("HarmonicBondForce"),
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_params
                ]
                for item in sublist
            ]
        ).get("HarmonicAngleForce"),
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_params
                ]
                for item in sublist
            ]
        ).get("PeriodicTorsionForce"),
        list_to_dict(
            [
                item
                for sublist in [
                    list(elem) for elem in prmtop_energy_decomposition_params
                ]
                for item in sublist
            ]
        ).get("NonbondedForce"),
    ]
    prmtop_energy_decomposition_params_list = [
        "HarmonicBondForce",
        "HarmonicAngleForce",
        "PeriodicTorsionForce",
        "NonbondedForce",
    ]
    df_energy_params = pd.DataFrame(
        list(
            zip(
                prmtop_energy_decomposition_params_list,
                prmtop_energy_decomposition_params_value,
            )
        ),
        columns=["Energy_term", "Energy_params"],
    )
    df_energy_params = df_energy_params.set_index("Energy_term")
    # print(df_energy_params)
    df_compare = pd.concat([df_energy_non_params, df_energy_params], axis=1)
    df_compare["Energy_difference"] = df_compare["Energy_params"].sub(
        df_compare["Energy_non_params"], axis=0
    )
    print(df_compare)


def rename_receptorligand_pdb(input_pdb):

    """
    The function restores the name of the initial PDB file.

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file.

    """
    command = "cp -r " + input_pdb + " " + input_pdb[:-4] + "_no_solvent.pdb"
    os.system(command)
    command = "cp -r " + input_pdb[:-4] + "_before_qmmm.pdb" + " " + input_pdb
    os.system(command)


def run_openmm_sim(input_pdb, forcefield_file, sim_steps, T):

    """
    Runs OpenMM simulation to test the validity of the
    reparameterized topology files with the existing
    PDB file.

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file, which is a copy of the PDB file
        for each anchor.

    forcefield_file: str
        User-defined topology file (prmtop/parm7 file), which is
        a copy of the topology file for each anchor.

    sim_steps: int
        Number of simulation steps for the OpenMM MD simulation.

    T: int
        OpenMM simulation temperature for NVT simulation.

    """
    prmtop = simtk.openmm.app.AmberPrmtopFile(forcefield_file)
    pdb = simtk.openmm.app.PDBFile(input_pdb)
    system = prmtop.createSystem(
        nonbondedCutoff=1 * simtk.unit.nanometer, constraints=simtk.openmm.app.HBonds
    )
    integrator = simtk.openmm.LangevinIntegrator(
        T * simtk.unit.kelvin, 1 / simtk.unit.picosecond, 0.002 * simtk.unit.picoseconds
    )
    simulation = simtk.openmm.app.Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
    simulation.minimizeEnergy(maxIterations=10000)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
    sim_output = input_pdb[:-4] + "_openmm_sim.pdb"
    simulation.reporters.append(
        simtk.openmm.app.PDBReporter(sim_output, int(sim_steps / 10))
    )
    simulation.reporters.append(
        simtk.openmm.app.StateDataReporter(
            stdout,
            int(sim_steps / 10),
            step=True,
            potentialEnergy=True,
            temperature=True,
        )
    )
    simulation.step(sim_steps)


def get_charge_diff_file(forcefield_file, ligand_pdb, ligand_charge_diff_file):

    """
    The function creates a new file where the old charges, new charges,
    and the charge difference for the ligand molecule are saved.

    Parameters
    ----------
    forcefield_file: str
        User-defined topology file (prmtop/parm7 file), which is
        a copy of the topology file for each anchor.

    ligand_pdb: str
        User-defined ligand PDB file.

    ligand_charge_diff_file: str
        User-defined charge difference file.

    """

    if forcefield_file[-6:] == "prmtop":
        non_qm_forcefield_file = forcefield_file[:-7] + "_before_qmmm.prmtop"
    if forcefield_file[-6:] == ".parm7":
        non_qm_forcefield_file = forcefield_file[:-6] + "_before_qmmm.parm7"
    command = "diff " + non_qm_forcefield_file + " " + forcefield_file + " > diff.txt"
    os.system(command)
    f = open("diff.txt", "r")
    lines = f.readlines()
    charges_list = []
    for i in lines:
        if len(extract_charges_str(i)) == 5:
            charges_list.append(extract_charges_str(i))
    list_1 = [
        item
        for sublist in charges_list[0 : int(len(charges_list) / 2)]
        for item in sublist
    ]
    list_2 = [
        item
        for sublist in charges_list[int(len(charges_list) / 2) :]
        for item in sublist
    ]
    substracted = list()
    for item1, item2 in zip(list_1, list_2):
        item = item1 - item2
        substracted.append(item)
    list1 = []
    list2 = []
    diff = []
    for i in range(len(substracted)):
        if substracted[i] != 0:
            list1.append(list_1[i])
            list2.append(list_2[i])
            diff.append(substracted[i])
    list1 = [i / 18.2223 for i in list1]
    list2 = [i / 18.2223 for i in list2]
    diff = [i / 18.2223 for i in diff]
    print("The sum of charges before paramterization: ", sum(list1))
    print("The sum of charges after paramterization: ", sum(list2))
    df1 = pd.read_csv(ligand_pdb, header=None, delimiter=r"\s+")
    df1 = df1[df1.columns[-1]]
    elements = df1.values.tolist()
    df = pd.DataFrame([elements, list1, list2, diff])
    df = df.transpose()
    df.columns = ["Element", "Before QMMM", "After QMMM", "Charge Diff."]
    df.to_csv(ligand_charge_diff_file, header=True, index=None, sep=",", mode="w")


def get_log_files(
    orca_pdb,
    orca_input_file,
    orca_out_file,
    receptor_pdb,
    ligand_pdb,
):

    """
    The function creates a new directory named,"log_files",
    where all the intermediate files and ORCA output files are
    stored.

    Parameters
    ----------
    orca_pdb: str
        User-defined ORCA PDB file.

    orca_input_file: str
        User-defined ORCA input file.

    orca_out_file: str
        User-defined ORCA output file.

    """
    os.system("rm -rf log_files")
    os.system("mkdir log_files")
    # TODO: consolidate
    command = (
        "mv "
        + orca_pdb
        + " "
        + orca_input_file
        + " "
        + orca_out_file
        + " "
        + receptor_pdb
        + " "
        + ligand_pdb
        + " log_files"
    )
    os.system(command)
    # TODO: this command is terrible! Perform all calcs within work directory
    #  so that this isn't necessary
    command = "mv *orca* *ORCAFF* *_before_qmmm* *_before_charge_replacement* *.txt* *_no_solvent* *_openmm_sim* log_files"
    os.system(command)
