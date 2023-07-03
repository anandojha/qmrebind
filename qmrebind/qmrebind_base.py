"""
qmrebind.py

Qmrebind reparametrizes the partial charges on atoms using a QMMM ONIOM 
calculation on the ligand in the bound state of the receptor.
"""

# Standard library imports
from itertools import zip_longest
import itertools
import shutil
import os

# Related third party imports
from biopandas.pdb import PandasPdb
import scipy.spatial as spatial
import pandas as pd
import numpy as np
import parmed
import simtk

# Local application-specific imports
import qmrebind.defaults as defaults

# TODO: put all this into a class in order to preserve useful intermediate
#  variables?

# TODO: set overwrite to False for production code
#def make_work_dir(work_dir=None, overwrite=False):
def make_work_dir(files_to_copy, work_dir=None, overwrite=True,
                  keep_old=False):
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
    
    if not keep_old:
        assert not os.path.exists(work_dir_abs), \
            f"Work directory {work_dir} already exists. Please remove "\
            "directory, choose a different work directory, or enable "\
            "overwrite flag in this function."        
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
    Create a list of all the directories
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
    Collect the data into fixed-length
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
    Convert an input list with mapped characters
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
    Extract the charges from the topology
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

def get_number_pdb_atoms(input_pdb):

    """
    Count the total number of atoms in the system,
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
    Return the atom indices of the QM region, i.e.,
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
    Return the lists of residues and its indices (beginning
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
    ligand_parmed = parmed.load_file(ligand_pdb)
    receptor_parmed = parmed.load_file(receptor_pdb)
    lig_coords = ligand_parmed.coordinates
    rec_coords = receptor_parmed.coordinates
    
    atom_list_within_dist = np.argwhere(
        spatial.distance.cdist(
            lig_coords, rec_coords).min(axis=0) < cut_off_distance)
    atom_list_within_dist = np.unique(atom_list_within_dist)
    atom_list_within_dist = map(int, atom_list_within_dist)
    receptor_parmed_list = []
    for i in atom_list_within_dist:
        receptor_parmed_list.append(receptor_parmed.atoms[i].residue)
    
    receptor_atom_index_set = set()
    receptor_index_set = set()
    for residue in receptor_parmed_list:
        receptor_index_set.add(residue.number)
        for atom in residue.atoms:
            receptor_atom_index_set.add(atom.idx)
    
    receptor_atom_index_list = list(receptor_atom_index_set)
    receptor_atom_index_list.sort()
    receptor_index_list = list(receptor_index_set)
    receptor_index_list.sort()
    return (receptor_index_list, receptor_atom_index_list)

def rename_receptorligand_pdb(input_pdb):
    """
    Restore the name of the initial PDB file.

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file.

    """
    input_pdb_base = os.path.splitext(input_pdb)[0]
    new_no_solvent_pdb_name = f"{input_pdb_base}_no_solvent.pdb"
    old_before_qmmm_pdb_name = f"{input_pdb_base}_before_qmmm.pdb"
    print("Renaming file:", input_pdb, "to:", new_no_solvent_pdb_name)
    os.rename(input_pdb, new_no_solvent_pdb_name)
    print("Renaming file:", old_before_qmmm_pdb_name, "to:", input_pdb)
    os.rename(old_before_qmmm_pdb_name, input_pdb)
    return

def run_check(my_check, skip_checks):
    """
    Run check and respond according to settings.
    """
    if not skip_checks:
        check_fail_str = "One or more fatal checks failed. It is highly "\
        "recommended that you address and correct each of these problems. "\
        "However, you can force Qmrebind to skip these checks by using "\
        "the --skip_checks (-x) argument."
        if not my_check:
            print(check_fail_str)
        assert my_check, "The Qmrebind calculation can not proceed due "\
            "to failed checks. Use argument '-x' to skip checks."
    return

def initialize_indices(indices):
    """
    Break down a string list of integers into a Python list of ints.
    """
    assert indices != ""
    integers = indices.split(",")
    for integer in integers:
        integer_int = int(integer)
    
    return list(map(int, integers))

def make_string_range(indices):
    """
    Make a string suitable for input to ORCA that defines ranges of
    atoms.
    """
    
    ranges = sum(
        (list(t) for t in zip(indices, indices[1:]) \
            if t[0] + 1 != t[1]), []
    )
    iranges = iter(indices[0:1] + ranges + indices[-1:])
    receptor_input_indices = " ".join([str(n) + ":" + str(next(iranges)) \
                                       for n in iranges])
    return receptor_input_indices