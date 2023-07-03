"""
Prepare non-ORCA files for QMRebind runs.
"""

import os
import shutil
import re

import numpy as np
import parmed

import qmrebind.qmrebind_base as base
import qmrebind.defaults as defaults

def prepare_pdb(input_pdb):
    """
    Use the pdb4amber module of AMBER to remove
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
    input_before_qmmm = f"{input_pdb_base}_before_qmmm.pdb"
    shutil.copyfile(input_pdb, input_before_qmmm)
    intermediate_file_I = f"{input_pdb_base}_intermediate_I.pdb"
    intermediate_file_II = f"{input_pdb_base}_intermediate_II.pdb"
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
        f"{intermediate_file_II_base}_nonprot.pdb",
        f"{intermediate_file_II_base}_renum.txt",
        f"{intermediate_file_II_base}_sslink",
        f"{intermediate_file_II_base}_water.pdb",
    )
    base.delete_files(to_delete)
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
    base.delete_files([intermediate_file_I, intermediate_file_II])
    return

def strip_topology(forcefield_file):

    """
    Remove the parameters for the solvent
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
    base.delete_files([stripped_parm_file, cpptraj_input_filename])
    return

def get_ligand_pdb(input_pdb, ligand_pdb, ligand_indices):

    """
    Read the PDB file, extracts the coordinate
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
    struct = parmed.load_file(input_pdb)
    new_struct = struct[np.array(ligand_indices)]
    print("Saving new structure:", ligand_pdb)
    new_struct.save(ligand_pdb, use_hetatoms=False, overwrite=True)
    
    """
    with open(input_pdb) as f1, open(ligand_pdb, "w") as f2:
        for line in f1:
            if ligand_resname in line:
                f2.write(line)
    """

def get_receptor_pdb(input_pdb, receptor_pdb, ligand_indices):
    """
    Read the PDB file, extracts the coordinate
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
    struct = parmed.load_file(input_pdb)
    receptor_indices = []
    for i, atom in enumerate(struct.atoms):
        if (atom.residue.name not in ["WAT, HOH"] + defaults.IONS) \
                and (i not in ligand_indices):
            receptor_indices.append(i)
        
    new_struct = struct[np.array(receptor_indices)]
    print("Saving new structure:", receptor_pdb)
    new_struct.save(receptor_pdb, use_hetatoms=False, overwrite=True)
        
    """
    with open(input_pdb) as f1, open(receptor_pdb, "w") as f2:
        for line in f1:
            if not re.search(f"{ligand_resname}|CRYST|WAT|HOH", line) \
                    and not any(ion in line for ion in defaults.IONS):
                f2.write(line)
    """
    return
