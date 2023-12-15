"""
Given two parm7 files, this script will transfer the partial charges from
the ligand (or other set of atoms) of one file to the other.
"""

import os
import argparse
import tempfile

import parmed

import qmrebind.qmrebind_base as base
import qmrebind.postprocessing as postprocessing
""" # TODO :marked for removal
def write_ligand_charges(prmtop_file, ligand_indices):
    struct = parmed.load_file(prmtop_file)
    charges = []
    for i, atom in enumerate(struct.atoms):
        if i in ligand_indices:
            charges.append(atom.charge)
            
    qm_charges_1_filename = "/tmp/qm_1.txt"#tempfile.NamedTemporaryFile().name
    with open(qm_charges_1_filename, "w") as f:
        for charge in charges:
            f.write(f"{charge:.6f}\n")
            
    return qm_charges_1_filename
    
def get_qmrebind_parm(old_parm7_file, new_parm7_file, input_pdb, 
                      ff_charges_qm_fmt_file):
    with open(old_parm7_file, "r") as f1:
        ff_lines = f1.readlines()
    for i, line in enumerate(ff_lines):
        if "%FLAG CHARGE" in line:
            to_begin = 0
            to_end = int(i)
    parm_lines_a = ff_lines[0 : to_end + 2]
    num_atoms = base.get_number_pdb_atoms(input_pdb=input_pdb)
    if num_atoms % 5 == 0:
        ff_lines_to_select = int(num_atoms / 5)
    else:
        ff_lines_to_select = int(num_atoms // 5 + 1)
    parm_lines_b = ff_lines[to_end + ff_lines_to_select + 2 :]
    with open(ff_charges_qm_fmt_file, "r") as f2:
        qm_ff_lines = f2.readlines()
    with open(new_parm7_file, "w") as f:
        for line in parm_lines_a:
            f.write(line)
        for line in qm_ff_lines:
            f.write(line)
        for line in parm_lines_b:
            f.write(line)
            
    return

def transfer_charges(input_pdb1, prmtop_file1, input_pdb2, prmtop_file2, 
                     new_prmtop_file, ligand_resname):
    
    qm_region_atom_indices_prmtop1 = base.get_indices_qm_region(
            input_pdb=input_pdb1, ligand_resname=ligand_resname)
    qm_region_atom_indices_prmtop2 = base.get_indices_qm_region(
            input_pdb=input_pdb2, ligand_resname=ligand_resname)
    
    print("qm_region_atom_indices_prmtop1:", qm_region_atom_indices_prmtop1)
    print("qm_region_atom_indices_prmtop2:", qm_region_atom_indices_prmtop2)
    
    qm_charges_1_filename = write_ligand_charges(
        prmtop_file1, qm_region_atom_indices_prmtop1)
    
    ff_charges_2_filename = "/tmp/ff_2.txt" #tempfile.NamedTemporaryFile().name
    postprocessing.get_ff_charges(
        forcefield_file=prmtop_file2,
        ff_charges_file=ff_charges_2_filename,
        input_pdb=input_pdb2,
    )
    ff_charges_qm_2_filename = "/tmp/ff_qm_2.txt" #tempfile.NamedTemporaryFile().name
    postprocessing.get_ff_qm_charges(
        qm_charge_file=qm_charges_1_filename,
        ff_charges_file=ff_charges_2_filename,
        ff_charges_qm_fmt_file=ff_charges_qm_2_filename,
        input_pdb=input_pdb2,
        ligand_indices=qm_region_atom_indices_prmtop2,
    )
    get_qmrebind_parm(prmtop_file2, new_prmtop_file, input_pdb2, 
                      ff_charges_qm_2_filename)
    
    return
"""

def transfer_charges(prmtop_file1, prmtop_file2, new_prmtop_file, 
                     ligand_resname):
    atom_name_charge_dict = {}
    prmtop1 = parmed.load_file(prmtop_file1)
    for residue in prmtop1.residues:
        for atom in residue.atoms:
            if residue.name == ligand_resname:
                atom_name_charge_dict[atom.name] = atom.charge
                
    prmtop2 = parmed.load_file(prmtop_file2)
    for residue in prmtop2.residues:
        for atom in residue.atoms:
            if residue.name == ligand_resname:
                assert atom.name in atom_name_charge_dict, \
                    f"Atom name not found: {atom.name}"
                atom.charge = atom_name_charge_dict[atom.name]
                atom_name_charge_dict.pop(atom.name)
                
    assert len(atom_name_charge_dict) == 0, "Mismatch between names or number "\
        "of atoms between prmtop files."

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "prmtop_file1", metavar="PRMTOP_FILE1", type=str, 
        help="The name of the .prmtop or .parm7 forcefield file whose atomic "\
        "point charges will be transferred FROM")
    argparser.add_argument(
        "prmtop_file2", metavar="PRMTOP_FILE2", type=str, 
        help="The name of the .prmtop or .parm7 forcefield file whose atomic "\
        "point charges will be transferred TO")
    argparser.add_argument(
        "new_prmtop_file", metavar="NEW_PRMTOP_FILE", type=str, 
        help="The name of the .prmtop or .parm7 forcefield file to write the "\
        "new system to.")
    argparser.add_argument(
        "-L", "--ligand_resname", dest="ligand_resname", 
        metavar="LIGAND_RESNAME", type=str, default="",
        help="The residue name of the ligand molecule for automatic index "\
        "selection. At this time, '-L' argument must be included.")
    
    args = argparser.parse_args()
    args = vars(args)
    prmtop_file1 = args["prmtop_file1"]
    prmtop_file2 = args["prmtop_file2"]
    new_prmtop_file = args["new_prmtop_file"]
    ligand_resname = args["ligand_resname"]
    
    transfer_charges(prmtop_file1, prmtop_file2, 
                     new_prmtop_file, ligand_resname)