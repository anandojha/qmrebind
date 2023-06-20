"""
Extra code that may be useful in the future.
"""

import os

import qmrebind.qmrebind_base as qmrebind

def get_initial_files(input_pdb, forcefield_file):

    """
    The function iterates through the building folder of
    each of the anchor directories and creates a copy
    of the topology file, i.e., prmtop/parm7 file
    and the PDB file for every anchor.

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file, which is a copy of the PDB file
        for each anchor.

    forcefield_file: str
        User-defined topology file (prmtop/parm7 file), which is
        a copy of the topology file for each anchor.

    """
    current_pwd = os.getcwd()
    anchor_list = []
    for anchor in os.listdir(current_pwd):
        if anchor.startswith("anchor_"):
            anchor_list.append(anchor)
    for i in anchor_list:
        pwd = os.chdir(current_pwd + "/" + i + "/" + "building")
        for f in os.listdir(pwd):
            if f.endswith("parm7"):
                command = "cp -r " + f + " " + forcefield_file
                os.system(command)
            if f.endswith("prmtop"):
                command = "cp -r " + f + " " + forcefield_file
                os.system(command)
            if f.endswith("pdb"):
                command = "cp -r " + f + " " + input_pdb
                os.system(command)
    os.chdir(current_pwd)

def get_system_charge(forcefield_file, input_pdb):

    """
    The function returns the charge of the input PDB
    file.

    Parameters
    ----------
    forcefield_file: str
        User-defined topology file (prmtop/parm7 file)

    input_pdb: str
        User-defined PDB file.

    """
    with open(forcefield_file, "r") as f:
        lines = f.readlines()
    no_atoms = qmrebind.get_pdb_atoms(input_pdb=input_pdb)
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
    charge_list_value = []
    for i in charge_list:
        charge_list_value.append(float(i))
    system_charge = sum(charge_list_value) / 18.2223
    print("The total charge of the system is: " + str(system_charge))
    system_charge = int(round(system_charge))
    return system_charge
