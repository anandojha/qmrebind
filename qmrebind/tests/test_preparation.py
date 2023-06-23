"""
test_preparation.py

Test the preparation.py module.
"""

import os

from .utils import get_data_filename
import qmrebind.preparation as preparation

def test_prepare_pdb():
    pdb_file = get_data_filename("system_solvent_b.pdb") 
    preparation.prepare_pdb(input_pdb=pdb_file)
    pdb_file_ = get_data_filename("system_solvent_b_before_qmmm.pdb")            
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    with open(pdb_file_, "r") as f_:
        lines_ = f_.readlines()
    #assert len(lines) != len(lines_)
    for i in lines:
        assert "HETATM" not in i
        assert "CONECT" not in i
        assert "TER" not in i
    command = "rm -rf " + pdb_file
    os.system(command)
    command = "mv " +  pdb_file_ + " " + pdb_file
    os.system(command)

def test_strip_topology():
    forcefield_file = get_data_filename("system_topology_solvent_a.parm7")           
    preparation.strip_topology(
        forcefield_file=forcefield_file)
    forcefield_file_ = get_data_filename(
        "system_topology_solvent_a_before_qmmm.parm7")           
    with open(forcefield_file, "r") as f:
        lines = f.readlines()
    with open(forcefield_file_, "r") as f_:
        lines_ = f_.readlines()
    assert len(lines) != len(lines_)
    command = "rm -rf " + forcefield_file
    os.system(command)
    command = "mv " +  forcefield_file_ + " " + forcefield_file
    os.system(command)
    
def test_get_receptor_pdb():
    pdb_file = get_data_filename("system_solvent_c.pdb") 
    receptor_pdb = get_data_filename("host_system_a.pdb") 
    preparation.get_receptor_pdb(input_pdb=pdb_file, receptor_pdb=receptor_pdb, 
                                   ligand_resname="APN")
    with open(receptor_pdb, "r") as f:
        lines = f.readlines()
    assert len(lines) == 149
    for i in lines:
        assert "WAT" not in i
        assert "HOH" not in i
        assert "APN" not in i