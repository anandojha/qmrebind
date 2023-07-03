"""
test_preparation.py

Test the preparation.py module.
"""

import os

from .utils import get_data_filename
import qmrebind.preparation as preparation

def test_prepare_pdb(tmpdir):
    os.chdir(tmpdir)
    pdb_file = get_data_filename("system_solvent_b.pdb") 
    preparation.prepare_pdb(input_pdb=pdb_file)
    pdb_file_ = "system_solvent_b_before_qmmm.pdb"         
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    with open(pdb_file_, "r") as f_:
        lines_ = f_.readlines()
    #assert len(lines) != len(lines_)
    for i in lines:
        assert "HETATM" not in i
        assert "CONECT" not in i
        assert "TER" not in i

def test_strip_topology(tmpdir):
    os.chdir(tmpdir)
    forcefield_file = get_data_filename("system_topology_solvent_a.parm7")           
    preparation.strip_topology(
        forcefield_file=forcefield_file)
    forcefield_file_ = "system_topology_solvent_a_before_qmmm.parm7"      
    with open(forcefield_file, "r") as f:
        lines = f.readlines()
    with open(forcefield_file_, "r") as f_:
        lines_ = f_.readlines()
    assert len(lines) != len(lines_)
    
def test_get_receptor_pdb(tmpdir):
    os.chdir(tmpdir)
    pdb_file = get_data_filename("system_solvent_c.pdb")
    receptor_pdb = "host_system_a.pdb"
    preparation.get_receptor_pdb(input_pdb=pdb_file, receptor_pdb=receptor_pdb, 
                                   ligand_indices=list(range(147,162)))
    with open(receptor_pdb, "r") as f:
        lines = f.readlines()
    atom_lines = []
    for line in lines:
        if line.startswith("ATOM"):
            atom_lines.append(line)
    assert len(atom_lines) == 147
    for i in lines:
        assert "WAT" not in i
        assert "HOH" not in i
        assert "APN" not in i