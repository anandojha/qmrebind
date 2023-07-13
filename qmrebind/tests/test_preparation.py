"""
test_preparation.py

Test the preparation.py module.
"""

import os

from .utils import get_data_filename, copy_to_tmpdir
import qmrebind.preparation as preparation

#TODO: run this command within a temporary work dir
def test_prepare_pdb(tmpdir):
    pdb_file = get_data_filename("hostguest_no_solvent.pdb")
    new_file = copy_to_tmpdir(pdb_file, tmpdir)
    os.chdir(tmpdir)
    preparation.prepare_pdb(input_pdb=new_file)
    pdb_file_ = "hostguest_no_solvent_before_qmmm.pdb"
    with open(new_file, "r") as f:
        lines = f.readlines()
    with open(pdb_file_, "r") as f_:
        lines_ = f_.readlines()
    #assert len(lines) != len(lines_)
    for i in lines:
        assert "HETATM" not in i
        assert "CONECT" not in i
        assert "TER" not in i
    return

#TODO: run this command within a temporary work dir
def test_strip_topology(tmpdir):
    forcefield_file = get_data_filename("hostguest_solvent2.parm7")
    new_ff_file = copy_to_tmpdir(forcefield_file, tmpdir)
    os.chdir(tmpdir)
    preparation.strip_topology(forcefield_file=new_ff_file)
    forcefield_file_ = "hostguest_solvent2_before_qmmm.parm7"
    with open(new_ff_file, "r") as f:
        lines = f.readlines()
    with open(forcefield_file_, "r") as f_:
        lines_ = f_.readlines()
    assert len(lines) != len(lines_)
    return
    
def test_get_receptor_pdb(tmpdir):
    pdb_file = get_data_filename("hostguest_solvent.pdb")
    new_file = copy_to_tmpdir(pdb_file, tmpdir)
    os.chdir(tmpdir)
    receptor_pdb = "host.pdb"
    ligand_indices = list(range(147, 162))
    preparation.get_receptor_pdb(input_pdb=new_file, receptor_pdb=receptor_pdb, 
                                   ligand_indices=ligand_indices)
    with open(receptor_pdb, "r") as f:
        lines = f.readlines()
    #assert len(lines) == 149
    for i in lines:
        assert "WAT" not in i
        assert "HOH" not in i
        assert "APN" not in i
    return