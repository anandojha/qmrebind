"""
Unit tests for run_qmrebind_amber.py.
"""

import os

from .utils import get_data_filename

import qmrebind.run_qmrebind_amber as run_qmrebind_amber

def test_run_qmrebind_amber_hostguest_system_a(tmpdir):
    os.chdir(tmpdir)
    input_pdb = get_data_filename("hostguest_solvent2.pdb") 
    forcefield_file = get_data_filename("hostguest_solvent2.parm7")
    ligand_resname = "APN"
    new_parm7 = run_qmrebind_amber.run_qmrebind_amber(
        input_pdb, forcefield_file, None, ligand_resname)
    assert os.path.exists(new_parm7)