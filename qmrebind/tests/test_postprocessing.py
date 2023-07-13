"""
test_postprocessing.py

Test the postprocessing.py module.
"""

import os

from .utils import get_data_filename, copy_to_tmpdir
import qmrebind.postprocessing as postprocessing

def test_get_qm_charges(qm_charges_persistent):
    tempdir, qm_charge_file = qm_charges_persistent
    os.chdir(tempdir)
    with open(qm_charge_file, "r") as f:
        lines = f.readlines()
    assert len(lines) == 15
    return

def test_get_ff_charges(ff_charges_persistent):
    tempdir, ff_charges_file = ff_charges_persistent
    os.chdir(tempdir)
    with open(ff_charges_file, "r") as f:
        lines = f.readlines()
    assert len(lines) == 162
    return

def test_get_ff_qm_charges(qm_ff_charges_persistent):
    tempdir, ff_charges_qm_fmt_file = qm_ff_charges_persistent
    os.chdir(tempdir)
    with open(ff_charges_qm_fmt_file, "r") as f:
        lines = f.readlines()
    assert len(lines) == 33
    return

def test_get_qmrebind_parm(qm_ff_charges_persistent):
    tempdir, ff_charges_qm_fmt_file = qm_ff_charges_persistent
    forcefield_file_ = get_data_filename("hostguest_no_solvent2.parm7")
    forcefield_file = get_data_filename("guest_no_solvent.parm7")
    input_pdb = get_data_filename("hostguest_no_solvent2.pdb")
    forcefield_file_new_ = copy_to_tmpdir(forcefield_file_, tempdir)
    forcefield_file_new = copy_to_tmpdir(forcefield_file, tempdir)
    input_pdb_new = copy_to_tmpdir(input_pdb, tempdir)
    os.chdir(tempdir)
    ret = postprocessing.get_qmrebind_parm(
        forcefield_file=forcefield_file_new, input_pdb=input_pdb_new, 
        ff_charges_qm_fmt_file=ff_charges_qm_fmt_file)
    return