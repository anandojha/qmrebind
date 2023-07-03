"""
test_postprocessing.py

Test the postprocessing.py module.
"""

import os

from .utils import get_data_filename
import qmrebind.postprocessing as postprocessing

def test_get_qm_charges(tmpdir):
    os.chdir(tmpdir)
    orca_out_file = get_data_filename("qmmm_calc.out") 
    input_pdb = get_data_filename("system_a.pdb") 
    qm_charge_file = "qm_charges.txt"
    ret = postprocessing.get_qm_charges(
        orca_out_file = orca_out_file, qm_charge_file=qm_charge_file, 
        input_pdb=input_pdb, ligand_indices=list(range(147,162)), 
        qm_charge_scheme = "CHELPG")
    with open(qm_charge_file, "r") as f:
        lines = f.readlines()
    assert len(lines) == 15

def test_get_ff_charges():
    forcefield_file = get_data_filename("system_a.parm7") 
    input_pdb = get_data_filename("system_a.pdb") 
    ff_charges_file = get_data_filename("ff_charges.txt") 
    ret = postprocessing.get_ff_charges(
        forcefield_file=forcefield_file, ff_charges_file=ff_charges_file, 
        input_pdb=input_pdb)
    with open(ff_charges_file, "r") as f:
        lines = f.readlines()
    assert len(lines) == 162

def test_get_ff_qm_charges():
    qm_charge_file = get_data_filename("qm_charges.txt") 
    ff_charges_file = get_data_filename("ff_charges.txt")
    ff_charges_qm_fmt_file = get_data_filename("ff_charge_qm_fmt.txt")
    input_pdb = get_data_filename("system_a.pdb")  
    ret = postprocessing.get_ff_qm_charges(
        qm_charge_file = qm_charge_file, ff_charges_file=ff_charges_file, 
        ff_charges_qm_fmt_file = ff_charges_qm_fmt_file, input_pdb=input_pdb,
        ligand_indices=list(range(147,162)))
    with open(ff_charges_qm_fmt_file, "r") as f:
        lines = f.readlines()
    assert len(lines) == 33

def test_get_qmrebind_parm():
    forcefield_file_ = get_data_filename("system_a.parm7") 
    forcefield_file = get_data_filename("system_b.parm7") 
    input_pdb = get_data_filename("system_a.pdb")  
    ff_charges_qm_fmt_file = get_data_filename("ff_charge_qm_fmt.txt")
    ret = postprocessing.get_qmrebind_parm(
        forcefield_file = forcefield_file, input_pdb = input_pdb, 
        ff_charges_qm_fmt_file = ff_charges_qm_fmt_file)