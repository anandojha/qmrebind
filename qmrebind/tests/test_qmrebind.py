"""
Unit and regression test for the qmrebind package.
"""
# Import package, test suite, and other packages as needed
import warnings
warnings.filterwarnings("ignore")
from .utils import get_data_filename
import qmrebind
import numpy as np
import pytest
import sys
import os
import re

def test_qmrebind_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "qmrebind" in sys.modules

def test_group_str():
    test_list = ["a", "b", "c", "d", "e", "f", "g", "h"]
    ret = list(qmrebind.qmrebind.group_str(test_list, 2))
    assert len(ret) == 4

def test_len_list_to_dict():
    test_list = ["a", "b", "c", "d", "e", "f", "g", "h"]
    ret = qmrebind.qmrebind.list_to_dict(test_list)
    assert len(ret) == 4

def test_strip_topology():
    forcefield_file = get_data_filename("system_topology_solvent_a.parm7")           
    qmrebind.qmrebind.strip_topology(
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

def test_get_system_charge():
    forcefield_file = get_data_filename("system_topology_solvent_a.parm7")  
    pdb_file = get_data_filename("system_solvent_a.pdb")  
    ret = qmrebind.qmrebind.get_system_charge(
        forcefield_file=forcefield_file, input_pdb = pdb_file)
    assert round(ret) == 0

def test_get_pdb_atoms():
    pdb_file = get_data_filename("system_solvent_a.pdb")  
    ret = qmrebind.qmrebind.get_pdb_atoms(input_pdb = pdb_file)
    assert ret == 5358

def test_prepare_pdb():
    pdb_file = get_data_filename("system_solvent_b.pdb") 
    ret = qmrebind.qmrebind.prepare_pdb(input_pdb = pdb_file)
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


def test_get_host_pdb():
    pdb_file = get_data_filename("system_solvent_c.pdb") 
    host_pdb = get_data_filename("host_system_a.pdb") 
    qmrebind.qmrebind.get_host_pdb(input_pdb = pdb_file, host_pdb = host_pdb, 
                                   guest_resname = "APN")
    with open(host_pdb, "r") as f:
        lines = f.readlines()
    assert len(lines) == 149
    for i in lines:
        assert "WAT" not in i
        assert "HOH" not in i
        assert "APN" not in i


def test_get_indices_qm_region():
    pdb_file = get_data_filename("system_a.pdb") 
    ret = qmrebind.qmrebind.get_indices_qm_region(input_pdb = pdb_file, 
                                                  guest_resname = "APN")
    assert len(ret) == 15

def test_get_indices_qm2_region():
    guest_file = get_data_filename("system_b.pdb") 
    host_file = get_data_filename("system_c.pdb") 
    ret = qmrebind.qmrebind.get_indices_qm2_region(
        guest_pdb = guest_file, host_pdb = host_file, cut_off_distance = 10)
    assert len(ret[0]) == 7
    assert len(ret[1]) == 147

def test_get_qm_charges():
    orca_out_file = get_data_filename("qmmm_calc.out") 
    input_pdb = get_data_filename("system_a.pdb") 
    qm_charge_file = get_data_filename("qm_charges.txt") 
    ret = qmrebind.qmrebind.get_qm_charges(
        orca_out_file = orca_out_file, qm_charge_file = qm_charge_file, 
        input_pdb = input_pdb, guest_resname = "APN", 
        qm_charge_scheme = "CHELPG")
    with open(qm_charge_file, "r") as f:
        lines = f.readlines()
    assert len(lines) == 15

def test_get_ff_charges():
    forcefield_file = get_data_filename("system_a.parm7") 
    input_pdb = get_data_filename("system_a.pdb") 
    ff_charges_file = get_data_filename("ff_charges.txt") 
    ret = qmrebind.qmrebind.get_ff_charges(
        forcefield_file = forcefield_file, ff_charges_file = ff_charges_file, 
        input_pdb = input_pdb)
    with open(ff_charges_file, "r") as f:
        lines = f.readlines()
    assert len(lines) == 162

def test_get_ff_qm_charges():
    qm_charge_file = get_data_filename("qm_charges.txt") 
    ff_charges_file = get_data_filename("ff_charges.txt")
    ff_charges_qm_fmt_file = get_data_filename("ff_charge_qm_fmt.txt")
    input_pdb = get_data_filename("system_a.pdb")  
    ret = qmrebind.qmrebind.get_ff_qm_charges(
        qm_charge_file = qm_charge_file, ff_charges_file = ff_charges_file, 
        ff_charges_qm_fmt_file = ff_charges_qm_fmt_file, input_pdb = input_pdb,
        guest_resname = "APN")
    with open(ff_charges_qm_fmt_file, "r") as f:
        lines = f.readlines()
    assert len(lines) == 33

def test_get_qmrebind_parm():
    forcefield_file_ = get_data_filename("system_a.parm7") 
    forcefield_file = get_data_filename("system_b.parm7") 
    input_pdb = get_data_filename("system_a.pdb")  
    ff_charges_qm_fmt_file = get_data_filename("ff_charge_qm_fmt.txt")
    ret = qmrebind.qmrebind.get_qmmmrebind_parm(
        forcefield_file = forcefield_file, input_pdb = input_pdb, 
        ff_charges_qm_fmt_file = ff_charges_qm_fmt_file)
            










