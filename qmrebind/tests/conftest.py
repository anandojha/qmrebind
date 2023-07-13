"""
conftest.py

configurations for qmrebind tests
"""

import os
import shutil
import pytest

import qmrebind.postprocessing as postprocessing
from .utils import get_data_filename, copy_to_tmpdir

@pytest.fixture(scope="session")
def qm_charges_persistent(tmpdir_factory):
    """
    Create a qm_charges file that is persistent across the tests.
    """
    tempdir = tmpdir_factory.mktemp("mytest1")
    os.chdir(tempdir)
    orca_out_file = get_data_filename("qmmm_calc.out")
    input_pdb = get_data_filename("hostguest_no_solvent2.pdb")
    qm_charge_file = "qm_charges.txt"
    ligand_indices = list(range(147, 162))
    ret = postprocessing.get_qm_charges(
        orca_out_file = orca_out_file, qm_charge_file=qm_charge_file, 
        input_pdb=input_pdb, ligand_indices=ligand_indices, 
        qm_charge_scheme = "CHELPG")
    return tempdir, qm_charge_file
    
        
@pytest.fixture(scope="session")
def ff_charges_persistent(tmpdir_factory):
    """
    Create ff_charges file that is persistent across the tests.
    """
    tempdir = tmpdir_factory.mktemp("mytest2")
    os.chdir(tempdir)
    forcefield_file = get_data_filename("hostguest_no_solvent2.parm7")
    input_pdb = get_data_filename("hostguest_no_solvent2.pdb")
    ff_charges_file = "ff_charges.txt"
    ret = postprocessing.get_ff_charges(
        forcefield_file=forcefield_file, ff_charges_file=ff_charges_file, 
        input_pdb=input_pdb)
    return tempdir, ff_charges_file

@pytest.fixture(scope="session")
def qm_ff_charges_persistent(tmpdir_factory):
    """
    Create ff_charges file that is persistent across the tests.
    """
    
    tempdir = tmpdir_factory.mktemp("mytest1")
    os.chdir(tempdir)
    orca_out_file = get_data_filename("qmmm_calc.out")
    input_pdb = get_data_filename("hostguest_no_solvent2.pdb")
    qm_charge_file = "qm_charges.txt"
    ligand_indices = list(range(147, 162))
    ret = postprocessing.get_qm_charges(
        orca_out_file = orca_out_file, qm_charge_file=qm_charge_file, 
        input_pdb=input_pdb, ligand_indices=ligand_indices, 
        qm_charge_scheme = "CHELPG")
    
    forcefield_file = get_data_filename("hostguest_no_solvent2.parm7")
    input_pdb = get_data_filename("hostguest_no_solvent2.pdb")
    ff_charges_file = "ff_charges.txt"
    ret = postprocessing.get_ff_charges(
        forcefield_file=forcefield_file, ff_charges_file=ff_charges_file, 
        input_pdb=input_pdb)
    
    ff_charges_qm_fmt_file = "ff_charge_qm_fmt.txt"
    input_pdb = get_data_filename("hostguest_no_solvent2.pdb")
    ligand_indices = list(range(147, 162))
    ret = postprocessing.get_ff_qm_charges(
        qm_charge_file=qm_charge_file, ff_charges_file=ff_charges_file, 
        ff_charges_qm_fmt_file=ff_charges_qm_fmt_file, input_pdb=input_pdb,
        ligand_indices=ligand_indices)
    
    return tempdir, ff_charges_qm_fmt_file