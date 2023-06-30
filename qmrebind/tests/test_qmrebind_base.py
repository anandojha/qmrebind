"""
test_qmrebind_base.py

Unit and regression test for the qmrebind_base module.
"""
# Import package, test suite, and other packages as needed
import warnings
warnings.filterwarnings("ignore")
from .utils import get_data_filename
import qmrebind.qmrebind_base as base

# TODO: improve file handling
def test_group_str():
    test_list = ["a", "b", "c", "d", "e", "f", "g", "h"]
    ret = list(base.group_str(test_list, 2))
    assert len(ret) == 4

def test_len_list_to_dict():
    test_list = ["a", "b", "c", "d", "e", "f", "g", "h"]
    ret = base.list_to_dict(test_list)
    assert len(ret) == 4

def test_get_pdb_atoms():
    pdb_file = get_data_filename("system_solvent_a.pdb")  
    ret = base.get_number_pdb_atoms(input_pdb=pdb_file)
    assert ret == 5358

def test_get_indices_qm_region():
    pdb_file = get_data_filename("system_a.pdb") 
    ret = base.get_indices_qm_region(input_pdb = pdb_file, 
                                                  ligand_resname = "APN")
    assert len(ret) == 15

def test_get_indices_qm2_region():
    ligand_file = get_data_filename("system_b.pdb") 
    receptor_file = get_data_filename("system_c.pdb") 
    ret = base.get_indices_qm2_region(
        ligand_pdb=ligand_file, receptor_pdb=receptor_file, 
        cut_off_distance=10)
    assert len(ret[0]) == 7
    assert len(ret[1]) == 147
