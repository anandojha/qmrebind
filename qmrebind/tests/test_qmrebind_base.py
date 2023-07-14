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
    assert ret == [('a', 'b'), ('c', 'd'), ('e', 'f'), ('g', 'h')]
    return

def test_len_list_to_dict():
    test_list = ["a", "b", "c", "d", "e", "f", "g", "h"]
    ret = base.list_to_dict(test_list)
    assert ret == {'a': 'b', 'c': 'd', 'e': 'f', 'g': 'h'}
    return

# TODO: make better names for these files
def test_get_pdb_atoms():
    pdb_file = get_data_filename("hostguest_solvent.pdb")  
    ret = base.get_number_pdb_atoms(input_pdb=pdb_file)
    assert ret == 5358

# TODO: make better names for these files
def test_get_indices_qm_region():
    pdb_file = get_data_filename("hostguest_solvent2.pdb") 
    ret = base.get_indices_qm_region(input_pdb=pdb_file, ligand_resname = "APN")
    assert len(ret) == 15

# TODO: make better names for these files
def test_get_indices_qm2_region():
    ligand_file = get_data_filename("guest_no_solvent.pdb") 
    receptor_file = get_data_filename("host_no_solvent.pdb") 
    ret = base.get_indices_qm2_region(
        ligand_pdb=ligand_file, receptor_pdb=receptor_file, 
        cut_off_distance=10)
    assert len(ret[0]) == 7
    assert len(ret[1]) == 147

def test_make_string_range():
    indices = [2917, 2918, 2919, 2920, 2921, 2922, 2923, 2924, 2925, 2926, 
               2927, 2928, 2929, 2930, 2931, 2932, 2933]
    result = base.make_string_range(indices)
    assert result == "2917:2933"
    
    indices2 = [605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 
                617, 618, 619, 620, 621, 798, 799, 800, 801, 802, 803, 804, 
                805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 
                817, 818]
    result2 = base.make_string_range(indices2)
    assert result2 == "605:621 798:818"