"""
test_extra.py

Test the functions in extra.py
"""

from .utils import get_data_filename
import qmrebind.extra as extra

def test_get_system_charge():
    forcefield_file = get_data_filename("system_topology_solvent_a.parm7")  
    pdb_file = get_data_filename("system_solvent_a.pdb")  
    ret = extra.get_system_charge(
        forcefield_file=forcefield_file, input_pdb = pdb_file)
    assert round(ret) == 0