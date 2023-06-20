"""
defaults.py

Defaults and constants to use for Qmrebind calculations.
"""

IONS = ["Na+", "Cs+", "K+", "Li+", "Rb+", "Cl-", "Br-", "F-", "I-"]

######## Default values can be used #########

# User-defined ligand PDB file.
ligand_pdb="ligand.pdb"

# User-defined receptor PDB file.
receptor_pdb="receptor.pdb"

# The string is set to "True" if the user chooses to
# optimize the geometry of the QM region. If not,
# the string is set to "False".
optimization="False"

# The string is set to "True" if the user chooses
# a frequency calculation for the QM region. If not,
# the string is set to "False".
frequency_calculation="False"

# User-defined ORCA PDB file.
orca_pdb="input_orca.pdb"

# User-defined ORCA input file.
orca_input_file="orca_qmmm.inp"

# User-defined ORCA output file.
orca_out_file="orca_qmmm.out"
 
# User-defined text file where the charges for the
# atoms of the QM region are saved.
qm_charge_file="ligand_qm_charges.txt"

# User-defined file where the charges defined in the
# topology file (prmtop/parm7 file) is stored in
# a text file.
ff_charges_file="ff_charges.txt"

# User-defined file where the replaced charges in the
# AMBER format is stored.
ff_charges_qm_fmt_file="ff_charge_qm_fmt.txt"

# User-defined charge difference file.
ligand_charge_diff_file="diff_charges.txt"

# Number of simulation steps for the OpenMM MD simulation. 
sim_steps=1000
   
# OpenMM simulation temperature for NVT simulation.
T=300
