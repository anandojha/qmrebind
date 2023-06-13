"""
TODO: add docstring
"""

######################################## Default values can be used #######################################

guest_pdb="guest.pdb"                           # User-defined guest PDB file.

host_pdb="host.pdb"                             # User-defined host PDB file.

optimization="False"                            # The string is set to "True" if the user chooses to
                                                # optimize the geometry of the QM region. If not,
                                                # the string is set to "False".

frequency_calculation="False"                   # The string is set to "True" if the user chooses
                                                # a frequency calculation for the QM region. If not,
                                                # the string is set to "False".

orca_pdb="input_orca.pdb"                       # User-defined ORCA PDB file.

orca_input_file="orca_qmmm.inp"                 # User-defined ORCA input file.

orca_out_file="orca_qmmm.out"                   # User-defined ORCA output file.
 
qm_charge_file="guest_qm_charges.txt"           # User-defined text file where the charges for the
                                                # atoms of the QM region are saved.


ff_charges_file="ff_charges.txt"                # User-defined file where the charges defined in the
                                                # topology file (prmtop/parm7 file) is stored in
                                                # a text file.

ff_charges_qm_fmt_file="ff_charge_qm_fmt.txt"   # User-defined file where the replaced charges in the
                                                # AMBER format is stored.

guest_charge_diff_file="diff_charges.txt"       # User-defined charge difference file.


sim_steps=1000                                  # Number of simulation steps for the OpenMM MD simulation. 
   

T=300                                           # OpenMM simulation temperature for NVT simulation.
