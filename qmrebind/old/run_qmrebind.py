import sys
import os

pwd_qmrebind = (
    "/home/aaojha/qmrebind/"  # pwd of qmrebind package
)
path_join = pwd_qmrebind + "qmrebind/"
module_path = os.path.abspath(os.path.join(path_join))
if module_path not in sys.path:
    sys.path.append(module_path)

# TODO: remove these global imports
from qmrebind import *
from parameters import *

# Getting started with the ORCA simulation using the modified intial PDB file

prepare_pdb(input_pdb=input_pdb)

strip_topology(forcefield_file=forcefield_file)

get_system_charge(forcefield_file=forcefield_file, input_pdb=input_pdb)

get_guest_pdb(input_pdb=input_pdb, guest_pdb=guest_pdb, 
              guest_resname=guest_resname)

get_host_pdb(input_pdb=input_pdb, host_pdb=host_pdb, 
             guest_resname=guest_resname)

print(
    "The indices for the atoms in the QM region are: "
    + str(get_indices_qm_region(input_pdb=input_pdb, 
                                guest_resname=guest_resname))
    + ", and the number of atoms are: "
    + str(len(get_indices_qm_region(input_pdb=input_pdb, 
                                    guest_resname=guest_resname)))
)

print(
    "The indices for atoms in the QM2 region are: "
    + str(
        get_indices_qm2_region(
            guest_pdb=guest_pdb, host_pdb=host_pdb, 
            cut_off_distance=cut_off_distance
        )[1]
    )
    + ", and the number of atoms are: "
    + str(
        len(
            get_indices_qm2_region(
                guest_pdb=guest_pdb,
                host_pdb=host_pdb,
                cut_off_distance=cut_off_distance,
            )[1]
        )
    )
)

print(
    "The indices for residues in the QM2 region are: "
    + str(
        get_indices_qm2_region(
            guest_pdb=guest_pdb, host_pdb=host_pdb, 
            cut_off_distance=cut_off_distance
        )[0]
    )
    + ", and the number of residues are: "
    + str(
        len(
            get_indices_qm2_region(
                guest_pdb=guest_pdb,
                host_pdb=host_pdb,
                cut_off_distance=cut_off_distance,
            )[0]
        )
    )
)

prepare_orca_pdb(
    input_pdb=input_pdb,
    guest_pdb=guest_pdb,
    orca_pdb=orca_pdb,
    guest_resname=guest_resname,
    host_pdb=host_pdb,
    cut_off_distance=cut_off_distance,
)

get_amber_to_orca_prms(forcefield_file=forcefield_file)

# ORCA simulation

get_orca_input(
    nprocs=nprocs,
    maxiter=maxiter,
    qm_method=qm_method,
    qm_basis_set=qm_basis_set,
    qm2_method=qm2_method,
    optimization=optimization,
    frequency_calculation=frequency_calculation,
    qm_charge_scheme=qm_charge_scheme,
    qm2_charge_scheme=qm2_charge_scheme,
    qm2_charge=qm2_charge,
    qm2_mult=qm2_mult,
    forcefield_file=forcefield_file,
    input_pdb=input_pdb,
    guest_resname=guest_resname,
    orca_pdb=orca_pdb,
    orca_input_file=orca_input_file,
    guest_pdb=guest_pdb,
    host_pdb=host_pdb,
    cut_off_distance=cut_off_distance,
    qm_charge=qm_charge,
    qm_mult=qm_mult,
)

"""
add_xtb_inputs(
    etemp=etemp,
    solvation=solvation,
    solvent=solvent,
    accuracy=accuracy,
    xtb_memory=xtb_memory,
    xtb_nprocs=xtb_nprocs,
    orca_input_file=orca_input_file,
    XTB_add_inputs=XTB_add_inputs,
)
"""
run_orca_qmmm(
    orca_dir_pwd=orca_dir_pwd,
    orca_input_file=orca_input_file,
    orca_out_file=orca_out_file,
)

# Post simulation

get_qm_charges(
    orca_out_file=orca_out_file,
    qm_charge_file=qm_charge_file,
    input_pdb=input_pdb,
    guest_resname=guest_resname,
    qm_charge_scheme=qm_charge_scheme,
)

get_ff_charges(
    forcefield_file=forcefield_file,
    ff_charges_file=ff_charges_file,
    input_pdb=input_pdb,
)

get_ff_qm_charges(
    qm_charge_file=qm_charge_file,
    ff_charges_file=ff_charges_file,
    ff_charges_qm_fmt_file=ff_charges_qm_fmt_file,
    input_pdb=input_pdb,
    guest_resname=guest_resname,
)

get_qmmmrebind_parm(
    forcefield_file=forcefield_file,
    input_pdb=input_pdb,
    ff_charges_qm_fmt_file=ff_charges_qm_fmt_file,
)

# Post Analysis

get_qmmmrebind_parm_solvent(
    input_pdb=input_pdb,
    forcefield_file=forcefield_file,
    ff_charges_file=ff_charges_file,
)

get_energy_diff_no_solvent(forcefield_file=forcefield_file, input_pdb=input_pdb)

get_energy_diff_solvent(forcefield_file=forcefield_file, input_pdb=input_pdb)

rename_hostguest_pdb(input_pdb=input_pdb)

run_openmm_sim(
    input_pdb=input_pdb, forcefield_file=forcefield_file, sim_steps=sim_steps, 
    T=T
)

get_charge_diff_file(
    forcefield_file=forcefield_file,
    guest_pdb=guest_pdb,
    guest_charge_diff_file=guest_charge_diff_file,
)

get_log_files(
    orca_pdb=orca_pdb,
    orca_input_file=orca_input_file,
    orca_out_file=orca_out_file,
    host_pdb=host_pdb,
    guest_pdb=guest_pdb,
)
