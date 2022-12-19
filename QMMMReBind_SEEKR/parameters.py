######################################## To be defined by the user ########################################

orca_dir_pwd="/home/aaojha/orca"                # PWD of the directory where ORCA is installed.

guest_resname="BEN"                             # Three-letter name for the guest residue.

cut_off_distance=2.20                           # Cut-off distance for the QM2 region within the 
                                                # vicinity of the QM region.
   
input_pdb="hostguest.pdb"                       # User-defined PDB file.

forcefield_file="hostguest.prmtop"              # User-defined topology file (prmtop/parm7 file)

nprocs=8                                        # Number of processors to be used for the MPI 
                                                # enabled calculation.

qm_basis_set="STO-3G"                           # Basis set for the QM level of theory. The options
                                                # include pople basis set, pople polarized basis set,
        	                                # pople polarized diffused basis set, correlation
        		                        # consistent basis set, DEF2 basis set, DEF2 diffused
                                                # basis set, auxiliary coulomb fitted basis_set, auxiliary
                                                # coulomb fitted and exchange basis set, and auxiliary
                                                # correlation consistent basis set.

qm_method="BLYP"                                # QM level of theory from a particular family of one of
                                                # the above-mentioned family of QM methods.

qm_charge_scheme="CHELPG"                       # Charge scheme for the calculation of QM charges for
                                                # the QM region. The options include Hirshfeld,
                                                # CHELPG, Mulliken, and the Loewdin charge calculation
        				        # methods.

qm2_method="XTB"                                # QM theory to be implemented for the lower layer,
                                                # i.e., the QM2 layer. The options for the QM2 method
                                                # include the semiempirical method, the tight binding DFT
                                                # method, and the composite method.

qm2_charge_scheme="CHELPG"                      # Charge scheme for the calculation of QM charges for
                                                # the QM2 region. The options include Hirshfeld,
                                                # CHELPG, Mulliken, and the Loewdin charge calculation
        				        # methods.

qm2_charge=1                                    # Charge of the QM2 (low) region. The user does not have
                                                # to provide the charge of the QM region since it is
                                                # calculated by ORCA itself. The user also does not
                                                # have to provide the charge of the MM region since it is
                                                # determined by the ORCA force field file.

qm2_mult=1                                      # The multiplicity of the QM2 (low) region. The user does
                                                # not have to provide the multiplicity of the QM region
                                                # since ORCA itself calculates it. The user also
                                                # does not have to provide the multiplicity of the MM
                                                # region since the ORCA forcefield file determines it.

qm_charge=1                                     # Charge of the QM region.

qm_mult=1                                       # Multiplicity of the QM region.

######################################## To be defined by the user ########################################

######################################## Default values can be used #######################################
guest_pdb="guest.pdb"                           # User-defined guest PDB file.

host_pdb="host.pdb"                             # User-defined host PDB file.

optimization="True"                             # The string is set to "True" if the user chooses to
                                                # optimize the geometry of the QM region. If not,
                                                # the string is set to "False".

frequency_calculation="False"                   # The string is set to "True" if the user chooses
                                                # a frequency calculation for the QM region. If not,
                                                # the string is set to "False".

orca_pdb="input_orca.pdb"                       # User-defined ORCA PDB file.

orca_input_file="orca_qmm.inp"                  # User-defined ORCA input file.

orca_out_file="orca_qmmm.out"                   # User-defined ORCA output file.
 
qm_charge_file="guest_qm_charges.txt"           # User-defined text file where the charges for the
                                                # atoms of the QM region are saved.


ff_charges_file="ff_charges.txt"                # User-defined file where the charges defined in the
                                                # topology file (prmtop/parm7 file) is stored in
                                                # a text file.

ff_charges_qm_fmt_file="ff_charge_qm_fmt.txt"   # User-defined file where the replaced charges in the
                                                # AMBER format is stored.

######################################## Default values can be used #######################################
