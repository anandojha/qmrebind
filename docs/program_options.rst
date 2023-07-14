Program Options
===============

The arguments and inputs of the various QMrebind programs are described here.

run_qmrebind_amber.py
---------------------

Run ``run_qmrebind_amber.py`` with the '-h' argument to see all available
arguments:

.. code-block:: bash

  python run_qmrebind_amber.py -h
  
usage::

  python run_qmrebind_amber.py [-h] [-l LIGAND_INDICES] [-L LIGAND_RESNAME]
                             [-o OUTPUT] [-c CUT_OFF_DISTANCE] [-n NPROCS]
                             [-i MAX_ITERATIONS] [-m QM_METHOD]
                             [-b QM_BASIS_SET] [-s QM_CHARGE_SCHEME]
                             [-q QM_CHARGE] [-u QM_MULTIPLICITY]
                             [-M QM2_METHOD] [-S QM2_CHARGE_SCHEME]
                             [-Q QM2_CHARGE] [-U QM2_MULTIPLICITY]
                             [-O ORCA_PATH] [-w WORK_DIR] [-x]
                             INPUT_PDB FORCEFIELD_FILE

  
Required Arguments
******************

**INPUT_PDB**
  The name of a PDB file of a molecular system. The PDB should contain box
  information in a CRYST line, and have receptor atoms, ligand atoms, and
  solvent atoms included as contiguous lines, in that order.
  
**FORCEFIELD_FILE**
  The name of a AMBER input .parm7 (or .prmtop) file, which contains
  existing system FF parameters, some of which will by overwritten (typically
  the partial charges of the ligand).

Optional Arguments
******************

  -h, --help
                        show help message and exit

  -l LIGAND_INDICES, --ligand_indices LIGAND_INDICES
                        A comma-separated list of integers defining site
                        within the ref_pdb structure. Ex: -l '1,2,0'. 
                        Either the '-l' or '-L' arguments must be included.
  -L LIGAND_RESNAME, --ligand_resname LIGAND_RESNAME
                        The residue name of the ligand molecule for automatic
                        index selection. Either the '-l' or '-L' arguments 
                        must be included.
  -o OUTPUT, --output OUTPUT
                        A path to an output file name for the forcefield file.
                        If left at the default of None, a file will be
                        generated in the work directory with the same name as
                        the input forcefield file. Default: None.
  -c CUT_OFF_DISTANCE, --cut_off_distance CUT_OFF_DISTANCE
                        The cut-off distance (in Angstroms) used to define the
                        QM2 region of the ONIOM calculation. Default: 3.0.
  -n NPROCS, --nprocs NPROCS
                        The number of processors to use for ORCA calculations. 
                        Default: 1.
  -i MAX_ITERATIONS, --max_iterations MAX_ITERATIONS
                        The maximum number of iterations for ORCA convergence.
                        Default: 2000.
  -m QM_METHOD, --qm_method QM_METHOD
                        The method to use for the QM region of the ONIOM
                        calculation. Please see the file
                        orca_methods_basis_sets.pdf for all possible options.
                        Default: B3LYP.
  -b QM_BASIS_SET, --qm_basis_set QM_BASIS_SET
                        The basis set to use for the QM region of the ONIOM
                        calculation. Please see the file
                        orca_methods_basis_sets.pdf for all possible options.
                        Default: 6-311G.
  -s QM_CHARGE_SCHEME, --qm_charge_scheme QM_CHARGE_SCHEME
                        The charge scheme to use for the QM region of the
                        ONIOM calculation. Please see the file
                        orca_methods_basis_sets.pdf for all possible options.
                        Default: CHELPG.
  -q QM_CHARGE, --qm_charge QM_CHARGE
                        The total charge of the QM region of the ONIOM
                        calculation. Default: 0.
  -u QM_MULTIPLICITY, --qm_multiplicity QM_MULTIPLICITY
                        The multiplicity of the QM region of the ONIOM
                        calculation. Default: 1.
  -M QM2_METHOD, --qm2_method QM2_METHOD
                        The method to use for the QM2 region of the ONIOM
                        calculation. Please see the file
                        orca_methods_basis_sets.pdf for all possible options.
                        Default: XTB.
  -S QM2_CHARGE_SCHEME, --qm2_charge_scheme QM2_CHARGE_SCHEME
                        The charge scheme to use for the QM2 region of the
                        ONIOM calculation. Please see the file
                        orca_methods_basis_sets.pdf for all possible options.
                        Default: CHELPG.
  -Q QM2_CHARGE, --qm2_charge QM2_CHARGE
                        The total charge of the QM2 region of the ONIOM
                        calculation. Default: 0.
  -U QM2_MULTIPLICITY, --qm2_multiplicity QM2_MULTIPLICITY
                        The multiplicity of the QM2 region of the ONIOM
                        calculation. Default: 1.
  -O ORCA_PATH, --orca_path ORCA_PATH
                        An absolute path to the ORCA program. If not
                        specified, ORCA will be found from the shutil.which()
                        command. Default: None.
  -w WORK_DIR, --work_dir WORK_DIR
                        A working directory for all temporary files, log
                        files, etc. If left at None, a work directory will be 
                        automatically generated in the current directory. 
                        Default: None.
  -x, --skip_checks     By default, checks will be run at various stages of
                        Qmrebind, and if the checks fail, the calculation will
                        not proceed. This argument bypasses those checks and
                        allows the calculation to proceed anyways. 
                        Default: False.