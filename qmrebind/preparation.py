"""
Prepare non-ORCA files for QMRebind runs.
"""

import os
import shutil

import numpy as np
import parmed

import qmrebind.qmrebind_base as base
import qmrebind.defaults as defaults

WATER_TO_TIP3P_PARMED_SCRIPT = \
"""
# This ParmEd script will eliminate all extra points on water and set the
# charges and vdW parameters for water to the TIP3P values. This will
# effectively convert TIP4P(ew) and TIP5P(ew) to TIP3P. All actions act on the
# current, 'active' prmtop. You should have a restart file loaded as well so
# that the coordinates of the EP-less system can be written out, too.

# Use this as <source water_to_tip3p.parmed>

# First strip the extra points
strip :WAT@%EP

# Next change the charges to match TIP3P
change charge @%OW -0.834
change charge @%HW  0.417

# Next change the LJ radii
changeLJSingleType @%HW 0.0 0.0
changeLJSingleType @%OW 1.7926 0.098

# Next change all of the ion parameters to the TIP3P parameters
changeLJSingleType @Li+      1.025     0.0279896
changeLJSingleType @Na+      1.369     0.0874393
changeLJSingleType @K+       1.705     0.1936829
changeLJSingleType @Rb+      1.813     0.3278219
changeLJSingleType @Cs+      1.976     0.4065394
changeLJSingleType @F-       2.303     0.0033640
changeLJSingleType @Cl-      2.513     0.0355910
changeLJSingleType @Br-      2.608     0.0586554
changeLJSingleType @I-       2.860     0.0536816
"""


def convert_to_TIP3P(input_parm7, input_pdb):
    """
    Load a filename and convert all waters to TIP3P
    """
    print(f"Converting waters in file: {input_parm7} to TIP3P.")
    pdb = parmed.load_file(input_pdb)
    parm = parmed.amber.AmberParm(input_parm7, xyz=pdb.coordinates)
    
    # Load PDB?
    rm_EP = parmed.tools.strip(parm, ":WAT@%EP")
    rm_EP.execute()
    
    # Modify charges
    ch_q_oxy = parmed.tools.change(parm, "CHARGE", "@%OW", "-0.834")
    ch_q_oxy.execute()
    ch_q_hyd = parmed.tools.change(parm, "CHARGE", "@%HW", "0.417")
    ch_q_hyd.execute()
    
    # Modify Water LJ params
    ch_lj_oxy = parmed.tools.changeLJSingleType(parm, "@%OW", "1.7926", "0.098")
    ch_lj_oxy.execute()
    ch_lj_hyd = parmed.tools.changeLJSingleType(parm, "@%HW", "0.0", "0.0")
    ch_lj_hyd.execute()
    
    # Modify ion LJ params
    ch_lj_Li = parmed.tools.changeLJSingleType(parm, "@Li+", "1.025", "0.0279896")
    ch_lj_Li.execute()
    ch_lj_Na = parmed.tools.changeLJSingleType(parm, "@Na+", "1.369", "0.0874393")
    ch_lj_Na.execute()
    ch_lj_K  = parmed.tools.changeLJSingleType(parm, "@K+",  "1.705", "0.1936829")
    ch_lj_K.execute()
    ch_lj_Rb = parmed.tools.changeLJSingleType(parm, "@Rb+", "1.813", "0.3278219")
    ch_lj_Rb.execute()
    ch_lj_Cs = parmed.tools.changeLJSingleType(parm, "@Cs+", "1.976", "0.4065394")
    ch_lj_Cs.execute()
    ch_lj_F  = parmed.tools.changeLJSingleType(parm, "@F-",  "2.303", "0.0033640")
    ch_lj_F.execute()
    ch_lj_Cl = parmed.tools.changeLJSingleType(parm, "@Cl-", "2.513", "0.0355910")
    ch_lj_Cl.execute()
    ch_lj_Br = parmed.tools.changeLJSingleType(parm, "@Br-", "2.608", "0.0586554")
    ch_lj_Br.execute()
    ch_lj_I  = parmed.tools.changeLJSingleType(parm, "@I-",  "2.860", "0.0536816")
    ch_lj_I.execute()
    
    parm.save(input_parm7, overwrite=True)
    parm.save(input_pdb, overwrite=True)
    
    return

def prepare_pdb(input_pdb, keep_solvent_molecules=False):
    """
    Use the pdb4amber module of AMBER to remove
    any "TER" and "CONECT" keyword in the PDB file. The function
    removes any solvent and ions present in the PDB file and
    replaces any "HETATM" keyword with the "ATOM" keyword to
    maintain homogeneity and backs up the original PDB file
    with an "_before_qmmm" extension.

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file.

    """
    input_pdb_base = os.path.splitext(input_pdb)[0]
    input_before_qmmm = f"{input_pdb_base}_before_qmmm.pdb"
    shutil.copyfile(input_pdb, input_before_qmmm)
    intermediate_file_I = f"{input_pdb_base}_intermediate_I.pdb"
    intermediate_file_II = f"{input_pdb_base}_intermediate_II.pdb"
    intermediate_file_II_base = input_pdb_base + "_intermediate_II"
    
    if keep_solvent_molecules:
        shutil.copyfile(input_pdb, intermediate_file_I)
        command = (
            f"pdb4amber -i {intermediate_file_I} -o {intermediate_file_II}" \
            " --noter --no-conect"
        )
    else:
        with open(input_pdb) as f1, open(intermediate_file_I, "w") as f2:
            for line in f1:
                if not any(ion in line for ion in defaults.IONS):
                    f2.write(line)
        command = (
            f"pdb4amber -i {intermediate_file_I} -o {intermediate_file_II}" \
            " --noter --no-conect --dry"
        )
        
    print("Running command:", command)
    os.system(command)
    to_delete = (
        f"{intermediate_file_II_base}_nonprot.pdb",
        f"{intermediate_file_II_base}_renum.txt",
        f"{intermediate_file_II_base}_sslink",
        f"{intermediate_file_II_base}_water.pdb",
    )
    base.delete_files(to_delete)
    with open(intermediate_file_II, "r") as f1:
        filedata = f1.readlines()
    filedata2 = []
    for line in filedata:
        line2 = line.replace("HETATM", "ATOM  ")
        if line.startswith("CONECT"):
            continue
        if "CL" in line2[12:16]:
            newname = line2[12:16].replace("CL", "Cl")
            line2 = line2[:12] + newname + line2[16:76] + "Cl" + line2[78:]
        
        filedata2.append(line2)
        
    with open(input_pdb, "w") as f2:
        f2.writelines(filedata2)
    base.delete_files([intermediate_file_I, intermediate_file_II])
    return

def strip_topology(forcefield_file, keep_solvent_molecules=False):

    """
    Remove the parameters for the solvent
    and the ions in the topology file. It does this using
    cpptraj

    Parameters
    ----------
    forcefield_file: str
        User-defined topology file (prmtop/parm7 file)

    """
    ff_base = os.path.splitext(forcefield_file)[0]
    ff_ext = os.path.splitext(forcefield_file)[1]
    ff_before_qmmm = ff_base  + f"_before_qmmm{ff_ext}"
    shutil.copyfile(forcefield_file, ff_before_qmmm)
    stripped_parm_file = f"{ff_base}_stripped{ff_ext}"
    cpptraj_input_filename = "strip_topology.cpptraj"
    if not keep_solvent_molecules:
        with open(cpptraj_input_filename, "w") as f:
            f.write(f"parm {forcefield_file}\n")
            f.write("parmstrip :WAT\n")
            for i in defaults.IONS:
                f.write(f"parmstrip :{i}\n")
            f.write(f"parmwrite out {stripped_parm_file}\n")
            f.write("run")
        command = f"cpptraj -i {cpptraj_input_filename}"
        print("Running command:", command)
        os.system(command)
        os.rename(stripped_parm_file, forcefield_file)
        base.delete_files([stripped_parm_file, cpptraj_input_filename])
        
    return

def get_ligand_pdb(input_pdb, ligand_pdb, ligand_indices):

    """
    Read the PDB file, extracts the coordinate
    information for the ligand molecule and saves it into a
    new PDB file.

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file.

    ligand_pdb: str
        User-defined ligand PDB file.

    ligand_resname: str
        Three-letter name for the ligand residue.

    """
    struct = parmed.load_file(input_pdb)
    new_struct = struct[np.array(ligand_indices)]
    print("Saving new structure:", ligand_pdb)
    new_struct.save(ligand_pdb, use_hetatoms=False, overwrite=True)
    return

def get_receptor_pdb(input_pdb, receptor_pdb, ligand_indices, 
                     keep_solvent_molecules=False):
    """
    Read the PDB file, extracts the coordinate
    information for the receptor and saves it into a PDB file.

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file.

    receptor_pdb: str
        User-defined receptor PDB coordinate file.

    ligand_resname: str
        Three-letter name for the ligand residue.

    """
    struct = parmed.load_file(input_pdb)
    receptor_indices = []
    for i, atom in enumerate(struct.atoms):
        if keep_solvent_molecules:
            if i not in ligand_indices:
                receptor_indices.append(i)
        else:
            if (atom.residue.name not in ["WAT", "HOH"] + defaults.IONS) \
                    and (i not in ligand_indices):
                receptor_indices.append(i)
        
    new_struct = struct[np.array(receptor_indices)]
    print("Saving new structure:", receptor_pdb)
    new_struct.save(receptor_pdb, use_hetatoms=False, overwrite=True)
    return
