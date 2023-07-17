Quickstart
==========

QMrebind needs an existing AMBER .parm7 with parameters defined for the
entire system. QMrebind will replace the partial charges for the relevant
atoms (typically of the ligand) based on the QMMM calculation, but keep all 
other parameters the same.

Make sure you have a .parm7 and .pdb file ready to reparametrize, if you
don't have any handy at this moment, then feel free to download these two:

:download:`hostguest.parm7 <media/hostguest.parm7>`

:download:`hostguest.pdb <media/hostguest.pdb>`

Input PDB file Requirements
---------------------------

QMrebind accepts the PDB input file with the following requirements:

* PDB file typically should have the box vector information.

* Ligand and the receptor must be assigned a residue name with the ligand following the receptor. 

Run QMrebind
------------

Make sure that the correct conda env is activated:

.. code-block:: bash

  conda activate QMMM

Assuming that one downloaded the files at the top of this page, QMrebind
can be run with the following command:

.. code-block:: bash

  python ~/qmrebind/qmrebind/run_qmrebind_amber.py hostguest.pdb hostguest.parm7 -L APN

As one can see, the PDB file and the PARM7 file must be provided, in that order,
then the resname defining the ligand is provided with the '-L' argument.
Alternatively, one may provide the exact atom indices of the ligand:

.. code-block:: bash

  python ~/qmrebind/qmrebind/run_qmrebind_amber.py hostguest.pdb hostguest.parm7 -l "147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161"

Many options exist as inputs to the ``run_qmrebind_amber.py``. One can see 
more of these options by visiting the 
:doc:`Program Options<program_options>` page or by running the 
``run_qmrebind_amber.py`` program with the '-h' argument.