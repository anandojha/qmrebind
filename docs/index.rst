QMrebind: Quantum Mechanical reparametrization at the receptor-ligand binding site
==================================================================================

:Release: |release|
:Date: |today|

Molecular mechanics (MM) forcefields (FF) seek to approximate the potential
energy surface (PES) of molecules in regions of phase space likely to be sampled
in a realistic situation (for example: near equilibrium, in a biological 
context, etc.) While the FFs of common features, such as amino acids, are
quite optimized, the parametrization of small molecules can be difficult and 
prone to errors. The QMrebind software automates the generation of
parameters for small molecules in a context of interest - within the binding 
site, for instance. Therefore, a highly accurate PES is generated for the
region of phase space where the ligand is in the bound state. The resulting
parameters may be then used as input to subsequent calculations, such as
SEEKR2.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   about
   installation
   quickstart
   tutorial
   program_options

Cite QMrebind
=============

If you wish to cite QMrebind, please cite the following paper:

* Ojha AA, Votapka LW, Amaro RE. QMrebind: incorporating quantum mechanical force field reparameterization at the ligand binding site for improved drug-target kinetics through milestoning simulations. Chem Sci. 2023 Oct 24;14(45):13159-13175. doi: 10.1039/d3sc04195f. PMID: 38023523; PMCID: PMC10664576.


Getting Involved
================

Please report **bugs** or **enhancement requests** through the `Issue 
Tracker`_. 

.. _Issue Tracker: https://github.com/seekrcentral/qmrebind/issues


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
