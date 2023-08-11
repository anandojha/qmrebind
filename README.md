qmrebind
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/anandojha/qmrebind/workflows/CI/badge.svg)](https://github.com/qmrebind/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/anandojha/qmrebind/branch/main/graph/badge.svg)](https://codecov.io/gh/anandojha/qmrebind/branch/main)

Quantum Mechanical - Molecular Mechanical reparameterization at the 
receptor-ligand Binding site (qmrebind) implemented in Simulation Enabled Estimation of 
Kinetic Rates (SEEKR) multiscale software package

Detailed installation and usage instructions can be found at 
https://qmrebind.readthedocs.io/en/latest/index.html.

## Software Requirements :
Make sure to install these packages before running qmrebind:

* ORCA
* SEEKR2 (optional)

One is likely to want to create a unique conda environment for qmrebind
calculations:

```bash
  conda create --name QMMM python=3.8
```

###Installing ORCA
********************** 
Follow the instructions on this page: https://qmrebind.readthedocs.io/en/latest

Alternatively, visit the official site for instructions to install ORCA: 
https://www.orcasoftware.de/tutorials_orca/first_steps/install.html

###Installing SEEKR2 (optional)
**********************
Use the instructions on this page to install SEEKR2, if desired:
https://seekr2.readthedocs.io/en/latest/installation.html

## Installation and Setup Instructions :

* Activate the previously created conda environment:
```bash
conda activate QMMM # activate the conda environment
conda install -c conda-forge ambertools biopandas pandas matplotlib parmed regex openmm
pip install PyPDF2
```
* Clone the *qmrebind* repository :

```bash
git clone https://github.com/seekrcentral/qmrebind.git
```
* Perform the following steps to get this package installed quickly on a local 
Linux machine (Installation in the home directory is recommended) : 

```bash
cd qmrebind
python setup.py install
python setup.py test # optional
```
Detailed documentation can be found at 
https://qmrebind.readthedocs.io/en/latest/installation.html.


Input PDB file Requirements.
**********************

qmrebind accepts the PDB input file with the following requirements:

* PDB file typically should have the box vector information.

* Ligand and the receptor must be assigned a residue name, with the ligand following the receptor. 

## Authors and Contributors
The following people have contributed directly to the coding and validation 
efforts of qmrebind (listed in alphabetical order of first name). 
The author would like to thank everyone who has helped or will help improve 
this project by providing feedback, bug reports, or other comments.

* Anupam Anand Ojha, UC San Diego (Author and Lead Developer)
* Eliseo Marin-Rimoldi, MoLSSI (Project Mentor and Collaborator)
* Lane W. Votapka, UC San Diego (Project Mentor and Developer)
* Rommie E. Amaro, UC San Diego (Principal Investigator)

### Citing qmrebind

If you use qmrebind, please cite the following paper:

* Ojha AA, Votapka L, Amaro R. QMrebind: Incorporating quantum mechanical 
force field reparameterization at the ligand binding site for improved 
drug-target kinetics through milestoning simulations. ChemRxiv. Cambridge: 
Cambridge Open Engage; 2023; This content is a preprint and has not been 
peer-reviewed.

### Copyright

Copyright (c) 2023, Anupam Anand Ojha


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
