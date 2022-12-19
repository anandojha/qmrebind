QMMMReBind_SEEKR
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/anandojha/QMMMReBind_SEEKR/workflows/CI/badge.svg)](https://github.com/anandojha/QMMMReBind_SEEKR/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/anandojha/QMMMReBind_SEEKR/branch/main/graph/badge.svg)](https://codecov.io/gh/anandojha/QMMMReBind_SEEKR/branch/main)

Quantum Mechanical - Molecular Mechanical Re-parameterization of the Receptor-Ligand Binding site implemented in Simulation Enabled Estimation of Kinetic Rates (SEEKR)


## Software Requirements :
Make sure to install these packages before running the QMMMReBind:

* ORCA

## Installation and Setup Instructions :

* Make sure [anaconda3](https://www.anaconda.com/) is installed on the local machine. Go to the  [download](https://www.anaconda.com/products/individual) page of anaconda3 and install the latest version of anaconda3.
* Create a new conda environment with python = 3.8 and install the package with the following commands in the terminal: 
```bash
conda create -n qmmmrebind_seekr python=3.8
conda activate qmmmrebind_seekr # activate the conda environment
conda install git # install git
conda install -c conda-forge ambertools biopandas pandas matplotlib parmed regex openmm # install the conda-forge dependencies for QMMMREBind_SEEKR
pip install PyPDF2 # install pip dependencies
```
* Clone the *QMMMReBind* repository :
```bash
git clone https://github.com/anandojha/QMMMReBind_SEEKR.git
```
* Perform the following steps to get this package installed quickly on a local linux machine (Installation in the home directory is recommended) : 
```bash
cd QMMMReBind_SEEKR # Enter the QMMMReBind_SEEKR directory
python install -e . # Install QMMMReBind_SEEKR
```
A detailed cocumentation can be found [here](https://qmmmrebind-seekr.readthedocs.io/en/latest/index.html).


## Authors and Contributors
The following people have contributed directly to the coding and validation efforts of QMMMReBind_SEEKR (listed in an alphabetical order of first name). 
The author would like to thank everyone who has helped or will help improve this project by providing feedback, bug reports, or other comments.

* Anupam Anand Ojha, UC San Diego (Author and Lead Developer)
* Eliseo Marin-Rimoldi, MoLSSI (Project Mentor and Collaborator)
* Lane Votapka, UC San Diego (Project Mentor and Collaborator)
* Rommie Amaro, UC San Diego (Principal Investigator)

### Copyright

Copyright (c) 2022, Anupam Anand Ojha


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
