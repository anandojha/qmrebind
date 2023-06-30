qmrebind
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/anandojha/qmrebind/workflows/CI/badge.svg)](https://github.com/qmrebind/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/anandojha/qmrebind/branch/main/graph/badge.svg)](https://codecov.io/gh/anandojha/qmrebind/branch/main)

Quantum Mechanical - Molecular Mechanical reparameterization at the 
receptor-ligand Binding site (qmrebind) implemented in Simulation Enabled Estimation of 
Kinetic Rates (SEEKR)


## Software Requirements :
Make sure to install these packages before running qmrebind:

* ORCA
* SEEKR2 (optional)

Section IA: Installing ORCA
********************** 
1. Go to https://orcaforum.kofo.mpg.de/ucp.php?mode=login. Create an account 
to log in with a username and a password. 

2. Go to Downloads.

3. Select ORCA 5.0.3

4. The ORCA tar files will be downloaded in three parts: Download ORCA 5.0.3, Linux, x86-64, .tar.xz Archive, Part 1/3, ORCA 5.0.3, Linux, x86-64, .tar.xz Archive, Part 2/3 and ORCA 5.0.3, Linux, x86-64, .tar.xz Archive, Part 3/3. These are separate downloaded tar files. 

5. Extract all three parts and copy all the contents into the folder named "orca."

6. Move the entire content to the home folder.

7. To assign the path variable and source it, open the bashrc file (vi ~/.bashrc) and add the following lines:

```bash
    #ORCA
    export PATH="$HOME/orca:$PATH"
    export LD_LIBRARY_PATH="$HOME/orca:$LD_LIBRARY_PATH"
```

8. Source the bashrc file:

```bash
    source ~/.bashrc
```

9. Run the orca using the following command by typing "orca" in the terminal. The expected outcome for a successful installation will be similar to the following:

```bash
    "This program requires the name of a parameter file as an argument 
    For example, ORCA TEST.INP"
```


Section IB: Installing OpenMPI
********************** 

1. Go to https://www.open-mpi.org/ and select Downloads.

2. Download the openmpi-4.1.1 release with this link: https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.bz2 (It *must* be this version)

3. Extract the file and rename the folder as "openmpi". Move this folder to the home directory. 

4. Go to the openmpi folder in the home directory. Open the terminal and execute the following command in the terminal:

```bash
    ./configure --prefix=$HOME/openmpi
    make all
    make install
```
8. To assign the path variable and source it, open the bashrc file (vi ~/.bashrc) and add the following lines:

```bash
    #OPENMPI
    export PATH=$HOME/openmpi/bin:$PATH
    export LD_LIBRARY_PATH="$HOME/openmpi/lib:$LD_LIBRARY_PATH"
```

9. Source the bashrc file:

```bash
    source ~/.bashrc
```


Section IC: Installing XTB
**********************

1. Go to https://github.com/grimme-lab/xtb/releases

2. Select the xtb version 6.5.1 or go to https://github.com/grimme-lab/xtb/releases/tag/v6.5.1 

3. Download the xtb tar file, xtb-6.5.1-linux-x86_64.tar.xz, and extract the file.

4. After extracting, the folder is named xtb-6.5.1-linux-x86_64.

5. Go to the folder, get into xtb-6.5.1/bin, copy the xtb executable to the orca folder in the home, and rename it as otool_xtb.


Section II: Installing SEEKR2 (optional)
**********************
1.  Make sure [anaconda3](https://www.anaconda.com/) is installed on the local machine. Go to the  [download](https://www.anaconda.com/products/individual) page of anaconda3 and install the latest version of anaconda3.

2. Create a new conda environment with the following commands in the terminal:

```bash
conda create -n SEEKR python=3.8 # Create a new conda environment
conda activate SEEKR # Activate the new environment
conda install git # Install git
pip install Cython # Install dependency for SEEKR2 installation
```
3. Once the conda environment is activated, install the SEEKR2 OpenMM Plugin:

```bash
conda install -c conda-forge seekr2_openmm_plugin cudatoolkit=10.2 openmm=7.7.0
```
4. Make sure you are in the home directory (for convenience), and install 
SEEKR2 with the following commands in the terminal: 

```bash
cd
git clone https://github.com/seekrcentral/seekr2.git
cd seekr2
python setup.py install
python setup.py test # optional tests
```
## Installation and Setup Instructions :

* Activate the previously created conda environment:
```bash
conda activate SEEKR # activate the conda environment
conda install -c conda-forge ambertools biopandas pandas matplotlib parmed regex openmm
pip install PyPDF2
```
* Clone the *qmrebind* repository :

```bash
git clone https://github.com/anandojha/qmrebind.git
```
* Perform the following steps to get this package installed quickly on a local 
Linux machine (Installation in the home directory is recommended) : 

```bash
cd qmrebind
python setup.py install
python setup.py test # optional
```
Detailed documentation can be found [here](https://qmmmrebind-seekr.readthedocs.io/en/latest/index.html).


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

### Copyright

Copyright (c) 2022, Anupam Anand Ojha


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
