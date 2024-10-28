Installation
============

At this time, QMrebind has only been tested on Linux systems. Therefore, all
installation instructions are for Linux only.

Install Conda
-------------

If you do not already have Conda, it can be easily installed by completing the
following steps:

Download Conda, run the script, and fill out the prompts:

.. code-block:: bash

  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh

Make sure Conda is installed by running:

.. code-block:: bash

  which conda

You will want to use Python 3.8, so you can see which version you are with
the command:

.. code-block:: bash

  python -V

If it says any other version besides Python 3.8, then enter:

.. code-block:: bash

  conda install python=3.8

If you want you can create a conda environment, 

.. code-block:: bash

  conda create --name QMMM python=3.8

but you can also just install all packages straight to the base environment
if you wish to. If using an environment, whenever you're installing or running 
anything involving QMrebind, make sure that you have activated your 
environment by running ``conda activate QMMM``.

Install Orca
------------

Make sure to install ORCA before running QMrebind:

Follow the instructions below or visit this page to install ORCA: 
https://www.orcasoftware.de/tutorials_orca/first_steps/install.html

1. Go to https://orcaforum.kofo.mpg.de/ucp.php?mode=login. Create an account 
to log in with a username and a password. 

2. Go to Downloads.

3. Select ORCA 5.0.3

4. The ORCA tar files will be downloaded in three parts: Download ORCA 5.0.3, 
Linux, x86-64, .tar.xz Archive, Part 1/3, ORCA 5.0.3, Linux, x86-64, .tar.xz 
Archive, Part 2/3 and ORCA 5.0.3, Linux, x86-64, .tar.xz Archive, Part 3/3. 
These are separate downloaded tar files. 

5. Extract all three parts and copy all the contents into the folder named 
"orca".

.. code-block:: bash

    tar -xf orca_5_0_3_linux_x86-64_openmpi411_part1.tar.xz
    tar -xf orca_5_0_3_linux_x86-64_openmpi411_part2.tar.xz
    tar -xf orca_5_0_3_linux_x86-64_openmpi411_part3.tar.xz
    mv orca_5_0_3_linux_x86-64_openmpi411_part1/* .
    mv orca_5_0_3_linux_x86-64_openmpi411_part2/* .
    mv orca_5_0_3_linux_x86-64_openmpi411_part3/* .


6. To assign the path variable and source it, open the bashrc file 
(vi ~/.bashrc) and add the following lines:

.. code-block:: bash

    export PATH="$HOME/orca:$PATH"
    export LD_LIBRARY_PATH="$HOME/orca:$LD_LIBRARY_PATH"

7. Source the bashrc file:

.. code-block:: bash

  source ~/.bashrc
  
The correct environment probably deactivated:

.. code-block:: bash

  conda activate QMMM

8. Run ORCA using the following command by typing "orca" in the terminal. 
The expected outcome for a successful installation will be similar to the 
following:

.. code-block::

    This program requires the name of a parameter file as an argument 
    For example, ORCA TEST.INP


Install OpenMPI
---------------

1. Go to https://www.open-mpi.org/ and select Downloads.

2. Download the openmpi-4.1.1 release with this link: 
https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.bz2 
(It *must* be this version)

3. Extract the file and rename the folder as "openmpi". Enter this directory.

.. code-block::

  tar -xzf openmpi-4.1.1.tar.gz
  mv openmpi-4.1.1 $HOME/openmpi
  cd $HOME/openmpi

4. Go into to the openmpi folder. Open the terminal and 
execute the following command in the terminal:

.. code-block::

  ./configure --prefix=$HOME/openmpi
  make all
  make install

8. To assign the path variable and source it, open the bashrc file (vi ~/.bashrc) and add the following lines:

.. code-block::

  export PATH=$HOME/openmpi/bin:$PATH
  export LD_LIBRARY_PATH="$HOME/openmpi/lib:$LD_LIBRARY_PATH"

9. Source the bashrc file:

.. code-block::

  source ~/.bashrc

The correct environment probably deactivated:

.. code-block:: bash

  conda activate QMMM

Install XTB
-----------

1. Go to https://github.com/grimme-lab/xtb/releases

2. Select the xtb version 6.5.1 or go to 
https://github.com/grimme-lab/xtb/releases/tag/v6.5.1 

3. Download the xtb tar file, xtb-6.5.1-linux-x86_64.tar.xz, and extract the 
file.

.. code-block:: bash

    tar -xf xtb-6.5.1-linux-x86_64.tar.xz

4. Copy the xtb executable to the 
orca folder in HOME, and rename it as otool_xtb.

.. code-block:: bash

    cp xtb-6.5.1/bin/xtb ~/orca/otool_xtb

Install QMrebind
----------------
1. Activate the previously created conda environment:

.. code-block:: bash

  conda activate QMMM # activate the conda environment
  conda install -c conda-forge ambertools 
  conda install -c conda-forge openmm
  pip install PyPDF2

2. Clone the *qmrebind* repository :

.. code-block:: bash

  git clone https://github.com/seekrcentral/qmrebind.git

3.  Perform the following steps to get this package installed quickly on a local 
Linux machine (Installation in the home directory is recommended) : 

.. code-block:: bash

  cd qmrebind
  python -m pip install .
  pytest # optional

