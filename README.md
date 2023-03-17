# Molecular docking with aminodeoxychorismate synthase component 1 (PabB), from *B. subtilis*

PabB is a bacterial, fungal and plant protein involved in the biosynthesis of folate. It has recently emerged as a potential drug target. In this set of experiments, the molecular structure of PabB from Bacillus subtilis 168 is used a receptor for high throughput screening using the molecular docking tool smina, a branch of AutoDock Vina. The input is a list of smiles, and the output is a list of smiles with their associated docking scores (minimized affinity). The aim is to find novel inhibitors of this enzyme, with a view to identifying potential antibiotic leads for preclinical testing and development.


## Conda and python

Whilst not essential, the use of conda environments is recommended. Conda helps with the installation and management of packages, and the use of environments allows specific packages to be run with specific versions of python.

To download miniconda (anaconda is also available) on Ubuntu, the following command can be used:

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```


Python will be automatically downloaded with conda.

Next, we will set up a conda environment, wihtin which we will download all of the packages necessary to perform the docking. The docking experiments have been tested with the latest version of python (3.11).

```
conda create -n docking-env python=3.11
```

To activate the environment for the installation of the packages in the next section;

```
conda activate docking-env
```


## Installation

The pipeline requires several packages to manipulate the molecules for docking. These are; RDkit, standardiser, and openbabel. For installation, ensure the conda environment is activated.

To install RDkit;

```
conda install -c conda-forge rdkit
```

To install openbabel on ubuntu (there will be a different installation process on Windows and Mac);

```
sudo apt-get install openbabel
```

To install standardiser;

```
conda install -c conda-forge standardiser
```

The repository, containing the smina executable file and the prepared protein, can be accessed by cloning directly from GitHub;

```
git clone https://github.com/ersilia-os/pabb-docking
```


## Usage

The docking function is designed to be run from the CLI with 3 inputs; path to scorer.py (the function), path to input.csv and path to output folder.


### CLI command

The function can be called from the command line with the following arguments as follows;

```
python PATH_TO_scorer.py PATH_TO_INPUT.csv PATH_TO_OUTPUT_FOLDER
```


### Input

The function is designed to work with .csv files containing a list of smiles. Whilst the smiles can be obtained from any source (databases, designed synthetically, from ML algorithms etc) it is important that the file contains 2 columns only, an 'ID' column (1) and a 'SMILES' column (2). Ensure 'SMILES' is capitalised. Ensure duplicates are removed and any rows containing an empty smiles cell are also removed. This can be easily done in excel.

Save the input file somewhere sensible.


### Output folder

The function will produce several files and store them in this folder. It is important to note that this folder must not already exist before calling the function - by defining the path to the folder the function will create it, which will subsequently be filled with the outputs of the docking.

Ensure this folder is created somewhere sensible, I suggest in a project folder along with the input.csv.


### Outputs

The scorer function returns 3 outputs; a docking-poses.sdf file, with the coordinates of the docked compound that can be viewed in pymol; a log.txt file; and an output.csv file, containing a list of the smiles and their associated ID and docking score. The top docking score only is returned for each input.


## Notes on smina and 'for developers'

Smina is a molecular docking tool adapted from AutoDock Vina, the most successful open-source docking tool with over 3000 citations in the literature. Smina incorporates numerous alterations to ADV, the most notable being a more adaptable scoring function and the ability to receive multi-molecule .sdf files as input. More information on smina can be found here https://sourceforge.net/projects/smina/ and here https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3726561/.

Within smina, we have adopted the vinardo scoring function, a force-field scoring function (vina default is an empirical scoring function) trained on newer datasets. More information on vinardo can be found here https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0155183. Exhaustiveness has been set to 8 (default).

Alteration of the docking parameters is possible and can be performed in the scorer.py file, although this is only recommended for those comfortable with docking experiments and python programming.

Another possible alteration may be the use of a different protein as the docking receptor. Again this is possible, if there is a protein of interest. The protein must be properly prepared using AutoDock Tools (either GUI or CLI - information here https://ccsb.scripps.edu/). The path to the protein can then be added in the scorer.py file, although again this is only recommended for those with a knowledge of python programming.

## About us

Learn about the [Ersilia Open Source Initiative](https://www.ersilia.io/)!

This function was developed in collaboration with Dr Lynden Rooms from the University of Bristol.