[![DOI](https://zenodo.org/badge/495005184.svg)](https://zenodo.org/badge/latestdoi/495005184)

# PePrMInt Dataset


This repository contains everything needed to generate the peprmint dataset. 
It is fully automated and should work on linux/mac and windows.

# Motivation
Membrane proteins (integral and peripheral) represent 20% to 30% of the sequenced genomes. Peripheral proteins are only one type of membrane proteins that only bind transiently to the surface of biological membranes where they accomplish their functions which can be a key step in signaling cascades, lipid transport or lipid processing to name a few. Membrane proteins are the target of 50% of drugs. Since the hydrophobic effect is important for drug binding and that the hydrophobic interactions are also for positioning proteins in the membrane, it is highly important to know the localization of the membrane-binding site. The current textbook model of peripheral membrane binding sites involves a combination of basic and hydrophobic amino acids but has recently been disputed. But unlike protein-protein or protein-DNA interfaces, interfacial binding sites of peripheral proteins are surprisingly poorly characterized.


# Objectives
Our main objective is to contribute to the update of the current model for peripheral membrane binding by revealing structural patterns and amino acid distribution at membrane-binding sites of peripheral proteins. We will develop a statistical framework that we will use to perform a bioinformatics survey of available structure and sequence datasets of peripheral proteins consisting of:
1. Establishing curated datasets of peripheral proteins structures and sequence to test additional features that will inform about other aspects of the binding mechanism
2. The development and the application of a model for hydrophobic protrusions. A model for hydrophobic protrusions was already developed (in the team) and shows that it can correctly identify the membrane-binding site of five prototypical membrane-binding (for a small number of domains). We will further develop this model adding new features relevant to membrane-binding mechanisms such as information about amino acids (e.g. lysines and arginines) neighboring the hydrophobic protrusions, investigate whether these improve discrimination between peripheral and control sets.
3. Analysis of amino acid propensities at the interfacial binding site of peripheral membrane protein. This step aims to map the repertoire of amino acids that Nature has used at the surface of peripheral proteins


# Methods

CATH (http://www.cathdb.info/) contains 3D protein domains classified into superfamilies.  
PROSITE (https://prosite.expasy.org/) also contains proteins domains classified into profiles, but mostly 1D (only the sequence). Moreover, we can find a multiple sequence alignment per domains with sequences that match with 3D structures (from CATH) and 1D sequences without structures.  
This is what we want: there is not a lot of structures, but we can enrich the dataset alphafold models.
Then we will be able to select a part of a protein (a loop, a helix, or just an amino acid) and we will be able to make statistics on more structures.

NB! All the CATH and AF models per domains are under `resources/structures`, while `resources/datasets` contains two datasets: 
  - `PePr2DS.csv.zip` that contains our dataset from CATH with all the features computes for each amino acids.
  - `fileS2.csv/xlsx` that contains the computed feature for CATH structure and AF models summed for each proteins.



The full process is divided into 3 phases: preparation, generation, and analysis.

## Preparation 
 1. All the necessary folders are first created
 2. All the alignment files from PROSITE are downloaded
 3. All the CATH PDBs for the selected domains are downloaded.
 
## Dataset generation
 1. All CATHpdbs are first cleaned to be more compatible with DSSP
 2. DSSP is used (called through the biopython package) to generate the secondary structures.
 3. The protrusions and the neighbor's list is then calculated
 4. For 3D Structures, the CATH cluster number is added
 5. Then a first mapping with UniProt is made to get the uniprot_acc and Uniprot_id which will allow us to map with the prosite alignment.
 6. A mapping of the 3D structure is made with the prosite alignment. Here it's a hard task because you have to map the residue number with the alignment position. This can be tricky since in the structure you can have some mutations.
 7. All other 1D sequences are added (all the sequences that do not have structures)
 8. A second mapping is done to get the uniref cluster's representative.


# Installation
The following assumes only that the current machine includes a recent 
installation of the Anaconda distribution (https://www.anaconda.com/).

## Environment
Conda is a well-known package and programming environment manager, allowing 
one to create a controlled setting with respect to installed libraries and 
their versions (a major goal is to avoid the "dependency hell"). Conda is 
currently being replaced by Mamba, especially for scientific packaging in 
Python. For all purposes of this project, consider mamba a thin, more efficient
layer on top of conda.

To download and install the latest version of mamba on a unix-like platform, 
open a terminal and run (we are asked to accept the license agreement, accept or 
suggest an install location, and accept to initialize the distribution)

```diff
 - ⚠️ ONLY IF YOU DON'T HAVE A PREVIOUS ANACONDA OR MINICONDA INSTALL! ⚠️
```
You can install it with the following command : `bash install_miniconda.sh`.  
also, if you have `conda` but not `mamba` installed you can install it with: `conda install -c conda-forge mamba -y`


When the process is over, close the current terminal and open a new one. Make 
sure that the installation was succesful by seeing a nice output in the terminal
after running

```bash
mamba info --envs
```

Optional: if you want your standard terminal to load the base conda/mamba environment
by default, enter
```bash
conda config --set auto_activate_base true
```


## Setup

We are now ready to setup the library dependencies to create the Peprmint dataset. First, enter the following to ensure that packages concerning the notebook support are in place.
```bash
mamba install -c conda-forge nb_conda_kernels nodejs ipywidgets jupyterlab
```

After downloading/cloning this repository (e.g. with `git clone https://github.com/reuter-group/tubiana_etal_2022.git`), one may create the adequate 
environment for this project either within jupyter (choose the file 
`environment.yml`) or directly in the terminal: navigate (`cd tubiana_etal_2022`) to the 
corresponding path and run

```bash
mamba env create -f environment.yml
```

Now activate the environement and install the package to make it visible everywhere in the system: 
```bash
mamba activate peprmint
pip install -e .
```


Have fun! (:


# Usage (under construction)

[Check again soon... currently preparing usage instructions!]

1. Adjust desired settings in the `peprmint_default.config` file (optional). 
Alternatively, deleting that file makes the execution fall back to factory 
settings, and populate datasets under folder `data` in the current directory.

2. Run `python main.py` in the terminal, or open `main-nb-example` for 
a notebook usage example.

The user should be able to choose which steps of the computation are actually 
executed, so that it is possible to interact with intermediate results to 
perform different analyses.

