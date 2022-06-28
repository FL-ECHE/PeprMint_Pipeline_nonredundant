# PePrMInt Dataset

```diff
- Readme still  Still WIP 
```

This repository contains everything needed to generate the peprmint dataset. 
It is fully automated and should work on linux/mac and windows.

# Motivation
Membrane proteins (integral and peripheral) represent 20% to 30% of the sequenced genomes. Peripheral proteins are only one type of membrane proteins that only bind transiently to the surface of biological membranes where they accomplish their functions which can be a key step in signaling cascades, lipid transport or lipid processing to name a few. Membrane proteins are the target of 50% of drugs. Since the hydrophobic effect is important for drug binding and that the hydrophobic interactions are also for positioning proteins in the membrane, it is highly important to know the localization of the membrane-binding site. The current textbook model of peripheral membrane binding sites involves a combination of basic and hydrophobic amino acids but has recently been disputed. But unlike protein-protein or protein-DNA interfaces, interfacial binding sites of peripheral proteins are surprisingly poorly characterized.


# Objectifs
Our main objective is to contribute to the update of the current model for peripheral membrane binding by revealing structural patterns and amino acid distribution at membrane-binding sites of peripheral proteins. We will develop a statistical framework that we will use to perform a bioinformatics survey of available structure and sequence datasets of peripheral proteins consisting of:
1. Establishing curated datasets of peripheral proteins structures and sequence to test additional features that will inform about other aspects of the binding mechanism
2. The development and the application of a model for hydrophobic protrusions. A model for hydrophobic protrusions was already developed (in the team) and shows that it can correctly identify the membrane-binding site of five prototypical membrane-binding (for a small number of domains). We will further develop this model adding new features relevant to membrane-binding mechanisms such as information about amino acids (e.g. lysines and arginines) neighboring the hydrophobic protrusions, investigate whether these improve discrimination between peripheral and control sets.
3. Analysis of amino acid propensities at the interfacial binding site of peripheral membrane protein. This step aims to map the repertoire of amino acids that Nature has used at the surface of peripheral proteins

## How-to?
CATH (http://www.cathdb.info/) contains 3D protein domains classified into superfamilies.  
PROSITE (https://prosite.expasy.org/) also contains proteins domains classified into profiles, but mostly 1D (only the sequence). Moreover, we can find a multiple sequence alignment per domains with sequences that match with 3D structures (from CATH) and 1D sequences without structures.  
This is what we want: there is not a lot of structures, but we can enrich the dataset alphafold models.
Then we will be able to select a part of a protein (a loop, a helix, or just an amino acid) and we will be able to make statistics on more structures.


# Methods
the full process can be divided into 3 phases, preparation, generation, analysis  

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
 

## 
# Dependencies
The easiest way to have a proper install of the dependencies is to create a new Conda environment with jupyter lab (or notebook) installed in the (base) environment.
Installation with an environment.

1. Install nb_conda_kernels to use your environment in conda (and only this because mkdssp will not be loaded in the PATH without this solution). Do that in "BASE"  
`conda deactivate
jupyter labextension install @jupyter-widgets/jupyterlab-manager
jupyter nbextension enable --py widgetsnbextension`

2. Create an environment for this PePrMInt dataset creation and use its kernel  
`conda create --name peprmint -c conda-forge -c salilab -c hydroid seaborn pandas biopandas biopython ipywidgets scipy numpy tqdm xmltodict pathlib dssp pytables termcolor requests bs4 lxml wordcloud nglview mdanalysis weasyprint openpyxl`

3. for progressbar usage (in base):  
```bash 
jupyter labextension install @jupyter-widgets/jupyterlab-manager 
jupyter nbextension enable --py widgetsnbextension```

3. Don't forget to choose the PePrMInt environment for the kernel "`Python [conda env:peprmint]`"

4. LINUX AND MAC ONLY:  
if you want some speedup sometimes, install pandarallel `pip install pandarallel`

# Installation
To install the package for the development version, use the command `pip install -e .` in the base folder.


# Usage
Some notebook start with a number, they are mandatory to run first to generate and download dataset.  
** THUS, it is mandatory to re-run 02 after 03 (downloading of alphafold structure) **
 - `00-SETUP.ipynb` is the configuration file. Change on this file the output folders, add some domains if you want.
   - You should only change `PEPRMINT_FOLDER`. other folders depend on this one.
   - _How do add new domains ?_
     - add the prosite ID of the domain in the `DOMAIN_PROSITE` dictionary. If you have one prositeID just add it like `"SH2": "PS50001"`, if you have a list:  "SH2": ["PS50001","PSXXXXX"]`
     - specify the cath domain in `DOMAIN_CATH` like `"PH": "2.30.29.30",`
 - `01-downloads_files.ipynb` will download all CATH pdbs and prosite alignment.
 - `02-Create_PePrMInt_Dataset.ipynb` will merge CATH and prosite data to create the PePrMInt dataset
 - `03-downloads_alphafold_structures.ipynb` is the notebook that will download alphafold structure
 - `04-04-GERENATE_FIGURE_FOR_ARTICLE.ipynb.ipynb` is the notebook used to generate figure for our main paper
 - `05-generate_CSV_fullraw_dataset.ipynb` is only to generate the full dataset in CSV foramt.

 # ressources

- In `Ressources/datasets` there are 2 datasets : 
  - `PePr2DS.csv.zip` that contains our dataset from CATH with all the features computes for each amino acids
  - `fileS2.csv/xlsx` that contains the computed feature for CATH structure and AF models summed for each proteins.

- `Ressources/structures` Contains all the CATH  and AF models per domains. : 





