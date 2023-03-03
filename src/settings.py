#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# MIT License
# 
# Copyright (c) 2022 Reuter Group
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


""" Global settings (previously on Notebook #0)

__author__ = ["Thibault Tubiana", "Phillippe Samer"]
__organization__ = "Computational Biology Unit, Universitetet i Bergen"
__copyright__ = "Copyright (c) 2022 Reuter Group"
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Phillippe Samer"
__email__ = "samer@uib.no"
__status__ = "Prototype"
"""

import os
from typing import Optional
import configparser

import pandas as pd
from tqdm.auto import trange, tqdm

from src.notebook_handle import NotebookHandle

class Settings:

    def __init__(self, path: Optional[str] = None):
        self._DEFAULT_CONFIG_FILE = '{0}/{1}'.format(os.getcwd(),
                                                     'peprmint_default.config')
        self.USING_NOTEBOOK = self._is_notebook()
        if self.USING_NOTEBOOK:
            self.NOTEBOOK_HANDLE = NotebookHandle()
        else:
            self.NOTEBOOK_HANDLE = None

        self._run_config(path)
        self._libs_setup()

        # create directory structure for peprmint
        cwd = self.config_file['GENERAL']['working_folder_path']
        try:
            print('Working directory: {0}'.format(cwd))
            os.makedirs(cwd)
            self.FORMER_WORKING_DIR = False
        except FileExistsError:
            self.FORMER_WORKING_DIR = True

        self.PEPRMINT_FOLDER = cwd
        self.SETUP = {}   # dictionary with ALL parameters
        self.define_folders()
        self.create_directories()
        self.map_cath_and_prosite()

    def _is_notebook(self) -> bool:
        try:
            shell = get_ipython().__class__.__name__
            if shell == 'ZMQInteractiveShell':
                return True   # jupyter notebook
            elif shell == 'TerminalInteractiveShell':
                return False  # terminal running IPython
            else:
                return False  # other type (?)
        except NameError:
            return False      # probably standard Python interpreter

    def _run_config(self, path: Optional[str] = None):
        """ Read configuration file
        First, try to read the user configurations, either at the standard file
        or at the location given in the optional constructor parameter.
        If neither of these work, genereate the original config file.
        """
        config_ready = False

        if path is not None:
            config_ready = self._read_config_file(path)
            if config_ready:
                print("Using configuration file '{0}'".format(path))
            else:
                print("Could not find configuration file '{0}'".format(path))

        if not config_ready:
            print("Reading standard configuration file... ", end = '')
            config_ready = self._read_config_file(self._DEFAULT_CONFIG_FILE)
            if config_ready:
                print("done")
            else:
                print("not found")
                print("Using factory configuration, saved locally at '{0}'".format(self._DEFAULT_CONFIG_FILE))
                self._write_default_config_file(self._DEFAULT_CONFIG_FILE)

    def _read_config_file(self, path: str) -> bool:
        if not os.path.isfile(path):
            return False

        # TO DO: maybe check if given file indeed has all of the expected sections and fields?
        self.config_file = configparser.ConfigParser()
        self.config_file.read(path)
        return True

    def _write_default_config_file(self, path: str):
        self.config_file = configparser.ConfigParser()

        self.config_file['GENERAL'] = {}
        # default folder: a new one in the current working directory
        self.config_file['GENERAL']['working_folder_path'] = '{0}/data'.format(os.getcwd())

        self.config_file['CATH'] = {}
        self.config_file['CATH']['version'] = 'v4_2_0'
        self.config_file['CATH']['domain_list_url'] = "http://download.cathdb.info/cath/releases/all-releases/{0}/cath-classification-data/cath-domain-list-{0}.txt".format(self.config_file['CATH']['version'])
        self.config_file['CATH']['fetch_pdb_url'] = "http://www.cathdb.info/version/{0}/api/rest/id/".format(self.config_file['CATH']['version'])

        self.config_file['UNIPROT'] = {}
        self.config_file['UNIPROT']['url'] = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz"

        self.config_file['PROSITE'] = {}
        self.config_file['PROSITE']['url'] = "ftp://ftp.expasy.org/databases/prosite/prosite_alignments.tar.gz"

        self.config_file['DATASET_MANAGER'] = {}
        self.config_file['DATASET_MANAGER']['recalculate'] = str(False)

        self.config_file['ALPHAFOLD_UTILS'] = {}
        self.config_file['ALPHAFOLD_UTILS']['rebuild'] = str(True)
        self.config_file['ALPHAFOLD_UTILS']['use_uniprot_boundaries'] = str(True)
        self.config_file['ALPHAFOLD_UTILS']['use_all_AFmodels'] = str(True)
        self.config_file['ALPHAFOLD_UTILS']['AF_pdbs_url_prefix'] = "https://alphafold.ebi.ac.uk/files/AF-"
        self.config_file['ALPHAFOLD_UTILS']['AF_pdbs_url_suffix'] = "-F1-model_v1.pdb"
        self.config_file['ALPHAFOLD_UTILS']['AF_interpro_url_prefix'] = "https://www.ebi.ac.uk/interpro/api/entry/"
        self.config_file['ALPHAFOLD_UTILS']['AF_interpro_url_middle'] = "/protein/reviewed/"
        self.config_file['ALPHAFOLD_UTILS']['AF_interpro_url_suffix'] = "/?page_size=200"
        
        with open(path, 'w') as configfile:
            self.config_file.write(configfile)

    def _libs_setup(self):
        # Place any additional settings for imported libraries here
        tqdm.pandas()   # activate tqdm progressbar for pandas

    def define_folders(self):
        self.WORKDIR = f"{self.PEPRMINT_FOLDER}/dataset/"
        self.CATHFOLDER = f"{self.PEPRMINT_FOLDER}/databases/cath/"
        self.ALPHAFOLDFOLDER = f"{self.PEPRMINT_FOLDER}/databases/alphafold/"
        self.PROSITEFOLDER = f"{self.PEPRMINT_FOLDER}/databases/prosite/"
        self.UNIPROTFOLDER = f"{self.PEPRMINT_FOLDER}/databases/uniprot/"
        self.FIGURESFOLDER = f"{self.PEPRMINT_FOLDER}/figures/"

        self.SETUP["PEPRMINT_FOLDER"] = self.PEPRMINT_FOLDER
        self.SETUP["WORKDIR"] = self.WORKDIR
        self.SETUP["CATHFOLDER"] = self.CATHFOLDER
        self.SETUP["PROSITEFOLDER"] = self.PROSITEFOLDER
        self.SETUP["ALPHAFOLDFOLDER"] = self.ALPHAFOLDFOLDER
        self.SETUP["UNIPROTFOLDER"] = self.UNIPROTFOLDER
        self.SETUP["FIGURESFOLDER"] = self.FIGURESFOLDER

        for k in self.SETUP:
            exec(f"self.{k}2 = self.SETUP['{k}']")

    def create_directories(self):
        if not os.path.exists(self.PEPRMINT_FOLDER):
            os.makedirs(self.PEPRMINT_FOLDER)
        if not os.path.exists(self.WORKDIR):
            os.makedirs(self.WORKDIR)
        if not os.path.exists(self.FIGURESFOLDER):
            os.makedirs(self.FIGURESFOLDER)
        if not os.path.exists(self.ALPHAFOLDFOLDER):
            os.makedirs(self.ALPHAFOLDFOLDER)
        if not os.path.exists(self.UNIPROTFOLDER):
            os.makedirs(self.UNIPROTFOLDER)
        if not os.path.exists(self.PROSITEFOLDER):
            #MSA will contains the alignments in "msa" format (FASTA)
            os.makedirs(self.PROSITEFOLDER)
        if not os.path.exists(self.CATHFOLDER):
            os.makedirs(self.CATHFOLDER)

    def map_cath_and_prosite(self):
        self.CATHVERSION = self.config_file['CATH']['version']
        self.DOMAIN_PROSITE = {
            "PH": "PS50003",
            "C2": ["PS50004","PS51547"],
            "C1": "PS50081",  # Note : no C1 prosite on SMART but 2 C1 ProSite on Interprot (PS50081,PS00479), I took PS50081 since the data in PS00479 are in PS50081.
            "PX": "PS50195",
            # "FYVE":"PS50178",
            "FYVE": ["PS50178",'PS50089', 'PS00518','PS50016','PS01359','PS50014','PS00633','PS50119'],  # FYVE CAN BE THIS ONE TOO....
            # "PPASE_MYOTUBULARIN":"PS51339",# no GRAM domain found on prosite. Has to do this manually. Go on http://smart.embl-heidelberg.de/smart/do_annotation.pl?DOMAIN=GRAM&BLAST=DUMMY
            "BAR": "PS51021",  # 1URU is missing on prosite
            # "GLA":"PS50963",
            "ENTH": "PS50942",
            "SH2": "PS50001",
            "SEC14": "PS50191",
            "START": "PS50848",
            "C2DIS":"PS50022",
            "GLA": "PS50998",
            "PLD":"PS50035",
            "PLA":"PS00118",
            "ANNEXIN":"PS00223",
        }
        # Invert keys and values to have PROSITEID ==> DOMAIN
        self.PROSITE_DOMAIN = {}
        for key, value in self.DOMAIN_PROSITE.items():
            if type(value) == type([]):
                for subvalues in value:
                    self.PROSITE_DOMAIN[subvalues] = key
            else:
                self.PROSITE_DOMAIN[value] = key
        # self.PROSITE_DOMAIN = {v: k for k, v in self.DOMAIN_PROSITE.items()}

        self.DOMAIN_CATH = {
            "PH": "2.30.29.30",
            "C2": "2.60.40.150",
            "C1": "3.30.60.20",
            "PX": "3.30.1520.10",
            "FYVE": "3.30.40.10",
            "BAR": "1.20.1270.60",
            "ENTH": "1.25.40.90",
            "SH2": "3.30.505.10",
            "SEC14": "3.40.525.10",
            "START": "3.30.530.20",
            "C2DIS": "2.60.120.260",
            "GLA":"2.40.20.10",
            "PLD":"3.20.20.190",
            "PLA":"1.20.90.10",
            "ANNEXIN":"1.10.220.10",
        }

        self.DOMAIN_INTERPRO = {
            "PH": "SSF50729",
            "C2": "SSF49562",
            "C1": None,
            "PX": "SSF64268",
            "FYVE": "SSF57903", #badly classified it looks like...
            "BAR": "SSF103657",
            "ENTH": "SSF48464",
            "SH2": "SSF55550",
            "SEC14": ["SSF52087","SSF46938"], #the CRAL TRIO domain is truncated in SSF.
            "START": "SSF55961",
            "C2DIS": "SSF49785",
            "GLA":None,
            "PLD":"SSF51695",
            "PLA":"G3DSA:1.20.90.10",
            "ANNEXIN":"SSF47874",
        }

        self.DOMAIN_INTERPRO_REFINE = {
            "PH": True,
            "C2": False,
            "C1": False,
            "PX": True,
            "FYVE": False,
            "BAR": False,
            "ENTH": False,
            "SH2": False,
            "SEC14": False,
            "START": True,
            "C2DIS": False,
            "GLA":False,
            "PLD":False,
            "PLA":True,
            "ANNEXIN":False,
        }

        # Invert keys and values to have CATHID ==> DOMAIN
        self.CATH_DOMAIN = {v: k for k, v in self.DOMAIN_CATH.items()}
        self.SUPERFAMILY = self.CATH_DOMAIN
        self.SETUP["DOMAIN_PROSITE"] = self.DOMAIN_PROSITE
        self.SETUP["PROSITE_DOMAIN"] = self.PROSITE_DOMAIN
        self.SETUP["DOMAIN_CATH"] = self.DOMAIN_CATH
        self.SETUP["CATH_DOMAIN"] = self.CATH_DOMAIN
        self.SETUP["SUPERFAMILY"] = self.SUPERFAMILY

    # NB! function 'selectUniquePerCluster(...)' from notebook #0
    # is redefined (and used) on notebook 02, so we are moving it
